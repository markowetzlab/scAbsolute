# Copyright 2022, Michael Schneider, All rights reserved.
import numpy as np
import os
import sys
import logging
import copy
import pickle
from datetime import datetime
#  startTime = datetime.now()
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
os.environ['TF_CPP_MAX_VLOG_LEVEL'] = '3'
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
NUM_THREADS=1
import tensorflow.compat.v2 as tf
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)
tf.config.threading.set_intra_op_parallelism_threads(NUM_THREADS)
tf.config.threading.set_inter_op_parallelism_threads(NUM_THREADS)
tf.enable_v2_behavior()
import tensorflow_probability as tfp
from tensorflow_probability import distributions as tfd
tfb = tfp.bijectors
tfd = tfp.distributions
logging.debug("Num GPUs Available: ", len(tf.config.list_physical_devices('GPU')))


# multiprocessing with better support for pickling objects
#  from pathos.multiprocessing import ProcessingPool as Pool
#  from multiprocessing import Pool
from functools import partial

# Set up logging.
#  stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
#  logdir = 'debug/func/%s' % stamp
#  writer = tf.summary.create_file_writer(logdir)

def simulate_data(true_rpc, true_zero, true_alpha, cov_info=None, n_size=1000):
  import scipy.stats
  #  covariate_info = np.ones(num_steps)

  state_means = [i*true_rpc for i in [1, 3, 4, 2, 4, 3, 4, 5]] + [true_zero]
  true_durations = [int(i*n_size) for i in [1, 1, 1, 2, 1, 1, 1, 1, 1]]
  true_states = np.repeat(state_means, true_durations)

  num_steps = sum(true_durations)
  #  if cov_info is None:
    #  covariate_info=np.exp(np.random.normal(loc=0.0, scale=0.2, size=num_steps))
  #  else:
  covariate_info = cov_info[0,:]

  def convert_params(mean, alpha):
      """
      Convert mean/dispersion parameterization of a negative binomial to the ones scipy supports

      See https://mathworld.wolfram.com/NegativeBinomialDistribution.html
      """
      #  alpha = 1./theta
      var = mean * (1. + alpha * mean)
      p = mean/var
      n = mean*p/(1.0 - p)
      return n, p

  true_means = [[i] * j for i,j in zip(state_means, true_durations)]
  true_means = [item for sublist in true_means for item in sublist]
  true_means = true_means * covariate_info

  true_means = np.where(true_means <= 0.1, 0.1, true_means)
  true_r, true_p = zip(*[convert_params(i, true_alpha) for i in true_means])

  observed_counts = np.concatenate([
    #  scipy.stats.nbinom(r, p).rvs(num_steps)
      #  for (r, p, num_steps) in zip(true_r, true_p, true_durations)
    scipy.stats.nbinom(r, p).rvs(1)
      for (r, p) in zip(true_r, true_p)
  ]).astype(np.float32)

  observed_counts = np.tile(observed_counts, (1, 1))

  return({"simulated_data" : observed_counts, "groundtruth" :true_states})


def run_hmm(observed_counts, rpc_prior_value, zero_prior, alpha_max_clip, alpha_prior_value, alpha_zero_prior_value, alpha_prior_scale,
            change_prob, chromosome_breakpoints, covariate_info, max_iterations, num_states, verbose, hmm_path,
            optimizer_learning_rate):

  nmax = np.max(observed_counts)
  epsilon = np.max([1e-12, 1e-2 / (nmax**2)])
  check_values = False
  logging.debug("Input data shape:\n{}".format(observed_counts.shape))
  logging.debug("Covariate info shape:\n{}".format(covariate_info.shape))
  logging.info("Epsilon value: {} (max: {:.2e})".format(epsilon, nmax))

  num_steps = len(observed_counts)
  # No info about initial states
  initial_state_logits = np.zeros([num_states], dtype=np.float32) # uniform distribution

  transition_probs = change_prob / (num_states-1) * np.ones(
      [num_states, num_states], dtype=np.float32)
  np.fill_diagonal(transition_probs,
                   1-change_prob)
  breakpoint_change_prob = 0.5
  transition_probs_breakpoint = breakpoint_change_prob / (num_states-1) * np.ones(
      [num_states, num_states], dtype=np.float32)
  np.fill_diagonal(transition_probs,
                   1-breakpoint_change_prob)

  time_varying_transition_distribution = np.repeat(transition_probs[np.newaxis, :, :], num_steps-1, axis=0)
  if(len(chromosome_breakpoints) >= 1):
    time_varying_transition_distribution[chromosome_breakpoints,:,:] = \
      np.repeat(transition_probs_breakpoint[np.newaxis, :, :], len(chromosome_breakpoints), axis=0)

  #  logging.debug("Input data shape:\n{}".format(observed_counts.shape))
  logging.debug("Initial parameters rpc {} alpha {} zero {}".format(rpc_prior_value, alpha_prior_value, zero_prior))
  logging.debug("Transition probs shape:\n{}".format(transition_probs.shape))
  logging.debug("Initial state logits:\n{}".format(initial_state_logits))
  logging.debug("Transition matrix:\n{}".format(transition_probs[0,:]))

  w = tf.constant([np.array(np.arange(1.0, num_states))], dtype=tf.float32)
  covariate_mean = tf.constant(covariate_info, dtype=tf.float32)

  def transform_rpc(rpc):
    # create rpc distance rpc_rate values (i.e. on grid)
    rpc_values = tf.math.multiply(rpc, w)
    ## add covariate correction (applied to mean - multiplicative!)
    # x = c * rpc * correction
    rpc_clean = rpc_values * covariate_mean[:, tf.newaxis]
    return(rpc_clean)

  def transform_zero(zero):
    ## add covariate correction (applied to mean - multiplicative!)
    # x = c * rpc * correction
    zero_clean = zero * covariate_mean
    return(zero_clean[:,tf.newaxis])

  #  bij_alpha =tfb.Sigmoid(low=tf.cast(tf.identity(epsilon), tf.float64), high=tf.cast(0.5, tf.float64))
  # alternatives
  bij_alpha = tfb.Chain([tfb.Shift(tf.cast(epsilon, tf.float64)), tfb.Exp()])
  bij_alpha_zero = tfb.Chain([tfb.Shift(tf.cast(1e-16, tf.float64)), tfb.Exp()])
  #  print("DEBUG")
  #  print(bij_alpha.forward(tf.cast(0.001, tf.float64)))
  #  print("END DEBUG")
  #  sys.exit(3)

  def transform_alpha_zero(alpha_zero_init):
    alpha_zero_init = tf.cast(alpha_zero_init, tf.float64)
    alpha_zero_init = tf.cast(bij_alpha_zero.forward(alpha_zero_init), dtype=tf.float32)
    alpha_zero = tf.broadcast_to(alpha_zero_init, [num_steps, 1])
    return(alpha_zero)

  def transform_alpha(alpha_init):
    # possible bug in tensorflow - when applying bijector
    alpha_init = tf.cast(alpha_init, tf.float64)
    alpha_init = tf.cast(bij_alpha.forward(alpha_init), dtype=tf.float32)
    # casting back and forth seems to resolve the issue

    alpha_init = tf.math.maximum(alpha_init, epsilon)
    alpha = tf.broadcast_to(alpha_init, [num_steps, num_states-1])
    return(alpha)

  #  bij_rpc = tfb.Chain(bijectors=[tfb.Softplus(hinge_softness=1e-3)])
  #  bij_alpha_zero = tfb.Chain([tfb.Shift(shift=epsilon), tfb.Softplus(hinge_softness=epsilon)])
  #  bij_alpha = tfb.Chain([tfb.Shift(shift=epsilon), tfb.Softplus(hinge_softness=epsilon)])

  trainable_rpc_zero = tf.Variable(zero_prior, name="zero", dtype=tf.float32,
                              constraint=lambda x: \
                                   tf.clip_by_value(x, clip_value_min=epsilon, clip_value_max=0.5))
  trainable_alpha_zero = tf.Variable(np.log(alpha_zero_prior_value), name="alpha_zero", dtype=tf.float32)
  #  trainable_alpha_zero = tf.Variable(alpha_zero_prior_value, name="alpha_zero", dtype=tf.float32,
                              #  constraint=lambda x: \
                                   #  tf.clip_by_value(x, clip_value_min=epsilon, clip_value_max=alpha_max_clip))
  trainable_rpc = tf.Variable(rpc_prior_value, name="rpc", dtype=tf.float32)
  trainable_alpha = tf.Variable(np.log(alpha_prior_value), name="alpha", dtype=tf.float32)
  #  trainable_alpha = tf.Variable(alpha_prior_value, name="alpha", dtype=tf.float32,
                              #  constraint=lambda x: \
                                   #  tf.clip_by_value(x, clip_value_min=epsilon, clip_value_max=alpha_max_clip))

  rpcs_full = tfp.util.DeferredTensor(trainable_rpc, transform_rpc, shape=(num_steps, num_states-1), dtype=tf.float32)
  rpcs_zero = tfp.util.DeferredTensor(trainable_rpc_zero, transform_zero, shape=(num_steps, 1), dtype=tf.float32)
  alphas_full = tfp.util.DeferredTensor(trainable_alpha, transform_alpha, shape=(num_steps, num_states-1), dtype=tf.float32)
  alphas_zero = tfp.util.DeferredTensor(trainable_alpha_zero, transform_alpha_zero, shape=(num_steps, 1), dtype=tf.float32)

  def transform_join(rpc1, rpc2, alpha1, alpha2):
    return dict(mean=tf.concat([rpc1, rpc2], axis=1),
                dispersion=tf.concat([alpha1, alpha2], axis=1))

  deferred_dist = tfp.experimental.util.DeferredModule(
    tfd.NegativeBinomial.experimental_from_mean_dispersion,
    args_fn=transform_join,
    rpc1=rpcs_zero,  # May be passed by position or by name.
    rpc2=rpcs_full,
    alpha1=alphas_zero,
    alpha2=alphas_full)

  hmm = tfd.HiddenMarkovModel(
    allow_nan_stats=(not check_values),
    validate_args=check_values,
    initial_distribution=tfd.Categorical(
        logits=initial_state_logits),
    transition_distribution=tfd.Categorical(probs=time_varying_transition_distribution),
    time_varying_transition_distribution=True,
    time_varying_observation_distribution=True,
    observation_distribution=deferred_dist,
    num_steps=num_steps)

  most_probable_states_pre = hmm.posterior_mode(observed_counts)
  bps = (np.where(most_probable_states_pre[:-1] != most_probable_states_pre[1:])[0])
  #  logging.debug("PRE: ", bps[0])

  rpc_prior = tfp.distributions.Normal(
                  loc=rpc_prior_value,
                  scale=0.15*rpc_prior_value, validate_args=check_values,
                  allow_nan_stats=(not check_values), name="rpc_prior")
  # Prior on dispersion parameter alpha
  # https://statmodeling.stat.columbia.edu/2018/04/03/justify-my-love/
  rpc_zero_prior = tfp.distributions.Exponential(rate=10.0, validate_args=check_values, allow_nan_stats=(not check_values))
  alpha_prior = tfp.distributions.Exponential(rate=alpha_prior_scale, validate_args=check_values, allow_nan_stats=(not check_values))
  zero_alpha_prior = tfp.distributions.Exponential(rate=alpha_prior_scale, validate_args=check_values,
                                                   allow_nan_stats=(not check_values))

  def log_prob():
    return ( \
          tf.reduce_sum(rpc_zero_prior.log_prob(trainable_rpc_zero)) +
          tf.reduce_sum(rpc_prior.log_prob(trainable_rpc)) +
          tf.reduce_sum(zero_alpha_prior.log_prob(tf.math.sqrt(
            tf.cast(bij_alpha_zero.forward(tf.cast(trainable_alpha_zero, tf.float64)), tf.float32)))) +
          tf.reduce_sum(alpha_prior.log_prob(tf.math.sqrt(
            tf.cast(bij_alpha.forward(tf.cast(trainable_alpha, tf.float64)), tf.float32)))) +
          hmm.log_prob(observed_counts))

  optimizer = tf.keras.optimizers.Adam(learning_rate=optimizer_learning_rate)

  @tf.function(autograph=True)
  def train_op():
    #  DEBUG tracing/execution -> need to re-create graph multiple times because of same initial alpha/rpc
    #  print("Tracing with lr = ", optimizer_learning_rate)
    #  tf.print("Executing with lr = ", optimizer_learning_rate)
    with tf.GradientTape() as tape:
      neg_log_prob = -log_prob()

    grads = tape.gradient(neg_log_prob, [trainable_rpc_zero,
                                         trainable_rpc,
                                         trainable_alpha_zero,
                                         trainable_alpha])

    optimizer.apply_gradients([(grads[0], trainable_rpc_zero),
                               (grads[1], trainable_rpc),
                               (grads[2], trainable_alpha_zero),
                               (grads[3], trainable_alpha)])

    return neg_log_prob, trainable_rpc_zero, trainable_rpc, \
      bij_alpha_zero.forward(tf.cast(trainable_alpha_zero,tf.float64)), bij_alpha.forward(tf.cast(trainable_alpha,tf.float64))


  # Create checkpoint for loading of model
  ckpt = tf.train.Checkpoint(trainable_alpha=trainable_alpha, trainable_rpc=trainable_rpc,
                             trainable_alpha_zero=trainable_alpha_zero, trainable_rpc_zero=trainable_rpc_zero,
                             optimizer=optimizer)
  logging.info("Checkpoint path: {}".format(os.path.expandvars(os.path.join(hmm_path, 'tf_ckpts'))))
  manager = tf.train.CheckpointManager(ckpt, os.path.expandvars(os.path.join(hmm_path, 'tf_ckpts')), max_to_keep=1)

  previous = -1. * np.inf
  best_loss = np.inf
  loss_tracker = []

  step1 = 0
  # TRAIN LOOP 1 - all parameters
  for step in range(max_iterations):
    step1 = step
    loss, rpc_zero_est, rpc_est, alpha_zero_est, alpha_est = [t.numpy() for t in train_op()]

    if step == 0:
      save_path = manager.save()
      logging.info("Saved initial checkpoint for step {} with loss {}: {}".format(int(step), loss, save_path))
      initial_loss = loss

    avg_loss = 1.0
    if step % 3 == 0 or verbose:
      logging.info("Round 1|step {} ({}): log prob {} lr {} rpc {:.2e} - {:.2e} alpha {:.2e} - {:.2e}".format(step, (-1.*loss/avg_loss), -loss,
        optimizer_learning_rate, rpc_zero_est, rpc_est, alpha_zero_est, alpha_est))

    loss_tracker.append(-1.*loss)
    if len(loss_tracker) > 20:
      loss_tracker = loss_tracker[1:]
      avg_loss = np.mean(loss_tracker[:-1])

    #  if step >= 9:
      #  logging.info("C1 : {}".format(loss > initial_loss))
      #  print(",".join([str(i) for i in loss_tracker]))
      #  logging.info("C2-content: {}".format((",".join([str(i) for i in loss_tracker]))))
      #  logging.info("C2-avg: {}".format(avg_loss))
      #  logging.info("C2-ratio : {}".format(-1.*loss/avg_loss))
      #  logging.info("C2 : {}".format((-1.*loss/avg_loss) > 0.99999))
      #  logging.info("C3 : {}".format(all(np.diff(loss_tracker) < 0)))

    if step >= 9 and \
            ( ( loss > initial_loss ) or all(np.diff(loss_tracker) < 0) ):
      logging.info("INCREASING LOSS - Loss {} / AVG Loss {}".format(-loss, avg_loss))
      loss = np.inf
      break

    if step >= 50 and \
            ((-1.*loss/avg_loss) > (1- 1e-5)):
      logging.info("Stop early - Loss {} / AVG Loss {} / ratio {}".format(-loss, avg_loss, -1.*loss/avg_loss))
      loss = np.inf
      break

    if np.isnan(loss):
      logging.info("INVALID VALUES - NAN {}".format(alpha_est))
      loss = np.inf
      break

    if loss < best_loss:
      save_path = manager.save()
      logging.info("Saved checkpoint for step {} with loss {}: {}".format(int(step), loss, save_path))
      best_loss = loss

  # restore best hmm so far
  status = ckpt.restore(manager.latest_checkpoint)

  # Restore variable values
  loss = best_loss
  alpha_est = bij_alpha.forward(tf.cast(trainable_alpha, tf.float64)).numpy()
  rpc_est = trainable_rpc.numpy()
  alpha_zero_est = bij_alpha_zero.forward(tf.cast(trainable_alpha_zero, tf.float64)).numpy()
  rpc_zero_est = trainable_rpc_zero.numpy()
  logging.info({"log_prob" : -loss, "learning_rate" : optimizer_learning_rate, "rpc_zero" : rpc_zero_est, "rpc" : rpc_est,
          "alpha_zero" : alpha_zero_est, "alpha" : alpha_est, "epoch1": step1, "num_states" : num_states, "hmm" : hmm})

  return({"log_prob" : -loss, "learning_rate" : optimizer_learning_rate, "rpc_zero" : rpc_zero_est, "rpc" : rpc_est,
          "alpha_zero" : alpha_zero_est, "alpha" : alpha_est, "epoch1": step1, "num_states" : num_states, "hmm" : hmm})


# here we use tensorflow definition of alpha var = mu + mu^2 * alpha
# that is alpha corresponds to reciprocal of size in rnbinom(mu, size)
def interface_hmm(data, rpc, alpha,
                  covariate_info,
                  transition_parameter, name,
                  chromosome_breakpoints=[],
                  max_iterations=201, num_states=10,
                  zero_prior = None, alpha_max_clip = 0.5, alpha_prior_scale = 1.0, alpha_zero_prior_value=None,
                  learning_rates=[0.1, 0.01, 0.001],
                  random_seed=2021, verbose=False,
                  hmm_path=""):

  if zero_prior is None:
    if len(data[data < 0.25 * rpc]) == 0:
      zero_prior = 1e-4
    else:
      zero_prior = np.mean(data[data < 0.25 * rpc]) + 1e-6

  if alpha_zero_prior_value is None:
    alpha_zero_prior_value = alpha

  logger = logging.getLogger()
  if(verbose):
      logger.setLevel(logging.INFO)
  else:
      logger.setLevel(logging.WARNING)

  data = tf.cast(data[:, 0], dtype=tf.float32)
  num_steps = len(data)
  #  covariate_info = tf.cast(covariate_info[0, :], dtype=tf.float32)
  covariate_info = covariate_info[0, :]
  breakpoint_dummy = data[:-1]

  ## debugging
  #  import pickle
  #  debug_dict = {"alpha" : alpha, "rpc" : rpc, "num_states" : num_states, "random_seed" : random_seed,
                #  "learning_rates" : learning_rates, "data" : data, "covariate_info" : covariate_info,
                #  "transition_parameter" : transition_parameter, "transition_parameter_breakpoint" :
                #  transition_parameter_breakpoint}
  #  pickle.dump(debug_dict, open(os.path.expandvars("$HOME/debug_dict.pkl"), "wb"))


  tf.random.set_seed(random_seed)
  np.random.seed(random_seed)

  rpc_prior_value = rpc
  alpha_prior_value = alpha
  #  change_prob_breakpoint = (1. - transition_parameter_breakpoint)

  results = {}
  f = partial(run_hmm, data, rpc, zero_prior, alpha_max_clip, alpha, alpha_zero_prior_value, alpha_prior_scale,
              transition_parameter, chromosome_breakpoints,
              covariate_info, max_iterations, num_states, verbose, hmm_path)

  if(len(learning_rates) == 1):
    results = {learning_rates[0] : f(learning_rates[0])}
  else:
    result_map = []
    for lr in learning_rates:
      result_map.append(f(lr))

    results = {learning_rates[i]: d for i, d in enumerate(result_map)}


  comparison = [(learning, loss["log_prob"], loss["epoch1"]) for learning, loss in results.items()]
  logprobs = list(list(zip(*comparison))[1])
  epochs = list(list(zip(*comparison))[2])

  # backup if all results are NaN
  if all([np.isnan(xx) for xx in logprobs]) or all([epo <= 10 for epo in epochs]):
    results = {1e-3 : f(1e-3)}
    comparison = [(learning, loss["log_prob"], loss["epoch1"]) for learning, loss in results.items()]
    logprobs = list(list(zip(*comparison))[1])
    epochs = list(list(zip(*comparison))[2])

  if all([np.isnan(xx) for xx in logprobs]) or all([epo <= 10 for epo in epochs]):
    results = {1e-4 : f(1e-4)}

  # issues with stability
  #  pool = Pool(n_workers)
  #  result_map = pool.map(lambda lr: run_hmm(data, rpc, zero_prior, alpha_max_clip, alpha, alpha_zero_prior_value, alpha_prior_scale,
            #  transition_parameter, covariate_info, max_iterations, num_states, verbose, lr), learning_rates)

  comparison = [(learning, loss["log_prob"]) for learning, loss in results.items()]
  logprobs = list(list(zip(*comparison))[1])
  logprobs = [-np.inf if np.isnan(xx) else xx for xx in logprobs]

  optimum_key = comparison[np.argmax(logprobs)][0]
  optimum = results[optimum_key]

  most_probable_states = optimum["hmm"].posterior_mode(data).numpy()
  bps = (np.where(most_probable_states[:-1] != most_probable_states[1:])[0])

  posterior_dists = optimum["hmm"].posterior_marginals(data)
  posterior_probs = posterior_dists.logits_parameter().numpy()

  class FinetuneSegmentation(tf.Module):
    def __init__(self, hmm, name):
      self._hmm = hmm
      self._name = tf.constant(name, shape=(), dtype=tf.string)

    @tf.function(input_signature=[])
    def get_name(self):
      return self._name

    @tf.function(input_signature=[tf.TensorSpec.from_tensor(data), tf.TensorSpec.from_tensor(breakpoint_dummy)])
    def update_segmentation(self, x, breakpoint_info):

      transition_probs = breakpoint_info[:, tf.newaxis, tf.newaxis] / (tf.constant(num_states, dtype=tf.float32)-1) * tf.ones(
          [num_states, num_states], dtype=tf.float32)
      diagonal = tf.ones((breakpoint_info.shape[0],num_states), dtype=tf.float32) - breakpoint_info[:, tf.newaxis]
      time_varying_transition_distribution = tf.linalg.set_diag(transition_probs, diagonal)

      hmm_params = self._hmm.parameters
      hmm_post = tfd.HiddenMarkovModel(
        initial_distribution=hmm_params["initial_distribution"],
        transition_distribution=tfd.Categorical(probs=time_varying_transition_distribution),
        observation_distribution=hmm_params["observation_distribution"],
        num_steps=hmm_params["num_steps"],
        time_varying_transition_distribution=True,
        time_varying_observation_distribution=True,
        validate_args=True,
        allow_nan_stats=False)

      most_probable_states_post = hmm_post.posterior_mode(x)
      return(most_probable_states_post)

    @tf.function(input_signature=[tf.TensorSpec.from_tensor(data),\
                                  tf.TensorSpec(shape=(), dtype=tf.int32),\
                                  tf.TensorSpec(shape=(), dtype=tf.int32),\
                                  tf.TensorSpec(shape=(), dtype=tf.int32),\
                                  tf.TensorSpec(shape=(), dtype=tf.int32),\
                                  tf.TensorSpec(shape=(), dtype=tf.int32)])
    def compute_pvalues(self, x, h0, h1, start, end, n_samples):

      hmm_params = self._hmm.parameters
      #  tf.print(hmm_params["observation_distribution"])
      tvod_h0 = hmm_params["observation_distribution"][start:end, h0]
      tvod_h1 = hmm_params["observation_distribution"][start:end, h1]

      samples = tvod_h0.sample(n_samples)

      lr = tf.math.divide(tf.math.exp(tf.reduce_logsumexp(tvod_h0.log_prob(x[start:end]))),
                          tf.math.exp(tf.reduce_logsumexp(tvod_h1.log_prob(x[start:end]))))
      lr_background = tf.math.divide(tf.math.exp(tf.reduce_logsumexp(tvod_h0.log_prob(samples), axis=1)),
                                tf.math.exp(tf.reduce_logsumexp(tvod_h1.log_prob(samples), axis=1)))

      # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC379178/
      mat = tf.where(tf.less_equal(lr_background, lr),
                     tf.constant(1.0, dtype=tf.float32),
                     tf.constant(0.0, dtype=tf.float32))
      extrema = tf.math.reduce_sum(mat)

      pval = tf.math.divide(tf.math.add(tf.cast(extrema, dtype=tf.float32), tf.constant(1.0, dtype=tf.float32)),
             tf.math.add(tf.cast(n_samples, dtype=tf.float32), tf.constant(1.0, dtype=tf.float32)))

      return(pval)

  hmm = optimum["hmm"]

  if not np.isnan(optimum["log_prob"]) and hmm_path != "":
      tf.saved_model.save(FinetuneSegmentation(hmm, name), hmm_path)

  # add state vector info
  observation_distribution = hmm.parameters["observation_distribution"]
  ndim = observation_distribution.batch_shape[1]
  state_vector = np.argmax(observation_distribution.prob(tf.convert_to_tensor(np.transpose(np.tile(data, (ndim, 1))))).numpy(), axis=1)
  oscillations =  np.count_nonzero(np.diff(state_vector))
  optimum["oscillations"] = oscillations

  return({"opt" : optimum, "results": results, "bps" : bps,
          "state_posterior" : most_probable_states,
          "state_marginals" : posterior_probs,
          "state_vector" : state_vector})



def interface_update_hmm(data, breakpoint_info, name, hmm_path, verbose=False):
  #  print(data.shape)
  data = data[0, :]
  breakpoint_info = breakpoint_info[0, :]

  logger = logging.getLogger()
  if(verbose):
      logger.setLevel(logging.INFO)
  else:
      logger.setLevel(logging.ERROR)

  assert len(data) == len(breakpoint_info) + 1

  q = tf.saved_model.load(os.path.expandvars(hmm_path))

  logging.info("Name: {}".format(q.get_name()))
  assert q.get_name() == name
  posterior_mode = q.update_segmentation(data, breakpoint_info)
  bps = (np.where(posterior_mode[:-1] != posterior_mode[1:])[0])

  return({"bps" : bps, "copynumber" : posterior_mode})


def interface_compute_pvalues(data, name, hmm_path, h0, h1, start, end, n_samples=10000):
  pvalues = []
  data = data[0, :]

  if not isinstance(h0, list):
      h0 = [h0]
      h1 = [h1]
      start = [start]
      end = [end]

  n_events = len(h0)
  assert(len(h0) == len(h1))
  assert(len(start) == len(end))
  assert(len(h0) == len(start))

  q = tf.saved_model.load(os.path.expandvars(hmm_path))

  logging.info("Name: {}".format(q.get_name()))
  assert q.get_name() == name

  for index in range(n_events):
    pval = q.compute_pvalues(tf.constant(data, dtype=tf.float32),
                           tf.constant(int(h0[index]), dtype=tf.int32),
                           tf.constant(int(h1[index]), dtype=tf.int32),
                           tf.constant(int(start[index]), dtype=tf.int32),
                           tf.constant(int(end[index]), dtype=tf.int32),
                           tf.constant(n_samples, dtype=tf.int32))
    pvalues.append(pval.numpy())

  return(pvalues)


def test_inference():

  rpc=7.2
  alpha=0.001
  zero=0.01
  n = 400
  max_iterations = 5
  learning_rates = [0.1, 0.01]
  #  learning_rates = [0.1]

  true_rpc=7.4; true_zero=0.3; true_alpha=0.02
  #  covariate_info = np.exp(np.random.normal(loc=0.0, scale=0.2, size=10*n))[tf.newaxis, :]
  covariate_info = np.ones((10*n,))[tf.newaxis,:]#(np.random.normal(loc=0.0, scale=0.2, size=10*n))
  simulation = simulate_data(true_rpc=true_rpc,
                             true_zero=true_zero, true_alpha=true_alpha,
                             cov_info=covariate_info, n_size=1*n)

  transition_parameter = 1e-3
  #  transition_parameter_breakpoint = (1. - 1e-9)

  data = simulation["simulated_data"][0,:]
  true_states = simulation["groundtruth"]
  #  assert(len(data.shape)==1)
  num_steps = data.shape[0]
  num_states = 10

  #  breakpoint_info = np.random.normal(size=num_steps, scale=0.1)
  #  print(transition_probs.shape)
  #  print(diagonal.shape)
  #  diag_view = np.einsum('...ii->...i',transition_probs)
  #  print(transition_probs.shape)
  #  print(transition_probs[0, :, :])
  #  print("---")
  #  print(transition_probs[2, :, :])
  #  print("---")
  #  print(transition_probs[3, :, :])
  #  sys.exit(2)
  #  print(diag_view.shape)
  #  diag_view[:] = diagonal
  #  print(transition_probs.shape)
  #  print(diagonal.shape)
  #  print(breakpoint_info.shape)
  #  np.fill_diagonal(transition_probs,
                   #  diagonal)
  #  print(transition_probs[0, :, :])
  #  print(transition_probs[3, :, :])
  #  sys.exit(2)
  chromosome_breakpoints = [399]
  print(covariate_info.shape)
  print(data.shape)
  print("START")
  import string
  import random
  letters = string.ascii_lowercase
  #  random_name=''.join(random.choice(letters) for i in range(6))
  random_name = "mymodel"
  print(random_name)
  #  print(type(random_name))
  dd = interface_hmm(data[:, tf.newaxis], rpc, alpha,
                      covariate_info,
                      transition_parameter, random_name,
                      chromosome_breakpoints,
                      max_iterations, num_states=8,
                      learning_rates=learning_rates,
                      random_seed=2015,verbose=True,
                      hmm_path="/tmp/hmm.model_"+random_name)

  estimate = dd["opt"]
  hmm = estimate["hmm"]
  posterior_dists = hmm.posterior_marginals(data)
  posterior_probs = posterior_dists.probs_parameter().numpy()

  most_probable_states = hmm.posterior_mode(data)
  bps = (np.where(most_probable_states[:-1] != most_probable_states[1:])[0])

  print(estimate)
  print("- - - - - -")
  print("Start values: rpc {} alpha {} zero {}".format(rpc, alpha, zero))
  print("True values: rpc {} alpha {} zero {}".format(true_rpc, true_alpha, true_zero))
  print("Inferred values: rpc {} alpha {} zero {}".format(estimate["rpc"], estimate["alpha"], estimate["rpc_zero"]))
  print("- - - - - -")
  print("Inferred breakpoints: ", bps)
  print("True breakpoints    : ", np.where(true_states[:-1] != true_states[1:])[0])

  return({"dd" : dd, "data" : data})


def test_inference_realworld():

  rpc = 1.03852
  #  alpha = 0.01500171
  alpha = 0.0001050171
  zero=0.01
  n = 400
  max_iterations = 101
  learning_rates = [0.1, 0.01, 0.05] #, 0.005, 0.001]

  covariate_info = np.loadtxt("/home/schnei01/scAbsolute/debug-segmentation/DEBUG-covariate-UID-10X-BREAST-B_SLX-00000_000011_AAAGATGTCTTGCGGG-1.csv")
  data = np.loadtxt("/home/schnei01/scAbsolute/debug-segmentation/DEBUG-data-UID-10X-BREAST-B_SLX-00000_000011_AAAGATGTCTTGCGGG-1.csv")
  covariate_info = covariate_info[np.newaxis, :]

  transition_parameter = 1e-6

  num_steps = data.shape[0]
  num_states = 8

  chromosome_breakpoints = [7221,14937,21400,27591,33376,38907,43875,48559,52164,56378,60690,
 65001,68165,71053,73603,76026,78512,80967,82785,84746,85844,86927,90664]
  print(covariate_info.shape)
  print(data.shape)
  print("START")
  import string
  import random
  letters = string.ascii_lowercase
  #  random_name=''.join(random.choice(letters) for i in range(6))
  random_name = "mymodel"
  print(random_name)
  #  print(type(random_name))
  dd = interface_hmm(data[:, tf.newaxis], rpc, alpha,
                      covariate_info,
                      transition_parameter, random_name,
                      chromosome_breakpoints,
                      max_iterations, num_states=8,
                      learning_rates=learning_rates,
                      random_seed=2021,verbose=True,
                      hmm_path="/tmp/hmm.model_"+random_name)

  estimate = dd["opt"]
  hmm = estimate["hmm"]
  posterior_dists = hmm.posterior_marginals(tf.cast(data, tf.float32))
  posterior_probs = posterior_dists.probs_parameter().numpy()

  most_probable_states = hmm.posterior_mode(tf.cast(data, tf.float32))
  bps = (np.where(most_probable_states[:-1] != most_probable_states[1:])[0])

  print(estimate)
  print("- - - - - -")
  print("Start values: rpc {} alpha {} zero {}".format(rpc, alpha, zero))
  #  print("True values: rpc {} alpha {} zero {}".format(true_rpc, true_alpha, true_zero))
  print("Inferred values: rpc {} alpha {} zero {}".format(estimate["rpc"], estimate["alpha"], estimate["rpc_zero"]))
  print("- - - - - -")
  print("Inferred breakpoints: ", bps)
  #  print("True breakpoints    : ", np.where(true_states[:-1] != true_states[1:])[0])

  return({"dd" : dd, "data" : data})

#  #  # Code for manual debugging
#  tf.random.set_seed(2020)
#  dd = test_inference()
#  import pandas as pd
#  df = pd.DataFrame.from_dict(dd["dd"]["results"], columns=["log_prob", "learning_rate", "rpc", "alpha", "alpha_zero",
        #  "rpc_zero", "epoch1"], orient="index")
#
#  data = dd["data"]
#  breakpoints = 0.1 * np.ones((data.shape[0]-1))
#  breakpoints[799] = 0.9
#  breakpoints[3199] = 0.9
#  breakpoints[1500] = 0.9
#  breakpoints[1199] = 0.9
#
#  q = tf.saved_model.load("/tmp/hmm.model_mymodel")
#  #  print(q.get_name)
#  name = q.get_name()
#  #  print(type(name))
#  start = [0, 400, 800, 1200, 2000]
#  end = [399, 799, 1199, 1999, 2399]
#  h1 = [1, 3, 4, 2, 4]
#  h0 = [0, 2, 3, 1, 3]
#  print(data.shape)
#  pvals = interface_compute_pvalues(data[tf.newaxis,:], name, "/tmp/hmm.model_mymodel",
        #  h0, h1, start, end, n_samples=10000)
#  print(pvals)
#
#  hmm = dd["dd"]["opt"]["hmm"]
#  posterior_dists = hmm.posterior_marginals(data)
#  posterior_probs = posterior_dists.probs_parameter().numpy()
#
#  posterior_mode = hmm.posterior_mode(data)


# below not working, variables not available
#  print(posterior_mode[800:1205])
#  print(np.where(posterior_mode[:-1] != posterior_mode[1:])[0])
#  n = 400
#  true_rpc=7.4 # see test_inference
#  state_means = [i*true_rpc for i in [1, 2, 3, 2, 4, 3, 4, 5]] + [true_zero]
#  true_durations = [int(i*n_size) for i in [1, 1, 1, 2, 1, 1, 1, 1, 1]]
#  a = tfd.Poisson(log_rate=[1., 2., 3., 4.])
#
#  rpc_zero_est = dd["dd"]["opt"]["rpc_zero"]
#  rpc_est = dd["dd"]["opt"]["rpc"]
#  num_states = dd["dd"]["opt"]["num_states"]
#  levels = np.array([rpc_zero_est] + [rpc_est * i for i in np.arange(1., num_states)])
#  levels = levels[posterior_mode]
#  from matplotlib import pylab as plt
#  fig = plt.figure(figsize=(10, 4))
#  ax = fig.add_subplot(1, 1, 1)
#  ax.plot(levels, c='green', lw=3, label='inferred rate')
#  ax.plot(data, c='black', alpha=0.3, label='observed counts')
#  ax.set_ylabel("latent rate")
#  ax.set_xlabel("time")
#  ax.set_title("Inferred latent rate over time")
#  ax.legend(loc=4)
#  plt.savefig('state.png')

## Test real world example
#  tf.random.set_seed(2020)
#  dd = test_inference_realworld()

## uncomment top one too
##  print(datetime.now() - startTime)

tf.random.set_seed(2020)
dd = test_inference()
