# Copyright 2022, Michael Schneider, All rights reserved.
import numpy as np
import pandas as pd
import tensorflow as tf
import logging
import warnings
import math
import os
import sys

def set_tf_loglevel(level):
    import os
    if level >= logging.FATAL:
        os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
    if level >= logging.ERROR:
        os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
    if level >= logging.WARNING:
        os.environ['TF_CPP_MIN_LOG_LEVEL'] = '1'
    else:
        os.environ['TF_CPP_MIN_LOG_LEVEL'] = '0'
    logging.getLogger('tensorflow').setLevel(level)

# pylint: disable-msg=C0103
# pylint: disable=too-many-instance-attributes
# pylint: disable-msg=too-many-arguments

class dpgmm():
    """Variational Bayesian estimation of a scaling parameter using a Gaussian mixture model.

    This class allows to infer an approximate posterior distribution over the
    parameters of a Gaussian mixture distribution, where the means are fixed
    except for a scaling parameter. The effective number of
    components can be inferred from the data using a truncated stick breaking
    representation. Note that input data is assumed to be 1-dimensional.

    References
    ----------
    1. Blei, David M. and Michael I. Jordan. (2006). "Variational
       inference for Dirichlet process mixtures". Bayesian analysis 1.1
    2. Bishop, C. M. (2006). Pattern recognition and machine learning.
       (New York: Springer).

    Parameters
    ----------
    T : int, defaults to 64.
        An upper bound for the number of mixture components.
    means_prior : array-like, shape (self.p_features,) | None
        The prior on the mean distribution (Gaussian).
        If it is None, it's set to 0 for all components.
    means_precision_prior : array-like, shape (self.p_features,) | None
        The prior on the mean distribution (Gaussian).
        If it is None, it's set to 1 for all components.
    a_prior : array-like, shape (self.p_features,) | None
    b_prior : array-like, shape (self.p_features,) | None
    verbose : int, default to 0.
        Enable verbose output. If 0, only warnings are printed. For 1, also
        output from the information level is given. For values greater than 1,
        debugging output is provided.
    verbose_interval : int, default to 20.
        Interval of steps before printing step counter, if the value is zero, nothing is printed

    """
    def __init__(self, T=64,
                 means_prior=None, means_precision_prior=None,
                 a_prior=None, b_prior=None,
                 verbose=0, verbose_interval=50, dtype=tf.float32):

        # Truncation parameter of infinite sum
        self.T = T
        self.DTYPE = dtype
        self.verbose = verbose
        self.verbose_interval = verbose_interval

        # prior values (just given for a single instance)
        # priors are identical across instances
        # here we define the prior values for a single instance
        self.means_prior = means_prior # m_k
        self.means_precision_prior = means_precision_prior  # \beta_k
        self.a_prior = a_prior
        self.b_prior = b_prior
        self.s_1_prior = 1.
        self.s_2_prior = .1

        # state variables
        # weights
        self.gamma_t1 = None
        self.gamma_t2 = None
        # mu
        self.means = None
        self.means_precision = None
        self.xi = None
        # Lambda
        self.a = None
        self.b = None

        self.covariance = None
        self.weights = None

        self.inf_params = None
        self.vlbs = None
        self.delta_vlbs = None

        self.isFitted = False

        # Setup logging
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s-%(levelname)s-%(message)s',
                            datefmt='%m/%d/%Y_%H:%M:%S')

        #  self.logger= logging.getLogger('tensorflow')
        #  self.logger = logging.getLogger('absVI')
        if self.verbose < 0:
            set_tf_loglevel(logging.DEBUG)
        elif self.verbose == 1:
            set_tf_loglevel(logging.INFO)
        elif self.verbose == 2:
            set_tf_loglevel(logging.WARN)
        else:
            set_tf_loglevel(logging.ERROR)


    def compute_weights(self, gamma_t1, gamma_t2):
        """Compute component weights"""
        weight_dirichlet_sum = (gamma_t1 + gamma_t2)
        tmp = gamma_t2 / weight_dirichlet_sum
        weights = (gamma_t1 / weight_dirichlet_sum *
                   tf.math.cumprod(tmp, exclusive=True, axis=1))

        weights = tf.cast(weights, dtype=tf.float64)
        weights = weights / tf.reduce_sum(input_tensor=weights, axis=1)[:, tf.newaxis]

        return weights

    def compute_covariance(self, a, b):
        covariance = b / a
        covariance = tf.linalg.diag(covariance)

        return covariance

    def check_convergence(self, threshold=1e-3):
        assert(self.isFitted)
        n_no_convergence = np.sum(np.sum(self.delta_vlbs <= threshold, axis=0) > 0)
        n_partial_convergence = np.sum(np.sum(self.delta_vlbs > threshold, axis=0)) - n_no_convergence
        n_all_convergence = self.delta_vlbs.shape[1] - n_no_convergence

        converged_items = np.sum(self.delta_vlbs <= threshold, axis=0) > 0

        return converged_items

    def summary(self):
        n_cells = self.means.shape[0]
        n_samples = self.means.shape[1]
        cell_idx = np.tile(range(1, n_cells+1), self.T)
        T_idx = np.repeat(range(1, self.T+1), n_cells, axis=0)

        d = {'means': np.reshape(self.means, -1, order='F'),
             'weights': np.reshape(self.weights, -1, order='F'),
             'xi': np.reshape(np.repeat(self.xi, self.T, axis=1), -1, order='F'),
             'cell' : cell_idx,
             'T': T_idx,
            }
        result = pd.DataFrame(d)

        return result

    #################### UPDATE EQUATIONS ####################

    def compute_eta_z(self, gamma_t1, gamma_t2):
        """"
        Returns
        -------
        eta_z : float

        see Blei, p.129 and supplementary materials
        E_{q} [\log V_t] + \sum_{i=1}^{t-1} E_{q} [\log (1-V_i)]...
        """
        digamma_sum = tf.math.digamma(gamma_t1 + gamma_t2)
        digamma_a = tf.math.digamma(gamma_t1)
        digamma_b = tf.math.digamma(gamma_t2)

        eta_z = digamma_a - digamma_sum + \
            tf.cumsum((digamma_b - digamma_sum), axis=1, exclusive=True, reverse=False)

        return eta_z

    def compute_eta_x(self, X, means, means_precision, a, b):
        precision = self.compute_covariance(b, a)
        diff = tf.subtract(means[:, tf.newaxis, :, :], X[:, :, tf.newaxis, :])
        cov_tiled = tf.tile(precision[:, tf.newaxis, :, :, :], [1, diff.shape[1], 1, 1, 1])
        s = tf.matmul(diff[: , :, :, tf.newaxis, :], cov_tiled)
        qf = tf.matmul(s, diff[: , :, :, :, tf.newaxis])

        eta_x = tf.constant(-0.5, dtype=self.DTYPE) * (
                  tf.math.log(2*tf.constant(math.pi, dtype=self.DTYPE))
                - tf.math.digamma(a)[:, tf.newaxis, :, :]
                + tf.math.log(tf.norm(tensor=tf.linalg.diag(b), ord="euclidean", axis=[2, 3]))[:, tf.newaxis, :, tf.newaxis]
                + (1./means_precision)[:, tf.newaxis, :, :]
                + tf.squeeze(qf, axis=4) )

        return eta_x

    def e_step(self, X, gamma_t1, gamma_t2, means, means_precision, a, b):
        """ update phi / equivalent to E step."""
        eta_z = self.compute_eta_z(gamma_t1, gamma_t2)
        eta_x = self.compute_eta_x(X, means, means_precision, a, b)

        phi_nk = (eta_z[:, tf.newaxis, :, tf.newaxis] + eta_x - 1)

        #  # solve numerical issues using the log-sum-exp-trick:
        #  # https://stats.stackexchange.com/questions/105602/example-of-how-the-log-sum-exp-trick-works-in-naive-bayes

        #  # normalize phi_nk to probability
        log_sum = tf.reduce_logsumexp(input_tensor=phi_nk, axis=2)
        log_resp = tf.squeeze(phi_nk - log_sum[:, :, tf.newaxis, :], axis=3)

        return log_resp

    def update_V(self, resp, w_1, w_2):
        """ update V.

        Parameters
        ----------
        resp : array-like, shape (n_samples, p_features)
            posterior probabilities (or responsibilities) of
            the point of each sample in X.

        """
        nk = tf.reduce_sum(input_tensor=resp, axis=1)

        gamma_t1 = tf.add(tf.constant(1.0, dtype=self.DTYPE), nk, name="gamma_t1")
        gamma_t2 = tf.add((w_1 / w_2), tf.cumsum(nk, axis=1, exclusive=True, reverse=True), name="gamma_t2")

        # Dirichlet process weight_concentration is a tuple
        # containing the two parameters of the beta distribution
        # see Blei, (Eq. 18, 19) on p. 129
        return gamma_t1, gamma_t2

    def update_alpha(self, gamma_t1, gamma_t2):
        dg_sum = tf.math.digamma(gamma_t1 + gamma_t2)
        dg0 = tf.math.digamma(gamma_t1) - dg_sum
        dg1 = tf.math.digamma(gamma_t2) - dg_sum

        w_1 = tf.add(tf.constant(self.s_1_prior, dtype=self.DTYPE), tf.constant(float(self.T), dtype=self.DTYPE) -
                     tf.constant(1.0, dtype=self.DTYPE), name = "w_1")
        w_2 = tf.subtract(tf.constant(self.s_2_prior, dtype=self.DTYPE), tf.reduce_sum(input_tensor=tf.reverse(dg1 - dg_sum, [1])),
                          name="w_2")

        return w_1, w_2

    def update_mu_dpgmm(self, X, resp, a, b):
        """ update \mu for DPGMM model - not used in the scAbsolute approach"""
        nk = tf.reduce_sum(input_tensor=resp, axis=1)
        means_precision = tf.add(self.means_precision_prior_, nk[:, :, tf.newaxis], name="means_precision")

        # Update mu
        means = tf.add(self.means_prior_ * self.means_precision_prior_, tf.einsum('aji,ajk->aik', resp, X),
                       name="means")
        means /= (self.means_precision_prior_ + nk[:, :, tf.newaxis])

        xi = tf.ones((tf.shape(input=resp)[0], 1), dtype=self.DTYPE, name="xi")

        return means, means_precision, xi

    def update_mu(self, X, resp, a, b):
        nk = tf.reduce_sum(input_tensor=resp, axis=1)
        means_precision = tf.add(self.means_precision_prior_, nk[:, :, tf.newaxis], name="means_precision")

        #Update Xi
        xi = tf.ones((tf.shape(input=resp)[0], 1), dtype=self.DTYPE, name="xi")
        k_vec = tf.tile(tf.range(1, self.T+1, 1.0)[tf.newaxis, :], [tf.shape(input=resp)[0], 1])

        precision = tf.squeeze(self.compute_covariance(b, a), axis=3)
        xi_nom = tf.reduce_sum(input_tensor=resp[:, :, :, tf.newaxis] * k_vec[:, tf.newaxis, :, tf.newaxis] *
                               precision[:, tf.newaxis, :, :] * X[:, :, tf.newaxis, :],
                               axis=[1, 2])
        xi_nom += tf.reduce_sum(input_tensor=k_vec[:, :, tf.newaxis] * self.means_prior_, axis=1)

        xi_denom = tf.reduce_sum(input_tensor=tf.square(k_vec[:, :, tf.newaxis]), axis=1) + \
                   tf.reduce_sum(input_tensor=resp[:, :, :, tf.newaxis] * tf.square(k_vec[:, tf.newaxis, :, tf.newaxis]) *
                                 precision[:, tf.newaxis, :, :], axis = [1, 2])

        xi = xi_nom / xi_denom

        # Means update
        raw = tf.tile(tf.range(1.0, self.T+1, delta=1.0, dtype=self.DTYPE)[tf.newaxis, :, tf.newaxis],
                            [tf.shape(input=xi)[0], 1, 1])
        means = tf.multiply(xi[:, :, tf.newaxis], raw, name="means")

        return means, means_precision, xi

    def update_Lambda(self, X, resp, means):
        """ update Lambda. """
        a = tf.constant(1.0, dtype=self.DTYPE) + (tf.constant(0.5, dtype=self.DTYPE) * tf.reduce_sum(input_tensor=resp, axis=1))
        a = a[:, :, tf.newaxis]

        b = tf.constant(1.0, dtype=self.DTYPE) + \
            (tf.constant(.5, dtype=self.DTYPE) * tf.reduce_sum(input_tensor=resp[:, :, :, tf.newaxis] *
                        tf.square(means[:, tf.newaxis, :, :] - X[:, :, :, tf.newaxis]), axis = 1)) + \
            (tf.constant(.5, dtype=self.DTYPE) * self.means_precision_prior_ * tf.square(means - self.means_prior_))

        return a, b

    def m_step(self, X, log_resp, w_1, w_2, a, b):
        """M step"""
        resp = tf.exp(log_resp)

        dummy = tf.ones_like(resp) * tf.constant(np.finfo(np.float32).eps, dtype=self.DTYPE)
        resp = tf.compat.v1.where(tf.equal(resp, 0.0), dummy, resp)

        gamma_t1, gamma_t2 = self.update_V(resp, w_1, w_2)
        w_1, w_2 = self.update_alpha(gamma_t1, gamma_t2)
        means, means_precision, xi = self.update_mu(X, resp, a, b)
        a, b = self.update_Lambda(X, resp, means)

        return w_1, w_2, gamma_t1, gamma_t2, means, means_precision, xi, a, b
    #################### END UPDATE EQUATIONS ####################

    #################### VARIATIONAL LOWER BOUND ####################
    def compute_variational_bound(self, X, log_resp, hyp_params, inf_params):
        """Estimate the lower bound of the model.

        The lower bound on the likelihood (of the training data with respect to
        the model) is used to detect the convergence and has to decrease at
        each iteration.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
        log_resp : array, shape (n_samples, T)
            Logarithm of the posterior probabilities (or responsibilities) of
            the point of each sample in X.

        Returns
        -------
        lower_bound : float
            value of the variational lower bound
        """
        lower_bound = tf.zeros((tf.shape(input=X)[0]), dtype=self.DTYPE)
        lp_Z = tf.zeros(tf.shape(input=X)[0], dtype=self.DTYPE)

        # unpack parameters
        w_1, w_2 = hyp_params
        gamma_t1, gamma_t2, means, means_precision, xi, a, b = inf_params
        resp = tf.exp(log_resp)
        covariance = self.compute_covariance(a, b)

        # precompute some quantities
        dg_sum = tf.math.digamma(gamma_t1 + gamma_t2)
        dg0 = tf.math.digamma(gamma_t1) - dg_sum
        dg1 = tf.math.digamma(gamma_t2) - dg_sum

        # \alpha ##
        # E_q[log p(alpha | s1, s2)]
        lp_alpha = tf.constant(-1, dtype=self.DTYPE) * (w_1 / w_2)
        lq_alpha = (-w_1 + tf.math.log(w_2) - tf.math.lgamma(w_1) -
                    ((tf.constant(1, dtype=self.DTYPE) - w_1) * tf.math.digamma(w_1)))
        #  tf.print("LB(Z    ): ", lp_Z)
        lower_bound += lp_alpha - lq_alpha

        # V ##
        # E_q[log p(V | 1, alpha)]
        alpha = w_1 / w_2
        lp_V = tf.constant(self.T, dtype=self.DTYPE) * (tf.math.lgamma(1 + alpha) - tf.math.lgamma(alpha)) \
            + (alpha - 1) * tf.reduce_sum(input_tensor=dg1, axis=1)
        # E_q[log q(V | gamma1, gamma2)]
        lq_V = tf.reduce_sum(input_tensor=tf.math.lgamma(gamma_t1 + gamma_t2)
                      - tf.math.lgamma(gamma_t1)
                      - tf.math.lgamma(gamma_t2)
                      + (gamma_t1 - 1) * dg0
                      + (gamma_t2 - 1) * dg1, axis=1)
        #  self.logger.debug("LB(V     ): %.5f (%.5f, %.5f)", lp_V - lq_V, lp_V, lq_V)
        lower_bound += lp_V - lq_V

        # \mu, \Lambda ##
        lpq_ml = tf.constant(-.5, dtype=self.DTYPE) * tf.reduce_sum(input_tensor=tf.square(means - self.means_prior_), axis=1)
        lpq_ml += tf.reduce_sum(input_tensor=a - tf.math.log(b) + tf.math.lgamma(a) + (tf.constant(1., dtype=self.DTYPE) - a) *
                                tf.math.digamma(a) - (a/b), axis = 1)
        #  self.logger.debug("LB(normal): %.5f", lpq_ml)
        lower_bound += tf.squeeze(lpq_ml, axis=1)

        # Z ##
        # Blei, p.129, E_q [\log (z_n | V)]
        # E_q[log p(Z | V)]
        dg_cumsum = tf.cumsum(dg1, exclusive=True, reverse=False)
        lp_Z = tf.reduce_sum(input_tensor=resp * dg_cumsum[:, tf.newaxis, :], axis=[1,2])
        lp_Z += tf.reduce_sum(input_tensor=resp * dg0[:, tf.newaxis, :], axis=[1,2])

        # E_q[log q(Z)]
        # neg. Entropy: \sum_{i} \phi_{n,k} \log(\phi_{n,k}
        lq_Z = tf.reduce_sum(input_tensor=log_resp * resp)

        #  tf.print(lp_Z, [lp_Z], message="LB(Z    )")
        lower_bound += (lp_Z - lq_Z)

        precision = self.compute_covariance(b, a)
        diff = tf.subtract(means[:, tf.newaxis, :, :], X[:, :, tf.newaxis, :])
        cov_tiled = tf.tile(precision[:, tf.newaxis, :, :, :], [1, diff.shape[1], 1, 1, 1])
        s = tf.matmul(diff[: , :, :, tf.newaxis, :], cov_tiled)
        qf = tf.matmul(s, diff[: , :, :, :, tf.newaxis])

        lp_X = tf.reduce_sum(input_tensor=-.5 * (tf.math.log(tf.constant(2*math.pi, dtype=self.DTYPE)) - tf.math.digamma(a)[:, tf.newaxis, :, :]
                    + tf.math.log(tf.norm(tensor=tf.linalg.diag(b), ord="euclidean", axis=[2, 3]))[:, tf.newaxis, :, tf.newaxis]
                    + (1./means_precision)[:, tf.newaxis, :, :]
                    + tf.squeeze(qf, axis=4) ), axis=[1,2,3])

        #  lp_X = tf.print(lp_X, [lp_X], message="LB(X    ): ")
        lower_bound += lp_X

        return lower_bound

    #################### END VARIATIONAL LOWER BOUND ####################

    def fit(self, data, batch_size=16, inner_batch_size=512, subset_ratio=0.8, n_steps=200, n_iterations=3,
            random_state=2018, tf_device="/cpu:0", n_cores=0):
        """Run variational EM algorithm."""
        m_instances, n_samples, p_features = np.shape(data)
        print("Instances %d\tSamples %d\tFeatures %d" % (m_instances, n_samples, p_features))
        # data format: (n_cells, n_bins, n_features)

        with tf.Graph().as_default():
            with tf.device(tf_device):

                # properly initialize random state
                tf.compat.v1.set_random_seed(random_state)

                # prepare input data and subset
                n_subset = tf.cast(math.floor(subset_ratio * n_samples), dtype=tf.int32)
                data_placeholder = tf.compat.v1.placeholder(self.DTYPE, (m_instances, n_samples, p_features))
                dataset = tf.data.Dataset.from_tensor_slices(data_placeholder)
                X = dataset.batch(batch_size)
                iterator = tf.compat.v1.data.make_initializable_iterator(X)
                batch = iterator.get_next()
                bsize = tf.shape(input=batch)[0]

                # optimal parameters -> return values
                opt_gamma_t1 = tf.zeros((m_instances, self.T), dtype=self.DTYPE, name="opt_gamma_t1")
                opt_gamma_t2 = tf.zeros((m_instances, self.T), dtype=self.DTYPE, name="opt_gamma_t2")
                opt_means = tf.zeros((m_instances, self.T, p_features), dtype=self.DTYPE, name="opt_means")
                opt_means_precision = tf.zeros((m_instances, self.T, p_features), dtype=self.DTYPE,
                                               name="opt_means_precision")
                opt_xi = tf.zeros((m_instances, 1), dtype=self.DTYPE, name="opt_xi")
                opt_a = tf.zeros((m_instances, self.T, p_features), dtype=self.DTYPE, name="opt_a")
                opt_b = tf.zeros((m_instances, self.T, p_features), dtype=self.DTYPE, name="opt_b")

                log_resp = tf.zeros((bsize, n_samples, self.T), dtype=self.DTYPE, name="log_resp")

                def body_initialize():
                    with tf.compat.v1.name_scope('initialization'):

                        # subsample for initialization
                        indices = tf.range(n_samples, name="indices")
                        tf.random.shuffle(indices)
                        subset_indices = tf.slice(indices, [0], [n_subset])
                        subset = tf.gather(batch, subset_indices, axis=1)

                        resp = tf.random.uniform([bsize, n_subset, self.T], dtype=self.DTYPE)

                        resp = resp + tf.constant(np.finfo(np.float32).eps, dtype=self.DTYPE) # probably necessary and fine
                        resp = tf.divide(resp, tf.reduce_sum(input_tensor=resp, axis=1, keepdims=True))

                        if self.means_prior is None:
                            #  self.logger.info("Mean prior set automatically")
                            self.means_prior_ = tf.zeros((bsize, self.T, p_features), dtype=self.DTYPE,
                                                         name="means_prior")
                        else:
                            self.means_prior_ = tf.cast(tf.tile(self.means_prior[tf.newaxis, :, :],
                                                                [bsize, 1, 1]),
                                                        self.DTYPE, name="means_prior")

                        if self.means_precision_prior is None:
                            #  self.logger.info("Mean precision prior set automatically")
                            self.means_precision_prior_ = tf.constant(1.0, dtype=self.DTYPE) * tf.ones((bsize, self.T, p_features),
                                                                   dtype=self.DTYPE, name="means_precision_prior")
                        else:
                            self.means_precision_prior_ = tf.cast(tf.tile(self.means_precision_prior[tf.newaxis, :, :],
                                                                          [bsize, 1, 1]),
                                                                  self.DTYPE, name="means_precision_prior")

                        if self.a_prior is None:
                            #  self.logger.info("Covariance gamma shape prior set automatically")
                            self.a_prior_ = tf.ones((bsize, self.T, p_features),
                                                     dtype=self.DTYPE, name="a_prior")
                        else:
                            self.a_prior_ = tf.cast(tf.tile(self.a_prior[tf.newaxis, :, :], [bsize, 1, 1]),
                                                    self.DTYPE, name="a_prior")

                        if self.b_prior is None:
                            #  self.logger.info("Covariance gamma rate prior set automatically")
                            self.b_prior_ = tf.ones((bsize, self.T, p_features), dtype=self.DTYPE,
                                                    name="b_prior")
                        else:
                            self.b_prior_ = tf.cast(tf.tile(self.b_prior[tf.newaxis, :, :], [bsize, 1, 1]),
                                                    self.DTYPE, name="b_prior")

                        # hyp_params
                        w_1 = tf.ones((bsize, 1), dtype=self.DTYPE, name="w_1")
                        w_2 = tf.ones((bsize, 1), dtype=self.DTYPE, name="w_2")
                        a = tf.ones((bsize, self.T, 1), dtype=self.DTYPE, name="a")
                        b = tf.ones((bsize, self.T, 1), dtype=self.DTYPE, name="b")

                        # initialize parameters
                        w_1, w_2, gamma_t1, gamma_t2, \
                            means, means_precision, xi, a, b = \
                            self.m_step(subset, tf.math.log(resp), w_1, w_2, a, b)

                        inf_params = gamma_t1, gamma_t2, means, means_precision, xi, a, b
                        hyp_params = w_1, w_2
                        return hyp_params, inf_params

                def condition(step, batch, hyp_params, inf_params, log_resp):
                    with tf.compat.v1.name_scope('condition'):
                        return step < n_steps

                def body_em(step, batch, hyp_params, inf_params, log_resp):
                    """
                    Update equations
                    """
                    with tf.compat.v1.name_scope('EM'):

                        def debug(step):
                            #  tf.print("Step - ", step)
                            new_step = tf.add(step, tf.constant(1))
                            return(new_step)

                        # update step counter
                        if self.verbose_interval > 0:
                            step = tf.cond(pred=tf.equal(tf.math.mod(step, self.verbose_interval), 0),
                                           true_fn=lambda: debug(step), #tf.add(step, tf.constant(1))
                                                       #tf.print("Step - ",x);x, # tf.add(step, tf.constant(1))),#, [step], message="Step - ",
                                                            #name="stepcounter"),
                                           false_fn=lambda: tf.add(step, tf.constant(1)))
                        else:
                            step = tf.add(step, tf.constant(1))

                        # inner batch loop - SVI
                        # see https://stackoverflow.com/questions/50235287/nested-while-loop-in-tensorflow

                        inner_n_batch = math.floor(int(batch.shape[1])/inner_batch_size)
                        splits = [inner_batch_size for _ in range(inner_n_batch)] + \
                                    [int(batch.shape[1]) - (inner_n_batch * inner_batch_size)]
                        # shuffle this here
                        batch = tf.transpose(a=batch, perm=[1, 0, 2])
                        batch = tf.random.shuffle(batch, name="random_shuffle")
                        batch = tf.transpose(a=batch, perm=[1, 0, 2])

                        items = tf.split(batch, num_or_size_splits=splits, axis=1)
                        items = items[:-1]
                        items = tf.stack(items, axis=1)

                        t_step = tf.constant(0)
                        inner_condition = lambda t_step, a, b, c : tf.less(t_step, tf.cast(inner_n_batch-1,
                                                                                           dtype=tf.int32))
                        self.tau = 1.0
                        self.kappa = 0.6
                        tau = tf.constant(self.tau)
                        kappa = tf.constant(self.kappa)

                        def t_steps(t_step, items, hyp_params, inf_params):

                            # unpack parameters
                            gamma_t1, gamma_t2, means, means_precision, xi, a, b = inf_params
                            w_1, w_2 = hyp_params

                            inner_batch = items[:, t_step, :, :] #tf.gather(items, t_step, axis=0)
                            log_resp = self.e_step(inner_batch, gamma_t1, gamma_t2, means, means_precision, a, b)

                            local_w_1, local_w_2, local_gamma_t1, local_gamma_t2, \
                                local_means, local_means_precision, local_xi, local_a, local_b = \
                                self.m_step(inner_batch, log_resp, w_1, w_2, a, b)

                            # weighted update
                            learning_rate = tf.math.pow(tf.add(tf.cast(t_step, dtype=tf.float32), tau), tf.constant(-1.0) * kappa)

                            w_1 = local_w_1
                            w_2 = local_w_2
                            gamma_t1 = (1.0 - learning_rate) * gamma_t1 + learning_rate * local_gamma_t1
                            gamma_t2 = (1.0 - learning_rate) * gamma_t2 + learning_rate * local_gamma_t2
                            means = (1.0 - learning_rate) * means + learning_rate * local_means
                            means_precision = (1.0 - learning_rate) * means_precision + learning_rate * local_means_precision
                            a = (1.0 - learning_rate) * a + learning_rate * local_a
                            b = (1.0 - learning_rate) * b + learning_rate * local_b

                            # collect results
                            inf_params = gamma_t1, gamma_t2, means, means_precision, xi, a, b
                            hyp_params = w_1, w_2

                            #  with tf.control_dependencies([print_op1, print_op2, print_op3]):
                                #  t_step = tf.add(t_step, tf.constant(1))
                            t_step = tf.add(t_step, tf.constant(1))

                            return (t_step, items, hyp_params, inf_params)

                        # test
                        _, test, hyp_params, inf_params = tf.while_loop(cond=inner_condition, body=t_steps,
                                                                     loop_vars=[t_step,items,hyp_params,inf_params],
                                                                     name="EM-updates-batch", parallel_iterations=1)

                        return (step, batch, hyp_params, inf_params, log_resp)

                def body_vlb(X, hyp_params, inf_params, log_resp):
                    """Compute the weighted log probabilities for each sample.

                    Parameters
                    ----------
                    X : array-like, shape (n_samples, self.p_features)
                        List of self.p_features-dimensional data points. Each row
                        corresponds to a single data point.

                    Returns
                    -------
                    log_prob : array, shape (n_samples,)
                        Log probabilities of each data point in X.
                    """
                    """
                    Variational lower bound
                    """
                    with tf.compat.v1.name_scope("VLB"):
                        # unpack parameters
                        gamma_t1, gamma_t2, means, means_precision, xi, a, b = inf_params
                        w_1, w_2 = hyp_params

                        # compute variational lower bound
                        vlb1 = self.compute_variational_bound(batch, log_resp, hyp_params, inf_params)
                        hyp_params, inf_params, log_resp = (None, None, None)

                        with tf.compat.v1.name_scope("vlb_update"):
                            # apply EM update equations -> single step
                            log_resp = self.e_step(batch, gamma_t1, gamma_t2, means, means_precision, a, b)

                            w_1, w_2, gamma_t1, gamma_t2, \
                                means, means_precision, xi, a, b_ = \
                                self.m_step(batch, log_resp, w_1, w_2, a, b)

                            # collect results
                            inf_params = gamma_t1, gamma_t2, means, means_precision, xi, a, b
                            hyp_params = w_1, w_2

                        vlb2 = self.compute_variational_bound(batch, log_resp, hyp_params, inf_params)
                        delta_vlb = vlb2 - vlb1

                        return vlb1, delta_vlb

                # setup computing graph
                n_batch = int(math.ceil(m_instances/batch_size))
                step = tf.Variable(0, name="step")
                #  iteration = tf.Variable(0, name="iteration")
                init_hyp_params, init_inf_params = body_initialize()

                step, batch, hyp_params, inf_params, \
                    log_resp = tf.while_loop(cond=condition, body=body_em, loop_vars=[step, batch,
                                                                  init_hyp_params,
                                                                  init_inf_params,
                                                                  log_resp], name="EM-updates-batch",
                                                                  parallel_iterations=1)
                vlb, delta_vlb = body_vlb(batch, hyp_params, inf_params, log_resp)

                # variables to store variational lower bound values across iterations
                vlbs = np.zeros((n_iterations, m_instances), dtype=np.float32)
                delta_vlbs = np.zeros((n_iterations, m_instances))
                opt_vlbs = np.ones((m_instances), dtype=np.float32) * np.inf * tf.constant(-1.0, dtype=self.DTYPE)

                metadata = tf.compat.v1.RunMetadata()
                #  run_options = tf.RunOptions(trace_level=tf.RunOptions.FULL_TRACE)
                #  ops = tf.GraphOptions(build_cost_model=5)
                session_conf = tf.compat.v1.ConfigProto(
                                        allow_soft_placement=True,
                                        log_device_placement=False,
                                        #  graph_options=ops,
                                        intra_op_parallelism_threads=n_cores,
                                        inter_op_parallelism_threads=n_cores)

                with tf.compat.v1.Session(config=session_conf) as session:

                    # global initialization
                    # remove graph writer -> creates a lot of files
                    # graph_writer = tf.summary.FileWriter(os.path.expandvars("$HOME/project1/graphs"), session.graph)
                    tf.compat.v1.global_variables_initializer().run()

                    # iterate over all iterations
                    for i_iter in range(n_iterations):
                        print("Iteration %d" % (i_iter+1))

                        # reinitialize dataset iterator for each iteration
                        session.run(iterator.initializer, feed_dict={data_placeholder: data})
                        #  session.run(self.inner_iterator.initializer)
                        results = [[], [], [], [], [], [], []]

                        try:
                            for p in range(n_batch):
                                print("BATCH %d/%d" % (p+1, int(math.ceil(m_instances/batch_size))))

                                #  run VI algorithm
                                s, hyp_params_run, inf_params_run, l_resp, vlb_run, delta_vlb_run = session.run([step,
                                                                                                         hyp_params,
                                                                                                         inf_params,
                                                                                                         log_resp,
                                                                                                         vlb,
                                                                                                         delta_vlb],
                                                                                                        {step: 0}) #,
                                                                                                        #  run_metadata=metadata)
                                                                                                        #  options=run_options
                                gamma_t1, gamma_t2, means, means_precision, xi, a, b = inf_params_run

                                results[0].append(gamma_t1)
                                results[1].append(gamma_t2)
                                results[2].append(means)
                                results[3].append(means_precision)
                                results[4].append(xi)
                                results[5].append(a)
                                results[6].append(b)

                                # compute variational lower bound for inferred parameters
                                vlbs[i_iter, p*batch_size:min((p+1)*batch_size, m_instances)] = vlb_run
                                delta_vlbs[i_iter, p*batch_size:min((p+1)*batch_size, m_instances)] = delta_vlb_run

                                # debugging artefact
                                #  graph_writer.close()


                        except tf.errors.OutOfRangeError:
                            # should be error
                            print("End of dataset")  # ==> "End of dataset"
                            sys.exit(2)

                        gamma_t1, gamma_t2, means, means_precision, xi, a, b  = [np.vstack(results[i]) for i in range(7)]

                        opt_gamma_t1 = tf.compat.v1.where(vlbs[i_iter, :] > opt_vlbs, gamma_t1, opt_gamma_t1)
                        opt_gamma_t2 = tf.compat.v1.where(vlbs[i_iter, :] > opt_vlbs, gamma_t2, opt_gamma_t2)
                        opt_means = tf.compat.v1.where(vlbs[i_iter, :] > opt_vlbs, means, opt_means)

                        opt_means_precision = tf.compat.v1.where(vlbs[i_iter, :] > opt_vlbs, means_precision,
                                                       opt_means_precision)
                        opt_xi = tf.compat.v1.where(vlbs[i_iter, :] > opt_vlbs, xi, opt_xi)
                        opt_a = tf.compat.v1.where(vlbs[i_iter, :] > opt_vlbs, a, opt_a)
                        opt_b = tf.compat.v1.where(vlbs[i_iter, :] > opt_vlbs, b, opt_b)
                        # do not forget to also update the optimal value of vlbs
                        opt_vlbs = tf.reduce_max(input_tensor=vlbs, axis=0)

                    # compute weights and covariance from other variables
                    opt_weights = self.compute_weights(opt_gamma_t1, opt_gamma_t2)
                    opt_covariance = self.compute_covariance(opt_a, opt_b)

                    # store results in object
                    self.gamma_t1 = opt_gamma_t1.eval()
                    self.gamma_t2 = opt_gamma_t2.eval()
                    self.means = opt_means.eval()
                    self.means_precision = opt_means_precision.eval()
                    self.xi = opt_xi.eval()
                    self.a = opt_a.eval()
                    self.b = opt_b.eval()

                    self.weights = opt_weights.eval()
                    self.covariance = opt_covariance.eval()

                    self.inf_params = (self.weights, self.means, self.covariance, self.xi)
                    self.vlbs = opt_vlbs.eval()
                    self.delta_vlbs = delta_vlbs
                    self.isFitted = True

                    session.close()

                return ((self.vlbs, self.delta_vlbs), self.inf_params)

    #  def sample(self, n_samples=1, random_state=None, p_features=1):
        #  """Generate random samples from the fitted Gaussian distribution.
#
        #  Parameters
        #  ----------
        #  n_samples : int, optional
            #  Number of samples to generate. Defaults to 1.
#
        #  Returns
        #  -------
        #  X : array, shape (n_samples, self.p_features)
            #  Randomly generated sample
        #  y : array, shape (nsamples,)
            #  Component labels
        #  """
        #  # everything is implemented for 1D case only
        #  assert(p_features == 1)
        #  assert(n_samples >= 1)
        #  assert(self.isFitted)
#
        #  # DEAL with this numpy issue concerning np.random.multinomial
        #  # https://github.com/numpy/numpy/issues/8317
        #  weights = self.weights
        #  weights = weights / weights.sum(axis=1, keepdims=True).astype("float64")
#
        #  n_samples_comp = np.vstack([np.random.multinomial(n_samples, weights[i, :]) for i in
                                         #  range(weights.shape[0])])
        #  X = np.vstack([mean + np.random.normal(size=(sample, p_features)) * np.sqrt(covariance)
                       #  for i in range(self.means.shape[0])
                       #  for (mean, covariance, sample) in
                       #  zip(self.means[i], self.covariance[i], n_samples_comp[i])])
        #  X = np.reshape(X, (self.means.shape[0], -1))# , order="F")
        #  y = np.concatenate([j * np.ones(sample, dtype=int)
                            #  for i in range(self.means.shape[0])
                            #  for j, sample in enumerate(n_samples_comp[i])])
        #  y = np.reshape(y, (self.means.shape[0], -1)) #, order="F")
#
        #  return (X, y)

def run_scAbsolute(data, truncation=128, batch_size=1,
                   means_precision_prior=None,
                   random_state=2018, verbose=1, verbose_interval=50,
                   n_init=3, n_steps=20):
    """
    Interface to call scAbsolute algorithms

    Parameters
    ----------
    see R code

    """
    if len(data.shape) == 2:
        data = data[:, :, np.newaxis]

    # it is absolutely necessary to rescale data to range of truncation -> locations are fixed to truncation!
    maxima = np.max(data, axis=1)
    prescale = float(truncation-1) / maxima
    data = np.matmul(np.diag(prescale[:, 0]), data[:, :, 0])
    data = data[:, :, np.newaxis]

    tf.compat.v1.set_random_seed(random_state)
    np.random.seed(random_state)

    means_prior = np.arange(1, truncation+1)[:, np.newaxis]

    if means_precision_prior is None:
        means_precision_prior = np.ones((truncation,1)) * 0.1

    model = dpgmm(T=truncation, means_prior=means_prior,
                       means_precision_prior=means_precision_prior,
                       verbose=verbose, verbose_interval=verbose_interval)

    (vlbs, delta_vlbs), inf_params = model.fit(data, batch_size=batch_size, random_state=random_state,
                                               n_iterations=n_init, n_steps=n_steps, tf_device="/cpu:0")
    sys.stdout.flush()

    # debugging
    #  import pickle
    #  params = [data, model.weights, model.means, model.xi, model.covariance]
    #  with open(os.path.expandvars('$HOME/scAbsolute/debug/params.pkl'), 'wb') as output:
      #  pickle.dump(params, output, pickle.HIGHEST_PROTOCOL)

    weights, means, covariance, xi = inf_params

    scale = model.xi[:, 0] / prescale
    solution = pd.DataFrame({'scale': scale[:, 0]})
    result = {'xi': xi[:, 0],
              'prescale': prescale[:, 0],
              'scale': scale[:, 0],
              'weights': weights[:, :],
              'means': means[:, :, 0],
              'covariance': covariance[:, :, 0, 0]}

    return (solution, result)

