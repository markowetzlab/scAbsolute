import logging
import sys

import torch
from torch.distributions import constraints

import pyro
import pyro.distributions as dist
from pyro import poutine
from pyro.infer.autoguide import AutoDelta
from pyro.infer import SVI, JitTraceEnum_ELBO, JitTrace_ELBO, TraceEnum_ELBO, TraceTMC_ELBO, infer_discrete, config_enumerate
from pyro.ops.indexing import Vindex
from pyro.optim import Adam
from pyro.util import ignore_jit_warnings, optional

import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture
from scipy.stats import skew
from scipy.spatial.distance import cityblock

from scdna_replication_tools.compute_consensus_clone_profiles import compute_consensus_clone_profiles, add_cell_ploidies, filter_ploidies
from scdna_replication_tools.normalize_by_cell import compute_cell_corrs

logging.basicConfig(format='%(relativeCreated) 9d %(message)s', level=logging.DEBUG)

# Add another handler for logging debugging events (e.g. for profiling)
# in a separate stream that can be captured.
log = logging.getLogger()
debug_handler = logging.StreamHandler(sys.stdout)
debug_handler.setLevel(logging.DEBUG)
debug_handler.addFilter(filter=lambda record: record.levelno <= logging.DEBUG)
log.addHandler(debug_handler)


class pert_infer_scRT():
    def __init__(self, cn_s, cn_g1, input_col='reads', gc_col='gc', rt_prior_col='mcf7rt',
                 clone_col='clone_id', cell_col='cell_id', library_col='library_id', 
                 chr_col='chr', start_col='start', cn_state_col='state', assign_col='copy',
                 rs_col='rt_state', frac_rt_col='frac_rt', cn_prior_method='g1_composite',
                 learning_rate=0.05, max_iter=2000, min_iter=100, rel_tol=1e-6,
                 max_iter_step1=None, min_iter_step1=None, max_iter_step3=None, min_iter_step3=None,
                 cuda=False, seed=0, P=13, K=4, upsilon=6, run_step3=True):
        '''
        initialise the pert_infer_scRT object
        :param cn_s: long-form dataframe containing copy number and read count information from S-phase cells. (pandas.DataFrame)
        :param cn_g1: long-form dataframe containing copy number and read count information from G1-phase cells. (pandas.DataFrame)
        :param input_col: column containing read count input. (str)
        :param gc_col: column for gc content of each bin. (str)
        :param rt_prior_col: column RepliSeq-determined replication timing values to be used as a prior. (str)
        :param clone_col: column for clone ID of each cell. (str)
        :param cell_col: column for cell ID of each cell. (str)
        :param library_col: column for library ID of each cell. (str)
        :param chr_col: column for chromosome of each bin. (str)
        :param start_col: column for start position of each bin. (str)
        :param cn_state_col: column for the HMMcopy-determined somatic copy number state of each bin; only needs to be present in cn_g1. (str)
        :param assign_col: column used for assigning S-phase cells to their closest matching G1-phase cells when constructing the cn prior concentrations. (str)
        :param rs_col: output column added containing the replication state of each bin for S-phase cells. (str)
        :param frac_rt_col: column added containing the fraction of replciated bins for each S-phase cell. (str)
        :param cn_prior_method: Method for building the cn prior concentrations. Options are 'hmmcopy', 'g1_cells', 'g1_clones', 'g1_composite', 'diploid'. (str)
        :param learning_rate: learning rate of Adam optimizer. (float)
        :param max_iter: max number of iterations of elbo optimization during Step 2 inference. (int)
        :param min_iter: min number of iterations of elbo optimization during Step 2 inference. (int)
        :param max_iter_step1: max number of iterations of elbo optimization during Step 1 inference. (int)
        :param min_iter_step1: min number of iterations of elbo optimization during Step 1 inference. (int)
        :param max_iter_step3: max number of iterations of elbo optimization during Step 3 inference. (int)
        :param min_iter_step3: min number of iterations of elbo optimization during Step 3 inference. (int)
        :param rel_tol: when the relative change in elbo drops to rel_tol, stop inference. (float)
        :param cuda: use cuda tensor type. (bool)
        :param seed: random number generator seed. (int)
        :param P: number of integer copy number states to include in the model, values range from 0 to P-1. (int)
        :param L: number of libraries. (int)
        :param K: maximum polynomial degree allowed for the GC bias curve. (int)
        :param upsilon: value that alpha and beta terms should sum to when creating a beta distribution prior for tau. (int)
        :param run_step3: Should the pretrained S-phase model be run on low variance cells to look for misclassifications. (bool)
        '''
        self.cn_s = cn_s
        self.cn_g1 = cn_g1

        self.input_col = input_col
        self.gc_col = gc_col
        self.rt_prior_col = rt_prior_col
        self.clone_col = clone_col
        self.cell_col = cell_col
        self.library_col = library_col
        self.chr_col = chr_col
        self.start_col = start_col
        self.cn_state_col = cn_state_col
        self.assign_col = assign_col

        self.rs_col = rs_col
        self.frac_rt_col = frac_rt_col

        self.learning_rate = learning_rate
        self.max_iter = max_iter
        self.min_iter = min_iter
        self.rel_tol = rel_tol
        self.cuda = cuda
        self.seed = seed

        # set max/min iter for step 1 and 3 to be half as those for step 2 if None
        if max_iter_step1 is None:
            self.max_iter_step1 = int(self.max_iter/2)
        else:
            self.max_iter_step1 = max_iter_step1
        if min_iter_step1 is None:
            self.min_iter_step1 = int(self.min_iter/2)
        else:
            self.min_iter_step1 = min_iter_step1
        if max_iter_step3 is None:
            self.max_iter_step3 = int(self.max_iter/2)
        else:
            self.max_iter_step3 = max_iter_step3
        if min_iter_step3 is None:
            self.min_iter_step3 = int(self.min_iter/2)
        else:
            self.min_iter_step3 = min_iter_step3
        
        self.cn_prior_method = cn_prior_method

        self.P = P
        self.L = None
        self.K = K

        self.upsilon = upsilon
        self.run_step3 = run_step3


    def process_input_data(self):
        # sort rows by correct genomic ordering
        self.cn_g1 = self.sort_by_cell_and_loci(self.cn_g1)
        self.cn_s = self.sort_by_cell_and_loci(self.cn_s)

        # drop any row where read count input is NaN
        self.cn_g1 = self.cn_g1[self.cn_g1[self.input_col].notna()]
        self.cn_s = self.cn_s[self.cn_s[self.input_col].notna()]

        # pivot to 2D matrix where each row is a unique cell, columns are loci
        cn_g1_reads_df = self.cn_g1.pivot_table(index=self.cell_col, columns=[self.chr_col, self.start_col], values=self.input_col)
        cn_g1_states_df = self.cn_g1.pivot_table(index=self.cell_col, columns=[self.chr_col, self.start_col], values=self.cn_state_col)
        cn_s_reads_df = self.cn_s.pivot_table(index=self.cell_col, columns=[self.chr_col, self.start_col], values=self.input_col)
        cn_s_states_df = self.cn_s.pivot_table(index=self.cell_col, columns=[self.chr_col, self.start_col], values=self.cn_state_col)

        cn_g1_reads_df = cn_g1_reads_df.dropna(axis=1)
        cn_g1_states_df = cn_g1_states_df.dropna(axis=1)
        cn_s_reads_df = cn_s_reads_df.dropna(axis=1)
        cn_s_states_df = cn_s_states_df.dropna(axis=1)

        assert cn_g1_states_df.shape == cn_g1_reads_df.shape
        print("DEBUG")
        print(cn_s_reads_df.columns)
        print("---")
        print(cn_g1_reads_df.columns)
        ## FIX assertion bug
        # Get the common index values using intersection
        common_index = cn_s_reads_df.columns.intersection(cn_g1_reads_df.columns)
        print('- - - - - ')
        print(common_index)
        cn_s_reads_df = cn_s_reads_df.iloc[:, cn_s_reads_df.columns.isin(common_index)]
        cn_g1_reads_df = cn_g1_reads_df.iloc[:, cn_g1_reads_df.columns.isin(common_index)]
        cn_s_states_df = cn_s_states_df.iloc[:, cn_s_states_df.columns.isin(common_index)]
        cn_g1_states_df = cn_g1_states_df.iloc[:, cn_g1_states_df.columns.isin(common_index)]
        print("Post fix")
        print(cn_s_reads_df.columns)
        print("---")
        print(cn_g1_reads_df.columns)
        assert cn_s_reads_df.shape[1] == cn_g1_reads_df.shape[1]
        assert cn_s_states_df.shape == cn_s_reads_df.shape

        # clean cn_s and cn_g1
        #  import pickle
        #  with open("/home/schnei01/scdna_replication_tools/debug1.pkl", "wb") as output_file:
          #  pickle.dump(cn_g1_reads_df, output_file)
        #  with open("/home/schnei01/scdna_replication_tools/debug2.pkl", "wb") as output_file:
          #  pickle.dump(cn_s_reads_df, output_file)
        #  #  with open("/home/schnei01/scdna_replication_tools/debug3.pkl", "wb") as output_file:
          #  #  pickle.dump(clone_cn_profiles, output_file)
        #  with open("/home/schnei01/scdna_replication_tools/debug0.pkl", "wb") as output_file:
            #  pickle.dump(self.cn_s, output_file)
        #print(self.cn_s)

        chr_values = cn_s_reads_df.columns.get_level_values("chr")
        start_values = cn_s_reads_df.columns.get_level_values("start")
        check_list = [(chr_values[i], start_values[i]) for i,_ in enumerate(list(chr_values))]
        #target_columns = ("chr", "start")
        #  print("ISSUE")
        #  print(check_list)
        #  print("---- ----- -----")
        print(self.cn_s.shape)
        print(self.cn_g1.shape)
        #  print(type(self.cn_s))
        self.cn_s = self.cn_s[self.cn_s[['chr', 'start']].apply(tuple, axis=1).isin(check_list)]
        self.cn_g1 = self.cn_g1[self.cn_g1[['chr', 'start']].apply(tuple, axis=1).isin(check_list)]
        #print(' - - - - - ')
        #self.cn_s = matches
        print(self.cn_s.shape)
        print(self.cn_g1.shape)
        #print(type(self.cn_s))

        #sys.exit(3)
        #print("PASS")

        cn_g1_reads_df = cn_g1_reads_df.T
        cn_g1_states_df = cn_g1_states_df.T
        cn_s_reads_df = cn_s_reads_df.T
        cn_s_states_df = cn_s_states_df.T

        # convert to tensor and unsqueeze the data dimension
        # convert to int64 before float32 to ensure that all values are rounded to the nearest int
        cn_g1_reads = torch.tensor(cn_g1_reads_df.values).to(torch.int64).to(torch.float32)
        cn_g1_states = torch.tensor(cn_g1_states_df.values).to(torch.int64).to(torch.float32)
        cn_s_reads = torch.tensor(cn_s_reads_df.values).to(torch.int64).to(torch.float32)
        cn_s_states = torch.tensor(cn_s_states_df.values).to(torch.int64).to(torch.float32)

        # get tensor of library_id index
        # need this because each library will have unique gc params
        libs_s, libs_g1 = self.get_libraries_tensor(self.cn_s, self.cn_g1)

        # make sure there's one library index per cell
        assert libs_s.shape[0] == cn_s_reads.shape[1]
        assert libs_g1.shape[0] == cn_g1_reads.shape[1]

        # get tensor for GC profile
        gammas = self.cn_s[[self.chr_col, self.start_col, self.gc_col]].drop_duplicates()
        gammas = gammas.dropna()
        gammas = gammas.loc[list(map(lambda x: x in check_list, zip(gammas['chr'], gammas['start']))),:]
        gammas = torch.tensor(gammas[self.gc_col].values).to(torch.float32)

        # get tensor for rt prior if provided
        if (self.rt_prior_col is not None) and (self.rt_prior_col in self.cn_s.columns):
            rt_prior_profile = self.cn_s[[self.chr_col, self.start_col, self.rt_prior_col]].drop_duplicates()
            rt_prior_profile = rt_prior_profile.dropna()
            rt_prior_profile = torch.tensor(rt_prior_profile[self.rt_prior_col].values).unsqueeze(-1).to(torch.float32)
            rt_prior_profile = self.convert_rt_prior_units(rt_prior_profile)
            assert cn_s_reads.shape[0] == gammas.shape[0] == rt_prior_profile.shape[0]
        else:
            rt_prior_profile = None

        return cn_g1_reads_df, cn_g1_states_df, cn_s_reads_df, cn_s_states_df, cn_g1_reads, cn_g1_states, cn_s_reads, cn_s_states, gammas, rt_prior_profile, libs_g1, libs_s


    def sort_by_cell_and_loci(self, cn):
        """ Sort long-form dataframe so each cell follows correct genomic ordering """
        cn[self.chr_col] = cn[self.chr_col].astype(str)
        cn[self.chr_col] = cn[self.chr_col].astype('category')
        chr_order = [str(i+1) for i in range(22)]
        chr_order.append('X')
        chr_order.append('Y')
        cn[self.chr_col] = cn[self.chr_col].cat.set_categories(chr_order)
        cn = cn.sort_values(by=[self.cell_col, self.chr_col, self.start_col])
        return cn


    def get_libraries_tensor(self, cn_s, cn_g1):
        """ Create a tensor of integers representing the unique library_id of each cell. """
        libs_s = cn_s[[self.cell_col, self.library_col]].drop_duplicates()
        libs_g1 = cn_g1[[self.cell_col, self.library_col]].drop_duplicates()

        # get all unique library ids found across cells of both cell cycle phases
        all_library_ids = pd.concat([libs_s, libs_g1])[self.library_col].unique()

        self.L = int(len(all_library_ids))
        
        # replace library_id strings with integer values
        for i, library_id in enumerate(all_library_ids):
            libs_s[self.library_col].replace(library_id, i, inplace=True)
            libs_g1[self.library_col].replace(library_id, i, inplace=True)
        
        # convert to tensors of type int (ints needed to index other tensors)
        libs_s = torch.tensor(libs_s[self.library_col].values).to(torch.int64)
        libs_g1 = torch.tensor(libs_g1[self.library_col].values).to(torch.int64)

        return libs_s, libs_g1


    def convert_rt_prior_units(self, rt_prior_profile):
        """ Make sure that early RT regions are close to 1, late RT regions are close to 0 """
        rt_prior_profile = rt_prior_profile / max(rt_prior_profile)
        return rt_prior_profile


    def build_trans_mat(self, cn):
        """ Use the frequency of state transitions in cn to build a new transition matrix. """
        trans_mat = torch.eye(self.P, self.P) + 1
        num_loci, num_cells = cn.shape
        for i in range(num_cells):
            for j in range(1, num_loci):
                cur_state = int(cn[j, i])
                prev_state = int(cn[j-1, i])
                trans_mat[prev_state, cur_state] += 1
        return trans_mat


    def build_cn_prior(self, cn, weight=1e6):
        """ Build a matrix with the cn prior concentration (eta) for each bin's cn state based on its value in cn. """
        num_loci, num_cells = cn.shape
        etas = torch.ones(num_loci, num_cells, self.P)
        for i in range(num_loci):
            for n in range(num_cells):
                state = int(cn[i, n].numpy())
                etas[i, n, state] = weight
        return etas


    def build_clone_cn_prior(self, cn, cn_df, cn_tensor, clone_cn_profiles, weight=1e6):
        """ Use the pseudobulk clone profiles as the prior concentration for copy number """
        cn_prior_input = torch.zeros(cn_tensor.shape)
        for i, cell_id in enumerate(cn_df.columns):
            cell_cn = cn.loc[cn[self.cell_col]==cell_id]  # get full cn data for this cell
            cell_clone = cell_cn[self.clone_col].values[0]  # get clone id
            cn_prior_input[:, i] = torch.tensor(clone_cn_profiles[cell_clone].values).to(torch.int64).to(torch.float32)  # assign consensus clone cn profile for this cell
        
        # build a proper prior over num_states using the consensus clone cn calls for each cell
        etas = self.build_cn_prior(cn_prior_input, weight=weight)

        return etas


    def build_composite_cn_prior(self, cn, clone_cn_profiles, weight=1e5, J=5):
        """ 
        Build a cn prior concentration matrix that uses both the consensus G1 clone profile and the closest G1-phase cell profiles.
        J represents the number of matching G1-phase cells to use for each S-phase cell.
        """

        temp_cn_g1 = self.cn_g1.copy()

        # remove G1 cells from certain clones that don't belong to the majority ploidy
        # i.e. remove tetraploid cells if clone is 90% diploid
        ploidy_col = 'ploidy'
        if ploidy_col not in temp_cn_g1.columns:
            temp_cn_g1 = add_cell_ploidies(temp_cn_g1, cell_col=self.cell_col, cn_state_col=self.cn_state_col, ploidy_col=ploidy_col)
        temp_cn_g1 = filter_ploidies(temp_cn_g1, clone_col=self.clone_col, ploidy_col=ploidy_col)

        num_loci, num_cells = cn.shape
        etas = torch.ones(num_loci, num_cells, self.P)

        for n, cell_id in enumerate(cn.columns):
            cell_cn = self.cn_s.loc[self.cn_s[self.cell_col]==cell_id]  # get full cn data for this cell
            cell_clone = cell_cn[self.clone_col].values[0]  # get clone id
            clone_cn_profile = torch.tensor(clone_cn_profiles[cell_clone].values).to(torch.int64).to(torch.float32)  # assign consensus clone cn profile for this cell
        
            # shrink set of G1 cells to those in the matching clone if clone_id has already been assigned
            # for this S-phase cell
            if self.clone_col is not None:
                clone_id = cell_cn[self.clone_col].values[0]
                clone_cn_g1 = temp_cn_g1.loc[temp_cn_g1[self.clone_col]==clone_id]
            else:
                clone_cn_g1 = temp_cn_g1
            
            # compute pearson correlations between this S-phase cell and all G1-phase cells in the same clone
            psi_mat = compute_cell_corrs(cell_cn, clone_cn_g1, cell_id, col=self.input_col,
                                         cell_col=self.cell_col, chr_col=self.chr_col, start_col=self.start_col)

            # get data from the top 5 G1 cells that best match the S-phase cell
            g1_cell_cns = np.zeros((num_loci, J))
            for j in range(J):
                g1_cell_id = psi_mat.iloc[j].g1_cell_id  # find the cell_id corresponding to this ranked match
                g1_cell_cn = clone_cn_g1.loc[clone_cn_g1[self.cell_col]==g1_cell_id]  # save the cn profile of this G1-phase cell
                temp_cn_profile = g1_cell_cn[self.cn_state_col].values
                g1_cell_cns[:, j] = temp_cn_profile

            # loop through all the loci for this cell and add the appropriate concentrations
            # based on the clone and cell cn profiles
            for i in range(num_loci):
                # add weight for consensus clone profile at this position
                temp_clone_state = int(clone_cn_profile[i].numpy())
                temp_weight = weight  * J * 2
                etas[i, n, temp_clone_state] += temp_weight

                # loop through top 5 matching G1-phase cells and add weight for those states too
                for j in range(J):
                    temp_cell_state = int(g1_cell_cns[i, j])
                    temp_weight = weight * (J - j)
                    etas[i, n, temp_cell_state] += temp_weight

        return etas


    def manhattan_binarization(self, X, MEAN_GAP_THRESH=0.7, EARLY_S_SKEW_THRESH=0.2, LATE_S_SKEW_THRESH=-0.2):
        """ Binarize X into binary replicated vs unreplicated states by drawing an optimal threshold through X that minimizes the manhattan distance of all points. """
        # center and scale the data
        X = (X - np.mean(X)) / np.std(X)
        
        # fit a 2-state GMM to the data
        gm = GaussianMixture(n_components=2, random_state=0)
        states = gm.fit_predict(X)
        
        # use GMM means to assign binary values for thresholding
        mean_0 = gm.means_[0][0]
        mean_1 = gm.means_[1][0]

        # find the distance between the two means for each state
        mean_gap = abs(mean_0 - mean_1)

        # assume means denote binary values
        binary_0 = min(mean_0, mean_1)
        binary_1 = max(mean_0, mean_1)
        
        X = X.flatten()
        
        # use skew to define the binary values if means are close together
        if mean_gap < MEAN_GAP_THRESH:
            cell_skew = skew(X)
            # positive skew indicates early S-phase
            if cell_skew > EARLY_S_SKEW_THRESH:
                binary_0 = np.percentile(X, 50)
                binary_1 = np.percentile(X, 95)
            # negative skew indicates late S-phase
            elif cell_skew < LATE_S_SKEW_THRESH:
                binary_0 = np.percentile(X, 5)
                binary_1 = np.percentile(X, 50)
            # assume mid-S when skew is neutral
            else:
                binary_0 = np.percentile(X, 25)
                binary_1 = np.percentile(X, 75)

        # now that binary values are selected, I must compute the Manhattan distance
        # between binarized data and X for 100 different thresholds
        threshs = np.linspace(binary_0, binary_1, 100)
        lowest_dist = np.inf
        best_t = None
        manhattan_dists = []
        for t in threshs:
            # set values to binary_1 when above t, to binary_0 when below t
            B = np.where(X>t, binary_1, binary_0)
            # compute Manhattan distance between two vectors
            dist = cityblock(X, B)
            manhattan_dists.append(dist)
            if dist < lowest_dist:
                lowest_dist = dist
                best_t = t

        # binarize X based on the best threshold
        cell_rt = np.where(X>best_t, 1, 0)
        # compute fraction of replicated bins (cell's time within s-phase)
        frac_rt = sum(cell_rt) / len(cell_rt)
        
        return cell_rt, frac_rt


    def guess_times(self, cn_s_reads, etas):
        """ 
        Come up with a guess for what each cell's time in S-phase should be by
        binarizing the cn-normalized read count.
        """
        num_loci, num_cells = cn_s_reads.shape
        t_init = torch.zeros(num_cells)
        t_alpha_prior = torch.zeros(num_cells)
        t_beta_prior = torch.zeros(num_cells)

        # normalize raw read count by whatever state has the highest probability in the cn prior
        # assume cn=0.5 in regions where there is a homozygous deletion
        ones = (torch.ones(cn_s_reads.shape) * 0.5).type(torch.float32)
        cn_states = torch.argmax(etas, dim=2).type(torch.float32)
        reads_norm_by_cn = cn_s_reads / torch.where(cn_states > 0.0, cn_states, ones)

        for i in range(num_cells):
            cell_profile = reads_norm_by_cn[:, i]
            
            X = cell_profile.numpy().reshape(-1, 1)
            y_pred2, t_guess = self.manhattan_binarization(X)
            
            t_init[i] = t_guess
            
            # use t_guess as the mean of a beta distribution parameterized by alpha and beta
            # where alpha and beta must sum to upsilon
            alpha = t_guess * self.upsilon
            beta = self.upsilon - alpha
            t_alpha_prior[i] = alpha
            t_beta_prior[i] = beta

        return t_init, t_alpha_prior, t_beta_prior


    def make_gc_features(self, x):
        """Builds features i.e. a matrix with columns [x, x^2, x^3, x^4]."""
        x = x.unsqueeze(1)
        return torch.cat([x ** i for i in reversed(range(0, self.K+1))], 1)


    def package_s_output(self, cn_s, trace_s, cn_s_reads_df, lambda_fit, losses_g, losses_s):
        """ Use model trace to extract all the fitted parameters and store them in DataFrames to be saved. """

        cn_s = cn_s.copy()

        # extract fitted parameters
        u_fit_s = trace_s.nodes['expose_u']['value']
        rho_fit_s = trace_s.nodes['expose_rho']['value']
        a_fit_s = trace_s.nodes['expose_a']['value']
        tau_fit_s = trace_s.nodes['expose_tau']['value']
        model_rep = trace_s.nodes['rep']['value']
        model_cn = trace_s.nodes['cn']['value']

        # add inferred CN and Rep states to the S-phase output df
        model_cn_df = pd.DataFrame(model_cn.detach().numpy(), index=cn_s_reads_df.index, columns=cn_s_reads_df.columns)
        model_rep_df = pd.DataFrame(model_rep.detach().numpy(), index=cn_s_reads_df.index, columns=cn_s_reads_df.columns)
        model_cn_df = model_cn_df.melt(ignore_index=False, value_name='model_cn_state').reset_index()
        model_rep_df = model_rep_df.melt(ignore_index=False, value_name='model_rep_state').reset_index()
        cn_s_out = pd.merge(cn_s, model_cn_df)
        cn_s_out = pd.merge(cn_s_out, model_rep_df)

        # add other inferred parameters to cn_s_out
        taus_out = pd.DataFrame(
            tau_fit_s.detach().numpy(),
            index=cn_s_reads_df.columns, 
            columns=['model_tau']).reset_index()
        Us_out = pd.DataFrame(
            u_fit_s.detach().numpy(),
            index=cn_s_reads_df.columns, 
            columns=['model_u']).reset_index()
        rho_out = pd.DataFrame(
            rho_fit_s.detach().numpy(),
            index=cn_s_reads_df.index,
            columns=['model_rho']).reset_index()
        cn_s_out = pd.merge(cn_s_out, taus_out)
        cn_s_out = pd.merge(cn_s_out, Us_out)
        cn_s_out = pd.merge(cn_s_out, rho_out)

        # create a separate df containing library- and sample-level params
        supp_out_df = []

        # add global lambda parameter
        supp_out_df.append(pd.DataFrame({
            'param': ['model_lambda'],
            'level': ['all'],
            'value': [lambda_fit.detach().numpy()[0]]
        }))

        # add global alpha parameter
        supp_out_df.append(pd.DataFrame({
            'param': ['model_a'],
            'level': ['all'],
            'value': [a_fit_s.detach().numpy()[0]]
        }))

        # add loss for each step in the G1/2-phase model
        supp_out_df.append(pd.DataFrame({
            'param': ['loss_g']*len(losses_g),
            'level': np.arange(len(losses_g)),
            'value': losses_g
        }))

        # add loss for each step in the S-phase model
        supp_out_df.append(pd.DataFrame({
            'param': ['loss_s']*len(losses_s),
            'level': np.arange(len(losses_s)),
            'value': losses_s
        }))

        # concatentate into one dataframe
        supp_out_df = pd.concat(supp_out_df, ignore_index=True)

        return cn_s_out, supp_out_df


    @config_enumerate
    def model_s(self, gammas, libs, cn0=None, rho0=None, num_cells=None, num_loci=None, data=None, etas=None, lamb=None, lambda_init=1e-1, t_alpha_prior=None, t_beta_prior=None, t_init=None):
        with ignore_jit_warnings():
            if data is not None:
                num_loci, num_cells = data.shape
            elif cn0 is not None:
                num_loci, num_cells = cn0.shape
            assert num_cells is not None
            assert num_loci is not None
            assert data is not None

        # controls the consistency of replicating on time
        a = pyro.sample('expose_a', dist.Gamma(torch.tensor([2.]), torch.tensor([0.2])))
        
        # variance of negative binomial distribution is governed by the success probability of each trial
        if lamb is None:
            lamb = pyro.param('expose_lambda', torch.tensor([lambda_init]), constraint=constraints.interval(0.001, 0.999))

        # gc bias params
        beta_means = pyro.sample('expose_beta_means', dist.Normal(0., 1.).expand([self.L, self.K+1]).to_event(2))
        beta_stds = pyro.param('expose_beta_stds', torch.logspace(start=0, end=-self.K, steps=(self.K+1)).reshape(1, -1).expand([self.L, self.K+1]),
                               constraint=constraints.positive)
        
        # define cell and loci plates
        loci_plate = pyro.plate('num_loci', num_loci, dim=-2)
        cell_plate = pyro.plate('num_cells', num_cells, dim=-1)

        if rho0 is not None:
            # fix replication timing as constant when input into model
            rho = rho0
        else:
            with loci_plate:
                # bulk replication timing profile
                rho = pyro.sample('expose_rho', dist.Beta(torch.tensor([1.]), torch.tensor([1.])))

        with cell_plate:

            # per cell time in S-phase (tau)
            # draw from prior if provided
            if (t_alpha_prior is not None) and (t_beta_prior is not None):
                tau = pyro.sample('expose_tau', dist.Beta(t_alpha_prior, t_beta_prior))
            elif t_init is not None:
                tau = pyro.param('expose_tau', t_init, constraint=constraints.unit_interval)
            else:
                tau = pyro.sample('expose_tau', dist.Beta(torch.tensor([1.5]), torch.tensor([1.5])))
            
            # per cell reads per copy per bin
            # u should be inversely related to tau and ploidy, positively related to reads
            if cn0 is not None:
                cell_ploidies = torch.mean(cn0.type(torch.float32), dim=0)
            elif etas is not None:
                temp_cn0 = torch.argmax(etas, dim=2).type(torch.float32)
                cell_ploidies = torch.mean(temp_cn0, dim=0)
            else:
                cell_ploidies = torch.ones(num_cells) * 2.

            u_guess = torch.mean(data.type(torch.float32), dim=0) / ((1 + tau) * cell_ploidies)
            u_stdev = u_guess / 10.
        
            u = pyro.sample('expose_u', dist.Normal(u_guess, u_stdev))

            # sample beta params for each cell based on which library the cell belongs to
            betas = pyro.sample('expose_betas', dist.Normal(beta_means[libs], beta_stds[libs]).to_event(1))
            
            with loci_plate:

                if cn0 is None:
                    if etas is None:
                        etas = torch.ones(num_loci, num_cells, self.P)
                    # sample cn probabilities of each bin from Dirichlet
                    pi = pyro.sample('expose_pi', dist.Dirichlet(etas))
                    # sample cn state from categorical based on cn_prob
                    cn = pyro.sample('cn', dist.Categorical(pi), infer={"enumerate": "parallel"})

                # per cell per bin late or early 
                t_diff = tau.reshape(-1, num_cells) - rho.reshape(num_loci, -1)

                # probability of having been replicated
                phi = 1 / (1 + torch.exp(-a * t_diff))

                # force phi to remain on the domain of 0.001 to 0.999
                phi[phi<0.001] = 0.001
                phi[phi>0.999] = 0.999

                # binary replicated indicator
                print("DEBUG AGAIN")
                print(phi)
                print(phi.shape)
                import pickle
                with open("/home/schnei01/scdna_replication_tools/debugphi.pkl", "wb") as output_file:
                    pickle.dump(phi, output_file)

                rep = pyro.sample('rep', dist.Bernoulli(phi), infer={"enumerate": "parallel"})
                #sys.exit(2)

                # total copy number accounting for replication
                chi = cn * (1. + rep)

                # find the gc bias rate of each bin using betas
                gc_features = self.make_gc_features(gammas).reshape(num_loci, 1, self.K+1)
                omega = torch.exp(torch.sum(torch.mul(betas, gc_features), 2))

                # expected reads per bin per cell
                theta = u * chi * omega

                # use lambda and the expected read count to define the number of trials (delta)
                # that should be drawn for each bin
                delta = theta * (1 - lamb) / lamb

                # replace all delta<1 values with 1 since delta should be >0
                # this avoids NaN errors when theta=0 at a given bin
                delta[delta<1] = 1
                
                reads = pyro.sample('reads', dist.NegativeBinomial(delta, probs=lamb), obs=data)


    def model_g1(self, gammas, libs, cn=None, num_cells=None, num_loci=None, lambda_init=1e-1, data=None):
        with ignore_jit_warnings():
            if data is not None:
                num_loci, num_cells = data.shape
            elif cn is not None:
                num_loci, num_cells = cn.shape
            assert num_cells is not None
            assert num_loci is not None
            assert (data is not None) and (cn is not None)
        
        # variance of negative binomial distribution is governed by the success probability of each trial
        lamb = pyro.param('expose_lambda', torch.tensor([lambda_init]), constraint=constraints.interval(0.001, 0.999))

        # gc bias params
        beta_means = pyro.sample('expose_beta_means', dist.Normal(0., 1.).expand([self.L, self.K+1]).to_event(2))
        beta_stds = pyro.param('expose_beta_stds', torch.logspace(start=0, end=-self.K, steps=(self.K+1)).reshape(1, -1).expand([self.L, self.K+1]),
                               constraint=constraints.positive)

        with pyro.plate('num_cells', num_cells):

            # per cell reads per copy per bin
            # u should be solved for when cn and read count are both observed
            cell_ploidies = torch.mean(cn.type(torch.float32), dim=0)
            u = torch.mean(data.type(torch.float32), dim=0) / cell_ploidies

            # sample beta params for each cell based on which library the cell belongs to
            betas = pyro.sample('expose_betas', dist.Normal(beta_means[libs], beta_stds[libs]).to_event(1))

            with pyro.plate('num_loci', num_loci):

                # copy number accounting for gc bias
                gc_features = self.make_gc_features(gammas).reshape(num_loci, 1, self.K+1)
                omega = torch.exp(torch.sum(torch.mul(betas, gc_features), 2))  # compute the gc 'rate' of each bin

                # expected reads per bin per cell
                theta = u * cn * omega

                # use lambda and the expected read count to define the number of trials (delta)
                # that should be drawn for each bin
                delta = theta * (1 - lamb) / lamb

                # replace all delta<1 values with 1 since delta should be >0
                # this avoids NaN errors when theta=0 at a given bin
                delta[delta<1] = 1

                reads = pyro.sample('reads', dist.NegativeBinomial(delta, probs=lamb), obs=data)


    def run_pert_model(self):
        if self.cuda:
            torch.set_default_tensor_type('torch.cuda.FloatTensor')

        logging.info('Loading data')

        logging.info('-' * 40)
        model_s = self.model_s
        model_g1 = self.model_g1

        cn_g1_reads_df, cn_g1_states_df, cn_s_reads_df, cn_s_states_df, \
            cn_g1_reads, cn_g1_states, cn_s_reads, cn_s_states, \
            gammas, rt_prior_profile, libs_g1, libs_s = self.process_input_data()

        # build transition matrix and cn prior for S-phase cells
        trans_mat = self.build_trans_mat(cn_g1_states)

        # compute consensus clone profiles for cn state
        clone_cn_profiles = compute_consensus_clone_profiles(
            self.cn_g1, self.cn_state_col, clone_col=self.clone_col, cell_col=self.cell_col, chr_col=self.chr_col,
            start_col=self.start_col, cn_state_col=self.cn_state_col
        )
        #import pickle
        #with open("/home/schnei01/scdna_replication_tools/debug1.pkl", "wb") as output_file:
        #    pickle.dump(cn_g1_reads_df, output_file)
        #with open("/home/schnei01/scdna_replication_tools/debug2.pkl", "wb") as output_file:
        #    pickle.dump(cn_s_reads_df, output_file)
        #with open("/home/schnei01/scdna_replication_tools/debug3.pkl", "wb") as output_file:
        #    pickle.dump(clone_cn_profiles, output_file)
        common_index = cn_s_reads_df.index.intersection(cn_g1_reads_df.index)
        temp = clone_cn_profiles["1"]
        a = temp.iloc[temp.index.isin(common_index)]
        clone_cn_profiles = pd.DataFrame(a, columns=pd.Index(data=["1"], name="clone_id"))

        if self.cn_prior_method == 'hmmcopy':
            # use hmmcopy states for the S-phase cells to build the prior
            etas = self.build_cn_prior(cn_s_states)
        elif self.cn_prior_method == 'g1_cells':
            # use G1-phase cell that has highest correlation to each S-phase cell as prior
            # raise ValueError("g1_cells method not implemented yet")

            cn_prior_input = torch.zeros(cn_s_states.shape)
            # loop through all S-phase cells
            for i, cell_id in enumerate(cn_s_reads_df.columns):
                cell_cn = self.cn_s.loc[self.cn_s[self.cell_col]==cell_id]  # get full cn data for this cell

                # shrink set of G1 cells to those in the matching clone if clone_id has already been assigned
                # for this S-phase cell
                if self.clone_col is not None:
                    clone_id = cell_cn[self.clone_col].values[0]
                    clone_cn_g1 = self.cn_g1.loc[self.cn_g1[self.clone_col]==clone_id]
                else:
                    clone_cn_g1 = self.cn_g1
                cell_cn = cell_cn[[self.chr_col, self.start_col, self.cell_col, self.input_col, self.cn_state_col]]
                
                # compute pearson correlations between this S-phase cell and all G1-phase cells in the same clone
                cell_corrs = compute_cell_corrs(cell_cn, clone_cn_g1, cell_id, col=self.input_col,
                                                cell_col=self.cell_col, chr_col=self.chr_col, start_col=self.start_col)
                
                # get data from the G1 cell that matches best
                g1_cell_id = cell_corrs.iloc[0].g1_cell_id
                g1_cell_cn = clone_cn_g1.loc[clone_cn_g1[self.cell_col]==g1_cell_id]

                # extract the cn_state profile from this best matching G1 cell
                cn_prior_input[:, i] = torch.tensor(g1_cell_cn[self.cn_state_col].values).to(torch.int64).to(torch.float32)

            # build a proper prior over num_states using the consensus clone cn calls for each cell
            etas = self.build_cn_prior(cn_prior_input)
        elif self.cn_prior_method == 'g1_clones':
            # use G1-phase clone that has highest correlation to each S-phase cell as prior
            etas = self.build_clone_cn_prior(self.cn_s, cn_s_reads_df, cn_s_states, clone_cn_profiles, weight=1e6)
        elif self.cn_prior_method == 'g1_composite':
            # use a composite of the principles in g1_clones and g1_cells to construct the prior
            etas = self.build_composite_cn_prior(cn_s_reads_df, clone_cn_profiles)
        elif self.cn_prior_method == 'diploid':
            # assume that every S-phase cell has a diploid prior
            num_loci, num_cells = cn_s_states.shape
            cn_s_diploid = torch.ones(num_loci, num_cells, self.P) * 2
            etas = self.build_cn_prior(cn_s_diploid)
        else:
            # assume uniform prior otherwise
            num_loci, num_cells = cn_s_states.shape
            etas = torch.ones(num_loci, num_cells, self.P) / self.P

        # fit GC params using G1-phase cells    
        guide_g1 = AutoDelta(poutine.block(model_g1, expose_fn=lambda msg: msg["name"].startswith("expose_")))

        optim = pyro.optim.Adam({'lr': self.learning_rate, 'betas': [0.8, 0.99]})
        elbo = JitTrace_ELBO(max_plate_nesting=2)

        svi = SVI(model_g1, guide_g1, optim, loss=elbo)

        # start inference
        logging.info('STEP 1: Learning reads to CN bias from low variance cells.')
        losses_g = []
        for i in range(self.max_iter_step1):
            loss = svi.step(gammas, libs_g1, cn=cn_g1_states, data=cn_g1_reads)

            losses_g.append(loss)
            logging.info('step: {}, loss: {}'.format(i, loss))

            # fancy convergence check that sees if the past 10 iterations have plateaued
            if i >= self.min_iter_step1:
                loss_diff = abs(max(losses_g[-10:-1]) - min(losses_g[-10:-1])) / abs(losses_g[0] - losses_g[-1])
                if loss_diff < self.rel_tol:
                    print('ELBO converged at iteration ' + str(i))
                    break
            
            # stop training if loss is NaN
            if np.isnan(loss):
                print('ELBO is NaN at iteration ' + str(i))
                break


        # replay model
        guide_trace_g1 = poutine.trace(guide_g1).get_trace(gammas, libs_g1, cn=cn_g1_states, data=cn_g1_reads)
        trained_model_g1 = poutine.replay(model_g1, trace=guide_trace_g1)

        # infer discrete sites and get model trace
        inferred_model_g1 = infer_discrete(
            trained_model_g1, temperature=0,
            first_available_dim=-3)
        trace_g1 = poutine.trace(inferred_model_g1).get_trace(gammas, libs_g1, cn=cn_g1_states, data=cn_g1_reads)

        # extract fitted parameters
        lambda_fit = trace_g1.nodes['expose_lambda']['value'].detach()
        betas_fit = trace_g1.nodes['expose_betas']['value'].detach()
        beta_means_fit = trace_g1.nodes['expose_beta_means']['value'].detach()
        beta_stds_fit = trace_g1.nodes['expose_beta_stds']['value'].detach()

        num_observations = float(cn_s_reads.shape[0] * cn_s_reads.shape[1])
        pyro.set_rng_seed(self.seed)
        pyro.clear_param_store()
        pyro.enable_validation(False)

        # condition gc betas of S-phase model using fitted results from G1-phase model
        model_s = poutine.condition(
            model_s,
            data={
                'expose_beta_means': beta_means_fit,
                'expose_beta_stds': beta_stds_fit,
            })

        # use manhattan binarization method to come up with an initial guess for each cell's time in S-phase
        t_init, t_alpha_prior, t_beta_prior = self.guess_times(cn_s_reads, etas)

        guide_s = AutoDelta(poutine.block(model_s, expose_fn=lambda msg: msg["name"].startswith("expose_")))
        optim_s = pyro.optim.Adam({'lr': self.learning_rate, 'betas': [0.8, 0.99]})
        elbo_s = JitTraceEnum_ELBO(max_plate_nesting=2)
        svi_s = SVI(model_s, guide_s, optim_s, loss=elbo_s)

        # start inference
        logging.info('STEP 2: Jointly infer replication and CN states in high variance cells.')
        losses_s = []
        for i in range(self.max_iter):
            loss = svi_s.step(gammas, libs_s, data=cn_s_reads, etas=etas, lamb=lambda_fit, t_init=t_init)

            losses_s.append(loss)
            logging.info('step: {}, loss: {}'.format(i, loss))

            # fancy convergence check that sees if the past 10 iterations have plateaued
            if i >= self.min_iter:
                loss_diff = abs(max(losses_s[-10:-1]) - min(losses_s[-10:-1])) / abs(losses_s[0] - losses_s[-1])
                if loss_diff < self.rel_tol:
                    print('ELBO converged at iteration ' + str(i))
                    break
            
            # stop training if loss is NaN
            if np.isnan(loss):
                print('ELBO is NaN at iteration ' + str(i))
                break


        # replay model
        guide_trace_s = poutine.trace(guide_s).get_trace(gammas, libs_s, data=cn_s_reads, etas=etas, lamb=lambda_fit, t_init=t_init)
        trained_model_s = poutine.replay(model_s, trace=guide_trace_s)

        # infer discrete sites and get model trace
        inferred_model_s = infer_discrete(
            trained_model_s, temperature=0,
            first_available_dim=-3)
        trace_s = poutine.trace(inferred_model_s).get_trace(gammas, libs_s, data=cn_s_reads, etas=etas, lamb=lambda_fit, t_init=t_init)

        # get output dataframes based on learned latent parameters and states
        cn_s_out, supp_s_out_df = self.package_s_output(self.cn_s, trace_s, cn_s_reads_df, lambda_fit, losses_g, losses_s)

        # run pre-trained S-phase model on cells labeled as G1/2-phase to see if
        # any of them are actually in S-phase
        if self.run_step3:
            # extract parameters learned in S-phase model that need to be conditioned
            rho_fit_s = trace_s.nodes['expose_rho']['value'].detach()
            a_fit_s = trace_s.nodes['expose_a']['value'].detach()
            
            pyro.clear_param_store()
            pyro.enable_validation(False)

            # condition gc betas of S-phase model using fitted results from G1-phase model
            model_s2 = poutine.condition(
                self.model_s,
                data={
                    'expose_rho': rho_fit_s,
                    'expose_a': a_fit_s,
                    'expose_beta_means': beta_means_fit,
                    'expose_beta_stds': beta_stds_fit,
                })

            # use clone pseudobulk profiles to set the CN prior of each cell
            etas2 = self.build_clone_cn_prior(self.cn_g1, cn_g1_reads_df, cn_g1_states, clone_cn_profiles, weight=1e6)

            # use manhattan binarization method to come up with an initial guess for each cell's time in S-phase
            t_init2, t_alpha_prior2, t_beta_prior2 = self.guess_times(cn_g1_reads, etas2)

            guide_s2 = AutoDelta(poutine.block(model_s2, expose_fn=lambda msg: msg["name"].startswith("expose_")))
            optim_s2 = pyro.optim.Adam({'lr': self.learning_rate, 'betas': [0.8, 0.99]})
            elbo_s2 = JitTraceEnum_ELBO(max_plate_nesting=2)
            svi_s2 = SVI(model_s2, guide_s2, optim_s2, loss=elbo_s2)

            # start inference
            logging.info('STEP 3: Running pre-trained S-phase model on low variance cells.')
            losses_s2 = []
            for i in range(self.max_iter_step3):
                loss = svi_s2.step(gammas, libs_g1, data=cn_g1_reads, etas=etas2, lamb=lambda_fit, t_init=t_init2)

                losses_s2.append(loss)
                logging.info('step: {}, loss: {}'.format(i, loss))

                # fancy convergence check that sees if the past 10 iterations have plateaued
                if i >= self.min_iter_step3:
                    loss_diff = abs(max(losses_s2[-10:-1]) - min(losses_s2[-10:-1])) / abs(losses_s2[0] - losses_s2[-1])
                    if loss_diff < self.rel_tol:
                        print('ELBO converged at iteration ' + str(i))
                        break
            
                # stop training if loss is NaN
                if np.isnan(loss):
                    print('ELBO is NaN at iteration ' + str(i))
                    break

            # replay model
            guide_trace_s2 = poutine.trace(guide_s2).get_trace(gammas, libs_g1, data=cn_g1_reads, etas=etas2, lamb=lambda_fit, t_init=t_init2)
            trained_model_s2 = poutine.replay(model_s2, trace=guide_trace_s2)

            # infer discrete sites and get model trace
            inferred_model_s2 = infer_discrete(
                trained_model_s2, temperature=0,
                first_available_dim=-3)
            trace_s2 = poutine.trace(inferred_model_s2).get_trace(gammas, libs_g1, data=cn_g1_reads, etas=etas2, lamb=lambda_fit, t_init=t_init2)

            # save results for these cells
            cn_g1_out, supp_g1_out_df = self.package_s_output(self.cn_g1, trace_s2, cn_g1_reads_df, lambda_fit, losses_g, losses_s2)
        else:
            cn_g1_out = None
            supp_g1_out_df = None 
        
        return cn_s_out, supp_s_out_df, cn_g1_out, supp_g1_out_df
