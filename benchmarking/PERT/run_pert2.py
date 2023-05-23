import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scdna_replication_tools.infer_scRT import scRT
from scdna_replication_tools.plot_pert_output import plot_model_results
from scdna_replication_tools.compute_ccc_features import compute_ccc_features
from scdna_replication_tools.predict_cycle_phase import predict_cycle_phase
from sklearn.metrics import confusion_matrix
import random
import re
import sys
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
random.seed(2023)

cn = pd.read_csv('UID-FML-PEO1-FUCCI-DOWNSAMPLED-EQUAL-FORDLP_reads.csv.gz', dtype={'chr': str})
cn['clone_id'] = '1'
cn['true_cellcycle'] = [re.search('193-(.*)-Plate',i).group(1) for i in cn.cell_id.values]
cn['library_id'] = ["SLX-"+re.search('SLX-(.*)_000',i).group(1) for i in cn.cell_id.values]
cn = cn[cn['valid']]
cn = cn[cn['chr'] != "Y"]
cn = cn[cn['chr'] != "X"]
n_train_samples = 200
n_test_samples = 85

# create groundtruth
g1_rows = cn[cn['true_cellcycle'] == 'G1']
s_rows = cn[cn['true_cellcycle'] == 'S']

cellnames = list(cn.cell_id.unique())
cellcycle = [re.search('193-(.*)-Plate',i).group(1) for i in cellnames]
numbers = pd.Series(cellcycle).value_counts()
print(numbers)

train_g1_set = random.sample([cellnames[i] for i,_ in enumerate(cellcycle) if cellcycle[i] == "G1"], n_train_samples)
test_g1_set = random.sample([cellnames[i] for i,_ in enumerate(cellcycle) if cellcycle[i] == "G1" and cellnames[i] not in train_g1_set], n_test_samples)
test_s_set = random.sample([cellnames[i] for i,_ in enumerate(cellcycle) if cellcycle[i] == "S"], n_test_samples)

train_g1 = cn[cn["cell_id"].isin(train_g1_set)]
test_g1 = cn[cn["cell_id"].isin(test_g1_set)]
test_s = cn[cn["cell_id"].isin(test_s_set)]
#test_g1_set = g1_rows[g1_rows['true_cellcycle'] == 'G1'][~g1_rows[g1_rows['true_cellcycle'] == 'G1'].index.isin(train_g1_set.index)]

cn = pd.concat([test_g1, test_s])

# temporarily remove columns that don't get used by PERT
temp_cn = cn[['cell_id', 'chr', 'start', 'end', 'gc', 'state', 'clone_id', 'library_id', 'copy']]
temp_cn_g1 = train_g1[['cell_id', 'chr', 'start', 'end', 'gc', 'clone_id', 'state', 'library_id', 'copy']]
temp_cn['chr_col'] = temp_cn['chr']
temp_cn_g1['chr_col'] = temp_cn_g1['chr']
# add the replication columns for the G1-phase cells
temp_cn_g1['true_rep'] = 0.0
temp_cn_g1['true_p_rep'] = 0.0
temp_cn_g1['true_t'] = 1.0

print('number of loci in S-phase: ', len(temp_cn[['chr', 'start']].drop_duplicates()))
print('number of loci in G1-phase: ', len(temp_cn_g1[['chr', 'start']].drop_duplicates()))

with open(r"/home/schnei01/scdna_replication_tools/checklist.pkl", "rb") as input_file:
    import pickle
    checklist = pickle.load(input_file)

temp_cn = temp_cn.loc[list(map(lambda x: x in checklist, zip(temp_cn['chr'], temp_cn['start']))),:]
temp_cn_g1 = temp_cn_g1.loc[list(map(lambda x: x in checklist, zip(temp_cn_g1['chr'], temp_cn_g1['start']))),:]

print('number of cells in S-phase: ', len(temp_cn['cell_id'].unique()))
print('number of cells in G1-phase: ', len(temp_cn_g1['cell_id'].unique()))
print('number of loci in S-phase: ', len(temp_cn[['chr', 'start']].drop_duplicates()))
print('number of loci in G1-phase: ', len(temp_cn_g1[['chr', 'start']].drop_duplicates()))

# Combine the three lists into one column
membership = train_g1_set + test_g1_set + test_s_set
list_indicator = ['Train'] * len(train_g1_set) + ['Test_G1'] * len(test_g1_set) + ['Test_S'] * len(test_s_set)
data = {'Name': membership, 'Set': list_indicator}
df = pd.DataFrame(data)
df.to_csv('membership_200.tsv', sep='\t', index=False)

#NOTE: not possible to run max_iter > 500 for step 2, at 0.01 learning rate -> program crashes with NaN everywhere
scrt = scRT(temp_cn, temp_cn_g1, input_col='copy', clone_col='clone_id', assign_col='copy', rt_prior_col=None,
            cn_state_col='state', gc_col='gc', cn_prior_method='g1_clones', max_iter=500,
            rel_tol=1e-6, learning_rate=0.01,
            max_iter_step1=5000, max_iter_step3=5000)

#  scrt = scRT(cn, g1_rows, input_col='reads', clone_col='clone_id', assign_col='reads', rt_prior_col=None,
            #  cn_state_col='state', gc_col='gc', cn_prior_method='g1_composite', max_iter=300)

# run inference using PERT
cn_s_with_scrt, supp_s_output, cn_g_with_scrt, supp_g_output = scrt.infer(level='pyro')


# concatenate the two dataframes containing PERT output
cn_out = pd.concat([cn_s_with_scrt, cn_g_with_scrt], ignore_index=True)

# predict the cycle phase for each clone based on the PERT output with default parameters
cn_s_out, cn_g_out, cn_lq_out = predict_cycle_phase(
    cn_out, frac_rt_col='cell_frac_rep', rep_state_col='model_rep_state',
    cn_state_col='model_cn_state', rpm_col='copy')

cn_all = pd.concat([cn_s_out, cn_g_out, cn_lq_out], ignore_index=True)
result = cn_all[["cell_id", "PERT_phase"]].drop_duplicates()
result['true_cellcycle'] = [re.search('193-(.*)-Plate',i).group(1) for i in result.cell_id.values]

import pickle
with open("/home/schnei01/scdna_replication_tools/result_200.pkl", "wb") as output_file:
    pickle.dump(result, output_file)

confusion_matrix = pd.crosstab(result['true_cellcycle'], result['PERT_phase'], rownames=['Actual'], colnames=['Predicted'])
print(confusion_matrix)
