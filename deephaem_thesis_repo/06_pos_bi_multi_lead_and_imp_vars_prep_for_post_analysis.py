import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
import pybedtools
import scanpy as sc
import seaborn as sns
import h5py
from matplotlib import pyplot as plt
from sklearn import preprocessing
from sklearn.model_selection import train_test_split

from tensorflow.keras.layers import Dense, Flatten, Dropout, Input, Conv1D, MaxPooling1D, BatchNormalization
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint, CSVLogger
from tensorflow.keras.initializers import GlorotUniform
from tensorflow.keras.models import Sequential
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.regularizers import l2
from tensorflow.keras import backend

import tensorflow as tf
import keras

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from IPython.display import Image
from IPython.core.display import HTML 

##pos controls 10 vars

with_vars_present_pos = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_shortlist_vars_prediction/with_vars_present_new_41_pos_control_prediction.csv', sep='\t') #was: with_vars_present_pos_control_prediction.csv
with_vars_present_pos.drop(['Unnamed: 0'], axis=1, inplace=True)
without_vars_present_pos = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_shortlist_vars_prediction/without_vars_present_new_41_pos_control_lead_and_imp_prediction.csv', sep='\t') #was: without_vars_present_pos_control_prediction.csv
without_vars_present_pos.drop(['Unnamed: 0'], axis=1, inplace=True)

# subtract region1_forward from region1_reverse 
cols_without_var_pos = with_vars_present_pos.columns.tolist()
cols_without_var_pos.remove('variant')
prediction_df_pos = with_vars_present_pos[cols_without_var_pos].subtract(without_vars_present_pos[cols_without_var_pos]) 
prediction_df_pos['variant'] = with_vars_present_pos['variant']
#prediction_df.to_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_shortlist_vars_prediction/pos_control_vars_damage_scores.csv', index=False)

long_prediction_df_pos =  pd.melt(prediction_df_pos, id_vars = ['variant'], value_vars=cols_without_var_pos)

long_prediction_df_pos = long_prediction_df_pos[['variant', 'variable', 'value']]
long_prediction_df_pos['extremes'] = ''
long_prediction_df_pos['extremes'][long_prediction_df_pos['value']>=0.1] = 'score_above_+0.1'
long_prediction_df_pos['extremes'][long_prediction_df_pos['value']<=-0.1] = 'score_below_-0.1'
long_prediction_df_pos['extremes'][long_prediction_df_pos['extremes']==''] = 'not_empirically_extreme'
long_prediction_df_pos['var_type'] = 'pos_control'

#############################################################################################################################

##bi-allelic lead imp 

with_vars_present_bi = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_shortlist_vars_prediction/with_vars_present_bi_allelic_3548_lead_and_imp_prediction.csv', sep='\t')
with_vars_present_bi.drop(['Unnamed: 0'], axis=1, inplace=True)
without_vars_present_bi = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_shortlist_vars_prediction/without_vars_present_bi_allelic_3548_lead_and_imp_prediction.csv', sep='\t')
without_vars_present_bi.drop(['Unnamed: 0'], axis=1, inplace=True)

with_vars_present_bi['var_type'] = 'bi'
without_vars_present_bi['var_type'] = 'bi'

# subtract region1_forward from region1_reverse 
cols_without_var_bi = with_vars_present_bi.columns.tolist()
cols_without_var_bi.remove('variant')
cols_without_var_bi.remove('var_type')
cols_without_var_bi.remove('REF')
cols_without_var_bi.remove('ALT')
prediction_df_bi = with_vars_present_bi[cols_without_var_bi].subtract(without_vars_present_bi[cols_without_var_bi]) 
prediction_df_bi['variant'] = with_vars_present_bi['variant']
#prediction_df.to_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_shortlist_vars_prediction/with_vars_present_bi_allelic_3548_lead_and_imp_vars_damage_scores.csv', index=False)

long_prediction_df_bi =  pd.melt(prediction_df_bi, id_vars = ['variant'], value_vars=cols_without_var_bi)

long_prediction_df_bi = long_prediction_df_bi[['variant', 'variable', 'value']]
long_prediction_df_bi['extremes'] = ''
long_prediction_df_bi['extremes'][long_prediction_df_bi['value']>=0.1] = 'score_above_+0.1'
long_prediction_df_bi['extremes'][long_prediction_df_bi['value']<=-0.1] = 'score_below_-0.1'
long_prediction_df_bi['extremes'][long_prediction_df_bi['extremes']==''] = 'not_empirically_extreme'
long_prediction_df_bi['var_type'] = 'lead_and_imp'

#############################################################################################################################

##multi-allelic lead imp 

with_vars_present_multi = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_shortlist_vars_prediction/with_vars_present_multi_allelic_847_lead_and_imp_prediction.csv', sep='\t')
with_vars_present_multi.drop(['Unnamed: 0'], axis=1, inplace=True)
without_vars_present_multi = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_shortlist_vars_prediction/without_vars_present_multi_allelic_847_lead_and_imp_prediction.csv', sep='\t')
without_vars_present_multi.drop(['Unnamed: 0'], axis=1, inplace=True)

# subtract region1_forward from region1_reverse 
cols_without_var_multi = with_vars_present_multi.columns.tolist()
cols_without_var_multi.remove('variant')
cols_without_var_multi.remove('var_type')
cols_without_var_multi.remove('REF')
cols_without_var_multi.remove('ALT')
prediction_df_multi = with_vars_present_multi[cols_without_var_multi].subtract(without_vars_present_multi[cols_without_var_multi]) 
prediction_df_multi['variant'] = with_vars_present_multi['variant']
#prediction_df.to_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_shortlist_vars_prediction/with_vars_present_multi_allelic_847_lead_and_imp_vars_damage_scores.csv', index=False)

long_prediction_df_multi =  pd.melt(prediction_df_multi, id_vars = ['variant'], value_vars=cols_without_var_multi)

long_prediction_df_multi = long_prediction_df_multi[['variant', 'variable', 'value']]
long_prediction_df_multi['extremes'] = ''
long_prediction_df_multi['extremes'][long_prediction_df_multi['value']>=0.1] = 'score_above_+0.1'
long_prediction_df_multi['extremes'][long_prediction_df_multi['value']<=-0.1] = 'score_below_-0.1'
long_prediction_df_multi['extremes'][long_prediction_df_multi['extremes']==''] = 'not_empirically_extreme'
long_prediction_df_multi['var_type'] = 'lead_and_imp'

#############################################################################################################################

##going on with them all concatenated together

prediction_df_bi['type'] = 'lead_and_imp'
prediction_df_multi['type'] = 'lead_and_imp'
prediction_df_pos['type'] = 'pos_control'
prediction_df_all = pd.concat([prediction_df_pos, prediction_df_bi, prediction_df_multi])
#prediction_df_all.to_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_shortlist_vars_prediction/with_vars_present_new_41_pos_controls_and_4395_lead_imp_vars_damage_scores.csv', index=False) #was: with_vars_present_pos_controls_and_4395_lead_imp_vars_damage_scores

long_prediction_df = pd.concat([long_prediction_df_pos, long_prediction_df_bi, long_prediction_df_multi])

#previously when 10 regulomedb pos controls inclusion:
# long_prediction_df['value'].mean()
# long_prediction_df['value'].std()
# #Mean is 6.684756802476176e-06
# #Std deviation is 0.0066753582243595806
# #2SD = 0.013350716448719161
# #3SD = 0.020026074673078743
# #Mean - SD = -0.006668673467557104
# #Mean + SD = 0.006682042981162057
# #Mean – 2SD = -0.013344031691916685
# #Mean + 2SD = 0.013357401205521638
# #Mean – 3SD = -0.02001938991627627
# #Mean + 3SD = 0.020032759429881218

#now with 41 new literature pos controls inclusion:
long_prediction_df['value'].mean()
long_prediction_df['value'].std()
#Mean is 2.0969435536983298e-05
#Std deviation is 0.006576753851799106
#2SD = 0.013153507703598213
#3SD = 0.019730261555397317
#Mean - SD = -0.006555784416262123
#Mean + SD = 0.00659772328733609
#Mean – 2SD = -0.013132538268061229
#Mean + 2SD = 0.013174477139135197
#Mean – 3SD = -0.019709292119860333
#Mean + 3SD = 0.0197512309909343

long_prediction_df['mean'] = long_prediction_df['value'].mean()
long_prediction_df['std'] = long_prediction_df['value'].std()
long_prediction_df['mean-std'] = long_prediction_df['mean']-long_prediction_df['std']
long_prediction_df['mean+std'] = long_prediction_df['mean']+long_prediction_df['std']
long_prediction_df['mean-2std'] = long_prediction_df['mean']-2*(long_prediction_df['std'])
long_prediction_df['mean+2std'] = long_prediction_df['mean']+2*(long_prediction_df['std'])
long_prediction_df['mean-3std'] = long_prediction_df['mean']-3*(long_prediction_df['std'])
long_prediction_df['mean+3std'] = long_prediction_df['mean']+3*(long_prediction_df['std'])
long_prediction_df['empirical_loss'] = -0.1
long_prediction_df['empirical_gain'] = 0.1

long_prediction_df['<=mean-std'] = ''
long_prediction_df['>=mean+std'] = ''
long_prediction_df['<=mean-2std'] = ''
long_prediction_df['>=mean+2std'] = ''
long_prediction_df['<=mean-3std'] = ''
long_prediction_df['>=mean+3std'] = ''
long_prediction_df['<=empirical_loss'] = ''
long_prediction_df['>=empirical_gain'] = ''

long_prediction_df['<=mean-std'][long_prediction_df['value']<=long_prediction_df['mean-std']] = 'Y'
long_prediction_df['>=mean+std'][long_prediction_df['value']>=long_prediction_df['mean+std']] = 'Y'
long_prediction_df['<=mean-2std'][long_prediction_df['value']<=long_prediction_df['mean-2std']] = 'Y'
long_prediction_df['>=mean+2std'][long_prediction_df['value']>=long_prediction_df['mean+2std']] = 'Y'
long_prediction_df['<=mean-3std'][long_prediction_df['value']<=long_prediction_df['mean-3std']] = 'Y'
long_prediction_df['>=mean+3std'][long_prediction_df['value']>=long_prediction_df['mean+3std']] = 'Y'
long_prediction_df['<=empirical_loss'][long_prediction_df['value']<=long_prediction_df['empirical_loss']] = 'Y'
long_prediction_df['>=empirical_gain'][long_prediction_df['value']>=long_prediction_df['empirical_gain']] = 'Y'

long_prediction_df.rename({'variant':'rsid'}, axis=1, inplace=True)

####

lead_imp_bi_snp_dataset = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/bi_allelic_3548_lead_and_imp_snp_dataset.csv')

lead_imp_multi_snp_dataset = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/multi_allelic_847_lead_and_imp_snp_dataset.csv')

lead_imp_snp_dataset = pd.concat([lead_imp_bi_snp_dataset, lead_imp_multi_snp_dataset])

merged = long_prediction_df.merge(lead_imp_snp_dataset[['rsid','chr','grch38_POS', 'REF', 'ALT']], how='left', on='rsid')

null_part = merged[merged['REF'].isnull()]
nonull_part = merged[merged['REF'].notnull()]
null_part[['chr', 'grch38_POS', 'REF', 'ALT']] = null_part['rsid'].str.split('_', expand=True)
new_merged = pd.concat([nonull_part, null_part])

new_merged['chr'] = new_merged['chr'].astype(int)
new_merged['grch38_POS'] = new_merged['grch38_POS'].astype(int)

new_merged_filt = new_merged[['rsid', 'chr', 'grch38_POS', 'REF', 'ALT', 'variable', 'value', '<=mean-std', '>=mean+std', '<=mean-2std', '>=mean+2std', '<=mean-3std', '>=mean+3std', '<=empirical_loss', '>=empirical_gain']]

new_merged_filt['strength'] = ''

new_rows = []
for index, row in new_merged_filt.iterrows():
    if row['<=empirical_loss']=='Y' or row['>=empirical_gain']=='Y':
        row['strength'] = 'empirical'
    elif row['<=mean-3std']=='Y' or row['>=mean+3std']=='Y':
        row['strength'] = '99.7%'
    elif row['<=mean-2std']=='Y' or row['>=mean+2std']=='Y':
        row['strength'] = '95%'
    elif row['<=mean-std']=='Y' or row['>=mean+std']=='Y':
        row['strength'] = '68%'
    else:
        row['strength'] = 'n/a'
    new_rows.append(row.values)
dff = pd.DataFrame(new_rows, columns=new_merged_filt.columns).reset_index()

dff.to_csv('/gpfs3/well/PROCARDIS/domwest/deephaem_prep/post_analysis/various_strength_post_analysis_pos_and_lead_and_imp_vars_damage_scores.csv', index=False)

##################################################################################################

#then wrangled this more to get list of rsids for input to opentarget in rescomp:

import pandas as pd
import os
import numpy as np
df = pd.read_csv('../various_strength_post_analysis_pos_and_lead_and_imp_vars_damage_scores.csv')
df_subset = df[df['strength'].isin(['95%','99.7%','empirical'])]
df_subset_filt = df_subset[df_subset['rsid'].str.startswith('rs')]
len(df_subset_filt['rsid'].unique())
rsids_of_interest = df_subset_filt['rsid'].unique().tolist()
#1019 unique snps
with open('snps_input_to_opentarget.txt', mode='wt', encoding='utf-8') as myfile:
    myfile.write('\n'.join(rsids_of_interest))


##ran /well/PROCARDIS/domwest/deephaem_prep/post_analysis/opentarget/testing_dom_runSNPsearchv3.pl

#now read in input files:

import glob
import pandas as pd
import os
import numpy as np

top_rows = []
for file in glob.glob("opentarget_v2g_*.out.txt"):
    print(file)
    df = pd.read_csv(file, sep='\t')
    df.sort_values(by='Overall V2G', ascending=False, inplace=True)
    df_top = df.iloc[0]
    top_rows.append(df_top.values)
top_df = pd.DataFrame(top_rows, columns=df.columns).reset_index()
top_df.drop(['index'], axis=1, inplace=True)
top_df.to_csv('/gpfs3/well/PROCARDIS/domwest/deephaem_prep/post_analysis/opentarget/opentarget_snp_to_gene.csv', index=False)

#in this top_df there are 85 unique genes:
['KYAT3','CRHR1','MLIP','NMB','PCMT1','PLEKHM1','ADAMTS7','MLST8','MAPT','RAPGEF1','GINM1','SLC23A1','SIPA1L1','HSPB7','NUP43','KANSL1','FAM167A','COPB1','ADPRHL1','ROCK2','STRN','SEMA3G','COL4A2','MITF','TCF7L2','XPO7','SYNPO2L','NCOA7','SLC6A6','SPATS2L','RRAS2','SPATA24','PRKCA','CCDC141','CEP68','TTC8','LSM3','BLK','TBX3','DNAJC18','GLYCTK','PROX1','ACTBL2','ALPK3','DMPK','CDKN1A','SVIL','LATS1','TPGS2','PROB1','PDE3A','LYZL1','LRP11','ERO1A','E2F6','ARHGAP27','MTSS1','CCDC136','FDFT1','TTN','MFSD9','SRL','PGR','PLN','VPREB3','EML5','ADSL','CLCNKA','KLF6','MAPRE2','BAG3','SPATA7','C1QTNF4','FERMT2','SOX8','CCT8','DNAJB4','PKN2','PROM1','ZC3H14','SSPN','LMF1','FNDC3B','FBXO32','HEY2']