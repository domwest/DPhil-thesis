#building matrices with indiv base snp coverage, whether in peak or not, ref prediction score, whether ref prediction score is >=0.8 or <= 0.2, damage score, and strength of damage score...
#for eryth snps (enhancer and promoter) as well as hcm snps

#OG ie endothelial classifiers

import glob
import os 
import numpy as np 
import pandas as pd 

####################################

huvec_ref_and_damage_scores_gain = []
gain_dfs = glob.glob('/well/PROCARDIS/domwest/deepmind/huvec_gain_prediction_dfs/*')
for df in gain_dfs:
	tmp = pd.read_csv(df)
	tmp.drop(['Unnamed: 0'], axis=1, inplace=True)
	tmp.rename({'classifier':'variable', 'pos':'grch38_POS', 'ref':'REF', 'alt':'ALT', 'damage_score':'damage_score_value', 'ref_score':'ref_score_value', 'alt_score':'alt_score_value'}, axis=1, inplace=True)
	tmp['chr'] = tmp['chr'].str.replace('chr','')
	tmp['var_type'] = 'gain'
	huvec_ref_and_damage_scores_gain.append(tmp)
huvec_ref_and_damage_scores_gain = pd.concat(huvec_ref_and_damage_scores_gain)

huvec_ref_and_damage_scores_loss = []
loss_dfs = glob.glob('/well/PROCARDIS/domwest/deepmind/huvec_loss_prediction_dfs/*')
for df in loss_dfs:
	tmp = pd.read_csv(df)
	tmp.drop(['Unnamed: 0'], axis=1, inplace=True)
	tmp.rename({'classifier':'variable', 'pos':'grch38_POS', 'ref':'REF', 'alt':'ALT', 'damage_score':'damage_score_value', 'ref_score':'ref_score_value', 'alt_score':'alt_score_value'}, axis=1, inplace=True)
	tmp['chr'] = tmp['chr'].str.replace('chr','')
	tmp['var_type'] = 'loss'
	huvec_ref_and_damage_scores_loss.append(tmp)
huvec_ref_and_damage_scores_loss = pd.concat(huvec_ref_and_damage_scores_loss)

huvec_ref_and_damage_scores_neutral = []
neutral_dfs = glob.glob('/well/PROCARDIS/domwest/deepmind/huvec_neutral_prediction_dfs/*')
for df in neutral_dfs:
	tmp = pd.read_csv(df)
	tmp.drop(['Unnamed: 0'], axis=1, inplace=True)
	tmp.rename({'classifier':'variable', 'pos':'grch38_POS', 'ref':'REF', 'alt':'ALT', 'damage_score':'damage_score_value', 'ref_score':'ref_score_value', 'alt_score':'alt_score_value'}, axis=1, inplace=True)
	tmp['chr'] = tmp['chr'].str.replace('chr','')
	tmp['var_type'] = 'neutral'
	huvec_ref_and_damage_scores_neutral.append(tmp)
huvec_ref_and_damage_scores_neutral = pd.concat(huvec_ref_and_damage_scores_neutral)

huvec_ref_and_damage_scores = pd.concat([huvec_ref_and_damage_scores_gain, huvec_ref_and_damage_scores_loss, huvec_ref_and_damage_scores_neutral])

huvec_ref_and_damage_scores['ref_score_ann'] = ''
huvec_ref_and_damage_scores['ref_score_ann'][huvec_ref_and_damage_scores['ref_score_value']>=0.8] = 'likely_open'
huvec_ref_and_damage_scores['ref_score_ann'][huvec_ref_and_damage_scores['ref_score_value']<=0.2] = 'likely_closed'
huvec_ref_and_damage_scores['ref_score_ann'][huvec_ref_and_damage_scores['ref_score_value']==''] = 'unsure'

####################################

#####***MERGE RSIDS ONTO hcm_ref_and_damage_scores 
#READ IN INDEX FILES
huvec_gain_index_file = pd.read_csv('../huvec_gain_index_file.txt', sep='\t', header=None)
huvec_loss_index_file = pd.read_csv('../huvec_loss_index_file.txt', sep='\t', header=None)
huvec_neutral_index_file = pd.read_csv('../huvec_neutral_index_file.txt', sep='\t', header=None)
huvec_all_index_files = pd.concat([huvec_gain_index_file, huvec_loss_index_file, huvec_neutral_index_file])

huvec_all_index_files['snp_test_id'] = huvec_all_index_files[2].str.replace('chr','') + '_' + huvec_all_index_files[3].astype(str) + '_' + huvec_all_index_files[4] + '_' + huvec_all_index_files[5]
huvec_all_index_files.rename({1:'rsid'}, axis=1, inplace=True)
huvec_ref_and_damage_scores = huvec_ref_and_damage_scores.merge(huvec_all_index_files[['rsid', 'snp_test_id']], how='left', on='snp_test_id')

huvec_ref_and_damage_scores.rename({'variable':'classifier'}, axis=1, inplace=True)

huvec_ref_and_damage_scores.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/total_huvec.csv', index=False)

# for clas in hcm_ref_and_damage_scores['classifier'].unique():
# 	classif = clas.replace(':','_')
# 	classif = classif.replace(' ','_')
# 	tmp = hcm_ref_and_damage_scores[hcm_ref_and_damage_scores['classifier']==clas]
# 	tmp['var_score_value'] = tmp['damage_score_value'] + tmp['ref_score_value']
# 	# tmp.sort_values(by=['cov'], ascending=False)
# 	tmp.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/' + classif + '_total_hcm.csv', index=False)
# 	tmp_top_gain = tmp[(tmp['ref_score_value']>=0.8) & (tmp['var_score_value']>tmp['ref_score_value'])]
# 	tmp_top_gain.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/' + classif + '_total_hcm_top_gain.csv', index=False)
# 	tmp_top_loss = tmp[(tmp['ref_score_value']>=0.8) & (tmp['var_score_value']<tmp['ref_score_value'])]
# 	tmp.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/' + classif + '_total_hcm_top_loss.csv', index=False)
# 	tmp_bottom_gain = tmp[(tmp['ref_score_value']<=0.2) & (tmp['var_score_value']>tmp['ref_score_value'])]
# 	tmp.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/' + classif + '_total_hcm_bottom_gain.csv', index=False)

########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################

#new neutral set of huvec snps (more neutral snps than previously)

import glob
import os 
import numpy as np 
import pandas as pd 

####################################

huvec_ref_and_damage_scores_neutral = []
neutral_dfs = glob.glob('/well/PROCARDIS/domwest/deepmind/huvec_more_neutral_prediction_dfs/*')
for df in neutral_dfs:
	tmp = pd.read_csv(df)
	tmp.drop(['Unnamed: 0'], axis=1, inplace=True)
	tmp.rename({'classifier':'variable', 'pos':'grch38_POS', 'ref':'REF', 'alt':'ALT', 'damage_score':'damage_score_value', 'ref_score':'ref_score_value', 'alt_score':'alt_score_value'}, axis=1, inplace=True)
	tmp['chr'] = tmp['chr'].str.replace('chr','')
	tmp['var_type'] = 'neutral'
	huvec_ref_and_damage_scores_neutral.append(tmp)
huvec_ref_and_damage_scores_neutral = pd.concat(huvec_ref_and_damage_scores_neutral)

huvec_ref_and_damage_scores = huvec_ref_and_damage_scores_neutral

huvec_ref_and_damage_scores['ref_score_ann'] = ''
huvec_ref_and_damage_scores['ref_score_ann'][huvec_ref_and_damage_scores['ref_score_value']>=0.8] = 'likely_open'
huvec_ref_and_damage_scores['ref_score_ann'][huvec_ref_and_damage_scores['ref_score_value']<=0.2] = 'likely_closed'
huvec_ref_and_damage_scores['ref_score_ann'][huvec_ref_and_damage_scores['ref_score_value']==''] = 'unsure'

####################################

#####***MERGE RSIDS ONTO hcm_ref_and_damage_scores 
#READ IN INDEX FILES
huvec_all_index_files = pd.read_csv('huvec_more_neutral_index_file.txt', sep='\t', header=None)

huvec_all_index_files['snp_test_id'] = huvec_all_index_files[2].str.replace('chr','') + '_' + huvec_all_index_files[3].astype(str) + '_' + huvec_all_index_files[4] + '_' + huvec_all_index_files[5]
huvec_all_index_files.rename({1:'rsid'}, axis=1, inplace=True)
huvec_ref_and_damage_scores = huvec_ref_and_damage_scores.merge(huvec_all_index_files[['rsid', 'snp_test_id']], how='left', on='snp_test_id')

huvec_ref_and_damage_scores.rename({'variable':'classifier'}, axis=1, inplace=True)

huvec_ref_and_damage_scores.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/total_huvec_only_more_neutrals.csv', index=False)

# for clas in hcm_ref_and_damage_scores['classifier'].unique():
# 	classif = clas.replace(':','_')
# 	classif = classif.replace(' ','_')
# 	tmp = hcm_ref_and_damage_scores[hcm_ref_and_damage_scores['classifier']==clas]
# 	tmp['var_score_value'] = tmp['damage_score_value'] + tmp['ref_score_value']
# 	# tmp.sort_values(by=['cov'], ascending=False)
# 	tmp.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/' + classif + '_total_hcm.csv', index=False)
# 	tmp_top_gain = tmp[(tmp['ref_score_value']>=0.8) & (tmp['var_score_value']>tmp['ref_score_value'])]
# 	tmp_top_gain.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/' + classif + '_total_hcm_top_gain.csv', index=False)
# 	tmp_top_loss = tmp[(tmp['ref_score_value']>=0.8) & (tmp['var_score_value']<tmp['ref_score_value'])]
# 	tmp.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/' + classif + '_total_hcm_top_loss.csv', index=False)
# 	tmp_bottom_gain = tmp[(tmp['ref_score_value']<=0.2) & (tmp['var_score_value']>tmp['ref_score_value'])]
# 	tmp.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/' + classif + '_total_hcm_bottom_gain.csv', index=False)

########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################

#Different classifiers ie other celltypes

import glob
import os 
import numpy as np 
import pandas as pd 

####################################

huvec_ref_and_damage_scores_gain = []
gain_dfs = glob.glob('/well/PROCARDIS/domwest/deepmind/huvec_gain_prediction_dfs/*')
for df in gain_dfs:
	tmp = pd.read_csv(df)
	tmp.drop(['Unnamed: 0'], axis=1, inplace=True)
	tmp.rename({'classifier':'variable', 'pos':'grch38_POS', 'ref':'REF', 'alt':'ALT', 'damage_score':'damage_score_value', 'ref_score':'ref_score_value', 'alt_score':'alt_score_value'}, axis=1, inplace=True)
	tmp['chr'] = tmp['chr'].str.replace('chr','')
	tmp['var_type'] = 'gain'
	huvec_ref_and_damage_scores_gain.append(tmp)
huvec_ref_and_damage_scores_gain = pd.concat(huvec_ref_and_damage_scores_gain)

huvec_ref_and_damage_scores_loss = []
loss_dfs = glob.glob('/well/PROCARDIS/domwest/deepmind/huvec_loss_prediction_dfs/*')
for df in loss_dfs:
	tmp = pd.read_csv(df)
	tmp.drop(['Unnamed: 0'], axis=1, inplace=True)
	tmp.rename({'classifier':'variable', 'pos':'grch38_POS', 'ref':'REF', 'alt':'ALT', 'damage_score':'damage_score_value', 'ref_score':'ref_score_value', 'alt_score':'alt_score_value'}, axis=1, inplace=True)
	tmp['chr'] = tmp['chr'].str.replace('chr','')
	tmp['var_type'] = 'loss'
	huvec_ref_and_damage_scores_loss.append(tmp)
huvec_ref_and_damage_scores_loss = pd.concat(huvec_ref_and_damage_scores_loss)

huvec_ref_and_damage_scores_neutral = []
neutral_dfs = glob.glob('/well/PROCARDIS/domwest/deepmind/huvec_neutral_prediction_dfs/*')
for df in neutral_dfs:
	tmp = pd.read_csv(df)
	tmp.drop(['Unnamed: 0'], axis=1, inplace=True)
	tmp.rename({'classifier':'variable', 'pos':'grch38_POS', 'ref':'REF', 'alt':'ALT', 'damage_score':'damage_score_value', 'ref_score':'ref_score_value', 'alt_score':'alt_score_value'}, axis=1, inplace=True)
	tmp['chr'] = tmp['chr'].str.replace('chr','')
	tmp['var_type'] = 'neutral'
	huvec_ref_and_damage_scores_neutral.append(tmp)
huvec_ref_and_damage_scores_neutral = pd.concat(huvec_ref_and_damage_scores_neutral)

huvec_ref_and_damage_scores = pd.concat([huvec_ref_and_damage_scores_gain, huvec_ref_and_damage_scores_loss, huvec_ref_and_damage_scores_neutral])

subset = huvec_ref_and_damage_scores[huvec_ref_and_damage_scores['variable'].isin(['DNASE:stomach_female_adult_(53_years)', 'DNASE:epidermal_melanocyte', 'DNASE:frontal_cortex_male_adult_(27_years)_and_male_adult_(35_years)', 'DNASE:skeletal_muscle_cell'])]

subset['ref_score_ann'] = ''
subset['ref_score_ann'][subset['ref_score_value']>=0.8] = 'likely_open'
subset['ref_score_ann'][subset['ref_score_value']<=0.2] = 'likely_closed'
subset['ref_score_ann'][subset['ref_score_value']==''] = 'unsure'

####################################

#####***MERGE RSIDS ONTO hcm_ref_and_damage_scores 
#READ IN INDEX FILES
huvec_gain_index_file = pd.read_csv('huvec_gain_index_file.txt', sep='\t', header=None)
huvec_loss_index_file = pd.read_csv('huvec_loss_index_file.txt', sep='\t', header=None)
huvec_neutral_index_file = pd.read_csv('huvec_neutral_index_file.txt', sep='\t', header=None)
huvec_all_index_files = pd.concat([huvec_gain_index_file, huvec_loss_index_file, huvec_neutral_index_file])

huvec_all_index_files['snp_test_id'] = huvec_all_index_files[2].str.replace('chr','') + '_' + huvec_all_index_files[3].astype(str) + '_' + huvec_all_index_files[4] + '_' + huvec_all_index_files[5]
huvec_all_index_files.rename({1:'rsid'}, axis=1, inplace=True)
subset = subset.merge(huvec_all_index_files[['rsid', 'snp_test_id']], how='left', on='snp_test_id')

subset.rename({'variable':'classifier'}, axis=1, inplace=True)

subset.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/total_huvec_other_celltypes.csv', index=False)

# for clas in hcm_ref_and_damage_scores['classifier'].unique():
# 	classif = clas.replace(':','_')
# 	classif = classif.replace(' ','_')
# 	tmp = hcm_ref_and_damage_scores[hcm_ref_and_damage_scores['classifier']==clas]
# 	tmp['var_score_value'] = tmp['damage_score_value'] + tmp['ref_score_value']
# 	# tmp.sort_values(by=['cov'], ascending=False)
# 	tmp.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/' + classif + '_total_hcm.csv', index=False)
# 	tmp_top_gain = tmp[(tmp['ref_score_value']>=0.8) & (tmp['var_score_value']>tmp['ref_score_value'])]
# 	tmp_top_gain.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/' + classif + '_total_hcm_top_gain.csv', index=False)
# 	tmp_top_loss = tmp[(tmp['ref_score_value']>=0.8) & (tmp['var_score_value']<tmp['ref_score_value'])]
# 	tmp.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/' + classif + '_total_hcm_top_loss.csv', index=False)
# 	tmp_bottom_gain = tmp[(tmp['ref_score_value']<=0.2) & (tmp['var_score_value']>tmp['ref_score_value'])]
# 	tmp.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/' + classif + '_total_hcm_bottom_gain.csv', index=False)


########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################

#HUVEC full set snps

import glob
import os 
import numpy as np 
import pandas as pd 

####################################

huvec_ref_and_damage_scores_full_set = []
full_set_dfs = glob.glob('/well/PROCARDIS/domwest/deepmind/huvec_full_set_prediction_dfs/*')
for df in full_set_dfs:
	tmp = pd.read_csv(df)
	tmp.drop(['Unnamed: 0'], axis=1, inplace=True)
	tmp.rename({'classifier':'variable', 'pos':'grch38_POS', 'ref':'REF', 'alt':'ALT', 'damage_score':'damage_score_value', 'ref_score':'ref_score_value', 'alt_score':'alt_score_value'}, axis=1, inplace=True)
	tmp['chr'] = tmp['chr'].str.replace('chr','')
	tmp['var_type'] = 'full_set'
	huvec_ref_and_damage_scores_full_set.append(tmp)
huvec_ref_and_damage_scores_full_set = pd.concat(huvec_ref_and_damage_scores_full_set)

huvec_ref_and_damage_scores = huvec_ref_and_damage_scores_full_set

huvec_ref_and_damage_scores['ref_score_ann'] = ''
huvec_ref_and_damage_scores['ref_score_ann'][huvec_ref_and_damage_scores['ref_score_value']>=0.8] = 'likely_open'
huvec_ref_and_damage_scores['ref_score_ann'][huvec_ref_and_damage_scores['ref_score_value']<=0.2] = 'likely_closed'
huvec_ref_and_damage_scores['ref_score_ann'][huvec_ref_and_damage_scores['ref_score_value']==''] = 'unsure'

####################################

#####***MERGE RSIDS ONTO hcm_ref_and_damage_scores 
#READ IN INDEX FILES
huvec_full_set_index_file = pd.read_csv('huvec_full_set_index_file.txt', sep='\t', header=None)
huvec_all_index_files = huvec_full_set_index_file

huvec_all_index_files['snp_test_id'] = huvec_all_index_files[2].str.replace('chr','') + '_' + huvec_all_index_files[3].astype(str) + '_' + huvec_all_index_files[4] + '_' + huvec_all_index_files[5]
huvec_all_index_files.rename({1:'rsid'}, axis=1, inplace=True)
huvec_ref_and_damage_scores = huvec_ref_and_damage_scores.merge(huvec_all_index_files[['rsid', 'snp_test_id']], how='left', on='snp_test_id')

huvec_ref_and_damage_scores.rename({'variable':'classifier'}, axis=1, inplace=True)

huvec_ref_and_damage_scores.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/total_huvec_full_set.csv', index=False)


# for clas in hcm_ref_and_damage_scores['classifier'].unique():
# 	classif = clas.replace(':','_')
# 	classif = classif.replace(' ','_')
# 	tmp = hcm_ref_and_damage_scores[hcm_ref_and_damage_scores['classifier']==clas]
# 	tmp['var_score_value'] = tmp['damage_score_value'] + tmp['ref_score_value']
# 	# tmp.sort_values(by=['cov'], ascending=False)
# 	tmp.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/' + classif + '_total_hcm.csv', index=False)
# 	tmp_top_gain = tmp[(tmp['ref_score_value']>=0.8) & (tmp['var_score_value']>tmp['ref_score_value'])]
# 	tmp_top_gain.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/' + classif + '_total_hcm_top_gain.csv', index=False)
# 	tmp_top_loss = tmp[(tmp['ref_score_value']>=0.8) & (tmp['var_score_value']<tmp['ref_score_value'])]
# 	tmp.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/' + classif + '_total_hcm_top_loss.csv', index=False)
# 	tmp_bottom_gain = tmp[(tmp['ref_score_value']<=0.2) & (tmp['var_score_value']>tmp['ref_score_value'])]
# 	tmp.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/' + classif + '_total_hcm_bottom_gain.csv', index=False)

########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################

#1000g set of snps, technically not huvec (only ~2000 of the 1000G snps)

import glob
import os 
import numpy as np 
import pandas as pd 

####################################

ref_and_damage_scores_1000g = []

for_1000g_dfs = glob.glob('/well/PROCARDIS/domwest/deepmind/1000g_prediction_dfs/*')
for df in for_1000g_dfs:
	tmp = pd.read_csv(df)
	tmp.drop(['Unnamed: 0'], axis=1, inplace=True)
	tmp.rename({'classifier':'variable', 'pos':'grch38_POS', 'ref':'REF', 'alt':'ALT', 'damage_score':'damage_score_value', 'ref_score':'ref_score_value', 'alt_score':'alt_score_value'}, axis=1, inplace=True)
	tmp['chr'] = tmp['chr'].str.replace('chr','')
	tmp['var_type'] = '1000g'
	ref_and_damage_scores_1000g.append(tmp)
ref_and_damage_scores_1000g = pd.concat(ref_and_damage_scores_1000g)

ref_and_damage_scores = ref_and_damage_scores_1000g

ref_and_damage_scores['ref_score_ann'] = ''
ref_and_damage_scores['ref_score_ann'][ref_and_damage_scores['ref_score_value']>=0.8] = 'likely_open'
ref_and_damage_scores['ref_score_ann'][ref_and_damage_scores['ref_score_value']<=0.2] = 'likely_closed'
ref_and_damage_scores['ref_score_ann'][ref_and_damage_scores['ref_score_value']==''] = 'unsure'

####################################

#####***MERGE RSIDS ONTO hcm_ref_and_damage_scores 
#READ IN INDEX FILES
for_1000g_index_files = pd.read_csv('1000g_index_file.txt', sep='\t', header=None)

for_1000g_index_files['snp_test_id'] = for_1000g_index_files[2].str.replace('chr','') + '_' + for_1000g_index_files[3].astype(str) + '_' + for_1000g_index_files[4] + '_' + for_1000g_index_files[5]
for_1000g_index_files.rename({1:'rsid'}, axis=1, inplace=True)
ref_and_damage_scores = ref_and_damage_scores.merge(for_1000g_index_files[['rsid', 'snp_test_id']], how='left', on='snp_test_id')

ref_and_damage_scores.rename({'variable':'classifier'}, axis=1, inplace=True)

ref_and_damage_scores.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/total_1000g.csv', index=False)

# for clas in hcm_ref_and_damage_scores['classifier'].unique():
# 	classif = clas.replace(':','_')
# 	classif = classif.replace(' ','_')
# 	tmp = hcm_ref_and_damage_scores[hcm_ref_and_damage_scores['classifier']==clas]
# 	tmp['var_score_value'] = tmp['damage_score_value'] + tmp['ref_score_value']
# 	# tmp.sort_values(by=['cov'], ascending=False)
# 	tmp.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/' + classif + '_total_hcm.csv', index=False)
# 	tmp_top_gain = tmp[(tmp['ref_score_value']>=0.8) & (tmp['var_score_value']>tmp['ref_score_value'])]
# 	tmp_top_gain.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/' + classif + '_total_hcm_top_gain.csv', index=False)
# 	tmp_top_loss = tmp[(tmp['ref_score_value']>=0.8) & (tmp['var_score_value']<tmp['ref_score_value'])]
# 	tmp.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/' + classif + '_total_hcm_top_loss.csv', index=False)
# 	tmp_bottom_gain = tmp[(tmp['ref_score_value']<=0.2) & (tmp['var_score_value']>tmp['ref_score_value'])]
# 	tmp.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/' + classif + '_total_hcm_bottom_gain.csv', index=False)

########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################

#1000g set of snps, technically not huvec (10000 of the 1000G snps and not eaf_range_filt)

import glob
import os 
import numpy as np 
import pandas as pd 

####################################

ref_and_damage_scores_1000g = []

for_1000g_dfs = glob.glob('/well/PROCARDIS/domwest/deepmind/1000g_10000_snps_prediction_dfs/*')
for df in for_1000g_dfs:
	tmp = pd.read_csv(df)
	tmp.drop(['Unnamed: 0'], axis=1, inplace=True)
	tmp.rename({'classifier':'variable', 'pos':'grch38_POS', 'ref':'REF', 'alt':'ALT', 'damage_score':'damage_score_value', 'ref_score':'ref_score_value', 'alt_score':'alt_score_value'}, axis=1, inplace=True)
	tmp['chr'] = tmp['chr'].str.replace('chr','')
	tmp['var_type'] = '1000g'
	ref_and_damage_scores_1000g.append(tmp)
ref_and_damage_scores_1000g = pd.concat(ref_and_damage_scores_1000g)

ref_and_damage_scores = ref_and_damage_scores_1000g

ref_and_damage_scores['ref_score_ann'] = ''
ref_and_damage_scores['ref_score_ann'][ref_and_damage_scores['ref_score_value']>=0.8] = 'likely_open'
ref_and_damage_scores['ref_score_ann'][ref_and_damage_scores['ref_score_value']<=0.2] = 'likely_closed'
ref_and_damage_scores['ref_score_ann'][ref_and_damage_scores['ref_score_value']==''] = 'unsure'

####################################

#####***MERGE RSIDS ONTO hcm_ref_and_damage_scores 
#READ IN INDEX FILES
for_1000g_index_files = pd.read_csv('1000g_index_file_10000_snps.txt', sep='\t', header=None)

for_1000g_index_files['snp_test_id'] = for_1000g_index_files[2].str.replace('chr','') + '_' + for_1000g_index_files[3].astype(str) + '_' + for_1000g_index_files[4] + '_' + for_1000g_index_files[5]
for_1000g_index_files.rename({1:'rsid'}, axis=1, inplace=True)
ref_and_damage_scores = ref_and_damage_scores.merge(for_1000g_index_files[['rsid', 'snp_test_id']], how='left', on='snp_test_id')

ref_and_damage_scores.rename({'variable':'classifier'}, axis=1, inplace=True)

ref_and_damage_scores.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/total_1000g_10000_snps.csv', index=False)

# for clas in hcm_ref_and_damage_scores['classifier'].unique():
# 	classif = clas.replace(':','_')
# 	classif = classif.replace(' ','_')
# 	tmp = hcm_ref_and_damage_scores[hcm_ref_and_damage_scores['classifier']==clas]
# 	tmp['var_score_value'] = tmp['damage_score_value'] + tmp['ref_score_value']
# 	# tmp.sort_values(by=['cov'], ascending=False)
# 	tmp.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/' + classif + '_total_hcm.csv', index=False)
# 	tmp_top_gain = tmp[(tmp['ref_score_value']>=0.8) & (tmp['var_score_value']>tmp['ref_score_value'])]
# 	tmp_top_gain.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/' + classif + '_total_hcm_top_gain.csv', index=False)
# 	tmp_top_loss = tmp[(tmp['ref_score_value']>=0.8) & (tmp['var_score_value']<tmp['ref_score_value'])]
# 	tmp.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/' + classif + '_total_hcm_top_loss.csv', index=False)
# 	tmp_bottom_gain = tmp[(tmp['ref_score_value']<=0.2) & (tmp['var_score_value']>tmp['ref_score_value'])]
# 	tmp.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/' + classif + '_total_hcm_bottom_gain.csv', index=False)

########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
