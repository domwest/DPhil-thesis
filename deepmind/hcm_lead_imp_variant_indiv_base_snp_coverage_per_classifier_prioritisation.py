import pandas as pd 
import os
import numpy as np 
import glob
import matplotlib.pyplot as plt
import pyBigWig

#TELLS US WHICH SNPS ARE PRESENT IN PEAKS PER CLASSIFIER!

# #first read in all lead and proxy variants regardless of any prioritisation filtering based on damage score:

# df = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deephaem_prep/post_analysis/various_strength_post_analysis_lead_and_imp_vars_damage_scores_with_alleles.csv')

# #get it into format for bedtools

# df['chr'] = 'chr' + df['chr'].astype(str)
# df['grch38_POS2'] = df['grch38_POS'] + 1

# df = df[['chr', 'grch38_POS', 'grch38_POS2', 'rsid', 'strength']]

# df.to_csv('/gpfs3/well/PROCARDIS/domwest/deephaem_prep/post_analysis/multiome/various_strength_post_analysis_lead_and_imp_vars_for_bedtools.bed', index=False, sep='\t', header=None)

#HCM!!!

#######

snps = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deephaem_prep/post_analysis/multiome/various_strength_post_analysis_lead_and_imp_vars_for_bedtools.bed', sep='\t', header=None)
# snps[0] = snps[0].str.replace('chr','')
# snps[0] = snps[0].astype(str)

#beds and bigwigs per classifier:

# bed_files = glob.glob('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/all_chromatin_data_beds_bigwigs/beds/*') #57
bigwig_files = glob.glob('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/all_chromatin_data_beds_bigwigs/bigwigs/*') #57

snp_cov_df = pd.DataFrame()
snp_rsid = []
snp_chr = []
snp_pos = []
snp_cov = []
classifiers = []
for file in bigwig_files:
	classifier = file.split('/')[-1].split('.bw')[0]
	for index,row in snps.iterrows():
		classifiers.append(classifier)
		snp_chr.append(row[0])
		snp_rsid.append(row[3])
		snp_pos.append(row[1])
		bw = pyBigWig.open(file)
		cov = bw.values(row[0], row[1], row[2])
		snp_cov.append(cov[0])
	# 	break
	# break

snp_cov_df['classifier'] = classifiers
snp_cov_df['rsid'] = snp_rsid
snp_cov_df['chr'] = snp_chr
snp_cov_df['pos'] = snp_pos
snp_cov_df['cov'] = snp_cov

snp_cov_df.to_csv('/well/PROCARDIS/domwest/deepmind/cov/HCM_vars_snp_indiv_base_cov_df_per_classifier.csv', index=False) #need to rename this to hcm_vars_blah

###############################

# snp_cov_df_check = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/HCM_vars_snp_indiv_base_cov_df_per_classifier.csv') #need to rename this to hcm_vars_blah

