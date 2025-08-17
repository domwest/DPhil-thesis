# https://ftp.ensembl.org/pub/release-105/variation/vcf/homo_sapiens/
# wget https://ftp.ensembl.org/pub/release-105/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz

#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO

#index file format:
# 0       rs2501376       chr1    24009551        T       A

import pandas as pd

dff = pd.read_excel('/well/PROCARDIS/domwest/deephaem_prep/gathering_snps/HCM_MTAG_OT_FUMA.xlsx')
df1 = pd.crosstab(pd.cut(dff['eaf'], bins=10), dff['eaf'])
df2 = pd.DataFrame(df1.sum(axis=1)).reset_index()
df2['eaf'] = df2['eaf'].astype(str)
df2[['start', 'end']] = df2['eaf'].str.split(',', expand=True)
df2['start'] = df2['start'].str.replace('(','')
df2['end'] = df2['end'].str.replace(']','')
df2['%'] =(df2[0] / 68)*100

df = pd.read_csv('1000GENOMES-phase_3.vcf', sep='\t', comment='#', header=None)
# df = df[df[0].isin([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22])]
# >>> dff['eaf'].min()
# 0.129
# >>> dff['eaf'].max()
# 0.987
df['att1'] = df[7].str.split('EUR=', expand=True)[1]
df['eaf'] = df['att1'].str.split(';', expand=True)[0]
df['eaf'] = df['eaf'].str.split(',', expand=True)[0]
df['eaf'] = df['eaf'].astype(float)
df = df[(df['eaf']>=0.129)&(df['eaf']<=0.987)] #went down to 3 million rows

#split into bins
f = pd.crosstab(pd.cut(df['eaf'], bins=10), df['eaf'])
f2 = pd.DataFrame(f.sum(axis=1)).reset_index()
f2['eaf'] = f2['eaf'].astype(str)
f2[['start', 'end']] = f2['eaf'].str.split(',', expand=True)
f2['start'] = f2['start'].str.replace('(','')
f2['end'] = f2['end'].str.replace(']','')
f2['%'] =(f2[0] / 3630496)*100

# #to get sample= nr:
# for index, row in df2.iterrows():
# 	val = (row['%']/100)*3630496
# 	print(val)

# bin1 = df[(df['eaf']>=0.128)&(df['eaf']<=0.215)]
# bin1 = bin1.sample(n=213559)

# bin2 = df[(df['eaf']>=0.215)&(df['eaf']<=0.301)]
# bin2 = bin2.sample(n=266948) #7.35 / 100 then * total nr rows ie 3630496

# bin3 = df[(df['eaf']>=0.301)&(df['eaf']<=0.386)]
# bin3 = bin3.sample(n=106779) #half
# bin3_2 = bin3.iloc[:53390,:]

# bin4 = df[(df['eaf']>=0.386)&(df['eaf']<=0.472)]
# bin4 = bin4.sample(n=106779) #half
# bin4_2 = bin4.iloc[:53390,:]

# bin5 = df[(df['eaf']>=0.472)&(df['eaf']<=0.558)]
# bin5 = bin5.sample(n=314009) ##had to size down here

# bin6 = df[(df['eaf']>=0.558)&(df['eaf']<=0.644)]
# bin6 = bin6.sample(n=246901) ##had to size down here

# bin7 = df[(df['eaf']>=0.644)&(df['eaf']<=0.73)]
# bin7 = bin7.sample(n=205954) ##had to size down here

# bin8 = df[(df['eaf']>=0.73)&(df['eaf']<=0.815)]
# bin8 = bin8.sample(n=161767) ##had to size down here

# bin9 = df[(df['eaf']>=0.815)&(df['eaf']<=0.901)]
# bin9 = bin9.sample(n=133128) ##had to size down here 

# bin10 = df[(df['eaf']>=0.901)&(df['eaf']<=0.987)]
# bin10 = bin10.sample(n=114201) ##had to size down here #half
# bin10_2 = bin10.iloc[:57101,:]

# bins = pd.concat([bin1, bin2, bin3_2, bin4_2, bin5, bin6, bin7, bin8, bin9, bin10_2])

# #split into bins
# chec = pd.crosstab(pd.cut(bins['eaf'], bins=10), bins['eaf'])
# check = pd.DataFrame(chec.sum(axis=1)).reset_index()
# check['eaf'] = check['eaf'].astype(str)
# check[['start', 'end']] = check['eaf'].str.split(',', expand=True)
# check['start'] = check['start'].str.replace('(','')
# check['end'] = check['end'].str.replace(']','')
# check['%'] =(check[0] / 1706147)*100

# binn = bins.sample(n=2500)

# chec = pd.crosstab(pd.cut(binn['eaf'], bins=10), binn['eaf'])
# check = pd.DataFrame(chec.sum(axis=1)).reset_index()
# check['eaf'] = check['eaf'].astype(str)
# check[['start', 'end']] = check['eaf'].str.split(',', expand=True)
# check['start'] = check['start'].str.replace('(','')
# check['end'] = check['end'].str.replace(']','')
# check['%'] =(check[0] / 2500)*100


# df_final = binn[[2, 0, 1, 3, 4]]
# df_final[4] = df_final[4].str.split(',', expand=True)[0]
# df_final[0] = 'chr' + df_final[0].astype(str)
# df_final.to_csv('1000g_index_file.txt', header=None, sep='\t')

df_subset = df.sample(n=10000)

f = pd.crosstab(pd.cut(df_subset['eaf'], bins=10), df_subset['eaf'])
f2 = pd.DataFrame(f.sum(axis=1)).reset_index()
f2['eaf'] = f2['eaf'].astype(str)
f2[['start', 'end']] = f2['eaf'].str.split(',', expand=True)
f2['start'] = f2['start'].str.replace('(','')
f2['end'] = f2['end'].str.replace(']','')
f2['%'] =(f2[0] / 10000)*100

df_subset = df_subset[[2, 0, 1, 3, 4]].reset_index()
df_subset.drop(['index'], axis=1, inplace=True)
df_subset[0] = 'chr' + df_subset[0].astype(str)
df_subset[4] = df_subset[4].str.split(',', expand=True)[0]
df_subset.to_csv('1000g_index_file_eaf_range_filt_10000_snps.txt', header=None, sep='\t')

dff = pd.read_excel('/well/PROCARDIS/domwest/deephaem_prep/gathering_snps/HCM_MTAG_OT_FUMA.xlsx')
df1 = pd.crosstab(pd.cut(dff['eaf'], bins=10), dff['eaf'])
df2 = pd.DataFrame(df1.sum(axis=1)).reset_index()
df2['eaf'] = df2['eaf'].astype(str)
df2[['start', 'end']] = df2['eaf'].str.split(',', expand=True)
df2['start'] = df2['start'].str.replace('(','')
df2['end'] = df2['end'].str.replace(']','')
df2['%'] =(df2[0] / 68)*100

df = pd.read_csv('1000GENOMES-phase_3.vcf', sep='\t', comment='#', header=None)
# df = df[df[0].isin([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22])]
# >>> dff['eaf'].min()
# 0.129
# >>> dff['eaf'].max()
# 0.987
# df['att1'] = df[7].str.split('EUR=', expand=True)[1]
# df['eaf'] = df['att1'].str.split(';', expand=True)[0]
# df['eaf'] = df['eaf'].str.split(',', expand=True)[0]
# df['eaf'] = df['eaf'].astype(float)

df_subset = df.sample(n=10000)

df_subset = df_subset[[2, 0, 1, 3, 4]].reset_index()
df_subset.drop(['index'], axis=1, inplace=True)
df_subset[0] = 'chr' + df_subset[0].astype(str)
df_subset[4] = df_subset[4].str.split(',', expand=True)[0]
df_subset.to_csv('1000g_index_file_10000_snps.txt', header=None, sep='\t')

######################################################################################################################################################################

import pandas as pd
import os

# DNASE_cardiac_muscle_cell_1_1000g.csv
# DNASE_cardiac_muscle_cell_2_1000g.csv
# DNASE_heart_left_ventricle_female_adult_(53_years)_1000g.csv
# DNASE_heart_left_ventricle_female_embryo_(101_day)_and_female_embryo_(103_days)_1000g.csv
# DNASE_heart_left_ventricle_female_embryo_(136_days)_1000g.csv

df = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_1_1000g.csv')
df

results = df.damage_score_value.sort_values()[::-1]
results

for i, v in enumerate(results[0:150]):
        print(v, (i + 1 + 1) / (len(results) + 1))

#find damage score corresponding to p-value of 0.05

df = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_1_total_hcm.csv')
df[df.damage_score_value > 0.0755509436130523].var_type.value_counts()
# var_type
# bi_allelic       46
# multi_allelic    18

huvec_greater_than_gain_threshold = df[df.damage_score_value > 0.0755509436130523]
huvec_greater_than_gain_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_1_1000g_greater_than_gain_threshold.csv', index=False)

df = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_1_1000g.csv')
df

results = list(df.damage_score_value.sort_values())
results

for i, v in enumerate(results[0:150]):
        print(v, (i + 1 + 1) / (len(results) + 1))

df = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_1_total_hcm.csv')
df[df.damage_score_value < -0.0367254391312599].var_type.value_counts()
# var_type
# bi_allelic       69
# multi_allelic    46

huvec_less_than_loss_threshold = df[df.damage_score_value < -0.0367254391312599]
huvec_less_than_loss_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_1_1000g_less_than_loss_threshold.csv', index=False)

########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################


import pandas as pd
import os

df = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_2_1000g.csv')
df

results = df.damage_score_value.sort_values()[::-1]
results

for i, v in enumerate(results[0:150]):
        print(v, (i + 1 + 1) / (len(results) + 1))

#find damage score corresponding to p-value of 0.05

df = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_2_total_hcm.csv')
df[df.damage_score_value > 0.0328718870878219].var_type.value_counts()
# var_type
# bi_allelic       64
# multi_allelic    31

huvec_greater_than_gain_threshold = df[df.damage_score_value > 0.0328718870878219]
huvec_greater_than_gain_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_2_1000g_greater_than_gain_threshold.csv', index=False)

df = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_2_1000g.csv')
df

results = list(df.damage_score_value.sort_values())
results

for i, v in enumerate(results[0:150]):
        print(v, (i + 1 + 1) / (len(results) + 1))

df = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_2_total_hcm.csv')
df[df.damage_score_value < -0.0260984227061271].var_type.value_counts()
# var_type
# bi_allelic       72
# multi_allelic    34

huvec_less_than_loss_threshold = df[df.damage_score_value < -0.0260984227061271]
huvec_less_than_loss_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_2_1000g_less_than_loss_threshold.csv', index=False)