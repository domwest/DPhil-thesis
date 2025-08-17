#taking most extreme thresholds out of the 2 endothelial classifiers ie:
#0.294990062713623
#-0.2931210994720459

import pandas as pd
import numpy as np

df1 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_1_total_hcm.csv')

# >>> df1[df1.damage_score_value > 0.294990062713623].var_type.value_counts()
# var_type
# bi_allelic       8
# multi_allelic    4
# Name: count, dtype: int64
# >>>
# >>> df1[df1.damage_score_value < -0.2931210994720459].var_type.value_counts()
# var_type
# bi_allelic       12
# multi_allelic    11
# Name: count, dtype: int64

df2 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_2_total_hcm.csv')

# >>> df2[df2.damage_score_value > 0.294990062713623].var_type.value_counts()
# var_type
# bi_allelic       8
# multi_allelic    1
# Name: count, dtype: int64
# >>>
# >>> df2[df2.damage_score_value < -0.2931210994720459].var_type.value_counts()
# var_type
# bi_allelic       6
# multi_allelic    2
# Name: count, dtype: int64

df3 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_heart_left_ventricle_female_adult_(53_years)_total_hcm.csv')

df4 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_heart_left_ventricle_female_embryo_(101_day)_and_female_embryo_(103_days)_total_hcm.csv')

df5 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_heart_left_ventricle_female_embryo_(136_days)_total_hcm.csv')

hcm_less_than_loss_threshold = df1[df1.damage_score_value < -0.2931210994720459]
hcm_less_than_loss_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_1_hcm_less_than_loss_threshold_with_lead_and_proxy_info.csv', index=False)

hcm_greater_than_gain_threshold = df1[df1.damage_score_value > 0.294990062713623]
hcm_greater_than_gain_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_1_hcm_greater_than_gain_threshold_with_lead_and_proxy_info.csv', index=False)

hcm_less_than_loss_threshold = df2[df2.damage_score_value < -0.2931210994720459]
hcm_less_than_loss_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_2_hcm_less_than_loss_threshold_with_lead_and_proxy_info.csv', index=False)

hcm_greater_than_gain_threshold = df2[df2.damage_score_value > 0.294990062713623]
hcm_greater_than_gain_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_2_hcm_greater_than_gain_threshold_with_lead_and_proxy_info.csv', index=False)

hcm_less_than_loss_threshold = df3[df3.damage_score_value < -0.2931210994720459]
hcm_less_than_loss_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_heart_left_ventricle_female_adult_(53_years)_hcm_less_than_loss_threshold_with_lead_and_proxy_info.csv', index=False)

hcm_greater_than_gain_threshold = df3[df3.damage_score_value > 0.294990062713623]
hcm_greater_than_gain_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_heart_left_ventricle_female_adult_(53_years)_hcm_greater_than_gain_threshold_with_lead_and_proxy_info.csv', index=False)

hcm_less_than_loss_threshold = df4[df4.damage_score_value < -0.2931210994720459]
hcm_less_than_loss_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_heart_left_ventricle_female_embryo_(101_day)_and_female_embryo_(103_days)_hcm_less_than_loss_threshold_with_lead_and_proxy_info.csv', index=False)

hcm_greater_than_gain_threshold = df4[df4.damage_score_value > 0.294990062713623]
hcm_greater_than_gain_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_heart_left_ventricle_female_embryo_(101_day)_and_female_embryo_(103_days)_hcm_greater_than_gain_threshold_with_lead_and_proxy_info.csv', index=False)

hcm_less_than_loss_threshold = df5[df5.damage_score_value < -0.2931210994720459]
hcm_less_than_loss_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_heart_left_ventricle_female_embryo_(136_days)_hcm_less_than_loss_threshold_with_lead_and_proxy_info.csv', index=False)

hcm_greater_than_gain_threshold = df5[df5.damage_score_value > 0.294990062713623]
hcm_greater_than_gain_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_heart_left_ventricle_female_embryo_(136_days)_hcm_greater_than_gain_threshold_with_lead_and_proxy_info.csv', index=False)


###making tables to report:

class_tot = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/all_heart_classifiers_hcm_prioritised_with_lead_and_proxy_info.csv')

class_tot1 = class_tot[class_tot['classifier']=='DNASE:cardiac_muscle_cell_1']
class_tot1 = class_tot1[['classifier', 'rsid', 'snp_test_id', 'var_type', 'ref_score_ann', 'ref_score_value', 'alt_score_value', 'damage_score_value', 'cov', 'is_in_peak']]
class_tot1 = class_tot1.drop_duplicates()
class_tot1.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_1_hcm_prioritised_with_lead_and_proxy_info_table_to_report.csv', index=False)

class_tot2 = class_tot[class_tot['classifier']=='DNASE:cardiac_muscle_cell_2']
class_tot2 = class_tot2[['classifier', 'rsid', 'snp_test_id', 'var_type', 'ref_score_ann', 'ref_score_value', 'alt_score_value', 'damage_score_value', 'cov', 'is_in_peak']]
class_tot2 = class_tot2.drop_duplicates()
class_tot2.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_2_hcm_prioritised_with_lead_and_proxy_info_table_to_report.csv', index=False)

class_tot3 = class_tot[class_tot['classifier']=='DNASE:heart_left_ventricle_female_adult_(53_years)']
class_tot3 = class_tot3[['classifier', 'rsid', 'snp_test_id', 'var_type', 'ref_score_ann', 'ref_score_value', 'alt_score_value', 'damage_score_value', 'cov', 'is_in_peak']]
class_tot3 = class_tot3.drop_duplicates()
class_tot3.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_heart_left_ventricle_female_adult_(53_years)_hcm_prioritised_with_lead_and_proxy_info_table_to_report.csv', index=False)

class_tot4 = class_tot[class_tot['classifier']=='DNASE:heart_left_ventricle_female_embryo_(101_day)_and_female_embryo_(103_days)']
class_tot4 = class_tot4[['classifier', 'rsid', 'snp_test_id', 'var_type', 'ref_score_ann', 'ref_score_value', 'alt_score_value', 'damage_score_value', 'cov', 'is_in_peak']]
class_tot4 = class_tot4.drop_duplicates()
class_tot4.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_heart_left_ventricle_female_embryo_(101_day)_and_female_embryo_(103_days)_hcm_prioritised_with_lead_and_proxy_info_table_to_report.csv', index=False)

class_tot5 = class_tot[class_tot['classifier']=='DNASE:heart_left_ventricle_female_embryo_(136_days)']
class_tot5 = class_tot5[['classifier', 'rsid', 'snp_test_id', 'var_type', 'ref_score_ann', 'ref_score_value', 'alt_score_value', 'damage_score_value', 'cov', 'is_in_peak']]
class_tot5 = class_tot5.drop_duplicates()
class_tot5.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_heart_left_ventricle_female_embryo_(136_days)_hcm_prioritised_with_lead_and_proxy_info_table_to_report.csv', index=False)





























######################################################################################################

#FOR EXTRA TRDN LOCUS ALONE:

#taking most extreme thresholds out of the 2 endothelial classifiers ie:
#0.294990062713623
#-0.2931210994720459

import pandas as pd
import numpy as np

df1 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/matrices/DNASE_cardiac_muscle_cell_1_total_hcm.csv')

# >>> df1[df1.damage_score_value > 0.294990062713623].var_type.value_counts()
# var_type
# bi_allelic       8
# multi_allelic    4
# Name: count, dtype: int64
# >>>
# >>> df1[df1.damage_score_value < -0.2931210994720459].var_type.value_counts()
# var_type
# bi_allelic       12
# multi_allelic    11
# Name: count, dtype: int64

df2 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/matrices/DNASE_cardiac_muscle_cell_2_total_hcm.csv')

# >>> df2[df2.damage_score_value > 0.294990062713623].var_type.value_counts()
# var_type
# bi_allelic       8
# multi_allelic    1
# Name: count, dtype: int64
# >>>
# >>> df2[df2.damage_score_value < -0.2931210994720459].var_type.value_counts()
# var_type
# bi_allelic       6
# multi_allelic    2
# Name: count, dtype: int64

df3 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/matrices/DNASE_heart_left_ventricle_female_adult_(53_years)_total_hcm.csv')

df4 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/matrices/DNASE_heart_left_ventricle_female_embryo_(101_day)_and_female_embryo_(103_days)_total_hcm.csv')

df5 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/matrices/DNASE_heart_left_ventricle_female_embryo_(136_days)_total_hcm.csv')

hcm_less_than_loss_threshold = df1[df1.damage_score_value < -0.2931210994720459]
#hcm_less_than_loss_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/matrices/DNASE_cardiac_muscle_cell_1_hcm_less_than_loss_threshold_with_lead_and_proxy_info.csv', index=False)

hcm_greater_than_gain_threshold = df1[df1.damage_score_value > 0.294990062713623]
#hcm_greater_than_gain_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/matrices/DNASE_cardiac_muscle_cell_1_hcm_greater_than_gain_threshold_with_lead_and_proxy_info.csv', index=False)

hcm_less_than_loss_threshold = df2[df2.damage_score_value < -0.2931210994720459]
#hcm_less_than_loss_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/matrices/DNASE_cardiac_muscle_cell_2_hcm_less_than_loss_threshold_with_lead_and_proxy_info.csv', index=False)

hcm_greater_than_gain_threshold = df2[df2.damage_score_value > 0.294990062713623]
#hcm_greater_than_gain_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/matrices/DNASE_cardiac_muscle_cell_2_hcm_greater_than_gain_threshold_with_lead_and_proxy_info.csv', index=False)

hcm_less_than_loss_threshold = df3[df3.damage_score_value < -0.2931210994720459]
#hcm_less_than_loss_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/matrices/DNASE_heart_left_ventricle_female_adult_(53_years)_hcm_less_than_loss_threshold_with_lead_and_proxy_info.csv', index=False)

hcm_greater_than_gain_threshold = df3[df3.damage_score_value > 0.294990062713623]
#hcm_greater_than_gain_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/matrices/DNASE_heart_left_ventricle_female_adult_(53_years)_hcm_greater_than_gain_threshold_with_lead_and_proxy_info.csv', index=False)

hcm_less_than_loss_threshold = df4[df4.damage_score_value < -0.2931210994720459]
#hcm_less_than_loss_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/matrices/DNASE_heart_left_ventricle_female_embryo_(101_day)_and_female_embryo_(103_days)_hcm_less_than_loss_threshold_with_lead_and_proxy_info.csv', index=False)

hcm_greater_than_gain_threshold = df4[df4.damage_score_value > 0.294990062713623]
#hcm_greater_than_gain_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/matrices/DNASE_heart_left_ventricle_female_embryo_(101_day)_and_female_embryo_(103_days)_hcm_greater_than_gain_threshold_with_lead_and_proxy_info.csv', index=False)

hcm_less_than_loss_threshold = df5[df5.damage_score_value < -0.2931210994720459]
#hcm_less_than_loss_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/matrices/DNASE_heart_left_ventricle_female_embryo_(136_days)_hcm_less_than_loss_threshold_with_lead_and_proxy_info.csv', index=False)

hcm_greater_than_gain_threshold = df5[df5.damage_score_value > 0.294990062713623]
#hcm_greater_than_gain_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/matrices/DNASE_heart_left_ventricle_female_embryo_(136_days)_hcm_greater_than_gain_threshold_with_lead_and_proxy_info.csv', index=False)

##all empty so not saving them

















######################################################################################################

#FOR EXTRA HRC 1 LOCUS ALONE:

#taking most extreme thresholds out of the 2 endothelial classifiers ie:
#0.294990062713623
#-0.2931210994720459

import pandas as pd
import numpy as np

df1 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/hrc_matrices_1/DNASE_cardiac_muscle_cell_1_total_hcm.csv')

# >>> df1[df1.damage_score_value > 0.294990062713623].var_type.value_counts()
# var_type
# bi_allelic       8
# multi_allelic    4
# Name: count, dtype: int64
# >>>
# >>> df1[df1.damage_score_value < -0.2931210994720459].var_type.value_counts()
# var_type
# bi_allelic       12
# multi_allelic    11
# Name: count, dtype: int64

df2 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/hrc_matrices_1/DNASE_cardiac_muscle_cell_2_total_hcm.csv')

# >>> df2[df2.damage_score_value > 0.294990062713623].var_type.value_counts()
# var_type
# bi_allelic       8
# multi_allelic    1
# Name: count, dtype: int64
# >>>
# >>> df2[df2.damage_score_value < -0.2931210994720459].var_type.value_counts()
# var_type
# bi_allelic       6
# multi_allelic    2
# Name: count, dtype: int64

df3 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/hrc_matrices_1/DNASE_heart_left_ventricle_female_adult_(53_years)_total_hcm.csv')

df4 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/hrc_matrices_1/DNASE_heart_left_ventricle_female_embryo_(101_day)_and_female_embryo_(103_days)_total_hcm.csv')

df5 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/hrc_matrices_1/DNASE_heart_left_ventricle_female_embryo_(136_days)_total_hcm.csv')

hcm_less_than_loss_threshold = df1[df1.damage_score_value < -0.2931210994720459]
#hcm_less_than_loss_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/matrices/DNASE_cardiac_muscle_cell_1_hcm_less_than_loss_threshold_with_lead_and_proxy_info.csv', index=False)

hcm_greater_than_gain_threshold = df1[df1.damage_score_value > 0.294990062713623]
#        snp_test_id  damage_score_value
# 3  19_49157632_A_C            0.297294
hcm_greater_than_gain_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/hrc_matrices_1/DNASE_cardiac_muscle_cell_1_hcm_greater_than_gain_threshold_with_lead_and_proxy_info.csv', index=False)

hcm_less_than_loss_threshold = df2[df2.damage_score_value < -0.2931210994720459]
#hcm_less_than_loss_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/matrices/DNASE_cardiac_muscle_cell_2_hcm_less_than_loss_threshold_with_lead_and_proxy_info.csv', index=False)

hcm_greater_than_gain_threshold = df2[df2.damage_score_value > 0.294990062713623]
#        snp_test_id  damage_score_value
# 4  19_49157632_A_C            0.298685
hcm_greater_than_gain_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/hrc_matrices_1/DNASE_cardiac_muscle_cell_2_hcm_greater_than_gain_threshold_with_lead_and_proxy_info.csv', index=False)

hcm_less_than_loss_threshold = df3[df3.damage_score_value < -0.2931210994720459]
#hcm_less_than_loss_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/matrices/DNASE_heart_left_ventricle_female_adult_(53_years)_hcm_less_than_loss_threshold_with_lead_and_proxy_info.csv', index=False)

hcm_greater_than_gain_threshold = df3[df3.damage_score_value > 0.294990062713623]
#hcm_greater_than_gain_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/matrices/DNASE_heart_left_ventricle_female_adult_(53_years)_hcm_greater_than_gain_threshold_with_lead_and_proxy_info.csv', index=False)

hcm_less_than_loss_threshold = df4[df4.damage_score_value < -0.2931210994720459]
#hcm_less_than_loss_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/matrices/DNASE_heart_left_ventricle_female_embryo_(101_day)_and_female_embryo_(103_days)_hcm_less_than_loss_threshold_with_lead_and_proxy_info.csv', index=False)

hcm_greater_than_gain_threshold = df4[df4.damage_score_value > 0.294990062713623]
#hcm_greater_than_gain_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/matrices/DNASE_heart_left_ventricle_female_embryo_(101_day)_and_female_embryo_(103_days)_hcm_greater_than_gain_threshold_with_lead_and_proxy_info.csv', index=False)

hcm_less_than_loss_threshold = df5[df5.damage_score_value < -0.2931210994720459]
#hcm_less_than_loss_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/matrices/DNASE_heart_left_ventricle_female_embryo_(136_days)_hcm_less_than_loss_threshold_with_lead_and_proxy_info.csv', index=False)

hcm_greater_than_gain_threshold = df5[df5.damage_score_value > 0.294990062713623]
#hcm_greater_than_gain_threshold.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/matrices/DNASE_heart_left_ventricle_female_embryo_(136_days)_hcm_greater_than_gain_threshold_with_lead_and_proxy_info.csv', index=False)
