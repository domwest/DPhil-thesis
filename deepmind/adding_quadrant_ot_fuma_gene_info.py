import pandas as pd
import numpy as np

##SNP TABLE PART 1

quadr = pd.read_csv('/well/PROCARDIS/domwest/replicate_old/compare_bulk_and_sc_quadrant_approach/final_bulk_and_sc_prioritised_genes_per_lead_snp_refseq_genes_only.csv')
quadr.rename({'rsid':'proxy_rsid'}, axis=1, inplace=True)
quadr = quadr[['proxy_rsid', 'bulk_gene_symbol', 'sc_gene_symbol', 'same_bulk_and_sc_gene', 'sc_enriched_in_cm']]

#latest scatter plots are here: C:\Users\domwest\Documents\oxford_dphil\compare_new_quadrant_appr_ilo_tads_to_old\new_scatter_plots_rescomp

ot_fuma = pd.read_excel('/well/PROCARDIS/domwest/deephaem_prep/gathering_snps/HCM_MTAG_OT_FUMA.xlsx')
ot_fuma.rename({'rsid':'proxy_rsid', 'snptestid':'snp_test_id', 'likely gene(s) (updated 17/11/22)':'OT/FUMA'}, axis=1, inplace=True)
ot_fuma['snp_test_id'] = ot_fuma['snp_test_id'].str.replace(':','_')
ot_fuma['snp_test_id'] = ot_fuma['snp_test_id'].str.replace('chr','')
ot_fuma = ot_fuma[['proxy_rsid', 'OT/FUMA']]
ot_fuma[['OT', 'FUMA']] = ot_fuma['OT/FUMA'].str.split('/', expand=True)
ot_fuma.drop(['OT/FUMA'], axis=1, inplace=True)
ot_fuma['FUMA'] = ot_fuma['FUMA'].fillna(ot_fuma['OT'])

gene_df =  ot_fuma.merge(quadr, on='proxy_rsid', how='left')
gene_df['lead/proxy'] = 'lead'

deepmind = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/all_heart_classifiers_hcm_prioritised_with_lead_and_proxy_info.csv')
deepmind = deepmind.sort_values("damage_score_value", ascending=False)
group = deepmind.groupby(['rsid', 'snp_test_id'])
df2 = pd.DataFrame(group.apply(lambda x: x['classifier'].unique())).reset_index()
df2.rename({0:'classifiers'}, axis=1, inplace=True)
df2 = deepmind.merge(df2, how='left', on=['rsid', 'snp_test_id'])

df3 = pd.DataFrame(group.apply(lambda x: x['damage_score_value'].unique())).reset_index()
df3.rename({0:'damage_score_values'}, axis=1, inplace=True)
df3 = df2.merge(df3, how='left', on=['rsid', 'snp_test_id'])

df3['classifiers'] = df3['classifiers'].tolist()
df3['classifiers'] = df3['classifiers'].astype(str)

df4 = pd.DataFrame(df3.groupby(['rsid','snp_test_id','classifiers'])['cov'].agg(list)).reset_index()
df4.rename({'cov':'covs'}, axis=1, inplace=True)
df4 = df3.merge(df4, how='left', on=['rsid', 'snp_test_id', 'classifiers'])

df5 = pd.DataFrame(df3.groupby(['rsid','snp_test_id','classifiers'])['is_in_peak'].agg(list)).reset_index()
df5.rename({'is_in_peak':'is_in_peaks'}, axis=1, inplace=True)
df5 = df4.merge(df5, how='left', on=['rsid', 'snp_test_id', 'classifiers'])

dfs = df5
dfs['var_type'] = ''
dfs['var_type'][dfs['damage_score_values'].astype(str).str.contains('-')] = 'loss'
dfs['var_type'][~(dfs['damage_score_values'].astype(str).str.contains('-'))] = 'gain'
dfs.rename({'rsid':'proxy_rsid'}, axis=1, inplace=True)
dfs = dfs[['proxy_rsid', 'snp_test_id', 'classifiers', 'damage_score_values', 'covs', 'is_in_peaks', 'var_type']]

dfs['damage_score_values'] = dfs['damage_score_values'].tolist()
dfs['damage_score_values'] = dfs['damage_score_values'].astype(str)

dfs['covs'] = dfs['covs'].tolist()
dfs['covs'] = dfs['covs'].astype(str)

dfs['is_in_peaks'] = dfs['is_in_peaks'].tolist()
dfs['is_in_peaks'] = dfs['is_in_peaks'].astype(str)

dfs = dfs.drop_duplicates() #53

#merge onto gene_df to get lead/proxy info

snps_og = pd.read_excel('/well/PROCARDIS/domwest/deephaem_prep/gathering_snps/dwest_deephaem_lead_imp_snps.xlsx')
dfs = dfs.merge(snps_og[['proxy_rsid', 'lead_rsid']], how='left', on='proxy_rsid')
gene_df = gene_df.merge(snps_og[['proxy_rsid', 'lead_rsid']], how='left', on='proxy_rsid')
gene_df.drop(['proxy_rsid'], axis=1, inplace=True)

deepmind_gene_df = dfs.merge(gene_df, how='left', on='lead_rsid')
deepmind_gene_df['lead/proxy'] = ''
deepmind_gene_df['lead/proxy'][deepmind_gene_df['lead_rsid']==deepmind_gene_df['proxy_rsid']] = 'lead'
deepmind_gene_df['lead/proxy'][deepmind_gene_df['lead_rsid']!=deepmind_gene_df['proxy_rsid']] = 'proxy'

#mark extra lead snps
extras = ['rs67491807', 'rs764462761', 'rs3218719', 'rs11570041', 'rs117847273']
deepmind_gene_df['lead/proxy'][deepmind_gene_df['proxy_rsid'].isin(extras)] = 'lead'
deepmind_gene_df['lead/proxy'][deepmind_gene_df['lead/proxy'].isnull()] = 'proxy'

#Add REgulamentary info

#reran first part of geneticsIntersection.py from regulamentary scripts, and produced the following csv:
intersect_final = pd.read_csv('/well/PROCARDIS/domwest/deepmind/regulamentary_merged_hcm_lead_and_proxy_snps_on_re_results.csv')
intersect_final['REgulamentary_element_coords'] = intersect_final['chrom'].str.replace('chr','') + '_' + intersect_final['start'].astype(str) + '_' + intersect_final['end'].astype(str) 
intersect_final.rename({'SNPS':'proxy_rsid', 'RE':'REgulamentary_element'}, axis=1, inplace=True)
intersect_final = intersect_final[['proxy_rsid', 'REgulamentary_element', 'REgulamentary_element_coords']]

deepmind_gene_df2 = deepmind_gene_df.merge(intersect_final, how='left', on='proxy_rsid')

#Add ABC model info

intersect_fin = pd.read_csv('/well/PROCARDIS/domwest/ABC-Enhancer-Gene-Prediction/abc_merged_prioritised_hcm_lead_and_proxy_snps_on_abc_results.csv')
intersect_fin['ABC_enhancer_coords'] = intersect_fin['chr'].str.replace('chr','') + '_' + intersect_fin['start'].astype(str) + '_' + intersect_fin['end'].astype(str) 
#cols available:
# Index(['chr', 'start', 'end', 'name', 'class', 'activity_base',
#        'activity_base_enh', 'activity_base_squared_enh', 'normalized_dhs_enh',
#        'normalized_h3k27ac_enh', 'TargetGene', 'TargetGeneTSS',
#        'TargetGeneExpression', 'TargetGenePromoterActivityQuantile',
#        'TargetGeneIsExpressed', 'TargetGeneEnsembl_ID', 'normalized_dhs_prom',
#        'normalized_h3k27ac_prom', 'distance', 'isSelfPromoter',
#        'powerlaw_contact', 'powerlaw_contact_reference', 'hic_contact',
#        'hic_contact_pl_scaled', 'hic_pseudocount', 'hic_contact_pl_scaled_adj',
#        'ABC.Score.Numerator', 'ABC.Score', 'powerlaw.Score.Numerator',
#        'powerlaw.Score', 'CellType', 'hic_contact_squared',
#        'ABC_enhancer_coords']
intersect_fin.rename({'SNPS':'proxy_rsid', 'ABC.Score':'ABC_score', 'TargetGene':'ABC_target_gene', 'class':'ABC_target_class'}, axis=1, inplace=True)
intersect_fin = intersect_fin[['proxy_rsid', 'ABC_enhancer_coords', 'ABC_target_gene', 'ABC_target_class', 'ABC_score']]

deepmind_gene_df2 = deepmind_gene_df2.merge(intersect_fin, how='left', on='proxy_rsid')

deepmind_gene_df2 = deepmind_gene_df2.sort_values("ABC_score", ascending=False)
group = deepmind_gene_df2.groupby(['proxy_rsid', 'snp_test_id'])
df2 = pd.DataFrame(group.apply(lambda x: x['ABC_target_gene'].unique())).reset_index()
df2.rename({0:'ABC_target_genes'}, axis=1, inplace=True)
df2 = deepmind_gene_df2.merge(df2, how='left', on=['proxy_rsid', 'snp_test_id'])

df2['ABC_target_genes'] = df2['ABC_target_genes'].tolist()
df2['ABC_target_genes'] = df2['ABC_target_genes'].astype(str)

# df3 = pd.DataFrame(group.apply(lambda x: x['ABC_score'].unique())).reset_index()
# df3.rename({0:'ABC_scores'}, axis=1, inplace=True)
# df3 = df2.merge(df3, how='left', on=['proxy_rsid', 'snp_test_id'])

df3 = pd.DataFrame(df2.groupby(['proxy_rsid','snp_test_id','ABC_target_genes'])['ABC_score'].agg(list)).reset_index()
df3.rename({'ABC_score':'ABC_scores'}, axis=1, inplace=True)
df3 = df2.merge(df3, how='left', on=['proxy_rsid', 'snp_test_id', 'ABC_target_genes'])

df4 = pd.DataFrame(df3.groupby(['proxy_rsid','snp_test_id','ABC_target_genes'])['ABC_target_class'].agg(list)).reset_index()
df4.rename({'ABC_target_class':'ABC_target_classes'}, axis=1, inplace=True)
df4 = df3.merge(df4, how='left', on=['proxy_rsid', 'snp_test_id', 'ABC_target_genes'])

dfs = df4

dfs.drop(['ABC_target_gene', 'ABC_target_class', 'ABC_score'], axis=1, inplace=True)

dfs['ABC_scores'] = dfs['ABC_scores'].tolist()
dfs['ABC_scores'] = dfs['ABC_scores'].astype(str)

dfs['ABC_target_classes'] = dfs['ABC_target_classes'].tolist()
dfs['ABC_target_classes'] = dfs['ABC_target_classes'].astype(str)

dfs = dfs.drop_duplicates() #53

deepmind_gene_df3 = dfs

#***add median and min and max damage scores found across all the other cell types

#add arhcr info

archr = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/post_analysis/multiome/with_whole_pipeline_all_ventricular_cm_peak_to_gene_id_only_filtered_for_lead_and_proxy_snps_wrangled.csv')
#should we take 'top score' as the top correlation for those wherever FDR is <=0.05...
archr = archr[archr['FDR']<=0.05]

###

archr1 = archr
archr1.rename({'rsid':'proxy_rsid'}, axis=1, inplace=True)
archr1 = archr1.sort_values("Correlation", ascending=False)
archr1 = archr1.drop_duplicates(subset=["proxy_rsid"], keep="first")
archr1.rename({'geneName':'[SNP]_ArchR_top_gene_FDR<=0.05'}, axis=1, inplace=True)
archr1['[SNP]_ArchR_top_gene_FDR<=0.05_scores'] = 'Corr:' + archr1['Correlation'].round(2).astype(str) + ';FDR:' + archr1['FDR'].round(2).astype(str)
deepmind_gene_df3 = deepmind_gene_df3.merge(archr1[['proxy_rsid', '[SNP]_ArchR_top_gene_FDR<=0.05', '[SNP]_ArchR_top_gene_FDR<=0.05_scores']], how='left', on='proxy_rsid')

###

archr2 = archr
archr2 = archr2[archr2['geneName'].isin(deepmind_gene_df3['sc_gene_symbol'].unique().tolist())] #2 being captured out of the 9
archr2 = archr2.sort_values("Correlation", ascending=False)
archr2 = archr2.drop_duplicates(subset=["proxy_rsid"], keep="first")
archr2.rename({'geneName':'sc_gene_symbol'}, axis=1, inplace=True)
archr2['[Q.A_SC_GENE]_ArchR_top_rsid_scores'] = 'Corr:' + archr2['Correlation'].round(5).astype(str) + ';FDR:' + archr2['FDR'].round(5).astype(str)
archr2.drop(['strength'], axis=1, inplace=True)
archr2 = archr2.drop_duplicates()
#manually edit:
archr2['proxy_rsid'].loc[82618] = 'rs75010486&rs16940806&rs17652748'
archr2['[Q.A_SC_GENE]_ArchR_top_rsid_scores'].loc[82618] = 'Corr:0.16254;FDR:0.00125&Corr:0.16254;FDR:0.00125&Corr:0.16254;FDR:0.00125'
archr2.drop([113735, 45279], inplace=True)
archr2['proxy_rsid'].loc[4248] = 'rs10803407&rs9442217'
archr2['[Q.A_SC_GENE]_ArchR_top_rsid_scores'].loc[4248] = 'Corr:0.10988;FDR:0.03894&Corr:0.10988;FDR:0.03894'
archr2.drop([181969], inplace=True)
archr2['proxy_rsid'].loc[67476] = 'rs4633690&rs2879828'
archr2['[Q.A_SC_GENE]_ArchR_top_rsid_scores'].loc[67476] = 'Corr:-0.12239;FDR:0.01944&Corr:-0.12239;FDR:0.01944'
archr2.drop([71600], inplace=True)
archr2['proxy_rsid'].loc[11346] = 'rs11838776&rs9515201'
archr2['[Q.A_SC_GENE]_ArchR_top_rsid_scores'].loc[11346] = 'Corr:-0.25815;FDR:0.0&Corr:-0.25815;FDR:0.0'
archr2.drop([52820], inplace=True)
archr2.rename({'proxy_rsid':'[Q.A_SC_GENE]_ArchR_top_rsid'}, axis=1, inplace=True)

archr2 = archr2.groupby('sc_gene_symbol').head(1)

deepmind_gene_df3 = deepmind_gene_df3.merge(archr2[['sc_gene_symbol', '[Q.A_SC_GENE]_ArchR_top_rsid', '[Q.A_SC_GENE]_ArchR_top_rsid_scores']], how='left', on='sc_gene_symbol')

##add open target gene per prioritised variant

# import pandas as pd
# import numpy as np

# deepmind_gene_df2 = pd.read_csv('/well/PROCARDIS/domwest/deepmind/prioritised_hcm_snps_final_result_table_part1.csv')
# len(deepmind_gene_df2['proxy_rsid'].unique()) #gives 49 instead of 53 and thats because multi-allelic snps will show the same rsid... so will get the same gene per snp in these scenarios

deepmind_ot_output = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/prioritised_snp_exploration/open_target/opentarget_snp_to_gene.csv')

deepmind_ot_output.rename({'SNP':'proxy_rsid', 'Gene':'[SNP]_OT_top_gene', 'Overall V2G':'[SNP]_OT_top_gene_overall_v2g_score'}, axis=1, inplace=True)

deepmind_gene_df3 = deepmind_gene_df3.merge(deepmind_ot_output[['proxy_rsid', '[SNP]_OT_top_gene', '[SNP]_OT_top_gene_overall_v2g_score']], how='left', on='proxy_rsid')

##add homer top motif for each variant

de_novo_motifs = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/prioritised_snp_exploration/homer/indiv_running/de_novo_wrangled_output_for_all_snps.csv')
de_novo_motifs.rename({'snp':'snp_test_id', 'motif_name':'[SNP]_homer_de_novo_top_motif_name', 'log_odds_detection_threshold':'[SNP]_homer_de_novo_top_motif_log_odds_detection_threshold', 'final_enrich_p_value':'[SNP]_homer_de_novo_top_motif_final_enrich_p_value', 'target_seq_occurence_%':'[SNP]_homer_de_novo_top_motif_target_seq_occurence_%', 'bg_seq_occurence_%':'[SNP]_homer_de_novo_top_motif_bg_seq_occurence_%'}, axis=1, inplace=True)

deepmind_gene_df3 = deepmind_gene_df3.merge(de_novo_motifs, how='left', on='snp_test_id')

known_motifs = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/prioritised_snp_exploration/homer/indiv_running/known_wrangled_output_for_all_snps.csv')
known_motifs.rename({'snp':'snp_test_id', 'Motif Name':'[SNP]_homer_known_top_motif_name', 'Consensus':'[SNP]_homer_known_top_motif_consensus', 'P-value':'[SNP]_homer_known_top_motif_final_enrich_p_value', '% of Target Sequences with Motif':'[SNP]_homer_known_top_motif_target_seq_occurence_%', '% of Background Sequences with Motif':'[SNP]_homer_known_top_motif_bg_seq_occurence_%'}, axis=1, inplace=True)

deepmind_gene_df3 = deepmind_gene_df3.merge(known_motifs, how='left', on='snp_test_id')

cols_to_drop = deepmind_gene_df3.columns[deepmind_gene_df3.columns.str.contains('homer')]
deepmind_gene_df3.drop(cols_to_drop, axis=1, inplace=True)

#add mtag p val and or cutoffs

loci = pd.read_excel('/well/PROCARDIS/domwest/deepmind/HCM_MTAG_OT_FUMA.xlsx')
summary_data = pd.read_csv('/well/PROCARDIS/domwest/deepmind/hcm.mtag.2024.format.tsv', sep='\t')

ld_info = pd.read_excel('/well/PROCARDIS/domwest/deepmind/dwest_deephaem_lead_imp_snps_b37.xlsx')

snp_test_id_info = pd.read_csv('/well/PROCARDIS/domwest/deepmind/all_chunks_00-All_merged_on_og_snps_filt.csv') #4395

ld_info_filt = ld_info[['lead_rsid', 'proxy_rsid']].drop_duplicates()
ld_info_filt2 = ld_info_filt[ld_info_filt['proxy_rsid'].isin(snp_test_id_info['rsid'].unique().tolist())] #4395

df1 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_1_total_hcm.csv')
df2 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_2_total_hcm.csv')
df3 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_heart_left_ventricle_female_adult_(53_years)_total_hcm.csv')
df4 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_heart_left_ventricle_female_embryo_(101_day)_and_female_embryo_(103_days)_total_hcm.csv')
df5 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_heart_left_ventricle_female_embryo_(136_days)_total_hcm.csv')
scores = pd.concat([df1, df2, df3, df4, df5])
len(scores['rsid'].unique())
#4395

loci.rename({'rsid':'lead_rsid'}, axis=1, inplace=True)
summary_data.rename({'rs_id':'proxy_rsid'}, axis=1, inplace=True)

lead_proxy_pvals_summary = ld_info_filt2.merge(summary_data[['proxy_rsid', 'variant_id', 'p_value', 'beta']], how='left', on='proxy_rsid')
lead_proxy_pvals_summary_and_loci = lead_proxy_pvals_summary.merge(loci[['lead_rsid', 'mtag_pval', 'mtag_beta']], how='left', on='lead_rsid')
lead_proxy_pvals_summary_and_loci['p_value'] = np.where(lead_proxy_pvals_summary_and_loci['p_value'].isnull(),lead_proxy_pvals_summary_and_loci['mtag_pval'],lead_proxy_pvals_summary_and_loci['p_value'])
lead_proxy_pvals_summary_and_loci['beta'] = np.where(lead_proxy_pvals_summary_and_loci['beta'].isnull(),lead_proxy_pvals_summary_and_loci['mtag_beta'],lead_proxy_pvals_summary_and_loci['beta'])

scores.rename({'rsid':'proxy_rsid', 'classifier':'variable', 'damage_score_value':'value'}, axis=1, inplace=True)
lead_proxy_pvals_summary_and_loci_values = lead_proxy_pvals_summary_and_loci.merge(scores[['proxy_rsid', 'variable', 'value']], how='left', on='proxy_rsid')
lead_proxy_pvals_summary_and_loci_values_filt = lead_proxy_pvals_summary_and_loci_values[['proxy_rsid', 'variant_id', 'variable', 'value', 'p_value', 'beta']].drop_duplicates()
lead_proxy_pvals_summary_and_loci_values_filt['-log_10_p_value'] = -(np.log10(lead_proxy_pvals_summary_and_loci_values_filt['p_value']))

#what is -log10(0.00000005) ... 7.3
priorit_loss_vars = lead_proxy_pvals_summary_and_loci_values_filt[(lead_proxy_pvals_summary_and_loci_values_filt['-log_10_p_value']>=7.3)&(lead_proxy_pvals_summary_and_loci_values_filt['value']<=-0.2931210994720459)]
priorit_gain_vars = lead_proxy_pvals_summary_and_loci_values_filt[(lead_proxy_pvals_summary_and_loci_values_filt['-log_10_p_value']>=7.3)&(lead_proxy_pvals_summary_and_loci_values_filt['value']>=0.294990062713623)]

lead_proxy_pvals_summary_and_loci_values_filt['OR'] = lead_proxy_pvals_summary_and_loci_values_filt.apply(lambda row: math.exp(row.beta), axis = 1)

priorit_loss_vars_or = lead_proxy_pvals_summary_and_loci_values_filt[(lead_proxy_pvals_summary_and_loci_values_filt['OR']>=1)&(lead_proxy_pvals_summary_and_loci_values_filt['value']<=-0.2931210994720459)]
priorit_gain_vars_or = lead_proxy_pvals_summary_and_loci_values_filt[(lead_proxy_pvals_summary_and_loci_values_filt['OR']>=1)&(lead_proxy_pvals_summary_and_loci_values_filt['value']>=0.294990062713623)]

p_val = priorit_loss_vars['proxy_rsid'].unique().tolist() + priorit_gain_vars['proxy_rsid'].unique().tolist()
o_r = priorit_loss_vars_or['proxy_rsid'].unique().tolist() + priorit_gain_vars_or['proxy_rsid'].unique().tolist()

deepmind_gene_df3['mtag_p_val_genome_wide_signif'] = ''
deepmind_gene_df3['mtag_or_above_1'] = ''

deepmind_gene_df3['mtag_p_val_genome_wide_signif'][deepmind_gene_df3['proxy_rsid'].isin(p_val)] = 'Y'
deepmind_gene_df3['mtag_or_above_1'][deepmind_gene_df3['proxy_rsid'].isin(o_r)] = 'Y'

deepmind_gene_df3.to_csv('/well/PROCARDIS/domwest/deepmind/prioritised_hcm_snps_final_result_table_part1.csv', index=False)


##SNP TABLE PART 2

deepmind_gene_df2 = pd.read_csv('/well/PROCARDIS/domwest/deepmind/prioritised_hcm_snps_final_result_table_part1.csv')

snps_og = pd.read_excel('/well/PROCARDIS/domwest/deephaem_prep/gathering_snps/dwest_deephaem_lead_imp_snps.xlsx')
snps_og['lead/proxy'] = ''
snps_og['lead/proxy'][snps_og['lead_rsid']==snps_og['proxy_rsid']] = 'lead'
snps_og['lead/proxy'][snps_og['lead_rsid']!=snps_og['proxy_rsid']] = 'proxy'

deepmind_gene_df2 = deepmind_gene_df2.merge(snps_og[['proxy_rsid', 'LD_r2']], how='left', on='proxy_rsid')
deepmind_gene_df2[['chr', 'pos_b38', 'ref', 'alt']] = deepmind_gene_df2['snp_test_id'].str.split('_', expand=True)

# summary_data = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/gathering_snps/regulome_db_pos_controls/hcm.fix.250621.n2.cm.txt', sep='\t')
summary_data = pd.read_csv('/well/PROCARDIS/domwest/replicate_old/compare_bulk_and_sc_quadrant_approach/hcm.mtag.2024.format.copy.tsv', sep='\t')

len(np.intersect1d(summary_data['rs_id'].unique().tolist(), deepmind_gene_df2['proxy_rsid'].unique().tolist()))
#47 out of 49 (when using hcm.fix and now 38 out of 49 when using hcm.mtag.2024)
summary_data_filt = summary_data[summary_data['rs_id'].isin(deepmind_gene_df2['proxy_rsid'].unique().tolist())]
summary_data_filt = summary_data_filt[['rs_id', 'effect_allele', 'other_allele', 'effect_allele_frequency', 'beta', 'standard_error', 'p_value']]
summary_data_filt.rename({'rs_id':'proxy_rsid'}, axis=1, inplace=True)

summary_data_filt_merg = summary_data_filt.merge(deepmind_gene_df2[['proxy_rsid', 'lead_rsid', 'chr', 'pos_b38', 'ref', 'alt', 'snp_test_id']], on='proxy_rsid', how='left')

leftover_snps = np.setdiff1d(deepmind_gene_df2['proxy_rsid'].unique().tolist(), summary_data['rs_id'].unique().tolist()).tolist()
leftover_summary_data_filt = deepmind_gene_df2[deepmind_gene_df2['proxy_rsid'].isin(leftover_snps)]
leftover_summary_data_filtt = deepmind_gene_df2[deepmind_gene_df2['proxy_rsid'].isin(leftover_snps)]
leftover_summary_data_filtt = leftover_summary_data_filtt[['lead_rsid', 'proxy_rsid', 'snp_test_id', 'chr', 'pos_b38', 'ref', 'alt']]
leftover_summary_data_filt = leftover_summary_data_filt[['proxy_rsid', 'chr', 'pos_b38', 'ref', 'alt']]

#add vep annotations

#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
1       154983396       .       A       G       .       .       .       .

summary_data_filt_merg['ID'] = '.'
summary_data_filt_merg['QUAL'] = ''
summary_data_filt_merg['FILTER'] = '.'
summary_data_filt_merg['INFO'] = '.'
summary_data_filt_merg['FORMAT'] = '.'
summary_data_filt_merg = summary_data_filt_merg[['chr', 'pos_b38', 'ID', 'ref', 'alt', 'QUAL', 'FILTER', 'INFO', 'FORMAT']]
summary_data_filt_merg.to_csv('/well/PROCARDIS/domwest/deepmind/prioritised_hcm_snps_vep_input_part2.vcf', index=False, sep='\t', header=None)


leftover_summary_data_filt['ID'] = '.'
leftover_summary_data_filt['QUAL'] = ''
leftover_summary_data_filt['FILTER'] = '.'
leftover_summary_data_filt['INFO'] = '.'
leftover_summary_data_filt['FORMAT'] = '.'
leftover_summary_data_filt = leftover_summary_data_filt[['chr', 'pos_b38', 'ID', 'ref', 'alt', 'QUAL', 'FILTER', 'INFO', 'FORMAT']]
leftover_summary_data_filt['chr'] = leftover_summary_data_filt['chr'].astype(int)
leftover_summary_data_filt['pos_b38'] = leftover_summary_data_filt['pos_b38'].astype(int)
leftover_summary_data_filt.sort_values(by=['chr', 'pos_b38'], inplace=True)
leftover_summary_data_filt.to_csv('/well/PROCARDIS/domwest/deepmind/leftover_prioritised_hcm_snps_vep_input_part2.vcf', index=False, sep='\t', header=None)


module purge
module load BCFtools

bgzip /well/PROCARDIS/domwest/deepmind/prioritised_hcm_snps_vep_input_part2.vcf
tabix -p vcf /well/PROCARDIS/domwest/deepmind/prioritised_hcm_snps_vep_input_part2.vcf.gz

bgzip /well/PROCARDIS/domwest/deepmind/leftover_prioritised_hcm_snps_vep_input_part2.vcf
tabix -p vcf /well/PROCARDIS/domwest/deepmind/leftover_prioritised_hcm_snps_vep_input_part2.vcf.gz

module purge
module load java/1.8.0_latest
module add perl/5.16.3
module load htslib/1.8-gcc5.4.0
module load Anaconda3/5.1.0
module load samtools/1.8-gcc5.4.0

export PERL5LIB=/well/PROCARDIS/cgrace/bin/perl_installs/lib/perl5/x86_64-linux:/well/PROCARDIS/cgrace/bin/perl_installs/lib/site_perl/5.16.3:/well/PROCARDIS/cgrace/bin/loftee_b38_v2/loftee:

bin=/well/PROCARDIS/cgrace/bin

/well/PROCARDIS/cgrace/KT.burden/vep -i /well/PROCARDIS/domwest/deepmind/prioritised_hcm_snps_vep_input_part2.vcf.gz -o /well/PROCARDIS/domwest/deepmind/prioritised_hcm_snps_vep_output_part2.out --assembly GRCh38 --vcf -offline --canonical --mane --hgvs --af_gnomade --af_gnomadg --pick --regulatory --dir_cache $bin/ensembl-vep/.vep/ --fasta /well/PROCARDIS/cgrace/KT.burden/Homo_sapiens.GRCh38.dna.primary_assembly.fa --plugin LoF,loftee_path:$bin/loftee_b38_v2/loftee,gerp_bigwig:$bin/loftee_b38_v2/files/gerp_conservation_scores.homo_sapiens.GRCh38.bw.1,human_ancestor_fa:$bin/loftee_b38_v2/files/human_ancestor.fa.gz,conservation_file:$bin/loftee_b38_v2/files/loftee.sql --dir_plugins $bin/vep_plugins/VEP_plugins --force_overwrite # --pick_order length,tsl,rank

/well/PROCARDIS/cgrace/KT.burden/vep -i /well/PROCARDIS/domwest/deepmind/leftover_prioritised_hcm_snps_vep_input_part2.vcf.gz -o /well/PROCARDIS/domwest/deepmind/leftover_prioritised_hcm_snps_vep_output_part2.out --assembly GRCh38 --vcf -offline --canonical --mane --hgvs --af_gnomade --af_gnomadg --pick --regulatory --dir_cache $bin/ensembl-vep/.vep/ --fasta /well/PROCARDIS/cgrace/KT.burden/Homo_sapiens.GRCh38.dna.primary_assembly.fa --plugin LoF,loftee_path:$bin/loftee_b38_v2/loftee,gerp_bigwig:$bin/loftee_b38_v2/files/gerp_conservation_scores.homo_sapiens.GRCh38.bw.1,human_ancestor_fa:$bin/loftee_b38_v2/files/human_ancestor.fa.gz,conservation_file:$bin/loftee_b38_v2/files/loftee.sql --dir_plugins $bin/vep_plugins/VEP_plugins --force_overwrite # --pick_order length,tsl,rank

# module purge
# module load java/1.8.0_latest
# module add perl/5.16.3
# module load htslib/1.8-gcc5.4.0
# module load Anaconda3/5.1.0
# module load samtools/1.8-gcc5.4.0

# export PERL5LIB=/well/PROCARDIS/cgrace/bin/perl_installs/lib/perl5/x86_64-linux:/well/PROCARDIS/cgrace/bin/perl_installs/lib/site_perl/5.16.3:/well/PROCARDIS/cgrace/bin/loftee_b38_v2/loftee:

# bin=/well/PROCARDIS/cgrace/bin

# /well/PROCARDIS/cgrace/KT.burden/vep -i /well/PROCARDIS/domwest/deepmind/leftover_prioritised_hcm_snps_vep_input_part2.vcf.gz -o /well/PROCARDIS/domwest/deepmind/leftover_prioritised_hcm_snps_vep_output_part2.out --assembly GRCh38 --vcf -offline --canonical --mane --hgvs --af_gnomade --af_gnomadg --pick --regulatory --dir_cache $bin/ensembl-vep/.vep/ --fasta /well/PROCARDIS/cgrace/KT.burden/Homo_sapiens.GRCh38.dna.primary_assembly.fa --plugin LoF,loftee_path:$bin/loftee_b38_v2/loftee,gerp_bigwig:$bin/loftee_b38_v2/files/gerp_conservation_scores.homo_sapiens.GRCh38.bw.1,human_ancestor_fa:$bin/loftee_b38_v2/files/human_ancestor.fa.gz,conservation_file:$bin/loftee_b38_v2/files/loftee.sql --dir_plugins $bin/vep_plugins/VEP_plugins --force_overwrite # --pick_order length,tsl,rank

import pandas as pd
import numpy as np

vep_out1 = pd.read_csv('/well/PROCARDIS/domwest/deepmind/prioritised_hcm_snps_vep_output_part2.out', sep='\t', comment='#', header=None)
vep_out2 = pd.read_csv('/well/PROCARDIS/domwest/deepmind/leftover_prioritised_hcm_snps_vep_output_part2.out', sep='\t', comment='#', header=None)

vep_out = pd.concat([vep_out1, vep_out2])

# info_str = 'Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|HGVS_OFFSET|gnomADe_AF|gnomADe_AFR_AF|gnomADe_AMR_AF|gnomADe_ASJ_AF|gnomADe_EAS_AF|gnomADe_FIN_AF|gnomADe_NFE_AF|gnomADe_OTH_AF|gnomADe_SAS_AF|gnomADg_AF|gnomADg_AFR_AF|gnomADg_AMI_AF|gnomADg_AMR_AF|gnomADg_ASJ_AF|gnomADg_EAS_AF|gnomADg_FIN_AF|gnomADg_MID_AF|gnomADg_NFE_AF|gnomADg_OTH_AF|gnomADg_SAS_AF|CLIN_SIG|SOMATIC|PHENO|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|LoF|LoF_filter|LoF_flags|LoF_info'

info_str = 'Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|HGVS_OFFSET|gnomADe_AF|gnomADe_AFR_AF|gnomADe_AMR_AF|gnomADe_ASJ_AF|gnomADe_EAS_AF|gnomADe_FIN_AF|gnomADe_NFE_AF|gnomADe_OTH_AF|gnomADe_SAS_AF|gnomADg_AF|gnomADg_AFR_AF|gnomADg_AMI_AF|gnomADg_AMR_AF|gnomADg_ASJ_AF|gnomADg_EAS_AF|gnomADg_FIN_AF|gnomADg_MID_AF|gnomADg_NFE_AF|gnomADg_OTH_AF|gnomADg_SAS_AF|CLIN_SIG|SOMATIC|PHENO|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS'

vep_out['vep_consequence'] = vep_out[7].str.split('|', expand=True)[1]
vep_out['vep_gene_symbol'] = vep_out[7].str.split('|', expand=True)[3]
vep_out['vep_picked_transcript'] = vep_out[7].str.split('|', expand=True)[6]
vep_out['canonical'] = vep_out[7].str.split('|', expand=True)[23]
vep_out['mane_select'] = vep_out[7].str.split('|', expand=True)[24]
vep_out['gnomADg_NFE_AF'] = vep_out[7].str.split('|', expand=True)[44]
vep_out['HGVSc'] = vep_out[7].str.split('|', expand=True)[10]
vep_out['HGVSp'] = vep_out[7].str.split('|', expand=True)[11]
vep_out['motif_name'] = vep_out[7].str.split('|', expand=True)[50]
vep_out['motif_pos'] = vep_out[7].str.split('|', expand=True)[51]
vep_out['high_inf_pos'] = vep_out[7].str.split('|', expand=True)[52]
vep_out['motif_score_change'] = vep_out[7].str.split('|', expand=True)[53]

vep_out['snp_test_id'] = vep_out[0].astype(str) + '_' + vep_out[1].astype(str) + '_' + vep_out[3] + '_' + vep_out[4]
vep_out = vep_out[['snp_test_id', 'vep_consequence', 'vep_gene_symbol', 'vep_picked_transcript', 'canonical', 'mane_select', 'gnomADg_NFE_AF', 'HGVSc', 'HGVSp', 'motif_name', 'motif_pos', 'high_inf_pos', 'motif_score_change']]

fin_part2 = summary_data_filt_merg.merge(vep_out, on='snp_test_id', how='right')

nonnull_fin_part2 = fin_part2[fin_part2['proxy_rsid'].notnull()]

null_fin_part2 = fin_part2[fin_part2['proxy_rsid'].isnull()]
null_fin_part2.drop(['proxy_rsid', 'effect_allele', 'other_allele', 'effect_allele_frequency', 'beta', 'standard_error', 'p_value', 'lead_rsid', 'chr', 'pos_b38', 'ref', 'alt'], axis=1, inplace=True)
null_fin_part3 = null_fin_part2.merge(leftover_summary_data_filtt, how='left', on='snp_test_id')
# Index(['snp_test_id', 'vep_consequence', 'vep_gene_symbol',
#        'vep_picked_transcript', 'canonical', 'mane_select', 'gnomADg_NFE_AF',
#        'HGVSc', 'HGVSp', 'motif_name', 'motif_pos', 'high_inf_pos',
#        'motif_score_change', 'proxy_rsid', 'chr', 'pos_b38', 'ref', 'alt'],
#       dtype='object')
null_fin_part2['effect_allele'] = ''
null_fin_part2['other_allele'] = ''
null_fin_part2['effect_allele_frequency'] = ''
null_fin_part2['beta'] = ''
null_fin_part2['standard_error'] = ''
null_fin_part2['p_value'] = ''

tot = pd.concat([nonnull_fin_part2, null_fin_part3])

tot.to_csv('/well/PROCARDIS/domwest/deepmind/prioritised_hcm_snps_final_result_table_part2.csv', index=False)

#######################################################################################################################################################################
#######################################################################################################################################################################
#######################################################################################################################################################################
#######################################################################################################################################################################
#######################################################################################################################################################################
#######################################################################################################################################################################
#######################################################################################################################################################################
#######################################################################################################################################################################
#######################################################################################################################################################################
#######################################################################################################################################################################


#add OT annotations 

#gonna do it manually

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

# top_df.to_csv('/gpfs3/well/PROCARDIS/domwest/deephaem_prep/post_analysis/opentarget/opentarget_snp_to_gene.csv', index=False) #1362 rows

top_df = top_df.sort_values(by=['Gene'])
x = pd.DataFrame(top_df[['SNP', 'Gene', 'Overall V2G']].groupby(['Gene', 'Overall V2G']).count()).reset_index()
x.rename({'SNP':'variant_count'}, axis=1, inplace=True)
top_df_x = top_df.merge(x, how='left', on=['Gene', 'Overall V2G'])
top_df_x = top_df_x[['Gene', 'Overall V2G', 'variant_count']].drop_duplicates()
# top_df_x.to_csv('/gpfs3/well/PROCARDIS/domwest/deephaem_prep/post_analysis/opentarget/opentarget_snp_to_gene_for_plotting.csv', index=False) #1362 rows

library(ggplot2)
library(scales)
theme_set(theme_classic())

top_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deephaem_prep/post_analysis/opentarget/opentarget_snp_to_gene_for_plotting.csv') #was: with_vars_present_pos_controls_and_4395_lead_imp_vars_damage_scores

# Change dot plot colors by groups
p<-ggplot(top_df, aes(x=Gene, y=Overall.V2G, fill=Gene)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.2)  +
  theme(legend.position = "none", panel.grid.major.x = element_line(color = "red", size = 0.5, linetype = 2))
p

ggsave(filename="/gpfs3/well/PROCARDIS/domwest/deephaem_prep/post_analysis/opentarget/overall_v2g.tiff", width=18, height=15, dpi = 300) 

q<-ggplot(top_df, aes(x=Gene, y=Overall.V2G, fill=factor(variant_count))) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.grid.major.x = element_line(color = "red", size = 0.5, linetype = 2))
q

ggsave(filename="/gpfs3/well/PROCARDIS/domwest/deephaem_prep/post_analysis/opentarget/overall_v2g2.tiff", width=18, height=15, dpi = 300) 








###################################################################################################################################################################

#making damage score tables per locus for chpater 4:

import pandas as pd
import numpy as np
##SNP TABLE PART 1
quadr = pd.read_csv('/well/PROCARDIS/domwest/replicate_old/compare_bulk_and_sc_quadrant_approach/final_bulk_and_sc_prioritised_genes_per_lead_snp_refseq_genes_only.csv')
quadr.rename({'rsid':'proxy_rsid'}, axis=1, inplace=True)
quadr = quadr[['proxy_rsid', 'bulk_gene_symbol', 'sc_gene_symbol', 'same_bulk_and_sc_gene', 'sc_enriched_in_cm']]
#latest scatter plots are here: C:\Users\domwest\Documents\oxford_dphil\compare_new_quadrant_appr_ilo_tads_to_old\new_scatter_plots_rescomp
ot_fuma = pd.read_excel('/well/PROCARDIS/domwest/deephaem_prep/gathering_snps/HCM_MTAG_OT_FUMA.xlsx')
ot_fuma.rename({'rsid':'proxy_rsid', 'snptestid':'snp_test_id', 'likely gene(s) (updated 17/11/22)':'OT/FUMA'}, axis=1, inplace=True)
ot_fuma['snp_test_id'] = ot_fuma['snp_test_id'].str.replace(':','_')
ot_fuma['snp_test_id'] = ot_fuma['snp_test_id'].str.replace('chr','')
ot_fuma = ot_fuma[['proxy_rsid', 'OT/FUMA']]
ot_fuma[['OT', 'FUMA']] = ot_fuma['OT/FUMA'].str.split('/', expand=True)
ot_fuma.drop(['OT/FUMA'], axis=1, inplace=True)
ot_fuma['FUMA'] = ot_fuma['FUMA'].fillna(ot_fuma['OT'])
gene_df =  ot_fuma.merge(quadr, on='proxy_rsid', how='left')
gene_df['lead/proxy'] = 'lead'
deepmind = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/all_heart_classifiers_hcm_prioritised_with_lead_and_proxy_info.csv')
len(deepmind['rsid'].unique())
deepmind.rename({'rsid':'proxy_rsid'}, axis=1, inplace=True)


dfs = deepmind
snps_og = pd.read_excel('/well/PROCARDIS/domwest/deephaem_prep/gathering_snps/dwest_deephaem_lead_imp_snps.xlsx')
dfs = dfs.merge(snps_og[['proxy_rsid', 'lead_rsid']], how='left', on='proxy_rsid')
df = dfs[['proxy_rsid', 'snp_test_id', 'lead_rsid', 'classifier', 'damage_score_value']].drop_duplicates()
df['proxy_rsid:snp_test_id'] = df['proxy_rsid'] + ':' + df['snp_test_id']
df.drop(['proxy_rsid', 'snp_test_id'], axis=1, inplace=True)
df = df.sort_values(by=['lead_rsid'])
for lead in df['lead_rsid'].unique().tolist():
	tmp = df[df['lead_rsid']==lead]
	tmp.drop(['lead_rsid'], axis=1, inplace=True)
	tmp_df = pd.DataFrame(tmp.pivot(index='proxy_rsid:snp_test_id', columns='classifier', values='damage_score_value')).reset_index()
	tmp_df = tmp_df.fillna('/')
	tmp_df.to_csv(lead + '_wide_df_snps_passing_threshold.csv', index=False)
