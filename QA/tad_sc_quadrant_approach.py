from cmapPy.pandasGEXpress.parse_gct import parse
from cmapPy.pandasGEXpress.write_gct import write
import pandas as pd
import gffpandas.gffpandas as gffpd
import os
import numpy as np
import statistics

#doing this now

#run scrnaseq_expr_data_quadrant_appr.py

##############################################################################################################################################################################################################

##FETCHING THE LOCI GENE SET

def quadr_appr_for_diff_cells(cell_type):
   genes = pd.read_csv('/well/PROCARDIS/domwest/replicate_old/using_tads/for_68_quadrant_approach_input.csv') #the ensembl output files that made up this whole df (via for_68_process_ensembl_api_output) is here: /well/PROCARDIS/domwest/hic_defined_TADs/for_68_togeth
   #genes.drop(['gene_symbol'], axis=1, inplace=True)
   #
   means_all_dfs = pd.read_csv('/well/PROCARDIS/domwest/replicate_old/using_tads/scrnaseq/' + cell_type + '_and_all_cols.csv')
   means_all_dfs.rename({cell_type + '_means':'mean_expr'}, axis=1, inplace=True)
   # means_all_dfs.rename({cell_type + '_median_of_means':'median_of_means'}, axis=1, inplace=True)
   means_all_dfs.rename({'fold_change_' + cell_type:'fold_change'}, axis=1, inplace=True)
   means_all_dfs = means_all_dfs[['gene_symbol', 'mean_expr', 'median_of_means', 'fold_change']]
   means_all_dfs = means_all_dfs.drop_duplicates()
   #
   combined_lv = genes.merge(means_all_dfs, on='gene_symbol')
   #
   loci_expr_and_foldchange = combined_lv[['rsid', 'gene_id', 'gene_symbol', 'mean_expr', 'fold_change']]
   loci_expr_and_foldchange = loci_expr_and_foldchange[~(loci_expr_and_foldchange['mean_expr']==0)]
   len(loci_expr_and_foldchange['rsid'].unique())
   #46
   #Quadrant for 3 loci:
   #first determine the median of median_expr in the heart so that we know what separates the different quadrants eg we know that fold change >= 1.3 is considered high and that below would be considered low
   #getting the df one with infinity fold change values included:
   #we want to give the inf values a more accurate comparison and so wherever there is an inf value in the fold change column, take the highest fold change value in the df (besides inf) and give it that value instead
   #loci_expr_and_foldchange['fold_change_lv'][~np.isinf(loci_expr_and_foldchange['fold_change_lv'])].max() #this gives 5873.9785
   #loci_expr_and_foldchange['fold_change_lv'] = loci_expr_and_foldchange['fold_change_lv'].replace([np.inf], 5873.9785)
   #trying a different approach here where we take the maximum fold change value within that particular locus rather than across all the loci:
   # if cell_type == 'cm':
   #    testing = []
   #    for un in loci_expr_and_foldchange['rsid'].unique():
   #       if un != 'rs6566955':
   #          temp = loci_expr_and_foldchange[loci_expr_and_foldchange['rsid']==un]
   #          temp2 = temp.loc[temp['fold_change'] != np.inf, 'fold_change'].max()
   #          temp['fold_change'].replace(np.inf,temp2,inplace=True)
   #          testing.append(temp)
   #    final_testing_1 = pd.concat(testing)
   #    to_add = loci_expr_and_foldchange[loci_expr_and_foldchange['rsid']=='rs6566955']
   #    final_testing = pd.concat([final_testing_1, to_add])
   #    #
   #    testing2 = []
   #    for un in loci_expr_and_foldchange['rsid'].unique():
   #       if un != 'rs6566955':
   #          temp = loci_expr_and_foldchange[loci_expr_and_foldchange['rsid']==un]
   #          temp2 = temp.loc[temp['fold_change'] != np.inf, 'fold_change'].max()
   #          testing2.append(temp2)
   #    median_of_inf_max = statistics.median(testing2)
   #    final_testing['fold_change'][final_testing['rsid']=='rs6566955'] = median_of_inf_max
   # else:
   #    testing = []
   #    for un in loci_expr_and_foldchange['rsid'].unique():
   #       temp = loci_expr_and_foldchange[loci_expr_and_foldchange['rsid']==un]
   #       temp2 = temp.loc[temp['fold_change'] != np.inf, 'fold_change'].max()
   #       temp['fold_change'].replace(np.inf,temp2,inplace=True)
   #       testing.append(temp)
   #    final_testing = pd.concat(testing)
   #
   testing = []
   for un in loci_expr_and_foldchange['rsid'].unique():
      temp = loci_expr_and_foldchange[loci_expr_and_foldchange['rsid']==un]
      temp2 = temp.loc[temp['fold_change'] != np.inf, 'fold_change'].max()
      temp['fold_change'].replace(np.inf,temp2,inplace=True)
      testing.append(temp)
   final_testing = pd.concat(testing)
   #
   loci_expr_and_foldchange = final_testing
   #
   #should maybe consider doing this without including genes that had expression = 0 in the mix. Since we want to get an idea of the typical range of expression among those genes that are actually important in the heart in order to determine an expression cutoff in the quadrants. We also want to include all genes in the heart, so must go back to older df long_initial_merge_lv
   means_all_dfs = means_all_dfs[~(means_all_dfs['mean_expr']==0)]
   median_of_mean_expr =means_all_dfs["mean_expr"].median() #0.014407334643091
   #
   loci_expr_and_foldchange = loci_expr_and_foldchange.sort_values(['fold_change', 'mean_expr'], ascending = [False, False])
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #using top gene metric to prioritise loci
   #Sort genes within loci by BOTH fold change as well as median_expr
   #the way it works is to sort by fold change and then expression
   loci_expr_and_foldchange['gradient_per_gene'] = loci_expr_and_foldchange['fold_change']/loci_expr_and_foldchange['mean_expr']
   #
   #now determine which quadrant it is in:
   loci_expr_and_foldchange['quadrant'] = ''
   loci_expr_and_foldchange['quadrant'][(loci_expr_and_foldchange['fold_change']>=1.3) & (loci_expr_and_foldchange['mean_expr']<median_of_mean_expr)] = 1
   loci_expr_and_foldchange['quadrant'][(loci_expr_and_foldchange['fold_change']>=1.3) & (loci_expr_and_foldchange['mean_expr']>=median_of_mean_expr)] = 2
   loci_expr_and_foldchange['quadrant'][(loci_expr_and_foldchange['fold_change']<1.3) & (loci_expr_and_foldchange['mean_expr']<median_of_mean_expr)] = 3
   loci_expr_and_foldchange['quadrant'][(loci_expr_and_foldchange['fold_change']<1.3) & (loci_expr_and_foldchange['mean_expr']>=median_of_mean_expr)] = 4
   #adding diagonal distance (ie hypotenuse to fold change and median expr) to help sort in addition (and after) sorting by custom quadrant
   loci_expr_and_foldchange['diagonal_distance'] = np.sqrt((loci_expr_and_foldchange['fold_change'] ** 2) + (loci_expr_and_foldchange['mean_expr'] ** 2))
   loci_expr_and_foldchange.sort_values(by=['rsid', 'quadrant', 'diagonal_distance'], inplace=True, ascending=[True,True,False])
   #
   #for ranking loci:
   loci_expr_and_foldchange['locus_prioritised'] = ''
   locus_prioritisation = pd.DataFrame()
   for targ in loci_expr_and_foldchange['rsid'].unique():
      quad2 = loci_expr_and_foldchange[(loci_expr_and_foldchange['rsid']==targ) & (loci_expr_and_foldchange['quadrant']==2)]
      quad1 = loci_expr_and_foldchange[(loci_expr_and_foldchange['rsid']==targ) & (loci_expr_and_foldchange['quadrant']==1)]
      quad3 = loci_expr_and_foldchange[(loci_expr_and_foldchange['rsid']==targ) & (loci_expr_and_foldchange['quadrant']==3)]
      quad4 = loci_expr_and_foldchange[(loci_expr_and_foldchange['rsid']==targ) & (loci_expr_and_foldchange['quadrant']==4)]
      if quad2.shape[0] > 0:
         #get top one ie already sorted by diagonal_distance, so this gets the longest distance
         loci_expr_and_foldchange['locus_prioritised'] = np.where( ( (loci_expr_and_foldchange['rsid']==targ) & (loci_expr_and_foldchange['quadrant']==2) ), 'Y', 'N')
         prior_loci_2 = loci_expr_and_foldchange[loci_expr_and_foldchange['locus_prioritised']=='Y'].head(1)
         #print(prior_loci_2)
         locus_prioritisation = pd.concat([locus_prioritisation, prior_loci_2])
         #print(locus_prioritisation)
         #loci_expr_and_foldchange['locus_prioritised'][(loci_expr_and_foldchange['gene_openTarget']==targ) & (loci_expr_and_foldchange['quadrant']==2)].head(1) = 'Y'
         #loci_expr_and_foldchange['locus_prioritised'][~(loci_expr_and_foldchange['gene_openTarget']==targ) & (loci_expr_and_foldchange['quadrant']==2)] = 'N'
      elif quad1.shape[0] > 0:
         loci_expr_and_foldchange['locus_prioritised'] = np.where( ( (loci_expr_and_foldchange['rsid']==targ) & (loci_expr_and_foldchange['quadrant']==1) ), 'Y', 'N')
         prior_loci_1 = loci_expr_and_foldchange[loci_expr_and_foldchange['locus_prioritised']=='Y'].head(1)   
         #print(prior_loci_1)
         locus_prioritisation = pd.concat([locus_prioritisation, prior_loci_1])
         #print(locus_prioritisation)
      elif quad4.shape[0] > 0:
         loci_expr_and_foldchange['locus_prioritised'] = np.where( ( (loci_expr_and_foldchange['rsid']==targ) & (loci_expr_and_foldchange['quadrant']==4) ), 'Y', 'N')
         prior_loci_4 = loci_expr_and_foldchange[loci_expr_and_foldchange['locus_prioritised']=='Y'].head(1)
         #print(prior_loci_4)
         locus_prioritisation = pd.concat([locus_prioritisation, prior_loci_4])
         #print(locus_prioritisation)
      else:
         loci_expr_and_foldchange['locus_prioritised'][loci_expr_and_foldchange['rsid']==targ] = 'Cannot prioritise locus due no gradients in quadrant 2, 1, or 4'
         cannot_prior = loci_expr_and_foldchange[loci_expr_and_foldchange['locus_prioritised']=='Cannot prioritise locus due no gradients in quadrant 2, 1, or 4']
         print(cannot_prior)
         locus_prioritisation = pd.concat([locus_prioritisation, cannot_prior])
         print(locus_prioritisation)
         #
   locus_prioritisation.drop(['locus_prioritised'], axis=1, inplace=True)
   #
   locus_prioritisation['quadrant'] = pd.Categorical(locus_prioritisation['quadrant'], [2, 1, 4, 3])
   locus_prioritisation.sort_values("quadrant", inplace=True)
   #
   # locus_prioritisation['mean_expr'] = locus_prioritisation['mean_expr'].round(decimals = 3)
   # locus_prioritisation['fold_change'] = locus_prioritisation['fold_change'].round(decimals = 3)
   # locus_prioritisation['gradient_per_gene'] = locus_prioritisation['gradient_per_gene'].round(decimals = 3)
   # locus_prioritisation['diagonal_distance'] = locus_prioritisation['diagonal_distance'].round(decimals = 3)
   locus_prioritisation.to_csv('/well/PROCARDIS/domwest/replicate_old/using_tads/' + cell_type + '_loci_prioritisation_sorted_by_quadrant.csv', index=False)
   #
   #RERUN_START>>>
   loci_expr_and_foldchange = final_testing
   #
   #should maybe consider doing this without including genes that had expression = 0 in the mix. Since we want to get an idea of the typical range of expression among those genes that are actually important in the heart in order to determine an expression cutoff in the quadrants. We also want to include all genes in the heart, so must go back to older df long_initial_merge_lv
   means_all_dfs = means_all_dfs[~(means_all_dfs['mean_expr']==0)]
   median_of_mean_expr =means_all_dfs["mean_expr"].median() #0.014407334643091
   #
   loci_expr_and_foldchange = loci_expr_and_foldchange.sort_values(['fold_change', 'mean_expr'], ascending = [False, False])
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #using top gene metric to prioritise loci
   #Sort genes within loci by BOTH fold change as well as median_expr
   #the way it works is to sort by fold change and then expression
   loci_expr_and_foldchange['gradient_per_gene'] = loci_expr_and_foldchange['fold_change']/loci_expr_and_foldchange['mean_expr']
   #
   #now determine which quadrant it is in:
   loci_expr_and_foldchange['quadrant'] = ''
   loci_expr_and_foldchange['quadrant'][(loci_expr_and_foldchange['fold_change']>=1.3) & (loci_expr_and_foldchange['mean_expr']<median_of_mean_expr)] = 1
   loci_expr_and_foldchange['quadrant'][(loci_expr_and_foldchange['fold_change']>=1.3) & (loci_expr_and_foldchange['mean_expr']>=median_of_mean_expr)] = 2
   loci_expr_and_foldchange['quadrant'][(loci_expr_and_foldchange['fold_change']<1.3) & (loci_expr_and_foldchange['mean_expr']<median_of_mean_expr)] = 3
   loci_expr_and_foldchange['quadrant'][(loci_expr_and_foldchange['fold_change']<1.3) & (loci_expr_and_foldchange['mean_expr']>=median_of_mean_expr)] = 4
   #adding diagonal distance (ie hypotenuse to fold change and median expr) to help sort in addition (and after) sorting by custom quadrant
   loci_expr_and_foldchange['diagonal_distance'] = np.sqrt((loci_expr_and_foldchange['fold_change'] ** 2) + (loci_expr_and_foldchange['mean_expr'] ** 2))
   loci_expr_and_foldchange.sort_values(by=['rsid', 'quadrant', 'diagonal_distance'], inplace=True, ascending=[True,True,False])
   #RERUN_END<<<
   #
   gene_prioritisation = loci_expr_and_foldchange
   # gene_prioritisation.sort_values(["rsid", "quadrant", "diagonal_distance"], inplace=True, ascending=[True, True, False])
   sorter = [2, 1, 4, 3]
   x = pd.DataFrame({'quadrant': sorter})
   x.index = x.index.set_names('number')
   x = x.reset_index()
   gene_prioritisation = pd.merge(gene_prioritisation, x, how='left', on='quadrant')
   gene_prioritisation.sort_values(['rsid', 'number', 'diagonal_distance'], ascending = [True, True, False], inplace=True)
   gene_prioritisation.drop(['number'], axis=1, inplace=True)
   # gene_prioritisation['mean_expr'] = gene_prioritisation['mean_expr'].round(decimals = 3)
   # gene_prioritisation['fold_change'] = gene_prioritisation['fold_change'].round(decimals = 3)
   # gene_prioritisation['gradient_per_gene'] = gene_prioritisation['gradient_per_gene'].round(decimals = 3)
   # gene_prioritisation['diagonal_distance'] = gene_prioritisation['diagonal_distance'].round(decimals = 3)
   gene_prioritisation.to_csv('/well/PROCARDIS/domwest/replicate_old/using_tads/' + cell_type + '_gene_prioritisation.csv', index=False)

quadr_appr_for_diff_cells('cm')
quadr_appr_for_diff_cells('ad')
quadr_appr_for_diff_cells('ec')
quadr_appr_for_diff_cells('fb')
quadr_appr_for_diff_cells('lymphoid')
quadr_appr_for_diff_cells('mast')
quadr_appr_for_diff_cells('myeloid')

###################################################################################################################################################

#CM ONLY extra plots to run AFTER comparing bulk and SC and getting the new list of SC genes


gene_prioritisation = pd.read_csv('/well/PROCARDIS/domwest/replicate_old/using_tads/cm_gene_prioritisation.csv')
gene_prioritisation['cell_type'] = 'CM'
gene_prioritisation.to_csv('/well/PROCARDIS/domwest/replicate_old/using_tads/cm_gene_prioritisation_with_celltype_col.csv', index=False)

# sc_genes_of_interest = ['HSPB7', 'TPD52L1', 'PXN-AS1', 'PROM1', 'PROB1', 'PLN', 'SPON1-AS1', 'PPP1R13L', 'PLPP7', 'ANKRD34C-AS1', 'LRP11', 'PLK2', 'PRKCZ', 'MAPT', 'MFHAS1', 'MAIP1', 'LINC00964', 'CDC42EP3', 'PLD1', 'DTNA', 'CHCHD10']
sc_genes_of_interest = ['LRP11', 'HSPB7', 'TPD52L1', 'NEXN', 'PROB1', 'PRKCZ', 'PXN-AS1', 'MFHAS1', 'RGS6', 'PLK2', 'PLPP7', 'NAXD', 'SPON1-AS1', 'DTNA', 'CHCHD10', 'MLIP-AS1', 'PLD1', 'MAPT', 'ANKRD34C-AS1', 'TBX3', 'ALPK2', 'MAIP1', 'LINC00964', 'PROM1', 'PLN', 'ALPK3', 'CDC42EP3', 'PPP1R13L']

gene_prioritisation_priorit = gene_prioritisation[gene_prioritisation['gene_symbol'].isin(sc_genes_of_interest)]
gene_prioritisation_priorit = gene_prioritisation_priorit.drop_duplicates()
gene_prioritisation_priorit.to_csv("/well/PROCARDIS/domwest/replicate_old/using_tads/gene_prioritisation_used_priorit_sc_genes_for_top_53_vars.csv", index=False)


.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/well/PROCARDIS/domwest/replicate_old/using_tads/cm_gene_prioritisation_with_celltype_col.csv')

p <- ggplot(long_prediction_df, aes(x=cell_type, y=mean_expr)) +
geom_boxplot(outlier.colour="black", outlier.shape=1, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
scale_fill_brewer(palette="Dark2")+
scale_y_continuous(n.breaks = 10) +
theme(text=element_text(size=18), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = 'none') +
labs(x = "Mean transcript counts distribution across genes in Cardiomyocyte") 

ggsave(filename="/well/PROCARDIS/domwest/replicate_old/using_tads/cm_mean_transcript_counts_distribution_boxplot_used_nonzero_counts.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

long_prediction_df <- read.csv('/well/PROCARDIS/domwest/replicate_old/using_tads/gene_prioritisation_used_priorit_sc_genes_for_top_53_vars.csv')

ggplot(long_prediction_df, aes(x=mean_expr, y=fold_change)) +
# geom_hline(yintercept = 1.3, color="#99CCCC", size =0.5) +
# geom_vline(xintercept = 1.0395, color="#99CCCC", size =0.5) +
geom_point(colour = "#336666") + 
geom_text_repel(aes(label = gene_symbol)) +
labs(title = "Fold change vs mean transcript counts for prioritised SC CM genes identified across case studies",
   subtitle = NULL,
   tag = NULL, 
   x = "Mean transcript counts",
   y= "Fold change values compared to other heart cell types",
   color = NULL) +
theme(text=element_text(size=18))
ggsave(filename="/well/PROCARDIS/domwest/replicate_old/using_tads/scatter_fold_change_cm_means_used_priorit_sc_genes_for_top_53_vars.tiff",width=18, height=15, dpi = 300)























































cell_type = 'cm'
genes = pd.read_csv('/well/PROCARDIS/domwest/replicate_old/using_tads/for_68_quadrant_approach_input.csv') #the ensembl output files that made up this whole df (via for_68_process_ensembl_api_output) is here: /well/PROCARDIS/domwest/hic_defined_TADs/for_68_togeth
#genes.drop(['gene_symbol'], axis=1, inplace=True)
#
means_all_dfs = pd.read_csv('/well/PROCARDIS/domwest/replicate_old/using_tads/scrnaseq/' + cell_type + '_and_all_cols.csv')
means_all_dfs.rename({cell_type + '_means':'mean_expr'}, axis=1, inplace=True)
# means_all_dfs.rename({cell_type + '_median_of_means':'median_of_means'}, axis=1, inplace=True)
means_all_dfs.rename({'fold_change_' + cell_type:'fold_change'}, axis=1, inplace=True)
means_all_dfs = means_all_dfs[['gene_symbol', 'mean_expr', 'median_of_means', 'fold_change']]
means_all_dfs = means_all_dfs.drop_duplicates()
#
combined_lv = genes.merge(means_all_dfs, on='gene_symbol')
#
loci_expr_and_foldchange = combined_lv[['rsid', 'gene_id', 'gene_symbol', 'mean_expr', 'fold_change']]
loci_expr_and_foldchange = loci_expr_and_foldchange[~(loci_expr_and_foldchange['mean_expr']==0)]
len(loci_expr_and_foldchange['rsid'].unique())
#46
#Quadrant for 3 loci:
#first determine the median of median_expr in the heart so that we know what separates the different quadrants eg we know that fold change >= 1.3 is considered high and that below would be considered low
#getting the df one with infinity fold change values included:
#we want to give the inf values a more accurate comparison and so wherever there is an inf value in the fold change column, take the highest fold change value in the df (besides inf) and give it that value instead
#loci_expr_and_foldchange['fold_change_lv'][~np.isinf(loci_expr_and_foldchange['fold_change_lv'])].max() #this gives 5873.9785
#loci_expr_and_foldchange['fold_change_lv'] = loci_expr_and_foldchange['fold_change_lv'].replace([np.inf], 5873.9785)
#trying a different approach here where we take the maximum fold change value within that particular locus rather than across all the loci:
# if cell_type == 'cm':
#    testing = []
#    for un in loci_expr_and_foldchange['rsid'].unique():
#       if un != 'rs6566955':
#          temp = loci_expr_and_foldchange[loci_expr_and_foldchange['rsid']==un]
#          temp2 = temp.loc[temp['fold_change'] != np.inf, 'fold_change'].max()
#          temp['fold_change'].replace(np.inf,temp2,inplace=True)
#          testing.append(temp)
#    final_testing_1 = pd.concat(testing)
#    to_add = loci_expr_and_foldchange[loci_expr_and_foldchange['rsid']=='rs6566955']
#    final_testing = pd.concat([final_testing_1, to_add])
#    #
#    testing2 = []
#    for un in loci_expr_and_foldchange['rsid'].unique():
#       if un != 'rs6566955':
#          temp = loci_expr_and_foldchange[loci_expr_and_foldchange['rsid']==un]
#          temp2 = temp.loc[temp['fold_change'] != np.inf, 'fold_change'].max()
#          testing2.append(temp2)
#    median_of_inf_max = statistics.median(testing2)
#    final_testing['fold_change'][final_testing['rsid']=='rs6566955'] = median_of_inf_max
# else:
#    testing = []
#    for un in loci_expr_and_foldchange['rsid'].unique():
#       temp = loci_expr_and_foldchange[loci_expr_and_foldchange['rsid']==un]
#       temp2 = temp.loc[temp['fold_change'] != np.inf, 'fold_change'].max()
#       temp['fold_change'].replace(np.inf,temp2,inplace=True)
#       testing.append(temp)
#    final_testing = pd.concat(testing)
#
testing = []
for un in loci_expr_and_foldchange['rsid'].unique():
   temp = loci_expr_and_foldchange[loci_expr_and_foldchange['rsid']==un]
   temp2 = temp.loc[temp['fold_change'] != np.inf, 'fold_change'].max()
   temp['fold_change'].replace(np.inf,temp2,inplace=True)
   testing.append(temp)
final_testing = pd.concat(testing)
#
loci_expr_and_foldchange = final_testing
#
#should maybe consider doing this without including genes that had expression = 0 in the mix. Since we want to get an idea of the typical range of expression among those genes that are actually important in the heart in order to determine an expression cutoff in the quadrants. We also want to include all genes in the heart, so must go back to older df long_initial_merge_lv
means_all_dfs = means_all_dfs[~(means_all_dfs['mean_expr']==0)]
median_of_mean_expr =means_all_dfs["mean_expr"].median() #0.014407334643091
#
loci_expr_and_foldchange = loci_expr_and_foldchange.sort_values(['fold_change', 'mean_expr'], ascending = [False, False])
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#using top gene metric to prioritise loci
#Sort genes within loci by BOTH fold change as well as median_expr
#the way it works is to sort by fold change and then expression
loci_expr_and_foldchange['gradient_per_gene'] = loci_expr_and_foldchange['fold_change']/loci_expr_and_foldchange['mean_expr']
#
#now determine which quadrant it is in:
loci_expr_and_foldchange['quadrant'] = ''
loci_expr_and_foldchange['quadrant'][(loci_expr_and_foldchange['fold_change']>=1.3) & (loci_expr_and_foldchange['mean_expr']<median_of_mean_expr)] = 1
loci_expr_and_foldchange['quadrant'][(loci_expr_and_foldchange['fold_change']>=1.3) & (loci_expr_and_foldchange['mean_expr']>=median_of_mean_expr)] = 2
loci_expr_and_foldchange['quadrant'][(loci_expr_and_foldchange['fold_change']<1.3) & (loci_expr_and_foldchange['mean_expr']<median_of_mean_expr)] = 3
loci_expr_and_foldchange['quadrant'][(loci_expr_and_foldchange['fold_change']<1.3) & (loci_expr_and_foldchange['mean_expr']>=median_of_mean_expr)] = 4
#adding diagonal distance (ie hypotenuse to fold change and median expr) to help sort in addition (and after) sorting by custom quadrant
loci_expr_and_foldchange['diagonal_distance'] = np.sqrt((loci_expr_and_foldchange['fold_change'] ** 2) + (loci_expr_and_foldchange['mean_expr'] ** 2))
loci_expr_and_foldchange.sort_values(by=['rsid', 'quadrant', 'diagonal_distance'], inplace=True, ascending=[True,True,False])
#
#for ranking loci:
loci_expr_and_foldchange['locus_prioritised'] = ''
locus_prioritisation = pd.DataFrame()
for targ in loci_expr_and_foldchange['rsid'].unique():
   quad2 = loci_expr_and_foldchange[(loci_expr_and_foldchange['rsid']==targ) & (loci_expr_and_foldchange['quadrant']==2)]
   quad1 = loci_expr_and_foldchange[(loci_expr_and_foldchange['rsid']==targ) & (loci_expr_and_foldchange['quadrant']==1)]
   quad3 = loci_expr_and_foldchange[(loci_expr_and_foldchange['rsid']==targ) & (loci_expr_and_foldchange['quadrant']==3)]
   quad4 = loci_expr_and_foldchange[(loci_expr_and_foldchange['rsid']==targ) & (loci_expr_and_foldchange['quadrant']==4)]
   if quad2.shape[0] > 0:
      #get top one ie already sorted by diagonal_distance, so this gets the longest distance
      loci_expr_and_foldchange['locus_prioritised'] = np.where( ( (loci_expr_and_foldchange['rsid']==targ) & (loci_expr_and_foldchange['quadrant']==2) ), 'Y', 'N')
      prior_loci_2 = loci_expr_and_foldchange[loci_expr_and_foldchange['locus_prioritised']=='Y'].head(1)
      #print(prior_loci_2)
      locus_prioritisation = pd.concat([locus_prioritisation, prior_loci_2])
      #print(locus_prioritisation)
      #loci_expr_and_foldchange['locus_prioritised'][(loci_expr_and_foldchange['gene_openTarget']==targ) & (loci_expr_and_foldchange['quadrant']==2)].head(1) = 'Y'
      #loci_expr_and_foldchange['locus_prioritised'][~(loci_expr_and_foldchange['gene_openTarget']==targ) & (loci_expr_and_foldchange['quadrant']==2)] = 'N'
   elif quad1.shape[0] > 0:
      loci_expr_and_foldchange['locus_prioritised'] = np.where( ( (loci_expr_and_foldchange['rsid']==targ) & (loci_expr_and_foldchange['quadrant']==1) ), 'Y', 'N')
      prior_loci_1 = loci_expr_and_foldchange[loci_expr_and_foldchange['locus_prioritised']=='Y'].head(1)   
      #print(prior_loci_1)
      locus_prioritisation = pd.concat([locus_prioritisation, prior_loci_1])
      #print(locus_prioritisation)
   elif quad4.shape[0] > 0:
      loci_expr_and_foldchange['locus_prioritised'] = np.where( ( (loci_expr_and_foldchange['rsid']==targ) & (loci_expr_and_foldchange['quadrant']==4) ), 'Y', 'N')
      prior_loci_4 = loci_expr_and_foldchange[loci_expr_and_foldchange['locus_prioritised']=='Y'].head(1)
      #print(prior_loci_4)
      locus_prioritisation = pd.concat([locus_prioritisation, prior_loci_4])
      #print(locus_prioritisation)
   else:
      loci_expr_and_foldchange['locus_prioritised'][loci_expr_and_foldchange['rsid']==targ] = 'Cannot prioritise locus due no gradients in quadrant 2, 1, or 4'
      cannot_prior = loci_expr_and_foldchange[loci_expr_and_foldchange['locus_prioritised']=='Cannot prioritise locus due no gradients in quadrant 2, 1, or 4']
      print(cannot_prior)
      locus_prioritisation = pd.concat([locus_prioritisation, cannot_prior])
      print(locus_prioritisation)
      #
locus_prioritisation.drop(['locus_prioritised'], axis=1, inplace=True)
#
locus_prioritisation['quadrant'] = pd.Categorical(locus_prioritisation['quadrant'], [2, 1, 4, 3])
locus_prioritisation.sort_values("quadrant", inplace=True)
#
# locus_prioritisation['mean_expr'] = locus_prioritisation['mean_expr'].round(decimals = 3)
# locus_prioritisation['fold_change'] = locus_prioritisation['fold_change'].round(decimals = 3)
# locus_prioritisation['gradient_per_gene'] = locus_prioritisation['gradient_per_gene'].round(decimals = 3)
# locus_prioritisation['diagonal_distance'] = locus_prioritisation['diagonal_distance'].round(decimals = 3)
locus_prioritisation.to_csv('/well/PROCARDIS/domwest/replicate_old/using_tads/' + cell_type + '_loci_prioritisation_sorted_by_quadrant.csv', index=False)
#
#RERUN_START>>>
loci_expr_and_foldchange = final_testing
#
#should maybe consider doing this without including genes that had expression = 0 in the mix. Since we want to get an idea of the typical range of expression among those genes that are actually important in the heart in order to determine an expression cutoff in the quadrants. We also want to include all genes in the heart, so must go back to older df long_initial_merge_lv
means_all_dfs = means_all_dfs[~(means_all_dfs['mean_expr']==0)]
median_of_mean_expr =means_all_dfs["mean_expr"].median() #0.014407334643091
#
loci_expr_and_foldchange = loci_expr_and_foldchange.sort_values(['fold_change', 'mean_expr'], ascending = [False, False])
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#using top gene metric to prioritise loci
#Sort genes within loci by BOTH fold change as well as median_expr
#the way it works is to sort by fold change and then expression
loci_expr_and_foldchange['gradient_per_gene'] = loci_expr_and_foldchange['fold_change']/loci_expr_and_foldchange['mean_expr']
#
#now determine which quadrant it is in:
loci_expr_and_foldchange['quadrant'] = ''
loci_expr_and_foldchange['quadrant'][(loci_expr_and_foldchange['fold_change']>=1.3) & (loci_expr_and_foldchange['mean_expr']<median_of_mean_expr)] = 1
loci_expr_and_foldchange['quadrant'][(loci_expr_and_foldchange['fold_change']>=1.3) & (loci_expr_and_foldchange['mean_expr']>=median_of_mean_expr)] = 2
loci_expr_and_foldchange['quadrant'][(loci_expr_and_foldchange['fold_change']<1.3) & (loci_expr_and_foldchange['mean_expr']<median_of_mean_expr)] = 3
loci_expr_and_foldchange['quadrant'][(loci_expr_and_foldchange['fold_change']<1.3) & (loci_expr_and_foldchange['mean_expr']>=median_of_mean_expr)] = 4
#adding diagonal distance (ie hypotenuse to fold change and median expr) to help sort in addition (and after) sorting by custom quadrant
loci_expr_and_foldchange['diagonal_distance'] = np.sqrt((loci_expr_and_foldchange['fold_change'] ** 2) + (loci_expr_and_foldchange['mean_expr'] ** 2))
loci_expr_and_foldchange.sort_values(by=['rsid', 'quadrant', 'diagonal_distance'], inplace=True, ascending=[True,True,False])
#RERUN_END<<<
#
gene_prioritisation = loci_expr_and_foldchange
# gene_prioritisation.sort_values(["rsid", "quadrant", "diagonal_distance"], inplace=True, ascending=[True, True, False])
sorter = [2, 1, 4, 3]
x = pd.DataFrame({'quadrant': sorter})
x.index = x.index.set_names('number')
x = x.reset_index()
gene_prioritisation = pd.merge(gene_prioritisation, x, how='left', on='quadrant')
gene_prioritisation.sort_values(['rsid', 'number', 'diagonal_distance'], ascending = [True, True, False], inplace=True)
gene_prioritisation.drop(['number'], axis=1, inplace=True)
# gene_prioritisation['mean_expr'] = gene_prioritisation['mean_expr'].round(decimals = 3)
# gene_prioritisation['fold_change'] = gene_prioritisation['fold_change'].round(decimals = 3)
# gene_prioritisation['gradient_per_gene'] = gene_prioritisation['gradient_per_gene'].round(decimals = 3)
# gene_prioritisation['diagonal_distance'] = gene_prioritisation['diagonal_distance'].round(decimals = 3)
gene_prioritisation.to_csv('/well/PROCARDIS/domwest/replicate_old/using_tads/' + cell_type + '_gene_prioritisation.csv', index=False)
