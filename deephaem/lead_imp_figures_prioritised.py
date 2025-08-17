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

dff = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deephaem_prep/post_analysis/various_strength_post_analysis_lead_and_imp_vars_damage_scores_with_alleles.csv')

##################################################################################################

dff_subset = dff[dff['strength'].isin(['95%','99.7%','empirical'])]

#just prioritised lead and imp snps:
dff_subset['value'].mean()
dff_subset['value'].std()
#Mean is 0.00025048752850254243
#Std deviation is 0.029438962970534385
#2SD = 0.05887792594106877
#3SD = 0.08831688891160315
#Mean - SD = -0.029188475442031844
#Mean + SD = 0.029689450499036926
#Mean – 2SD = -0.05862743841256623
#Mean + 2SD = 0.05912841346957131
#Mean – 3SD = -0.08806640138310061
#Mean + 3SD = 0.08856737644010569

dff_subset.to_csv('/gpfs3/well/PROCARDIS/domwest/deephaem_prep/post_analysis/prioritised_post_analysis_lead_and_imp_vars_damage_scores_with_alleles.csv', index=False)

x = pd.DataFrame(dff_subset[['rsid', 'REF', 'ALT', 'variable']].drop_duplicates().groupby(['rsid', 'REF', 'ALT']).count()).reset_index()
x.drop(['REF', 'ALT'], axis=1, inplace=True)
x.rename({'variable':'count'}, axis=1, inplace=True)
dff_subset_counts = dff_subset.merge(x, how='left', on='rsid')
dff_subset_counts.to_csv('/gpfs3/well/PROCARDIS/domwest/deephaem_prep/post_analysis/prioritised_post_analysis_lead_and_imp_vars_damage_scores_and_classifier_counts_with_alleles.csv', index=False)

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deephaem_prep/post_analysis/prioritised_post_analysis_lead_and_imp_vars_damage_scores_and_classifier_counts_with_alleles.csv')
#1210 unique snp_test_ids

value_stats <- long_prediction_df |>
summarize(mean_value = mean(value))
value_stats

# c25 <- c(
#   "dodgerblue2", "#E31A1C", # red
#   "green4",
#   "#6A3D9A", # purple
#   "#FF7F00", # orange
#   "black", "gold1",
#   "skyblue2", "#FB9A99", # lt pink
#   "palegreen2",
#   "#CAB2D6", # lt purple
#   "#FDBF6F", # lt orange
#   "gray70", "khaki2",
#   "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
#   "darkturquoise", "green1", "yellow4", "yellow3",
#   "darkorange4", "brown"
# )

# e <- ggplot(long_prediction_df, aes(x=variable, y=value)) +
# geom_jitter(aes(size = count, color = variable, alpha = 0.5), position = position_jitter(0.2)) +
# #scale_size(range = c(5, 10)) +
# stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
# scale_colour_manual(values=c25) +
# theme(axis.text.x=element_text(angle = 90, hjust = 0.5, vjust = 0.5), axis.ticks.x=element_blank(), legend.position = 'none')

# ggsave(filename="/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/ds_scatter_plots_Jun4/prioritised_lead_and_imp_stripchart_grouped_by_classifier.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

e <- ggplot(long_prediction_df, aes(x=variable, y=value)) +
geom_jitter(aes(size = count, color = variable, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
# scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom")

ggsave(filename="/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/ds_scatter_plots_Jun4/prioritised_lead_and_imp_stripchart_grouped_by_classifier.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


#first in python:
dff_subset_wide = dff_subset.pivot(index=['rsid', 'chr', 'grch38_POS', 'REF', 'ALT'], columns='variable', values='value').reset_index()
dff_subset_wide.rename({'rsid':'variant'}, axis=1, inplace=True)
dff = pd.read_excel('/well/PROCARDIS/domwest/deephaem_prep/gathering_snps/HCM_MTAG_OT_FUMA.xlsx')
dff_filt = dff[['rsid']]
dff_filt['type'] = 'lead'
#df.drop(['type'], axis=1, inplace=True)
dff_filt.rename({'rsid':'variant'}, axis=1, inplace=True)
new_df = dff_subset_wide.merge(dff_filt, how='left', on='variant')
extras = ['rs67491807', 'rs764462761', 'rs3218719', 'rs11570041', 'rs117847273']
new_df['type'][new_df['variant'].isin(extras)] = 'lead'
len(new_df['variant'][new_df['type']=='lead'].unique()) #23
new_df['type'][new_df['type'].isnull()] = 'proxy'
new_df[['variant', 'REF', 'ALT']].drop_duplicates() #1568
new_df.to_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_shortlist_vars_prediction/with_vars_present_prioritised_1568_lead_imp_vars_damage_scores_type_added.csv', index=False)

#just prioritised lead and imp snps:
dff_subset['value'].mean()
dff_subset['value'].std()
#Mean is 0.00025048752850254243
#Std deviation is 0.029438962970534385
#2SD = 0.05887792594106877
#3SD = 0.08831688891160315
#Mean - SD = -0.029188475442031844
#Mean + SD = 0.029689450499036926
#Mean – 2SD = -0.05862743841256623
#Mean + 2SD = 0.05912841346957131
#Mean – 3SD = -0.08806640138310061
#Mean + 3SD = 0.08856737644010569

make_scatter_plots <- function(cell_type) {
    #
    prediction_df <- read.csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_shortlist_vars_prediction/with_vars_present_prioritised_1568_lead_imp_vars_damage_scores_type_added.csv') #was: with_vars_present_pos_controls_and_4395_lead_imp_vars_damage_scores
    #
    prediction_df_filt <- prediction_df[, c("variant", "type", cell_type)]
    #
    prediction_df_filt <- prediction_df_filt[order(prediction_df_filt[[cell_type]], decreasing=TRUE),]
    #
    basic <- 
    ggplot(prediction_df_filt, aes_string(x='variant', y=cell_type, colour="type")) + 
    geom_point() +
    geom_text(data=subset(prediction_df_filt, prediction_df_filt[[cell_type]] >= 0.1 | prediction_df_filt[[cell_type]] <= -0.1), aes(label=variant), position = position_dodge(width = 1), vjust = -0.5) +
    #geom_point(data=prediction_df_filt[prediction_df_filt[[cell_type]]<=0.1,], pch=21, fill=NA, size=4, colour="black", stroke=1) +
    geom_hline(yintercept=0.1, color="darkmagenta", size =0.5) +
    geom_hline(yintercept=-0.1, color="darkmagenta", size =0.5) +
    geom_hline(yintercept=0, color="black", size =0.3, linetype="dashed") +
    geom_hline(yintercept=-0.029188475442031844, color="cornflowerblue", size =0.5) + #1SD
    geom_hline(yintercept=0.029689450499036926, color="cornflowerblue", size =0.5) + #1SD
    geom_hline(yintercept=-0.05862743841256623, color="chartreuse4", size =0.5) + #2SD
    geom_hline(yintercept=0.05912841346957131, color="chartreuse4", size =0.5) + #2SD
    geom_hline(yintercept=-0.08806640138310061, color="brown2", size =0.5) + #3SD
    geom_hline(yintercept=0.08856737644010569, color="brown2", size =0.5) + #3SD
    labs(title = paste(cell_type," classifier's damage scores"),
        subtitle = NULL,
        tag = NULL, 
        x = "Variants",
        y= "Damage score",
        color = NULL) +
    theme(
        plot.title=element_text(size=18),
        axis.text.x = element_blank(), #axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 16))
    basic
    print(basic) #in rstudio will show it
    #
    ggsave(filename=paste("/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/ds_scatter_plots_Jun4/1568_prioritised_lead_imp_",cell_type,".tiff",sep=""), width=18, height=15, dpi = 300) } #was: scatter_plots/all/

make_scatter_plots("fibroblast_ATAC")
make_scatter_plots("Endocardial_ATAC")
make_scatter_plots("Pericyte_General_4_ATAC")
make_scatter_plots("Fetal_V_Cardiomyocyte_ATAC")
make_scatter_plots("Fetal_Mesothelial_ATAC")
make_scatter_plots("bulk_heart_LV_tissue_ATAC2")
make_scatter_plots("cardiac_muscle_cell_DNASE2")
make_scatter_plots("Mesothelial_ATAC")
make_scatter_plots("Pericyte_ATAC")
make_scatter_plots("Fetal_Skeletal_Myocyte_1_ATAC")
make_scatter_plots("Cardiac_Pericyte_2_ATAC")
make_scatter_plots("Endothelial_Myocardial_ATAC")
make_scatter_plots("SMC2_ATAC")
make_scatter_plots("smooth_muscle_ATAC")
make_scatter_plots("atrial_cm_ATAC")
make_scatter_plots("Fetal_Cardiac_Fibroblast_ATAC")
make_scatter_plots("Type_I_Skeletal_Myocyte_ATAC")
make_scatter_plots("SMC1_ATAC")
make_scatter_plots("bulk_heart_LV_tissue_DNASE2")
make_scatter_plots("Vasc_Sm_Muscle_1_ATAC")
make_scatter_plots("Pericyte_General_2_ATAC")
make_scatter_plots("Vasc_Sm_Muscle_2_ATAC")
make_scatter_plots("Fetal_Skeletal_Myocyte_3_ATAC")
make_scatter_plots("Ery_Don002_hg38_ATAC")
make_scatter_plots("Cardiac_Pericyte_1_ATAC")
make_scatter_plots("Adipocyte_ATAC")
make_scatter_plots("adipocyte_ATAC")
make_scatter_plots("ventricular_cm_ATAC")
make_scatter_plots("Ery_Don003_hg38_ATAC")
make_scatter_plots("Type_II_Skeletal_Myocyte_ATAC")
make_scatter_plots("Mast.1_ATAC")
make_scatter_plots("Endothelial_ATAC")
make_scatter_plots("endothelial_ATAC")
make_scatter_plots("A_Cardiomyocyte_ATAC")
make_scatter_plots("Cardiac_Fibroblast_ATAC")
make_scatter_plots("Ery_Don001_hg38_ATAC")
make_scatter_plots("lymphocyte_ATAC")
make_scatter_plots("Pericyte_General_3_ATAC")
make_scatter_plots("cardiac_muscle_cell_DNASE1")
make_scatter_plots("Cardiac_Pericyte_4_ATAC")
make_scatter_plots("cardiac_muscle_cell_CHIP_HISTONE_H3K4me3_single")
make_scatter_plots("cardiac_muscle_cell_CHIP_HISTONE_H3K27me3_single")
make_scatter_plots("cardiac_muscle_cell_DNASE_single")
make_scatter_plots("cardiac_muscle_cell_CHIP_HISTONE_H3K27ac_single")
make_scatter_plots("cardiac_muscle_cell_CHIP_TF_single")
make_scatter_plots("bulk_heart_LV_tissue_ATAC_paired")
make_scatter_plots("bulk_heart_LV_tissue_DNASE_paired")
make_scatter_plots("bulk_heart_LV_tissue_CHIP_TF_paired")
make_scatter_plots("bulk_heart_LV_tissue_CHIP_HISTONE_H3K4me3_single")
make_scatter_plots("bulk_heart_LV_tissue_CHIP_HISTONE_H3K27me3_single")
make_scatter_plots("bulk_heart_LV_tissue_CHIP_HISTONE_H3K4me1_single")
make_scatter_plots("bulk_heart_LV_tissue_CHIP_HISTONE_H3K27ac_single")
make_scatter_plots("bulk_heart_LV_tissue_CHIP_TF_single")
make_scatter_plots("cardiac_muscle_cell_DNASE_paired")
make_scatter_plots("Pericyte_General_1_ATAC")
make_scatter_plots("Cardiac_Pericyte_3_ATAC")
make_scatter_plots("Fetal_Endocardial_ATAC")
make_scatter_plots("Mast_ATAC")
make_scatter_plots("Macrophage_ATAC")
make_scatter_plots("macrophage_ATAC")
make_scatter_plots("V_Cardiomyocyte_ATAC")
make_scatter_plots("Fetal_A_Cardiomyocyte_ATAC")
make_scatter_plots("Pericyte_Muscularis_ATAC")
make_scatter_plots("Fetal_Skeletal_Myocyte_2_ATAC")
make_scatter_plots("Fibroblast_ATAC")