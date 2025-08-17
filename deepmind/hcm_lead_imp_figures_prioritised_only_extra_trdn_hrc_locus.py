import pandas as pd
import numpy as np

classes = ['DNASE_cardiac_muscle_cell_2', 'DNASE_heart_left_ventricle_female_adult_(53_years)', 'DNASE_cardiac_muscle_cell_1', 'DNASE_heart_left_ventricle_female_embryo_(101_day)_and_female_embryo_(103_days)', 'DNASE_heart_left_ventricle_female_embryo_(136_days)']

overarch = []
for clas in classes:
    tmp = pd.read_csv('/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/hrc_matrices_1/' + clas +'_total_hcm.csv')
    overarch.append(tmp)
overarch = pd.concat(overarch)

overarch.to_csv('/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/hrc_matrices_1/all_classifiers_total_hcm.csv', index=False)

###

classes = ['DNASE_cardiac_muscle_cell_2', 'DNASE_heart_left_ventricle_female_adult_(53_years)', 'DNASE_cardiac_muscle_cell_1', 'DNASE_heart_left_ventricle_female_embryo_(101_day)_and_female_embryo_(103_days)', 'DNASE_heart_left_ventricle_female_embryo_(136_days)']

overarch = []
for clas in classes:
    tmp = pd.read_csv('/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/hrc_matrices_2/' + clas +'_total_hcm.csv')
    overarch.append(tmp)
overarch = pd.concat(overarch)

overarch.to_csv('/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/hrc_matrices_2/all_classifiers_total_hcm.csv', index=False)

#########################################################################################################################################################

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/hrc_matrices_1/all_classifiers_total_hcm.csv')
#1210 unique snp_test_ids

# value_stats <- long_prediction_df |>
# summarize(mean_value = mean(value))
# value_stats

c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

# e <- ggplot(long_prediction_df, aes(x=variable, y=value)) +
# geom_jitter(aes(size = count, color = variable, alpha = 0.5), position = position_jitter(0.2)) +
# #scale_size(range = c(5, 10)) +
# stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
# scale_colour_manual(values=c25) +
# theme(axis.text.x=element_text(angle = 90, hjust = 0.5, vjust = 0.5), axis.ticks.x=element_blank(), legend.position = 'none')

# ggsave(filename="/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/ds_scatter_plots_Jun4/prioritised_lead_and_imp_stripchart_grouped_by_classifier.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

e <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value)) +
geom_jitter(aes(color = classifier), position = position_jitter(0.2)) + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
geom_hline(yintercept=0.294990062713623, color="darkmagenta", size =0.5) +
geom_hline(yintercept=-0.2931210994720459, color="darkmagenta", size =0.5) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/vars_all_classifiers_hcm_hrc_1_lead_and_imp_damage_score_stripchart.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

#########################################################################################################################################################

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/hrc_matrices_2/all_classifiers_total_hcm.csv')
#1210 unique snp_test_ids

# value_stats <- long_prediction_df |>
# summarize(mean_value = mean(value))
# value_stats

c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

# e <- ggplot(long_prediction_df, aes(x=variable, y=value)) +
# geom_jitter(aes(size = count, color = variable, alpha = 0.5), position = position_jitter(0.2)) +
# #scale_size(range = c(5, 10)) +
# stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
# scale_colour_manual(values=c25) +
# theme(axis.text.x=element_text(angle = 90, hjust = 0.5, vjust = 0.5), axis.ticks.x=element_blank(), legend.position = 'none')

# ggsave(filename="/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/ds_scatter_plots_Jun4/prioritised_lead_and_imp_stripchart_grouped_by_classifier.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

e <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value)) +
geom_jitter(aes(color = classifier), position = position_jitter(0.2)) + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
geom_hline(yintercept=0.294990062713623, color="darkmagenta", size =0.5) +
geom_hline(yintercept=-0.2931210994720459, color="darkmagenta", size =0.5) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/extra_trdn_hrc_locus/vars_all_classifiers_hcm_hrc_2_lead_and_imp_damage_score_stripchart.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/
