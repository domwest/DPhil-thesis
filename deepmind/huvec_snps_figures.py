import pandas as pd
import numpy as np

#OG ie endothelial classifiers

overarch = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/total_huvec.csv')

for clas in overarch['classifier'].unique():
    classif = clas.replace(':','_')
    tmp = overarch[overarch['classifier']==clas]
    tmp.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/' + classif + '_huvec.csv', index=False)

#############################################################################################

.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

# -rw-r--r-- 1 wiq135 watkins  71286 Aug  9 10:06 DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec.csv
# -rw-r--r-- 1 wiq135 watkins  71888 Aug  9 10:06 DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec.csv


long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec.csv')

long_prediction_df$var_type=factor(long_prediction_df$var_type,levels=c("gain","neutral","loss"))

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

e <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value)) +
geom_jitter(aes(color = var_type), position = position_jitter(0.2)) + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right") +
labs(x = "DNASE:endothelial_cell_of_umbilical_vein_newborn_1; index=29")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_damage_score_stripchart.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value, fill=var_type)) +
geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
scale_fill_manual(values=c25) +
scale_y_continuous(n.breaks = 10) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") +
labs(x = "DNASE:endothelial_cell_of_umbilical_vein_newborn_1; index=29") 

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_damage_score_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


###

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec.csv')

long_prediction_df$var_type=factor(long_prediction_df$var_type,levels=c("gain","neutral","loss"))

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

e <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value)) +
geom_jitter(aes(color = var_type), position = position_jitter(0.2)) + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right") +
labs(x = "DNASE:endothelial_cell_of_umbilical_vein_newborn_2; index=118")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_damage_score_stripchart.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value, fill=var_type)) +
geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
scale_fill_manual(values=c25) +
scale_y_continuous(n.breaks = 10) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") +
labs(x = "DNASE:endothelial_cell_of_umbilical_vein_newborn_2; index=118")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_damage_score_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

#############################################################################################

import pandas as pd
import numpy as np

gains = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/HUVEC.gains.tsv', sep='\t')
gains['var_type'] = 'gain'
losses = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/HUVEC.losses.tsv', sep='\t')
losses['var_type'] = 'loss'
neutral = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/HUVEC.neutral.tsv', sep='\t')
neutral['var_type'] = 'neutral'
skew_df = pd.concat([gains, losses, neutral])

class1 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec.csv')
skew_df.rename({'TEST.SNP.ID':'rsid'}, axis=1, inplace=True)
class1_new = class1.merge(skew_df[['rsid', 'SKEW.ATAC']], how='left', on='rsid')
class1_new.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec_with_skew.csv', index=False)

class2 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec.csv')
# skew_df.rename({'TEST.SNP.ID':'rsid'}, axis=1, inplace=True)
class2_new = class2.merge(skew_df[['rsid', 'SKEW.ATAC']], how='left', on='rsid')
class2_new.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec_with_skew.csv', index=False)

#############################################################################################


.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec_with_skew.csv')

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

e <- ggplot(long_prediction_df, aes(x=SKEW.ATAC, y=damage_score_value, shape=var_type, color=var_type)) +
geom_point() + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
#stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_damage_score_skew_scatterplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

#

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec_with_skew.csv')

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

e <- ggplot(long_prediction_df, aes(x=SKEW.ATAC, y=damage_score_value, shape=var_type, color=var_type)) +
geom_point() + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
#stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_damage_score_skew_scatterplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/



#############################################################################################

import pandas as pd
import numpy as np

gains = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/HUVEC.gains.tsv', sep='\t')
gains['var_type'] = 'gain'
losses = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/HUVEC.losses.tsv', sep='\t')
losses['var_type'] = 'loss'
neutral = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/HUVEC.neutral.tsv', sep='\t')
neutral['var_type'] = 'neutral'
skew_df = pd.concat([gains, losses, neutral])

class1 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec.csv')
skew_df.rename({'TEST.SNP.ID':'rsid'}, axis=1, inplace=True)
class1_new = class1.merge(skew_df[['rsid', 'SKEW.ATAC']], how='left', on='rsid')

class2 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec.csv')
# skew_df.rename({'TEST.SNP.ID':'rsid'}, axis=1, inplace=True)
class2_new = class2.merge(skew_df[['rsid', 'SKEW.ATAC']], how='left', on='rsid')

class_tot = pd.concat([class1_new, class2_new])

m_categories = ['gain', 'neutral', 'loss']
class_tot["var_type"] = pd.Categorical(class_tot["var_type"], categories = m_categories)
class_tot = class_tot.sort_values(by = "var_type")
class_tot.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_huvec_with_skew.csv', index=False)

#********************FILTER FOR NEUTRAL!!!!!!!!

class_tot_neutral = class_tot[class_tot['var_type']=='neutral']
class_tot_neutral.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_huvec_neutral_with_skew.csv', index=False)

#############################################################################################


.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_huvec_with_skew.csv')

p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value, fill=var_type)) +
geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
# scale_fill_manual(values=c26) +
scale_y_continuous(n.breaks = 10) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") 

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_damage_score_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


#

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_huvec_with_skew.csv')


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

e <- ggplot(long_prediction_df, aes(x=SKEW.ATAC, y=damage_score_value, shape=var_type, color=classifier)) +
geom_point() + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
#stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_damage_score_skew_scatterplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

################################################################################################################################################################

.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_huvec_neutral_with_skew.csv')

value_stats <- long_prediction_df |>
summarize(mean_value = mean(damage_score_value))
value_stats

# #probability density plot
# ggplot(long_prediction_df, aes(x = damage_score_value, y = after_stat(density))) +
#     geom_histogram(bins = 200) +
#     geom_vline(aes(xintercept = mean_value), value_stats, color = "darkolivegreen1", linewidth = 2) +
#     geom_density(color = "darkcyan", linewidth = 2) +
#     scale_x_continuous(n.breaks = 10)
# #default is binwidth 30

# ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_damage_score_prob_density_plot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


long_prediction_df <- long_prediction_df |>
    mutate(variable = factor(classifier))

ggplot(data = long_prediction_df, aes(x = damage_score_value, fill = classifier)) +
    geom_histogram(bins = 200) +
    scale_x_continuous(n.breaks = 10)

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_huvec_neutral_damage_score_histogram.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

##

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec_only_more_neutrals.csv')

value_stats <- long_prediction_df |>
summarize(mean_value = mean(damage_score_value))
value_stats

long_prediction_df <- long_prediction_df |>
    mutate(variable = factor(classifier))

ggplot(data = long_prediction_df, aes(x = damage_score_value)) +
    geom_histogram(bins = 200, fill = "darkturquoise", color = "white") +
    geom_vline(xintercept=0.1793844103813171, color = "darkslateblue") +
    geom_vline(xintercept=-0.2028707265853881, color = "darkslateblue") +
    scale_x_continuous(n.breaks = 10) +
    ggtitle("Extended HUVEC neutral set distribution: DNASE:endothelial_cell_of_umbilical_vein_newborn_1; index=29")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_more_neutrals_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec_damage_score_histogram.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

##

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec_only_more_neutrals.csv')

value_stats <- long_prediction_df |>
summarize(mean_value = mean(damage_score_value))
value_stats

long_prediction_df <- long_prediction_df |>
    mutate(variable = factor(classifier))

ggplot(data = long_prediction_df, aes(x = damage_score_value)) +
    geom_histogram(bins = 200, fill = "darkturquoise", color = "white") +
    geom_vline(xintercept=0.294990062713623, color = "darkslateblue") +
    geom_vline(xintercept=-0.2931210994720459, color = "darkslateblue") +
    scale_x_continuous(n.breaks = 10) +
    ggtitle("Extended HUVEC neutral set distribution: DNASE:endothelial_cell_of_umbilical_vein_newborn_2; index=118")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_more_neutrals_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec_damage_score_histogram.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

##

# long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec_with_skew.csv')

# value_stats <- long_prediction_df |>
# summarize(mean_value = mean(damage_score_value))
# value_stats

# # #probability density plot
# # ggplot(long_prediction_df, aes(x = damage_score_value, y = after_stat(density))) +
# #     geom_histogram(bins = 200) +
# #     geom_vline(aes(xintercept = mean_value), value_stats, color = "darkolivegreen1", linewidth = 2) +
# #     geom_density(color = "darkcyan", linewidth = 2) +
# #     scale_x_continuous(n.breaks = 10)
# # #default is binwidth 30

# # ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_damage_score_prob_density_plot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


# long_prediction_df <- long_prediction_df |>
#     mutate(variable = factor(classifier))

# ggplot(data = long_prediction_df, aes(x = damage_score_value, fill = classifier)) +
#     geom_histogram(bins = 200) +
#     scale_x_continuous(n.breaks = 10)

# ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_damage_score_histogram.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


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

import pandas as pd
import numpy as np

#new neutral set of huvec snps (more neutral snps than previously)

overarch = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/total_huvec_only_more_neutrals.csv')

for clas in overarch['classifier'].unique():
    classif = clas.replace(':','_')
    tmp = overarch[overarch['classifier']==clas]
    tmp.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/' + classif + '_huvec_only_more_neutrals.csv', index=False)

#############################################################################################

.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

# -rw-r--r-- 1 wiq135 watkins  71286 Aug  9 10:06 DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec.csv
# -rw-r--r-- 1 wiq135 watkins  71888 Aug  9 10:06 DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec.csv


long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec_only_more_neutrals.csv')

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

e <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value)) +
geom_jitter(aes(color = var_type), position = position_jitter(0.2)) + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_only_more_neutrals_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_damage_score_stripchart.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value, fill=var_type)) +
geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
scale_fill_manual(values=c25) +
scale_y_continuous(n.breaks = 10) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") 

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_only_more_neutrals_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_damage_score_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


###

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec_only_more_neutrals.csv')

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

e <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value)) +
geom_jitter(aes(color = var_type), position = position_jitter(0.2)) + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_only_more_neutrals_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_damage_score_stripchart.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value, fill=var_type)) +
geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
scale_fill_manual(values=c25) +
scale_y_continuous(n.breaks = 10) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") 

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_only_more_neutrals_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_damage_score_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

#############################################################################################

import pandas as pd
import numpy as np

neutral = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/Simone.HUVEC.neutral.tsv', sep='\t')
neutral['var_type'] = 'extended_neutrals'
skew_df = neutral

class1 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec_only_more_neutrals.csv')
skew_df.rename({'TEST.SNP.ID':'rsid'}, axis=1, inplace=True)
class1_new = class1.merge(skew_df[['rsid', 'SKEW.ATAC']], how='left', on='rsid')
class1_new.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec_only_more_neutrals_with_skew.csv', index=False)

class2 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec_only_more_neutrals.csv')
# skew_df.rename({'TEST.SNP.ID':'rsid'}, axis=1, inplace=True)
class2_new = class2.merge(skew_df[['rsid', 'SKEW.ATAC']], how='left', on='rsid')
class2_new.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec_only_more_neutrals_with_skew.csv', index=False)

#############################################################################################


.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec_only_more_neutrals_with_skew.csv')

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

e <- ggplot(long_prediction_df, aes(x=SKEW.ATAC, y=damage_score_value, shape=var_type, color=var_type)) +
geom_point() + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
#stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_only_more_neutrals_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_damage_score_skew_scatterplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

#

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec_only_more_neutrals_with_skew.csv')

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

e <- ggplot(long_prediction_df, aes(x=SKEW.ATAC, y=damage_score_value, shape=var_type, color=var_type)) +
geom_point() + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
#stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_only_more_neutrals_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_damage_score_skew_scatterplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/



#############################################################################################

import pandas as pd
import numpy as np

neutral = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/Simone.HUVEC.neutral.tsv', sep='\t')
neutral['var_type'] = 'extended_neutrals'
skew_df = neutral

class1 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec_only_more_neutrals.csv')
skew_df.rename({'TEST.SNP.ID':'rsid'}, axis=1, inplace=True)
class1_new = class1.merge(skew_df[['rsid', 'SKEW.ATAC']], how='left', on='rsid')

class2 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec_only_more_neutrals.csv')
# skew_df.rename({'TEST.SNP.ID':'rsid'}, axis=1, inplace=True)
class2_new = class2.merge(skew_df[['rsid', 'SKEW.ATAC']], how='left', on='rsid')

class_tot = pd.concat([class1_new, class2_new])

class_tot.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_huvec_only_more_neutrals_with_skew.csv', index=False)

#############################################################################################


.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_huvec_only_more_neutrals_with_skew.csv')

p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value, fill=var_type)) +
geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
# scale_fill_manual(values=c26) +
scale_y_continuous(n.breaks = 10) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") 

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_only_more_neutrals_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_damage_score_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


#

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_huvec_only_more_neutrals_with_skew.csv')


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

e <- ggplot(long_prediction_df, aes(x=SKEW.ATAC, y=damage_score_value, shape=var_type, color=classifier)) +
geom_point() + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
#stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_only_more_neutrals_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_damage_score_skew_scatterplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

################################################################################################################################################################

.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_huvec_only_more_neutrals_with_skew.csv')

value_stats <- long_prediction_df |>
summarize(mean_value = mean(damage_score_value))
value_stats

# #probability density plot
# ggplot(long_prediction_df, aes(x = damage_score_value, y = after_stat(density))) +
#     geom_histogram(bins = 200) +
#     geom_vline(aes(xintercept = mean_value), value_stats, color = "darkolivegreen1", linewidth = 2) +
#     geom_density(color = "darkcyan", linewidth = 2) +
#     scale_x_continuous(n.breaks = 10)
# #default is binwidth 30

# ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_damage_score_prob_density_plot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


long_prediction_df <- long_prediction_df |>
    mutate(variable = factor(classifier))

ggplot(data = long_prediction_df, aes(x = damage_score_value, fill = classifier)) +
    geom_histogram(bins = 200) +
    scale_x_continuous(n.breaks = 10)

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_huvec_only_more_neutrals_damage_score_histogram.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

##


# long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec_with_skew.csv')

# value_stats <- long_prediction_df |>
# summarize(mean_value = mean(damage_score_value))
# value_stats

# # #probability density plot
# # ggplot(long_prediction_df, aes(x = damage_score_value, y = after_stat(density))) +
# #     geom_histogram(bins = 200) +
# #     geom_vline(aes(xintercept = mean_value), value_stats, color = "darkolivegreen1", linewidth = 2) +
# #     geom_density(color = "darkcyan", linewidth = 2) +
# #     scale_x_continuous(n.breaks = 10)
# # #default is binwidth 30

# # ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_damage_score_prob_density_plot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


# long_prediction_df <- long_prediction_df |>
#     mutate(variable = factor(classifier))

# ggplot(data = long_prediction_df, aes(x = damage_score_value, fill = classifier)) +
#     geom_histogram(bins = 200) +
#     scale_x_continuous(n.breaks = 10)

# ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_damage_score_histogram.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


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

import pandas as pd
import numpy as np

#OG ie endothelial classifiers but only 'prioritised HUVEC snps' based on more neutrals used to define thresholds at 0.05 p value

#############################################################################################

import pandas as pd
import numpy as np

class1_1 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec_less_than_loss_threshold_with_skew.csv')
class1_2 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec_greater_than_gain_threshold_with_skew.csv')
class1_new = pd.concat([class1_1, class1_2])

class2_1 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec_less_than_loss_threshold_with_skew.csv')
class2_2 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec_greater_than_gain_threshold_with_skew.csv')
class2_new = pd.concat([class2_1, class2_2])

class_tot = pd.concat([class1_new, class2_new])

m_categories = ['gain', 'neutral', 'loss']
class_tot["var_type"] = pd.Categorical(class_tot["var_type"], categories = m_categories)
class_tot = class_tot.sort_values(by = "var_type")
class_tot.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_huvec_prioritised_with_skew.csv', index=False)

#############################################################################################


.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_huvec_prioritised_with_skew.csv')

p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value, fill=var_type)) +
geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
# scale_fill_manual(values=c26) +
scale_y_continuous(n.breaks = 10) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") 

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_prioritised_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_damage_score_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


#

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_huvec_prioritised_with_skew.csv')


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

e <- ggplot(long_prediction_df, aes(x=SKEW.ATAC, y=damage_score_value, shape=var_type, color=classifier)) +
geom_point() + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
#stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_prioritised_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_damage_score_skew_scatterplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

################################################################################################################################################################

.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec_with_skew.csv')

value_stats <- long_prediction_df |>
summarize(mean_value = mean(damage_score_value))
value_stats

long_prediction_df <- long_prediction_df |>
    mutate(variable = factor(classifier))

ggplot(data = long_prediction_df, aes(x = damage_score_value)) +
    geom_histogram(bins = 200, fill = "darkturquoise", color = "white") +
    geom_vline(xintercept=0.1793844103813171, color = "darkslateblue") +
    geom_vline(xintercept=-0.2028707265853881, color = "darkslateblue") +
    scale_x_continuous(n.breaks = 10) +
    ggtitle("DNASE:endothelial_cell_of_umbilical_vein_newborn_1; index=29")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_prioritised_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec_damage_score_histogram.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

##

.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec_with_skew.csv')

value_stats <- long_prediction_df |>
summarize(mean_value = mean(damage_score_value))
value_stats

long_prediction_df <- long_prediction_df |>
    mutate(variable = factor(classifier))

ggplot(data = long_prediction_df, aes(x = damage_score_value)) +
    geom_histogram(bins = 200, fill = "darkturquoise", color = "white") +
    geom_vline(xintercept=0.294990062713623, color = "darkslateblue") +
    geom_vline(xintercept=-0.2931210994720459, color = "darkslateblue") +
    scale_x_continuous(n.breaks = 10) +
    ggtitle("DNASE:endothelial_cell_of_umbilical_vein_newborn_2; index=118")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_prioritised_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec_damage_score_histogram.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


##


# long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec_with_skew.csv')

# value_stats <- long_prediction_df |>
# summarize(mean_value = mean(damage_score_value))
# value_stats

# # #probability density plot
# # ggplot(long_prediction_df, aes(x = damage_score_value, y = after_stat(density))) +
# #     geom_histogram(bins = 200) +
# #     geom_vline(aes(xintercept = mean_value), value_stats, color = "darkolivegreen1", linewidth = 2) +
# #     geom_density(color = "darkcyan", linewidth = 2) +
# #     scale_x_continuous(n.breaks = 10)
# # #default is binwidth 30

# # ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_damage_score_prob_density_plot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


# long_prediction_df <- long_prediction_df |>
#     mutate(variable = factor(classifier))

# ggplot(data = long_prediction_df, aes(x = damage_score_value, fill = classifier)) +
#     geom_histogram(bins = 200) +
#     scale_x_continuous(n.breaks = 10)

# ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_damage_score_histogram.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


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


import pandas as pd
import numpy as np

overarch = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/total_huvec.csv')

for clas in overarch['classifier'].unique():
    classif = clas.replace(':','_')
    tmp = overarch[overarch['classifier']==clas]
    tmp.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/' + classif + '_huvec.csv', index=False)

#############################################################################################

.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

# -rw-r--r-- 1 wiq135 watkins  71286 Aug  9 10:06 DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec.csv
# -rw-r--r-- 1 wiq135 watkins  71888 Aug  9 10:06 DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec.csv


long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec.csv')

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

e <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value)) +
geom_jitter(aes(color = var_type), position = position_jitter(0.2)) + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_damage_score_stripchart.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value, fill=var_type)) +
geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
scale_fill_manual(values=c25) +
scale_y_continuous(n.breaks = 10) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") 

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_damage_score_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


###

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec.csv')

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

e <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value)) +
geom_jitter(aes(color = var_type), position = position_jitter(0.2)) + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_damage_score_stripchart.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value, fill=var_type)) +
geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
scale_fill_manual(values=c25) +
scale_y_continuous(n.breaks = 10) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") 

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_damage_score_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

#############################################################################################

import pandas as pd
import numpy as np

gains = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/HUVEC.gains.tsv', sep='\t')
gains['var_type'] = 'gain'
losses = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/HUVEC.losses.tsv', sep='\t')
losses['var_type'] = 'loss'
neutral = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/HUVEC.neutral.tsv', sep='\t')
neutral['var_type'] = 'neutral'
skew_df = pd.concat([gains, losses, neutral])

class1 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec.csv')
skew_df.rename({'TEST.SNP.ID':'rsid'}, axis=1, inplace=True)
class1_new = class1.merge(skew_df[['rsid', 'SKEW.ATAC']], how='left', on='rsid')
class1_new.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec_with_skew.csv', index=False)

class2 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec.csv')
# skew_df.rename({'TEST.SNP.ID':'rsid'}, axis=1, inplace=True)
class2_new = class2.merge(skew_df[['rsid', 'SKEW.ATAC']], how='left', on='rsid')
class2_new.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec_with_skew.csv', index=False)

#############################################################################################


.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec_with_skew.csv')

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

e <- ggplot(long_prediction_df, aes(x=SKEW.ATAC, y=damage_score_value, shape=var_type, color=var_type)) +
geom_point() + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
#stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_damage_score_skew_scatterplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

#

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec_with_skew.csv')

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

e <- ggplot(long_prediction_df, aes(x=SKEW.ATAC, y=damage_score_value, shape=var_type, color=var_type)) +
geom_point() + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
#stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_damage_score_skew_scatterplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/



#############################################################################################

import pandas as pd
import numpy as np

gains = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/HUVEC.gains.tsv', sep='\t')
gains['var_type'] = 'gain'
losses = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/HUVEC.losses.tsv', sep='\t')
losses['var_type'] = 'loss'
neutral = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/HUVEC.neutral.tsv', sep='\t')
neutral['var_type'] = 'neutral'
skew_df = pd.concat([gains, losses, neutral])

class1 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec.csv')
skew_df.rename({'TEST.SNP.ID':'rsid'}, axis=1, inplace=True)
class1_new = class1.merge(skew_df[['rsid', 'SKEW.ATAC']], how='left', on='rsid')

class2 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec.csv')
# skew_df.rename({'TEST.SNP.ID':'rsid'}, axis=1, inplace=True)
class2_new = class2.merge(skew_df[['rsid', 'SKEW.ATAC']], how='left', on='rsid')

class_tot = pd.concat([class1_new, class2_new])

m_categories = ['gain', 'neutral', 'loss']
class_tot["var_type"] = pd.Categorical(class_tot["var_type"], categories = m_categories)
class_tot = class_tot.sort_values(by = "var_type")
class_tot.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_huvec_with_skew.csv', index=False)

#********************FILTER FOR NEUTRAL!!!!!!!!

class_tot_neutral = class_tot[class_tot['var_type']=='neutral']
class_tot_neutral.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_huvec_neutral_with_skew.csv', index=False)

#############################################################################################


.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_huvec_with_skew.csv')

p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value, fill=var_type)) +
geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
# scale_fill_manual(values=c26) +
scale_y_continuous(n.breaks = 10) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") 

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_damage_score_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


#

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_huvec_with_skew.csv')


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

e <- ggplot(long_prediction_df, aes(x=SKEW.ATAC, y=damage_score_value, shape=var_type, color=classifier)) +
geom_point() + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
#stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_damage_score_skew_scatterplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

################################################################################################################################################################

.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_huvec_neutral_with_skew.csv')

value_stats <- long_prediction_df |>
summarize(mean_value = mean(damage_score_value))
value_stats

# #probability density plot
# ggplot(long_prediction_df, aes(x = damage_score_value, y = after_stat(density))) +
#     geom_histogram(bins = 200) +
#     geom_vline(aes(xintercept = mean_value), value_stats, color = "darkolivegreen1", linewidth = 2) +
#     geom_density(color = "darkcyan", linewidth = 2) +
#     scale_x_continuous(n.breaks = 10)
# #default is binwidth 30

# ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_damage_score_prob_density_plot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


long_prediction_df <- long_prediction_df |>
    mutate(variable = factor(classifier))

ggplot(data = long_prediction_df, aes(x = damage_score_value, fill = classifier)) +
    geom_histogram(bins = 200) +
    scale_x_continuous(n.breaks = 10)

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_huvec_neutral_damage_score_histogram.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

##


# long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec_with_skew.csv')

# value_stats <- long_prediction_df |>
# summarize(mean_value = mean(damage_score_value))
# value_stats

# # #probability density plot
# # ggplot(long_prediction_df, aes(x = damage_score_value, y = after_stat(density))) +
# #     geom_histogram(bins = 200) +
# #     geom_vline(aes(xintercept = mean_value), value_stats, color = "darkolivegreen1", linewidth = 2) +
# #     geom_density(color = "darkcyan", linewidth = 2) +
# #     scale_x_continuous(n.breaks = 10)
# # #default is binwidth 30

# # ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_damage_score_prob_density_plot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


# long_prediction_df <- long_prediction_df |>
#     mutate(variable = factor(classifier))

# ggplot(data = long_prediction_df, aes(x = damage_score_value, fill = classifier)) +
#     geom_histogram(bins = 200) +
#     scale_x_continuous(n.breaks = 10)

# ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_damage_score_histogram.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


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


import pandas as pd
import numpy as np

#HUVEC gain, loss, and neutral (not extended) sets but predictions done in other celltypes randomly chosen -- 4 classifiers other than endothelial

overarch = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/total_huvec_other_celltypes.csv')

for clas in overarch['classifier'].unique():
    classif = clas.replace(':','_')
    tmp = overarch[overarch['classifier']==clas]
    tmp.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/' + classif + '_huvec.csv', index=False)

#############################################################################################

.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

# -rw-r--r-- 1 wiq135 watkins  71286 Aug  9 10:06 DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec.csv
# -rw-r--r-- 1 wiq135 watkins  71888 Aug  9 10:06 DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec.csv


# long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_frontal_cortex_male_adult_(27_years)_and_male_adult_(35_years)_huvec.csv')

# c25 <- c(
#   "gold1",
#   "skyblue2", "#FB9A99", # lt pink
#   "dodgerblue2", "#E31A1C", # red
#   "green4",
#   "#6A3D9A", # purple
#   "#FF7F00", # orange
#   "black", 
#   "palegreen2",
#   "#CAB2D6", # lt purple
#   "#FDBF6F", # lt orange
#   "gray70", "khaki2",
#   "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
#   "darkturquoise", "green1", "yellow4", "yellow3",
#   "darkorange4", "brown"
# )

# e <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value)) +
# geom_jitter(aes(color = var_type), position = position_jitter(0.2)) + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
# #scale_size(range = c(5, 10)) +
# stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
# scale_colour_manual(values=c25) +
# theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

# ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_frontal_cortex_male_adult_(27_years)_and_male_adult_(35_years)_damage_score_stripchart.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

# p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value, fill=var_type)) +
# geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
# stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
# scale_fill_manual(values=c25) +
# scale_y_continuous(n.breaks = 10) +
# theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") 

# ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_frontal_cortex_male_adult_(27_years)_and_male_adult_(35_years)_damage_score_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


# ###

# long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_skeletal_muscle_cell_huvec.csv')

# c25 <- c(
#   "gold1",
#   "skyblue2", "#FB9A99", # lt pink
#   "dodgerblue2", "#E31A1C", # red
#   "green4",
#   "#6A3D9A", # purple
#   "#FF7F00", # orange
#   "black", 
#   "palegreen2",
#   "#CAB2D6", # lt purple
#   "#FDBF6F", # lt orange
#   "gray70", "khaki2",
#   "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
#   "darkturquoise", "green1", "yellow4", "yellow3",
#   "darkorange4", "brown"
# )

# e <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value)) +
# geom_jitter(aes(color = var_type), position = position_jitter(0.2)) + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
# #scale_size(range = c(5, 10)) +
# stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
# scale_colour_manual(values=c25) +
# theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

# ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_skeletal_muscle_cell_damage_score_stripchart.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


# p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value, fill=var_type)) +
# geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
# stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
# scale_fill_manual(values=c25) +
# scale_y_continuous(n.breaks = 10) +
# theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") 

# ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_skeletal_muscle_cell_damage_score_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

#############################################################################################

# import pandas as pd
# import numpy as np

# gains = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/HUVEC.gains.tsv', sep='\t')
# gains['var_type'] = 'gain'
# losses = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/HUVEC.losses.tsv', sep='\t')
# losses['var_type'] = 'loss'
# neutral = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/HUVEC.neutral.tsv', sep='\t')
# neutral['var_type'] = 'neutral'
# skew_df = pd.concat([gains, losses, neutral])

# class1 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_frontal_cortex_male_adult_(27_years)_and_male_adult_(35_years)_huvec.csv')
# skew_df.rename({'TEST.SNP.ID':'rsid'}, axis=1, inplace=True)
# class1_new = class1.merge(skew_df[['rsid', 'SKEW.ATAC']], how='left', on='rsid')
# class1_new.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_frontal_cortex_male_adult_(27_years)_and_male_adult_(35_years)_huvec_with_skew.csv', index=False)

# class2 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_skeletal_muscle_cell_huvec.csv')
# # skew_df.rename({'TEST.SNP.ID':'rsid'}, axis=1, inplace=True)
# class2_new = class2.merge(skew_df[['rsid', 'SKEW.ATAC']], how='left', on='rsid')
# class2_new.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_skeletal_muscle_cell_huvec_with_skew.csv', index=False)

# #############################################################################################


# .libPaths('/well/PROCARDIS/domwest/R')

# library(tidyverse)
# #library(showtext)
# #library(ggtext)
# library(ggrepel)
# library(dplyr) 
# library(RColorBrewer)

# long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_frontal_cortex_male_adult_(27_years)_and_male_adult_(35_years)_huvec_with_skew.csv')

# c25 <- c(
#   "gold1",
#   "skyblue2", "#FB9A99", # lt pink
#   "dodgerblue2", "#E31A1C", # red
#   "green4",
#   "#6A3D9A", # purple
#   "#FF7F00", # orange
#   "black", 
#   "palegreen2",
#   "#CAB2D6", # lt purple
#   "#FDBF6F", # lt orange
#   "gray70", "khaki2",
#   "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
#   "darkturquoise", "green1", "yellow4", "yellow3",
#   "darkorange4", "brown"
# )

# e <- ggplot(long_prediction_df, aes(x=SKEW.ATAC, y=damage_score_value, shape=var_type, color=var_type)) +
# geom_point() + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
# #scale_size(range = c(5, 10)) +
# #stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
# scale_colour_manual(values=c25) +
# theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

# ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_frontal_cortex_male_adult_(27_years)_and_male_adult_(35_years)_damage_score_skew_scatterplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

# #

# long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_skeletal_muscle_cell_huvec_with_skew.csv')

# c25 <- c(
#   "gold1",
#   "skyblue2", "#FB9A99", # lt pink
#   "dodgerblue2", "#E31A1C", # red
#   "green4",
#   "#6A3D9A", # purple
#   "#FF7F00", # orange
#   "black", 
#   "palegreen2",
#   "#CAB2D6", # lt purple
#   "#FDBF6F", # lt orange
#   "gray70", "khaki2",
#   "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
#   "darkturquoise", "green1", "yellow4", "yellow3",
#   "darkorange4", "brown"
# )

# e <- ggplot(long_prediction_df, aes(x=SKEW.ATAC, y=damage_score_value, shape=var_type, color=var_type)) +
# geom_point() + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
# #scale_size(range = c(5, 10)) +
# #stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
# scale_colour_manual(values=c25) +
# theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

# ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_skeletal_muscle_cell_damage_score_skew_scatterplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/



#############################################################################################

import pandas as pd
import numpy as np

gains = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/HUVEC.gains.tsv', sep='\t')
gains['var_type'] = 'gain'
losses = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/HUVEC.losses.tsv', sep='\t')
losses['var_type'] = 'loss'
neutral = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/HUVEC.neutral.tsv', sep='\t')
neutral['var_type'] = 'neutral'
skew_df = pd.concat([gains, losses, neutral])

class1 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_frontal_cortex_male_adult_(27_years)_and_male_adult_(35_years)_huvec.csv')
skew_df.rename({'TEST.SNP.ID':'rsid'}, axis=1, inplace=True)
class1_new = class1.merge(skew_df[['rsid', 'SKEW.ATAC']], how='left', on='rsid')

class2 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_skeletal_muscle_cell_huvec.csv')
# skew_df.rename({'TEST.SNP.ID':'rsid'}, axis=1, inplace=True)
class2_new = class2.merge(skew_df[['rsid', 'SKEW.ATAC']], how='left', on='rsid')

class3 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_stomach_female_adult_(53_years)_huvec.csv')
# skew_df.rename({'TEST.SNP.ID':'rsid'}, axis=1, inplace=True)
class3_new = class3.merge(skew_df[['rsid', 'SKEW.ATAC']], how='left', on='rsid')

class4 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_epidermal_melanocyte_huvec.csv')
# skew_df.rename({'TEST.SNP.ID':'rsid'}, axis=1, inplace=True)
class4_new = class4.merge(skew_df[['rsid', 'SKEW.ATAC']], how='left', on='rsid')

class_tot = pd.concat([class1_new, class2_new, class3_new, class4_new])

m_categories = ['gain', 'neutral', 'loss']
class_tot["var_type"] = pd.Categorical(class_tot["var_type"], categories = m_categories)
class_tot = class_tot.sort_values(by = "var_type")

class_tot.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/4_other_celltypes_huvec_with_skew.csv', index=False)

#############################################################################################


.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/4_other_celltypes_huvec_with_skew.csv')

p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value, fill=var_type)) +
geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
# scale_fill_manual(values=c26) +
scale_y_continuous(n.breaks = 10) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") 

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_4_other_celltypes_damage_score_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


#

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/4_other_celltypes_huvec_with_skew.csv')


c25 <- c(
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", 
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)


e <- ggplot(long_prediction_df, aes(x=SKEW.ATAC, y=damage_score_value, shape=var_type, color=classifier)) +
geom_point() + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
#stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_4_other_celltypes_damage_score_skew_scatterplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


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


import pandas as pd
import numpy as np

#checking distribution of ~400 top huvec gain, loss, neutral snps but in other celltypes (4 other random celltypes)

overarch = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/total_huvec_other_celltypes.csv')

for clas in overarch['classifier'].unique():
    classif = clas.replace(':','_')
    tmp = overarch[overarch['classifier']==clas]
    tmp.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/' + classif + '_huvec.csv', index=False)

#############################################################################################

.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

# -rw-r--r-- 1 wiq135 watkins  71286 Aug  9 10:06 DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec.csv
# -rw-r--r-- 1 wiq135 watkins  71888 Aug  9 10:06 DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec.csv


long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_frontal_cortex_male_adult_(27_years)_and_male_adult_(35_years)_huvec.csv')

c25 <- c(
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", 
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

e <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value)) +
geom_jitter(aes(color = var_type), position = position_jitter(0.2)) + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_frontal_cortex_male_adult_(27_years)_and_male_adult_(35_years)_damage_score_stripchart.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value, fill=var_type)) +
geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
scale_fill_manual(values=c25) +
scale_y_continuous(n.breaks = 10) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") 

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_frontal_cortex_male_adult_(27_years)_and_male_adult_(35_years)_damage_score_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


###

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_skeletal_muscle_cell_huvec.csv')

c25 <- c(
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", 
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

e <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value)) +
geom_jitter(aes(color = var_type), position = position_jitter(0.2)) + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_skeletal_muscle_cell_damage_score_stripchart.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value, fill=var_type)) +
geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
scale_fill_manual(values=c25) +
scale_y_continuous(n.breaks = 10) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") 

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_skeletal_muscle_cell_damage_score_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


###

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_epidermal_melanocyte_huvec.csv')

c25 <- c(
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", 
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

e <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value)) +
geom_jitter(aes(color = var_type), position = position_jitter(0.2)) + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_epidermal_melanocyte_damage_score_stripchart.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value, fill=var_type)) +
geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
scale_fill_manual(values=c25) +
scale_y_continuous(n.breaks = 10) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") 

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_epidermal_melanocyte_damage_score_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

###

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_stomach_female_adult_(53_years)_huvec.csv')

c25 <- c(
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", 
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

e <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value)) +
geom_jitter(aes(color = var_type), position = position_jitter(0.2)) + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_stomach_female_adult_(53_years)_damage_score_stripchart.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value, fill=var_type)) +
geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
scale_fill_manual(values=c25) +
scale_y_continuous(n.breaks = 10) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") 

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_stomach_female_adult_(53_years)_damage_score_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

#############################################################################################

import pandas as pd
import numpy as np

gains = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/HUVEC.gains.tsv', sep='\t')
gains['var_type'] = 'gain'
losses = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/HUVEC.losses.tsv', sep='\t')
losses['var_type'] = 'loss'
neutral = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/HUVEC.neutral.tsv', sep='\t')
neutral['var_type'] = 'neutral'
skew_df = pd.concat([gains, losses, neutral])

class1 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_frontal_cortex_male_adult_(27_years)_and_male_adult_(35_years)_huvec.csv')
skew_df.rename({'TEST.SNP.ID':'rsid'}, axis=1, inplace=True)
class1_new = class1.merge(skew_df[['rsid', 'SKEW.ATAC']], how='left', on='rsid')
class1_new.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_frontal_cortex_male_adult_(27_years)_and_male_adult_(35_years)_huvec_with_skew.csv', index=False)

class2 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_skeletal_muscle_cell_huvec.csv')
# skew_df.rename({'TEST.SNP.ID':'rsid'}, axis=1, inplace=True)
class2_new = class2.merge(skew_df[['rsid', 'SKEW.ATAC']], how='left', on='rsid')
class2_new.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_skeletal_muscle_cell_huvec_with_skew.csv', index=False)

class3 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_epidermal_melanocyte_huvec.csv')
# skew_df.rename({'TEST.SNP.ID':'rsid'}, axis=1, inplace=True)
class3_new = class3.merge(skew_df[['rsid', 'SKEW.ATAC']], how='left', on='rsid')
class3_new.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_epidermal_melanocyte_huvec_with_skew.csv', index=False)

class4 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_stomach_female_adult_(53_years)_huvec.csv')
# skew_df.rename({'TEST.SNP.ID':'rsid'}, axis=1, inplace=True)
class4_new = class4.merge(skew_df[['rsid', 'SKEW.ATAC']], how='left', on='rsid')
class4_new.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_stomach_female_adult_(53_years)_huvec_with_skew.csv', index=False)

#############################################################################################


.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_frontal_cortex_male_adult_(27_years)_and_male_adult_(35_years)_huvec_with_skew.csv')

c25 <- c(
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", 
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

e <- ggplot(long_prediction_df, aes(x=SKEW.ATAC, y=damage_score_value, shape=var_type, color=var_type)) +
geom_point() + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
#stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_frontal_cortex_male_adult_(27_years)_and_male_adult_(35_years)_damage_score_skew_scatterplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

#

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_skeletal_muscle_cell_huvec_with_skew.csv')

c25 <- c(
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", 
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

e <- ggplot(long_prediction_df, aes(x=SKEW.ATAC, y=damage_score_value, shape=var_type, color=var_type)) +
geom_point() + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
#stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_skeletal_muscle_cell_damage_score_skew_scatterplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

#

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_epidermal_melanocyte_huvec_with_skew.csv')

c25 <- c(
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", 
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

e <- ggplot(long_prediction_df, aes(x=SKEW.ATAC, y=damage_score_value, shape=var_type, color=var_type)) +
geom_point() + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
#stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_epidermal_melanocyte_damage_score_skew_scatterplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

#

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_stomach_female_adult_(53_years)_huvec_with_skew.csv')

c25 <- c(
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", 
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

e <- ggplot(long_prediction_df, aes(x=SKEW.ATAC, y=damage_score_value, shape=var_type, color=var_type)) +
geom_point() + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
#stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_stomach_female_adult_(53_years)_damage_score_skew_scatterplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


#############################################################################################

import pandas as pd
import numpy as np

gains = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/HUVEC.gains.tsv', sep='\t')
gains['var_type'] = 'gain'
losses = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/HUVEC.losses.tsv', sep='\t')
losses['var_type'] = 'loss'
neutral = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/HUVEC.neutral.tsv', sep='\t')
neutral['var_type'] = 'neutral'
skew_df = pd.concat([gains, losses, neutral])

class1 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_frontal_cortex_male_adult_(27_years)_and_male_adult_(35_years)_huvec.csv')
skew_df.rename({'TEST.SNP.ID':'rsid'}, axis=1, inplace=True)
class1_new = class1.merge(skew_df[['rsid', 'SKEW.ATAC']], how='left', on='rsid')

class2 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_skeletal_muscle_cell_huvec.csv')
# skew_df.rename({'TEST.SNP.ID':'rsid'}, axis=1, inplace=True)
class2_new = class2.merge(skew_df[['rsid', 'SKEW.ATAC']], how='left', on='rsid')

class3 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_epidermal_melanocyte_huvec.csv')
# skew_df.rename({'TEST.SNP.ID':'rsid'}, axis=1, inplace=True)
class3_new = class3.merge(skew_df[['rsid', 'SKEW.ATAC']], how='left', on='rsid')

class4 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_stomach_female_adult_(53_years)_huvec.csv')
# skew_df.rename({'TEST.SNP.ID':'rsid'}, axis=1, inplace=True)
class4_new = class4.merge(skew_df[['rsid', 'SKEW.ATAC']], how='left', on='rsid')

class_tot = pd.concat([class1_new, class2_new, class3_new, class4_new])

m_categories = ['gain', 'neutral', 'loss']
class_tot["var_type"] = pd.Categorical(class_tot["var_type"], categories = m_categories)
class_tot = class_tot.sort_values(by = "var_type")

class_tot.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_frontal_cortex_male_adult_(27_years)_and_male_adult_(35_years)_and_DNASE_skeletal_muscle_cell_and_DNASE_epidermal_melanocyte_and_DNASE_stomach_female_adult_(53_years)_huvec_with_skew.csv', index=False)

#############################################################################################


.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_frontal_cortex_male_adult_(27_years)_and_male_adult_(35_years)_and_DNASE_skeletal_muscle_cell_and_DNASE_epidermal_melanocyte_and_DNASE_stomach_female_adult_(53_years)_huvec_with_skew.csv')

p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value, fill=var_type)) +
geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
# scale_fill_manual(values=c26) +
scale_y_continuous(n.breaks = 10) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") 

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_frontal_cortex_male_adult_(27_years)_and_male_adult_(35_years)_and_DNASE_skeletal_muscle_cell_and_DNASE_epidermal_melanocyte_and_DNASE_stomach_female_adult_(53_years)_damage_score_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


#

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_frontal_cortex_male_adult_(27_years)_and_male_adult_(35_years)_and_DNASE_skeletal_muscle_cell_and_DNASE_epidermal_melanocyte_and_DNASE_stomach_female_adult_(53_years)_huvec_with_skew.csv')


c25 <- c(
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", 
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)


e <- ggplot(long_prediction_df, aes(x=SKEW.ATAC, y=damage_score_value, shape=var_type, color=classifier)) +
geom_point() + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
#stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_frontal_cortex_male_adult_(27_years)_and_male_adult_(35_years)_and_DNASE_skeletal_muscle_cell_and_DNASE_epidermal_melanocyte_and_DNASE_stomach_female_adult_(53_years)_damage_score_skew_scatterplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

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

import pandas as pd
import numpy as np

#OG ie endothelial classifiers BUT FULL SET OF SNPS (not just 100 gain, 100 neutral, 100 loss)

overarch = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/total_huvec_full_set.csv')

for clas in overarch['classifier'].unique():
    classif = clas.replace(':','_')
    tmp = overarch[overarch['classifier']==clas]
    tmp.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/full_set_' + classif + '_huvec.csv', index=False)

#############################################################################################

.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

# -rw-r--r-- 1 wiq135 watkins  71286 Aug  9 10:06 DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec.csv
# -rw-r--r-- 1 wiq135 watkins  71888 Aug  9 10:06 DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec.csv


long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/full_set_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec.csv')

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

e <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value)) +
geom_jitter(aes(color = var_type), position = position_jitter(0.2)) + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_full_set_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_damage_score_stripchart.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value)) +
geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
scale_fill_manual(values=c25) +
scale_y_continuous(n.breaks = 10) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") 

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_full_set_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_damage_score_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


###

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/full_set_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec.csv')

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

e <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value)) +
geom_jitter(aes(color = var_type), position = position_jitter(0.2)) + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_full_set_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_damage_score_stripchart.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value)) +
geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
scale_fill_manual(values=c25) +
scale_y_continuous(n.breaks = 10) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") 

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_full_set_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_damage_score_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

#############################################################################################

import pandas as pd
import numpy as np

skew_df = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/Simone.HUVEC.full.tsv', sep='\t')

class1 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/full_set_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec.csv')
skew_df.rename({'TEST.SNP.ID':'rsid'}, axis=1, inplace=True)
class1_new = class1.merge(skew_df[['rsid', 'SKEW.ATAC']], how='left', on='rsid')
class1_new.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/full_set_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec_with_skew.csv', index=False)

class2 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/full_set_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec.csv')
# skew_df.rename({'TEST.SNP.ID':'rsid'}, axis=1, inplace=True)
class2_new = class2.merge(skew_df[['rsid', 'SKEW.ATAC']], how='left', on='rsid')
class2_new.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/full_set_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec_with_skew.csv', index=False)

#############################################################################################


.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/full_set_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec_with_skew.csv')

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

e <- ggplot(long_prediction_df, aes(x=SKEW.ATAC, y=damage_score_value)) +
geom_point() + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
#stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_full_set_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_damage_score_skew_scatterplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

#

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/full_set_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec_with_skew.csv')

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

e <- ggplot(long_prediction_df, aes(x=SKEW.ATAC, y=damage_score_value)) +
geom_point() + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
#stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_full_set_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_damage_score_skew_scatterplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

#############################################################################################

import pandas as pd
import numpy as np

skew_df = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/Simone.HUVEC.full.tsv', sep='\t')

class1 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/full_set_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec.csv')
skew_df.rename({'TEST.SNP.ID':'rsid'}, axis=1, inplace=True)
class1_new = class1.merge(skew_df[['rsid', 'SKEW.ATAC']], how='left', on='rsid')

class2 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/full_set_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec.csv')
# skew_df.rename({'TEST.SNP.ID':'rsid'}, axis=1, inplace=True)
class2_new = class2.merge(skew_df[['rsid', 'SKEW.ATAC']], how='left', on='rsid')

class_tot = pd.concat([class1_new, class2_new])
class_tot.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/full_set_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_huvec_with_skew.csv', index=False)

#############################################################################################


.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/full_set_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_huvec_with_skew.csv')

p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value)) +
geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
# scale_fill_manual(values=c26) +
scale_y_continuous(n.breaks = 10) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") 

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_full_set_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_damage_score_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


#

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/full_set_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_huvec_with_skew.csv')


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

e <- ggplot(long_prediction_df, aes(x=SKEW.ATAC, y=damage_score_value, color=classifier)) +
geom_point() + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
#stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_full_set_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_and_2_damage_score_skew_scatterplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


# long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec_with_skew.csv')

# value_stats <- long_prediction_df |>
# summarize(mean_value = mean(damage_score_value))
# value_stats

# # #probability density plot
# # ggplot(long_prediction_df, aes(x = damage_score_value, y = after_stat(density))) +
# #     geom_histogram(bins = 200) +
# #     geom_vline(aes(xintercept = mean_value), value_stats, color = "darkolivegreen1", linewidth = 2) +
# #     geom_density(color = "darkcyan", linewidth = 2) +
# #     scale_x_continuous(n.breaks = 10)
# # #default is binwidth 30

# # ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_damage_score_prob_density_plot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


# long_prediction_df <- long_prediction_df |>
#     mutate(variable = factor(classifier))

# ggplot(data = long_prediction_df, aes(x = damage_score_value, fill = classifier)) +
#     geom_histogram(bins = 200) +
#     scale_x_continuous(n.breaks = 10)

# ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_damage_score_histogram.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

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

import pandas as pd
import numpy as np

#1000g set, similar to approach for huvec more neutrals -- so technically these are not huvec - only ~2000 of the 1000G snps

overarch = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/total_1000g.csv')

for clas in overarch['classifier'].unique():
    classif = clas.replace(':','_')
    tmp = overarch[overarch['classifier']==clas]
    tmp.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/' + classif + '_1000g.csv', index=False)


######################################################################################################

# DNASE_cardiac_muscle_cell_1_1000g.csv
# DNASE_cardiac_muscle_cell_2_1000g.csv
# DNASE_heart_left_ventricle_female_adult_(53_years)_1000g.csv
# DNASE_heart_left_ventricle_female_embryo_(101_day)_and_female_embryo_(103_days)_1000g.csv
# DNASE_heart_left_ventricle_female_embryo_(136_days)_1000g.csv

.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_1_1000g.csv')

c25 <- c(
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown",
    "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2"
)

e <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value)) +
geom_jitter(aes(color = var_type), position = position_jitter(0.2)) + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/set_1000g_DNASE_cardiac_muscle_cell_1_damage_score_stripchart.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value, fill=var_type)) +
geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
scale_fill_manual(values=c25) +
scale_y_continuous(n.breaks = 10) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") 

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/set_1000g_DNASE_cardiac_muscle_cell_1_damage_score_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


###

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_2_1000g.csv')

c25 <- c(
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown",
    "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2"
)

e <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value)) +
geom_jitter(aes(color = var_type), position = position_jitter(0.2)) + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/set_1000g_DNASE_cardiac_muscle_cell_2_damage_score_stripchart.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value, fill=var_type)) +
geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
scale_fill_manual(values=c25) +
scale_y_continuous(n.breaks = 10) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") 

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/set_1000g_DNASE_cardiac_muscle_cell_2_damage_score_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

#############################################################################################

import pandas as pd
import numpy as np

class1_new = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_1_1000g.csv')

class2_new = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_2_1000g.csv')

class_tot = pd.concat([class1_new, class2_new])

class_tot.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_1_and_2_1000g.csv', index=False)

#############################################################################################


.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_1_and_2_1000g.csv')

p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value, fill=var_type)) +
geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
# scale_fill_manual(values=c26) +
scale_y_continuous(n.breaks = 10) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") 

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/set_1000g_DNASE_cardiac_muscle_cell_1_and_2_damage_score_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

################################################################################################################################################################

.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_1_and_2_1000g.csv')

value_stats <- long_prediction_df |>
summarize(mean_value = mean(damage_score_value))
value_stats

# #probability density plot
# ggplot(long_prediction_df, aes(x = damage_score_value, y = after_stat(density))) +
#     geom_histogram(bins = 200) +
#     geom_vline(aes(xintercept = mean_value), value_stats, color = "darkolivegreen1", linewidth = 2) +
#     geom_density(color = "darkcyan", linewidth = 2) +
#     scale_x_continuous(n.breaks = 10)
# #default is binwidth 30

# ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_damage_score_prob_density_plot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


long_prediction_df <- long_prediction_df |>
    mutate(variable = factor(classifier))

ggplot(data = long_prediction_df, aes(x = damage_score_value, fill = classifier)) +
    geom_histogram(bins = 200) +
    scale_x_continuous(n.breaks = 10)

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/set_1000g_DNASE_cardiac_muscle_cell_1_and_2_damage_score_histogram.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

##


# long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec_with_skew.csv')

# value_stats <- long_prediction_df |>
# summarize(mean_value = mean(damage_score_value))
# value_stats

# # #probability density plot
# # ggplot(long_prediction_df, aes(x = damage_score_value, y = after_stat(density))) +
# #     geom_histogram(bins = 200) +
# #     geom_vline(aes(xintercept = mean_value), value_stats, color = "darkolivegreen1", linewidth = 2) +
# #     geom_density(color = "darkcyan", linewidth = 2) +
# #     scale_x_continuous(n.breaks = 10)
# # #default is binwidth 30

# # ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_damage_score_prob_density_plot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


# long_prediction_df <- long_prediction_df |>
#     mutate(variable = factor(classifier))

# ggplot(data = long_prediction_df, aes(x = damage_score_value, fill = classifier)) +
#     geom_histogram(bins = 200) +
#     scale_x_continuous(n.breaks = 10)

# ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_damage_score_histogram.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


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

import pandas as pd
import numpy as np

#OG ie hcm classifiers but only 'prioritised snps' based on 1000g set above used to define thresholds at 0.05 p value 

#############################################################################################

import pandas as pd
import numpy as np

class1_1 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_1_1000g_less_than_loss_threshold.csv')
class1_2 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_1_1000g_greater_than_gain_threshold.csv')
class1_new = pd.concat([class1_1, class1_2])

class2_1 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_2_1000g_less_than_loss_threshold.csv')
class2_2 = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_2_1000g_greater_than_gain_threshold.csv')
class2_new = pd.concat([class2_1, class2_2])

class_tot = pd.concat([class1_new, class2_new])

class_tot.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_1_and_2_1000g_prioritised.csv', index=False)

#############################################################################################


.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_1_and_2_1000g_prioritised.csv')

p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value, fill=classifier)) +
geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
# scale_fill_manual(values=c26) +
scale_y_continuous(n.breaks = 10) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") 

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/set_1000g_prioritised_DNASE_cardiac_muscle_cell_1_and_2_damage_score_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

################################################################################################################################################################

.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_1_1000g.csv')

value_stats <- long_prediction_df |>
summarize(mean_value = mean(damage_score_value))
value_stats

long_prediction_df <- long_prediction_df |>
    mutate(variable = factor(classifier))

ggplot(data = long_prediction_df, aes(x = damage_score_value, fill = classifier)) +
    geom_histogram(bins = 200) +
    geom_vline(xintercept=0.0048274043947458) +
    geom_vline(xintercept=-0.0045122168958187) +
    scale_x_continuous(n.breaks = 10)

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/set_1000g_prioritised_DNASE_cardiac_muscle_cell_1_1000g_damage_score_histogram.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

##

.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_2_1000g.csv')

value_stats <- long_prediction_df |>
summarize(mean_value = mean(damage_score_value))
value_stats

long_prediction_df <- long_prediction_df |>
    mutate(variable = factor(classifier))

ggplot(data = long_prediction_df, aes(x = damage_score_value, fill = classifier)) +
    geom_histogram(bins = 200) +
    geom_vline(xintercept=0.00795772485435) +
    geom_vline(xintercept=-0.0056563951075077) +
    scale_x_continuous(n.breaks = 10)

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/set_1000g_prioritised_DNASE_cardiac_muscle_cell_2_1000g_damage_score_histogram.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


##


# long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec_with_skew.csv')

# value_stats <- long_prediction_df |>
# summarize(mean_value = mean(damage_score_value))
# value_stats

# # #probability density plot
# # ggplot(long_prediction_df, aes(x = damage_score_value, y = after_stat(density))) +
# #     geom_histogram(bins = 200) +
# #     geom_vline(aes(xintercept = mean_value), value_stats, color = "darkolivegreen1", linewidth = 2) +
# #     geom_density(color = "darkcyan", linewidth = 2) +
# #     scale_x_continuous(n.breaks = 10)
# # #default is binwidth 30

# # ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_damage_score_prob_density_plot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


# long_prediction_df <- long_prediction_df |>
#     mutate(variable = factor(classifier))

# ggplot(data = long_prediction_df, aes(x = damage_score_value, fill = classifier)) +
#     geom_histogram(bins = 200) +
#     scale_x_continuous(n.breaks = 10)

# ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_damage_score_histogram.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


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

import pandas as pd
import numpy as np

#1000g set, similar to approach for huvec more neutrals -- so technically these are not huvec (10000 of the 1000G snps - without eaf_range_filt)

overarch = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/total_1000g_10000_snps.csv')

for clas in overarch['classifier'].unique():
    classif = clas.replace(':','_')
    tmp = overarch[overarch['classifier']==clas]
    tmp.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/' + classif + '_1000g_10000_snps.csv', index=False)


######################################################################################################

# DNASE_cardiac_muscle_cell_1_1000g.csv
# DNASE_cardiac_muscle_cell_2_1000g.csv
# DNASE_heart_left_ventricle_female_adult_(53_years)_1000g.csv
# DNASE_heart_left_ventricle_female_embryo_(101_day)_and_female_embryo_(103_days)_1000g.csv
# DNASE_heart_left_ventricle_female_embryo_(136_days)_1000g.csv

.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_1_1000g_10000_snps.csv')

c25 <- c(
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown",
    "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2"
)

e <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value)) +
geom_jitter(aes(color = var_type), position = position_jitter(0.2)) + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/set_1000g_10000_snps_DNASE_cardiac_muscle_cell_1_damage_score_stripchart.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value, fill=var_type)) +
geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
scale_fill_manual(values=c25) +
scale_y_continuous(n.breaks = 10) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") 

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/set_1000g_10000_snps_DNASE_cardiac_muscle_cell_1_damage_score_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


###

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_2_1000g_10000_snps.csv')

c25 <- c(
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown",
    "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2"
)

e <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value)) +
geom_jitter(aes(color = var_type), position = position_jitter(0.2)) + #geom_jitter(aes(size = damage_score_value, color = classifier, alpha = 0.5), position = position_jitter(0.2)) +
#scale_size(range = c(5, 10)) +
stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
scale_colour_manual(values=c25) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/set_1000g_10000_snps_DNASE_cardiac_muscle_cell_2_damage_score_stripchart.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value, fill=var_type)) +
geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
scale_fill_manual(values=c25) +
scale_y_continuous(n.breaks = 10) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") 

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/set_1000g_10000_snps_DNASE_cardiac_muscle_cell_2_damage_score_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

#############################################################################################

import pandas as pd
import numpy as np

class1_new = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_1_1000g_10000_snps.csv')

class2_new = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_2_1000g_10000_snps.csv')

class_tot = pd.concat([class1_new, class2_new])

class_tot.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_1_and_2_1000g_10000_snps.csv', index=False)

#############################################################################################


.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_1_and_2_1000g_10000_snps.csv')

p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value, fill=var_type)) +
geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
# scale_fill_manual(values=c26) +
scale_y_continuous(n.breaks = 10) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") 

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/set_1000g_10000_snps_DNASE_cardiac_muscle_cell_1_and_2_damage_score_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

################################################################################################################################################################

.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_1_and_2_1000g_10000_snps.csv')

value_stats <- long_prediction_df |>
summarize(mean_value = mean(damage_score_value))
value_stats

# #probability density plot
# ggplot(long_prediction_df, aes(x = damage_score_value, y = after_stat(density))) +
#     geom_histogram(bins = 200) +
#     geom_vline(aes(xintercept = mean_value), value_stats, color = "darkolivegreen1", linewidth = 2) +
#     geom_density(color = "darkcyan", linewidth = 2) +
#     scale_x_continuous(n.breaks = 10)
# #default is binwidth 30

# ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_1_damage_score_prob_density_plot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


long_prediction_df <- long_prediction_df |>
    mutate(variable = factor(classifier))

ggplot(data = long_prediction_df, aes(x = damage_score_value, fill = classifier)) +
    geom_histogram(bins = 1000) +
    scale_x_continuous(n.breaks = 10)

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/set_1000g_10000_snps_DNASE_cardiac_muscle_cell_1_and_2_damage_score_histogram.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_1_1000g_10000_snps.csv')

value_stats <- long_prediction_df |>
summarize(mean_value = mean(damage_score_value))
value_stats

long_prediction_df <- long_prediction_df |>
    mutate(variable = factor(classifier))

ggplot(data = long_prediction_df, aes(x = damage_score_value)) +
    geom_histogram(bins = 1000, fill = "darkturquoise") +
    geom_vline(xintercept=0.0755509436130523, color = "darkslateblue") +
    geom_vline(xintercept=-0.0367254391312599, color = "darkslateblue") +
    scale_x_continuous(n.breaks = 10) +
    ggtitle("10000 snps from the 1000G set distribution: DNASE:cardiac_muscle_cell_1; index=85") #2 is 586

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/set_1000g_10000_snps_DNASE_cardiac_muscle_cell_1_damage_score_histogram.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices/DNASE_cardiac_muscle_cell_2_1000g_10000_snps.csv')

value_stats <- long_prediction_df |>
summarize(mean_value = mean(damage_score_value))
value_stats

long_prediction_df <- long_prediction_df |>
    mutate(variable = factor(classifier))

ggplot(data = long_prediction_df, aes(x = damage_score_value)) +
    geom_histogram(bins = 1000, fill = "darkturquoise") +
    geom_vline(xintercept=0.0328718870878219, color = "darkslateblue") +
    geom_vline(xintercept=-0.0260984227061271, color = "darkslateblue") +
    scale_x_continuous(n.breaks = 10) +
    ggtitle("10000 snps from the 1000G set distribution: DNASE:cardiac_muscle_cell_2; index=586") #2 is 586

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/set_1000g_10000_snps_DNASE_cardiac_muscle_cell_2_damage_score_histogram.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


##


# long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec_with_skew.csv')

# value_stats <- long_prediction_df |>
# summarize(mean_value = mean(damage_score_value))
# value_stats

# # #probability density plot
# # ggplot(long_prediction_df, aes(x = damage_score_value, y = after_stat(density))) +
# #     geom_histogram(bins = 200) +
# #     geom_vline(aes(xintercept = mean_value), value_stats, color = "darkolivegreen1", linewidth = 2) +
# #     geom_density(color = "darkcyan", linewidth = 2) +
# #     scale_x_continuous(n.breaks = 10)
# # #default is binwidth 30

# # ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_damage_score_prob_density_plot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


# long_prediction_df <- long_prediction_df |>
#     mutate(variable = factor(classifier))

# ggplot(data = long_prediction_df, aes(x = damage_score_value, fill = classifier)) +
#     geom_histogram(bins = 200) +
#     scale_x_continuous(n.breaks = 10)

# ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_DNASE_endothelial_cell_of_umbilical_vein_newborn_2_damage_score_histogram.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


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
