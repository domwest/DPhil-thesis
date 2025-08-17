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

overarch = pd.read_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/total_huvec.csv')

##################################################################################################

#first layer: gof or lof
overarch['damage_score_ann'] =''
overarch['damage_score_ann'][overarch['damage_score_value']<0] = 'LoF'
overarch['damage_score_ann'][overarch['damage_score_value']>0] = 'GoF'
overarch['damage_score_ann'][overarch['damage_score_value']==0] = 'No change'
#
inaccurate_gains = overarch[(overarch['var_type']=='gain')&(overarch['damage_score_ann']!='GoF')]
inaccurate_losses = overarch[(overarch['var_type']=='loss')&(overarch['damage_score_ann']!='LoF')]
inaccurate_neutrals = overarch[(overarch['var_type']=='neutral')&(overarch['damage_score_ann']!='No change')]
# overarch = overarch[~(overarch['damage_score_ann']=='No change')]
overarch.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/total_huvec_with_damage_score_ann.csv', index=False)

#2nd layer: top or bottom gain vs lof
overarch_open = overarch[overarch['ref_score_value']>=0.8]
overarch_open['damage_score_ann2'] =''
overarch_open['damage_score_ann2'][(overarch_open['var_score_value']>overarch_open['ref_score_value'])] = 'top_gain'
overarch_open['damage_score_ann2'][(overarch_open['var_score_value']<overarch_open['ref_score_value'])] = 'top_loss'
overarch_open['damage_score_ann2'][(overarch_open['var_score_value']==overarch_open['ref_score_value'])] = 'no_change'
# overarch_open = overarch_open[~(overarch_open['damage_score_ann2']=='no_change')]
overarch_open.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/total_huvec_open_with_damage_score_ann.csv', index=False)

overarch_closed = overarch[overarch['ref_score_value']<=0.2]
overarch_closed['damage_score_ann2'] =''
overarch_closed['damage_score_ann2'][(overarch_closed['var_score_value']>overarch_closed['ref_score_value'])] = 'bottom_gain'
overarch_closed['damage_score_ann2'][(overarch_closed['var_score_value']<overarch_closed['ref_score_value'])] = 'bottom_loss'
overarch_closed['damage_score_ann2'][(overarch_closed['var_score_value']==overarch_closed['ref_score_value'])] = 'no_change'
# overarch_closed = overarch_closed[~(overarch_closed['damage_score_ann2']=='no_change')]
overarch_closed.to_csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/total_huvec_closed_with_damage_score_ann.csv', index=False)

#############################################################################################################################

##in R

.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/total_huvec_with_damage_score_ann.csv')

value_stats <- long_prediction_df |>
summarize(mean_value = mean(damage_score_value))
value_stats

c26 <- c("#CAB2D6", "darkturquoise", "#CAB2D2")

p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value, fill=damage_score_ann)) +
geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
scale_fill_manual(values=c26) +
scale_y_continuous(n.breaks = 10) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") 

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_comparing_gof_lof_damage_score_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


####

long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/total_huvec_open_with_damage_score_ann.csv')

p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value, fill=damage_score_ann2)) +
geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
scale_fill_manual(values=c26) +
scale_y_continuous(n.breaks = 10) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") 

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/huvec_comparing_open_damage_score_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


long_prediction_df <- read.csv('/gpfs3/well/PROCARDIS/domwest/deepmind/matrices_huvec/all_classifiers_total_hcm_closed_with_damage_score_ann.csv')

p <- ggplot(long_prediction_df, aes(x=classifier, y=damage_score_value, fill=damage_score_ann2)) +
geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
scale_fill_manual(values=c26) +
scale_y_continuous(n.breaks = 10) +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") 

ggsave(filename="/well/PROCARDIS/domwest/deepmind/ds_scatter_plots/total_huvec_closed_with_damage_score_ann.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/
