#IN TERMINAL!!!

#!/usr/bin/env python
# coding: utf-8

##first in bash:
# module purge
# module load BEDTools
# module load python
# module load SciPy-bundle 
# module load pybedtools                                                                                           
# python3

import warnings
warnings.filterwarnings('ignore')

#from IPython.core.display import display, HTML
#display(HTML('<style>.container {width: 90% !important; }</style>'))

import pandas as pd
import glob
import os
import shutil
import pybedtools


# In[31]:


if not os.path.isdir("/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/tmp"):
    os.mkdir("tmp")
pybedtools.set_tempdir('tmp')


# In[32]:


beds1 = glob.glob("/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/successful_merged_runs/*/genetics/CATCH-UP/analysis/results/09_peak_calling/*_L-tron.bed") 
#contains: UpStreamPipeline_success_single_merged, UpStreamPipeline_success_paired_merged, UpStreamPipeline_success_single_merged_SC, UpStreamPipeline_success_paired_merged_SC

# beds2 = glob.glob("/well/PROCARDIS/domwest/upstr_processing/atac_sc/lt_web_find_and_score_peaks_output/*.tsv")
#already captured in the below

beds2 = glob.glob("/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/dh_input_data_220524/*")
#copied the files over from /well/PROCARDIS/domwest/celltype_enrichment/peaks/ and then renamed them to include the word 'ATAC' in files that didn't already have ATAC or DNase labelled in the name
#this is because Simone's files (starting with a capital letter) were the ones that didn't have a label specifying the methodology, and I want my classifiers to have this level of info in the name

print(len(beds1), "nr files") #files not cell_types, because multiple files for bulk lv tissue and likewise for cardiac_muscle_cells
#19
print(len(beds2), "nr files")
#51

#let's remove some as suggested by Jim: ie H3K36 and H3K9 data...
#len([ x for x in beds1 if "H3K9" not in x ])
#17
#len([ x for x in beds1 if "H3K36" not in x ])
#17
#so 2 H3K9 and 2 H3K36 files

beds1 = [ x for x in beds1 if "H3K9" not in x ]
beds1 = [ x for x in beds1 if "H3K36" not in x ] #now 15

#got blacklisted regions from here:https: //github.com/Boyle-Lab/Blacklist/tree/master
blackregions = pd.read_csv("/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/hg38-blacklist.v2.bed", sep='\t', header=None)
blackregions=blackregions[[0, 1, 2]]
blackregions.rename({0:"chrom", 1:"start", 2:"end"}, axis=1, inplace=True)
#636 rows
blackregions.head()
blackregions = pybedtools.BedTool.from_dataframe(blackregions)

#i already filtered all the peaks to be >0.5 peak score for this beds2 dir, so won't be chekcing for nr peaks before score filtering here...

beds2_dict = {}
for i in beds2:
    name_pre = i.split("/")[-1].replace("_L-tron.bed", "")
    df = pd.read_csv(i, sep='\t', header=None)
    print(df.head())
    # (1) removing blacklist regions
    col_names = df.columns.tolist()
    print('col_names')
    print(col_names)
    df = pybedtools.BedTool.from_dataframe(df) #creates bedtool from pandas df
    print('df.head')
    print(df.head())
    print(len(df))
    #blackregions = pybedtools.BedTool.from_dataframe(blackregions)
    print('blackregions.head')
    print(blackregions.head())
    print(len(blackregions))
    df = df.intersect(blackregions, v=True) 
    #not too sure what intersect does, but guessing v=True means it removes rows where they interesect ie black regions?
    #for 1 file eg, why does blacklist df = 636, original df = 35664, and intersected df = 32058
    print(len(df))
    print('intersected_df')
    print(df.head())
    df = df.to_dataframe(names=col_names)
    print('df_with_colnames')
    print(df.head())
    df.reset_index(drop=True, inplace=True)
    df["cell_type"] = name_pre
    beds2_dict[name_pre] = df
    print(beds2_dict)

keys = []
nr_peaks = []
for i in beds2_dict.keys():
    keys.append(i)
    nr_peaks.append(beds2_dict[i].shape[0])
    print("%s has %s peaks after filtering (.5 lanceotron threshold)"%(i, beds2_dict[i].shape[0]))
beds2_dict_df = pd.DataFrame()
beds2_dict_df['key'] = keys 
beds2_dict_df['nr_peaks_post_LT_0.5_filter'] = nr_peaks

beds1_dict = {}
for i in beds1:
    if 'paired' in i.split('/')[7]:
        seq = 'paired'
    else:
        seq = 'single'
    if 'SC' in i.split('/')[7]:
        name = 'cardiac_muscle_cell'
    else:
        name = 'bulk_heart_LV_tissue'
    print(seq)
    print(name)
    methodology = i.split("/")[-1].replace("_hg38_L-tron.bed", "")
    print(methodology)
    final_name = name + '_' + methodology + '_' + seq
    print(final_name)
    break
#single
#cardiac_muscle_cell
#CHIP_HISTONE_H3K4me3
#cardiac_muscle_cell_CHIP_HISTONE_H3K4me3_single

# beds1_dict_non = {}
# for i in beds1:
#     print(i)
#     if 'paired' in i.split('/')[7]:
#         seq = 'paired'
#     else:
#         seq = 'single'
#     if 'SC' in i.split('/')[7]:
#         name = 'cardiac_muscle_cell'
#     else:
#         name = 'bulk_heart_LV_tissue'
#     print(name)
#     methodology = i.split("/")[-1].replace("_hg38_L-tron.bed", "")
#     print(methodology)
#     final_name = name + '_' + methodology + '_' + seq
#     print(final_name)
#     df = pd.read_csv(i, sep='\t')
#     print(df.head())
#     # (1) removing blacklist regions
#     col_names = df.columns.tolist()
#     print('col_names')
#     print(col_names)
#     df = pybedtools.BedTool.from_dataframe(df) #creates bedtool from pandas df
#     print('df.head')
#     print(df.head())
#     print(len(df))
#     #blackregions = pybedtools.BedTool.from_dataframe(blackregions)
#     print('blackregions.head')
#     print(blackregions.head())
#     print(len(blackregions))
#     df = df.intersect(blackregions, v=True) 
#     #not too sure what intersect does, but guessing v=True means it removes rows where they interesect ie black regions?
#     #for 1 file eg, why does blacklist df = 636, original df = 35664, and intersected df = 32058
#     print(len(df))
#     print('intersected_df')
#     print(df.head())
#     df = df.to_dataframe(names=col_names)
#     print('df_with_colnames')
#     print(df.head())
#     # (2) filtering peaks
#     df.reset_index(drop=True, inplace=True)
#     df["cell_type"] = final_name
#     beds1_dict_non[final_name] = df
#     print(beds1_dict_non)

# for i in beds1_dict_non.keys():
#     print("%s has %s peaks without filtering (.5 lanceotron threshold)"%(i, beds1_dict_non[i].shape[0]))

beds1_dict = {}
for i in beds1:
    print(i)
    if 'paired' in i.split('/')[7]:
        seq = 'paired'
    else:
        seq = 'single'
    if 'SC' in i.split('/')[7]:
        name = 'cardiac_muscle_cell'
    else:
        name = 'bulk_heart_LV_tissue'
    print(name)
    methodology = i.split("/")[-1].replace("_hg38_L-tron.bed", "")
    print(methodology)
    final_name = name + '_' + methodology + '_' + seq
    print(final_name)
    df = pd.read_csv(i, sep='\t')
    print(df.head())
    # (1) removing blacklist regions
    col_names = df.columns.tolist()
    print('col_names')
    print(col_names)
    df = pybedtools.BedTool.from_dataframe(df) #creates bedtool from pandas df
    print('df.head')
    print(df.head())
    print(len(df))
    #blackregions = pybedtools.BedTool.from_dataframe(blackregions)
    print('blackregions.head')
    print(blackregions.head())
    print(len(blackregions))
    df = df.intersect(blackregions, v=True) 
    #not too sure what intersect does, but guessing v=True means it removes rows where they interesect ie black regions?
    #for 1 file eg, why does blacklist df = 636, original df = 35664, and intersected df = 32058
    print(len(df))
    print('intersected_df')
    print(df.head())
    df = df.to_dataframe(names=col_names)
    print('df_with_colnames')
    print(df.head())
    # (2) filtering peaks
    df = df[df['overall_peak_score'] >= 0.5]
    df.reset_index(drop=True, inplace=True)
    df["cell_type"] = final_name
    beds1_dict[final_name] = df
    print(beds1_dict)

keys = []
nr_peaks = []
for i in beds1_dict.keys():
    keys.append(i)
    nr_peaks.append(beds1_dict[i].shape[0])
    print("%s has %s peaks after filtering (.5 lanceotron threshold)"%(i, beds1_dict[i].shape[0]))
beds1_dict_df = pd.DataFrame()
beds1_dict_df['key'] = keys 
beds1_dict_df['nr_peaks_post_LT_0.5_filter'] = nr_peaks

#print(beds_dict.keys())
#keys are ['cardiac_muscle_cell', 'bulk_heart_LV_tissue']
#print(beds_dict.values())
#value are all the columns ie chrom, start, end, overall_peak_score, shape_score etc.

#print(beds_dict_extra.keys())
#keys are ['smooth_muscle', 'macrophage', 'atrial_cm', 'endothelial', 'fibroblast', 'lymphocyte', 'adipocyte', 'ventricular_cm']
#print(beds_dict_extra.values())

#should I compare this to the nr peaks before / does there need to be a checkpoint here of how many peaks are expected per file / per cell_type?

total_peaks_across_classifiers_post_LT_5_filter = pd.concat([beds1_dict_df, beds2_dict_df]).sort_values(by='nr_peaks_post_LT_0.5_filter', ascending=False)
total_peaks_across_classifiers_post_LT_5_filter.to_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/total_peaks_across_classifiers_post_LT_5_filter.csv', index=False)

# In[ ]:


##Intersecting bed 

#So this is turning the dictionary into a dataframe again - with all the peak info per file except this time just chrom, start, end, cell_type

dfs2 = []
for i in beds2_dict.keys():
    df = beds2_dict[i][[0,1,2,'cell_type']]
    print(df.head())
    dfs2.append(df) 
    df.rename({0:'Chr', 1:'Start', 2:'End'}, axis=1, inplace=True)
    #df.to_csv("/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/filtered_beds_post_L-tron_220524/%s_hg38_L-tron_filt.bed"%(i), sep='\t', index=False, header=False)
dfs2 = pd.concat(dfs2)
dfs2.rename({'Chr':'chrom', 'Start':'start', 'End':'end'}, axis=1, inplace=True)
dfs2.head()
dfs2.shape[0]
#2281676

dfs1 = []
for i in beds1_dict.keys():
    df = beds1_dict[i][['chrom','start', 'end', 'cell_type']]
    dfs1.append(df) 
    #df.to_csv("/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/filtered_beds_post_L-tron_220524/%s_hg38_L-tron_filt.bed"%(i), sep='\t', index=False, header=False)
dfs1 = pd.concat(dfs1)
dfs1.head()
dfs1.shape[0]
#406961

#removing datasets where there were so few peaks... ie less than 1000 peaks and so this dataset with 1 peak
dfs1 = dfs1[~(dfs1['cell_type'].isin(['cardiac_muscle_cell_CHIP_HISTONE_H3K4me1_paired']))]
dfs1.shape[0]
#406960

both_dfs = pd.concat([dfs1, dfs2])
both_dfs.shape[0] 
#2688636
both_dfs.head()
print(both_dfs['cell_type'].unique())
['cardiac_muscle_cell_CHIP_HISTONE_H3K4me3_single'
 'cardiac_muscle_cell_CHIP_HISTONE_H3K27me3_single'
 'cardiac_muscle_cell_DNASE_single'
 'cardiac_muscle_cell_CHIP_HISTONE_H3K27ac_single'
 'cardiac_muscle_cell_CHIP_TF_single' 'bulk_heart_LV_tissue_ATAC_paired'
 'bulk_heart_LV_tissue_DNASE_paired' 'bulk_heart_LV_tissue_CHIP_TF_paired'
 'bulk_heart_LV_tissue_CHIP_HISTONE_H3K4me3_single'
 'bulk_heart_LV_tissue_CHIP_HISTONE_H3K27me3_single'
 'bulk_heart_LV_tissue_CHIP_HISTONE_H3K4me1_single'
 'bulk_heart_LV_tissue_CHIP_HISTONE_H3K27ac_single'
 'bulk_heart_LV_tissue_CHIP_TF_single' 'cardiac_muscle_cell_DNASE_paired'
 'Pericyte_General_1_ATAC' 'Cardiac_Pericyte_3_ATAC'
 'Fetal_Endocardial_ATAC' 'Mast_ATAC' 'Macrophage_ATAC' 'macrophage_ATAC'
 'V_Cardiomyocyte_ATAC' 'Fetal_A_Cardiomyocyte_ATAC'
 'Pericyte_Muscularis_ATAC' 'Fetal_Skeletal_Myocyte_2_ATAC'
 'Fibroblast_ATAC' 'fibroblast_ATAC' 'Endocardial_ATAC'
 'Pericyte_General_4_ATAC' 'Fetal_V_Cardiomyocyte_ATAC'
 'Fetal_Mesothelial_ATAC' 'bulk_heart_LV_tissue_ATAC2'
 'cardiac_muscle_cell_DNASE2' 'Mesothelial_ATAC' 'Pericyte_ATAC'
 'Fetal_Skeletal_Myocyte_1_ATAC' 'Cardiac_Pericyte_2_ATAC'
 'Endothelial_Myocardial_ATAC' 'SMC2_ATAC' 'smooth_muscle_ATAC'
 'atrial_cm_ATAC' 'Fetal_Cardiac_Fibroblast_ATAC'
 'Type_I_Skeletal_Myocyte_ATAC' 'SMC1_ATAC' 'bulk_heart_LV_tissue_DNASE2'
 'Vasc_Sm_Muscle_1_ATAC' 'Pericyte_General_2_ATAC' 'Vasc_Sm_Muscle_2_ATAC'
 'Fetal_Skeletal_Myocyte_3_ATAC' 'Ery_Don002_hg38_ATAC'
 'Cardiac_Pericyte_1_ATAC' 'Adipocyte_ATAC' 'adipocyte_ATAC'
 'ventricular_cm_ATAC' 'Ery_Don003_hg38_ATAC'
 'Type_II_Skeletal_Myocyte_ATAC' 'Mast.1_ATAC' 'Endothelial_ATAC'
 'endothelial_ATAC' 'A_Cardiomyocyte_ATAC' 'Cardiac_Fibroblast_ATAC'
 'Ery_Don001_hg38_ATAC' 'lymphocyte_ATAC' 'Pericyte_General_3_ATAC'
 'cardiac_muscle_cell_DNASE1' 'Cardiac_Pericyte_4_ATAC']
#65 cell_types going forwards

both_dfs[both_dfs['chrom'].str.startswith('GL')]
#empty 
both_dfs[both_dfs['chrom'].str.startswith('KI')]
#empty

# both_dfs = both_dfs[~(both_dfs['chrom'].str.startswith('KI'))]
# both_dfs = both_dfs[~(both_dfs['chrom'].str.startswith('GL'))]

df = pybedtools.BedTool.from_dataframe(both_dfs)
df = df.sort()

df = df.merge(c=[4], o=['collapse'])
df = df.to_dataframe()
print(df.head())
print(df.shape[0])
#321280

if os.path.isdir("tmp"):
    shutil.rmtree("tmp") 

new_name = []
for i in df['name'].tolist():
    new_name.append(",".join(list(set(i.split(',')))))
#print('new_name')
#print(new_name)
df['name'] = new_name

pct_cellTypes = []
for i in df['name'].str.split(','):
    print(i)
    pct_cellTypes.append(round(len(i)/65, 2)) #denominator is 65 because print(len(both_dfs['cell_type'].unique())) = 65
    print(round(len(i)))
    print(round(len(i)/65, 2))
#print('pct_cellTypes')
#print(pct_cellTypes)
df['pct'] = pct_cellTypes

df.to_csv("/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/intersect.bed", sep='\t', index=False, header=False)

###########################################################

import pandas as pd
import glob
import os
import shutil
import pybedtools

df = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/intersect.bed', sep='\t', header=None)
df['nr_classifiers'] = round(df[4]*65)
df.rename({0:'chr', 1:'start', 2:'end'}, axis=1, inplace=True)
df.drop([3,4], axis=1, inplace=True)
df['nr_classifiers'] = df['nr_classifiers'].astype(int)
new_df = pd.DataFrame(df['nr_classifiers'].value_counts()).reset_index()
new_df.rename({'index':'nr_classifiers', 'nr_classifiers':'nr_classifier_counts'}, axis=1, inplace=True)
new_df.to_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/intersect_wrangled_piechart_final.csv', index=False)
# celltype = []
# counts = []
# for ty in both_dfs['cell_type'].unique():
#     celltype.append(ty)
#     counts.append(len(df[df[3].str.contains(ty)]))
# counts_df = pd.DataFrame()
# counts_df['cell_tissue'] = celltype #24
# counts_df['count'] = counts # divide ocunt by 2485794 (and then multiply by 100) to get % of how often that celltype falls within a peak
# counts_df.to_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/intersect_wrangled_piechart_final.csv', index=False)

# df = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/intersect.bed', sep='\t', header=None)
# df[5] = df[4]*100
# df = df.sort_values(by=[5])
# df[5] = df[5].astype(int)
# # df[5] = df[5].astype(str) + '%'
# counts_df = pd.DataFrame(df[5].value_counts()).reset_index()
# counts_df.rename({'index':5, 5:'count'}, axis=1, inplace=True)
# merged = df.merge(counts_df, on=5, how='left')
# merged = merged[[5, 'count']].drop_duplicates()
# merged.to_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/intersect_wrangled_histogram.csv', index=False)

###################################

# Load ggplot2
library(ggplot2)
library(dplyr)

# Create Data
# data <- read.csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/intersect_wrangled_histogram.csv')

# data$X5 <- factor(data$X5, levels = data$X5)

# p<-ggplot(data=data, aes(x=X5, y=count)) +
#   geom_bar(stat="identity", fill="lightblue", color="darkblue") +
#   scale_y_continuous(breaks = round(seq(min(data$count), max(data$count), by = 4000),1)) +
#   xlab('Percentage of classifiers sharing a peak', size=9) +
#   ylab('Number of peaks', size=9) +
#   theme(axis.text.x = element_text(angle = 90, size=7, vjust = 0.5, hjust=1), axis.text.y = element_text(size=7), axis.title.x = element_text(size=9), axis.title.y = element_text(size=9))

# ggsave("/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/intersect_histogram.png")

# library(scales)
# library(ggplot2)
# library(tidyverse)
# library(ggrepel)

# data <- read.csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/intersect_wrangled_piechart_final.csv')

# blank_theme <- theme_minimal()+
#   theme(
#   axis.title.x = element_blank(),
#   axis.title.y = element_blank(),
#   #panel.border = element_blank(),
#   #panel.grid=element_text(size=5),
#   #axis.ticks = element_blank(),
#   plot.title=element_text(size=14, face="bold")
#   )

# png(file="/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/intersect_piechart.png",width=3300 ,height=2000,res=150)

# bp<- ggplot(data, aes(x="", y=count, fill=cell_tissue))+
# geom_bar(width = 1, stat = "identity")
# bp

# pie <- bp + coord_polar("y", start=0)
# pie

# pie + scale_colour_manual(values = c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "black", "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise", "green1", "yellow4", "yellow3", "darkorange4", "brown"))

# pie + blank_theme

# pie + theme(axis.text.x=element_blank())

# pie + geom_text(aes(label = final_perc), position = position_stack(vjust = 0.5), size=2)

# dev.off()

library(scales)
library(ggplot2)
library(tidyverse)
library(ggrepel)

data <- read.csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/intersect_wrangled_piechart_final.csv')

png(file="/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/intersect_piechart.png",width=3300 ,height=2000,res=150)

pie <- ggplot(df, aes(x = "", y=nr_classifier_counts, fill = factor(nr_classifiers))) +
geom_bar(width = 1, stat = "identity") +
theme(axis.line = element_blank(), plot.title = element_text(hjust=0.5)) +
labs(fill="nr_classifiers", x=NULL, y=NULL, title="Proportion of classifiers with a given peak" +
coord_polar("y", start=0)

pie

dev.off() 


###########################################################

##Building dataset for ML
df = pd.read_csv("/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/intersect.bed", sep='\t', header=False)

ds = pd.DataFrame(df['chrom']+'-'+df['start'].astype(str)+'-'+df['end'].astype(str), columns=['region'])
ds.shape[0]
#321280

print(beds1_dict.keys())
print(beds2_dict.keys())
all_beds_dict_keys = list(beds1_dict.keys()) + list(beds2_dict.keys())
print(len(all_beds_dict_keys))
#66
all_beds_dict_keys.remove("cardiac_muscle_cell_CHIP_HISTONE_H3K4me1_paired")
print(len(all_beds_dict_keys))
#65

# It takes ~15 mins
from tqdm import tqdm

dfs = []
for idx, i in tqdm(zip(ds['region'].tolist(), df['name'].tolist())):
    print(i)
    values = i.split(',')
    print(values)
    df_tmp = pd.DataFrame(data=[[1.0]*len(values)], columns=values, index=[idx])
    print(df_tmp)
    missing_col = set(all_beds_dict_keys).difference(set(values))
    print(missing_col)
    if len(missing_col) > 0:
        for c in missing_col:
            df_tmp[c] = 0.0
    df_tmp = df_tmp[all_beds_dict_keys]
    print(df_tmp)
    df_tmp = df_tmp.loc[:,~df_tmp.columns.duplicated()]
    print(df_tmp)
    dfs.append(df_tmp)


# In[ ]:

df_concat = pd.concat(dfs, axis=0, ignore_index=False)
#321280

#here lets remove all those regions that we cant fetch sequences for ie chrom starting with either GL or KL

df_concat.to_csv("/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/concat_ML_hot.csv", sep='\t')

