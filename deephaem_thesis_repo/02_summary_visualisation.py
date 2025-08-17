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

# from IPython.core.display import display, HTML
# display(HTML('<style>.container {width: 90% !important; }</style>'))#

import pandas as pd
import numpy as np
import glob
import anndata
import scanpy as sc

##df_intersect

df_intersect = pd.read_csv("/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/intersect.bed", sep='\t', header=None)
print(df_intersect.head())
print(df_intersect[df_intersect[4] <= 0.1].shape)
print(df_intersect[df_intersect[4] >= 1].shape)

##df_concat

df_concat = pd.read_csv("/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/concat_ML_hot.csv", sep='\t', index_col=0)
print(df_concat.head())
print(df_concat.columns.tolist())
print(df_concat.head())

##Object creation

obj = anndata.AnnData(df_concat)
#print(obj.X) #The initial data we passed are accessible as a sparse matrix using adata.X
print(obj.obs_names)
# Index(['chr1-794319-795830', 'chr1-816174-816505', 'chr1-816640-818479',
#        'chr1-818565-819327', 'chr1-825291-825361', 'chr1-826401-828490',
#        'chr1-828850-829276', 'chr1-829380-829673', 'chr1-829944-830502',
#        'chr1-830719-831176',
#        ...
#        'chrY-26319013-26320651', 'chrY-26352859-26353378',
#        'chrY-26359669-26360992', 'chrY-26363313-26363851',
#        'chrY-26364735-26365295', 'chrY-26405101-26405510',
#        'chrY-26408684-26409642', 'chrY-26524968-26525472',
#        'chrY-26562949-26563399', 'chrY-26589728-26590646'],
#       dtype='object', length=321280)
print(obj.var_names)
# Index(['cardiac_muscle_cell_CHIP_HISTONE_H3K4me3_single',
#        'cardiac_muscle_cell_CHIP_HISTONE_H3K27me3_single',
#        'cardiac_muscle_cell_DNASE_single',
#        'cardiac_muscle_cell_CHIP_HISTONE_H3K27ac_single',
#        'cardiac_muscle_cell_CHIP_TF_single',
#        'bulk_heart_LV_tissue_ATAC_paired', 'bulk_heart_LV_tissue_DNASE_paired',
#        'bulk_heart_LV_tissue_CHIP_TF_paired',
#        'bulk_heart_LV_tissue_CHIP_HISTONE_H3K4me3_single',
#        'bulk_heart_LV_tissue_CHIP_HISTONE_H3K27me3_single',
#        'bulk_heart_LV_tissue_CHIP_HISTONE_H3K4me1_single',
#        'bulk_heart_LV_tissue_CHIP_HISTONE_H3K27ac_single',
#        'bulk_heart_LV_tissue_CHIP_TF_single',
#        'cardiac_muscle_cell_DNASE_paired', 'Pericyte_General_1_ATAC',
#        'Cardiac_Pericyte_3_ATAC', 'Fetal_Endocardial_ATAC', 'Mast_ATAC',
#        'Macrophage_ATAC', 'macrophage_ATAC', 'V_Cardiomyocyte_ATAC',
#        'Fetal_A_Cardiomyocyte_ATAC', 'Pericyte_Muscularis_ATAC',
#        'Fetal_Skeletal_Myocyte_2_ATAC', 'Fibroblast_ATAC', 'fibroblast_ATAC',
#        'Endocardial_ATAC', 'Pericyte_General_4_ATAC',
#        'Fetal_V_Cardiomyocyte_ATAC', 'Fetal_Mesothelial_ATAC',
#        'bulk_heart_LV_tissue_ATAC2', 'cardiac_muscle_cell_DNASE2',
#        'Mesothelial_ATAC', 'Pericyte_ATAC', 'Fetal_Skeletal_Myocyte_1_ATAC',
#        'Cardiac_Pericyte_2_ATAC', 'Endothelial_Myocardial_ATAC', 'SMC2_ATAC',
#        'smooth_muscle_ATAC', 'atrial_cm_ATAC', 'Fetal_Cardiac_Fibroblast_ATAC',
#        'Type_I_Skeletal_Myocyte_ATAC', 'SMC1_ATAC',
#        'bulk_heart_LV_tissue_DNASE2', 'Vasc_Sm_Muscle_1_ATAC',
#        'Pericyte_General_2_ATAC', 'Vasc_Sm_Muscle_2_ATAC',
#        'Fetal_Skeletal_Myocyte_3_ATAC', 'Ery_Don002_hg38_ATAC',
#        'Cardiac_Pericyte_1_ATAC', 'Adipocyte_ATAC', 'adipocyte_ATAC',
#        'ventricular_cm_ATAC', 'Ery_Don003_hg38_ATAC',
#        'Type_II_Skeletal_Myocyte_ATAC', 'Mast.1_ATAC', 'Endothelial_ATAC',
#        'endothelial_ATAC', 'A_Cardiomyocyte_ATAC', 'Cardiac_Fibroblast_ATAC',
#        'Ery_Don001_hg38_ATAC', 'lymphocyte_ATAC', 'Pericyte_General_3_ATAC',
#        'cardiac_muscle_cell_DNASE1', 'Cardiac_Pericyte_4_ATAC'],
#       dtype='object')
obj.obs['pct'] = df_intersect[4].tolist()
obj
# AnnData object with n_obs × n_vars = 321280 × 65
#     obs: 'pct'
#print(obj.obs)

import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors

# Inital setting for plot size
from matplotlib import rcParams
FIGSIZE=(10,10)
rcParams['figure.figsize']=FIGSIZE

import warnings
warnings.filterwarnings('ignore')

sc.tl.pca(obj, svd_solver='arpack')
sc.pl.pca_variance_ratio(obj, log=True, save='.png')

sc.pp.neighbors(obj, n_neighbors=10)
#I tried 2 n_neighbours. #Compute a neighborhood graph of observations

sc.tl.umap(obj, n_components=3)
#Embed the neighborhood graph using UMAP. UMAP (Uniform Manifold Approximation and Projection) is a manifold learning technique suitable for visualizing high-dimensional data. Besides tending to be faster than tSNE, it optimizes the embedding such that it best reflects the topology of the data, which we represent throughout Scanpy using a neighborhood graph.

sc.tl.leiden(obj, resolution=0.4) 
#I tried 0.25. #Cluster cells into subgroups using the Leiden algorithm.
#244332 x 24

obj.write('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/cm.h5ad', compression='gzip')

print(obj.obs['leiden'])

df_intersect['UMAP1'] = obj.obsm['X_umap'][:,0]
df_intersect['UMAP2'] = obj.obsm['X_umap'][:,1]
df_intersect['UMAP3'] = obj.obsm['X_umap'][:,2]
df_intersect['clusters'] = obj.obs['leiden'].tolist()
df_intersect.columns = ['chromosome', 'start', 'end', 'names', 'pct', 'UMAP1', 'UMAP2', 'UMAP3', 'clusters']
df_intersect.to_csv("/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/intersect_MDV.bed", sep='\t', index=False, header=True)

##Visualisation

sc.pl.umap(obj, save='.png')
sc.pl.umap(obj, color=['pct'], palette= 'crest', projection='3d', save='.png')

obj = sc.read_h5ad("/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/cm.h5ad")
sc.pl.umap(obj, color=['leiden'], projection='3d', palette='crest', save='.png') #had the following option as well: , projection='3d'

# for c in df_concat.columns.tolist():
#     sc.pl.umap(obj, color=c, projection='2d', vmin=0.0, vmax=1.0, components='all')

import matplotlib.pyplot as plt

fig, (ax1) = plt.subplots(1, 1, figsize=(40,40), gridspec_kw={'wspace':0.9}, layout="constrained")
sc.pl.matrixplot(obj, obj.var_names, groupby='leiden', show=False, cmap='plasma', ax=ax1)
ax = sc.pl.heatmap(obj, obj.var_names, groupby='leiden', cmap='plasma', dendrogram=False)
fig.savefig('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/02_summary_vis_plot.png')
