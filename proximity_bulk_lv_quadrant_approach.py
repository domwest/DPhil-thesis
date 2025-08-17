from cmapPy.pandasGEXpress.parse_gct import parse
from cmapPy.pandasGEXpress.write_gct import write
import pandas as pd
import gffpandas.gffpandas as gffpd
import os
import numpy as np

#testing my approach on genes that have already been marked valid during pathogenicity assessment of their involvement in HCM
#will run the code from the very beginning and summarise the relevant parts here... 
#approach will differ a bit in the sense that the loci will not be determined by 500 mb up and downstream of the variant; instead it will use the midpoint of the gene of interest (here the ones with evidence for HCM) to define these loci

#took info in the order of: generate_fold_changes_using_median.py (INCLUDE), from_fold_changes_to_top_expr_genes.py (EXCLUDE), top_expr_genes_lv_analysis.py (INCLUDE), labeled_scatter_plot_script.R

##GETTING EXPR DATA FROM GTEX

#can do:
#data.data_df shows the actual median TPM values
#data.col_metadata_df explains what the columns represent (ie which tissues are represented by which cols)
#data.row_metadata_df explains which gene symbols are depicted by the ensembl ID in the index of each row 

#54 cols so 54 tissues
#need to use these tissue medians to determine whether genes in my list are heart-specific or heart-enhanced genes just like this paper did:
#Integrative Analysis Revealing Human Heart-Specific Genes and Consolidating Heart-Related Phenotypes: https://www.frontiersin.org/articles/10.3389/fgene.2020.00777/full
#Their approach for this is as follows:
#***Human protein-coding and non-coding genes were filtered by general thresholds of at least 20% of samples having TPM > 0.1 and median TPM > 0.5 in the heart AA and LV
#...Evidence has yet to be explored to link between genes expressed exclusively or abundantly in the heart (e.g., heart-specific genes and heart-enhanced genes described in the current study) and heart traits or diseases...
#Heart-specific genes were defined as genes with fold changes of the median expression levels (FCMs) higher than 5.0 in the heart versus all other tissues. 
#Heart-enhanced genes were defined as genes having FCMs higher than 5.0 in the heart versus all other tissues, with an exception of 1 or 2 other tissues, and a higher median expression in the heart in the case of exceptions.
#data = parse('/home/domwest/Documents/oxford_dphil/meta_analysis_variants/exploring_ensembl/downloaded_GTEx_datasets/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct')
data = parse('GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct')

#get df with all tissues in desirable format first:
all_tissues = data.data_df
all_tissues['gene_symbol'] = data.row_metadata_df
all_tissues = all_tissues.reset_index()
#process before merging
all_tissues.rename(columns={"rid": "gene_id"}, inplace=True)
all_tissues['gene_id'] = all_tissues['gene_id'].str.split('.').str[0]

############################################################################################################################################################################################################

##GENERATING MEDIAN_OF_MEDIANS AND THEN FOLD_CHANGES, SO A BIT OF WRANGLING BEFORE DOING THIS

#AA:

#***Not going to do this for now... initial_merge_aa = all_tissues[all_tissues['gene_symbol'].isin(genes_thresh_aa)]
#***ERROR HAPPENING: where a gene symbol might have > 0.2, this exact pair combination  of gene_id and gene_symbol might be different to that in genes df... like the gene_symbol is the same but the gene_id captured in genes_df happens to be a combo that had 0.00 median expr in the all_tissues df
initial_merge_aa = all_tissues
long_initial_merge_aa=pd.melt(initial_merge_aa,id_vars=['gene_id', 'gene_symbol'],var_name='tissue', value_name='median_expr')
long_initial_merge_aa = long_initial_merge_aa[(long_initial_merge_aa['tissue']=='Heart - Atrial Appendage')]
#for merge 2
#***Not going to do this for now... wide_all_tissues_aa = all_tissues[all_tissues['gene_symbol'].isin(genes_thresh_aa)]
wide_all_tissues_aa = all_tissues
cols = wide_all_tissues_aa.columns.tolist()
cols.remove('gene_id')
cols.remove('gene_symbol')
cols.remove('Heart - Atrial Appendage')
#need to remove LV as well since LV should not be in  denominator when calculating median_of_medians for atrial appendage and vice versa
cols.remove('Heart - Left Ventricle')
wide_all_tissues_aa['median_of_medians'] = wide_all_tissues_aa[cols].median(axis=1) #wide_all_tissues_aa['mean_of_medians'] = wide_all_tissues_aa[cols].mean(axis=1) 
#calc fold change by doing gene / mean_of_medians for both heart aa and heart lv 
wide_all_tissues_aa['fold_change_aa'] = wide_all_tissues_aa['Heart - Atrial Appendage']/wide_all_tissues_aa['median_of_medians']

wide_all_tissues_aa_copy = wide_all_tissues_aa
wide_all_tissues_aa_copy.rename({'Heart - Atrial Appendage':'Heart_AA_median_TPM'}, axis=1, inplace=True)
wide_all_tissues_aa_copy = wide_all_tissues_aa_copy[['gene_id', 'Heart_AA_median_TPM', 'fold_change_aa']]
wide_all_tissues_aa_copy['tissue'] = 'Heart_AA'

wide_all_tissues_aa_copy.to_csv("fold_change_aa_medians_used.csv", index=False)

#LV:

#***Not going to do this for now... initial_merge_aa = all_tissues[all_tissues['gene_symbol'].isin(genes_thresh_aa)]
#***ERROR HAPPENING: where a gene symbol might have > 0.2, this exact pair combination  of gene_id and gene_symbol might be different to that in genes df... like the gene_symbol is the same but the gene_id captured in genes_df happens to be a combo that had 0.00 median expr in the all_tissues df
initial_merge_lv = all_tissues
long_initial_merge_lv=pd.melt(initial_merge_lv,id_vars=['gene_id', 'gene_symbol'],var_name='tissue', value_name='median_expr')
long_initial_merge_lv = long_initial_merge_lv[(long_initial_merge_lv['tissue']=='Heart - Left Ventricle')]
#for merge 2
#***Not going to do this for now... wide_all_tissues_aa = all_tissues[all_tissues['gene_symbol'].isin(genes_thresh_aa)]
wide_all_tissues_lv = all_tissues
cols = wide_all_tissues_lv.columns.tolist()
cols.remove('gene_id')
cols.remove('gene_symbol')
cols.remove('Heart - Left Ventricle')
#need to remove AA as well since AA should not be in  denominator when calculating median_of_medians for left ventricle and vice versa
cols.remove('Heart - Atrial Appendage')
wide_all_tissues_lv['median_of_medians'] = wide_all_tissues_lv[cols].median(axis=1) 
#calc fold change by doing gene / mean_of_medians for both heart aa and heart lv 
wide_all_tissues_lv['fold_change_lv'] = wide_all_tissues_lv['Heart - Left Ventricle']/wide_all_tissues_lv['median_of_medians']

wide_all_tissues_lv_copy = wide_all_tissues_lv
wide_all_tissues_lv_copy.rename({'Heart - Left Ventricle':'Heart_LV_median_TPM'}, axis=1, inplace=True)
wide_all_tissues_lv_copy = wide_all_tissues_lv_copy[['gene_id', 'Heart_LV_median_TPM', 'fold_change_lv']]
wide_all_tissues_lv_copy['tissue'] = 'Heart_LV'

wide_all_tissues_lv_copy.to_csv("fold_change_lv_medians_used.csv", index=False)

wide_all_tissues_lv_copy = pd.read_csv("fold_change_lv_medians_used.csv")

# >>> wide_all_tissues_lv_copy['Heart_LV_median_TPM'].quantile(0.95)
# 17.33
# >>> wide_all_tissues_lv_copy['Heart_LV_median_TPM'].describe()
# count    56202.000000
# mean        17.054625
# std        785.932716
# min          0.000000
# 25%          0.000000
# 50%          0.000000
# 75%          0.641400
# max      82390.000000

>>> len(wide_all_tissues_lv_copy['gene_id'].unique())
56202
>>> len(wide_all_tissues_lv_copy['gene_id'][wide_all_tissues_lv_copy['Heart_LV_median_TPM']!=0].unique())
25202
>>> wide_all_tissues_lv_copy_nonzero_expr = wide_all_tissues_lv_copy[wide_all_tissues_lv_copy['Heart_LV_median_TPM']!=0]

>>> wide_all_tissues_lv_copy_nonzero_expr['Heart_LV_median_TPM'].quantile(0.95)
37.52899999999994
>>> wide_all_tissues_lv_copy_nonzero_expr['Heart_LV_median_TPM'].describe()
count    25202.000000
mean        38.032857
std       1173.336789
min          0.000816
25%          0.137025
50%          1.039500
75%          6.284000
max      82390.000000
Name: Heart_LV_median_TPM, dtype: float64

wide_all_tissues_lv_copy_nonzero_expr.to_csv("fold_change_lv_medians_used_nonzero_median_tpm.csv", index=False)

wide_all_tissues_lv_copy = pd.read_csv("fold_change_lv_medians_used.csv")

bulk_genes_of_interest = ['HSPB7', 'HEY2', 'ACADS', 'FGFBP2', 'PROB1', 'PLN', 'PSMA1', 'CKM', 'PLPP7', 'DNAJA4', 'PPP1R14C', 'GPBP1', 'GNB1', 'PLCD3', 'MFHAS1', 'MAIP1', 'NDUFB9', 'VIT', 'TNFSF10', 'DTNA', 'CHCHD10', 'ALPK2', 'NEXN', 'TBX3', 'PCNX1', 'ALPK3', 'MLIP', 'COL4A2']

id_symbol = pd.read_csv('ids_symbols.csv')

id_symbol[id_symbol['symbol'].isin(bulk_genes_of_interest)].shape[0]
#25
len(bulk_genes_of_interest)
#28
filt_id_symbol = id_symbol[id_symbol['symbol'].isin(bulk_genes_of_interest)]
captured = filt_id_symbol['symbol'].unique().tolist()
import numpy as np
np.setdiff1d(bulk_genes_of_interest, captured)
#array(['GNB1', 'HEY2', 'PROB1'], dtype='<U8')
captured_gene_ids = filt_id_symbol['gene'].unique().tolist()
to_add = ['ENSG00000078369', 'ENSG00000135547', 'ENSG00000228672']
all_gene_ids = captured_gene_ids + to_add
len(all_gene_ids)
#28

wide_all_tissues_lv_copy_priorit = wide_all_tissues_lv_copy[wide_all_tissues_lv_copy['gene_id'].isin(all_gene_ids)]
id_symbol.rename({'gene':'gene_id'}, axis=1, inplace=True)
wide_all_tissues_lv_copy_priorit = wide_all_tissues_lv_copy_priorit.merge(id_symbol, how='left', on='gene_id')
wide_all_tissues_lv_copy_priorit['symbol'][wide_all_tissues_lv_copy_priorit.index==0] = 'GNB1'
wide_all_tissues_lv_copy_priorit['symbol'][wide_all_tissues_lv_copy_priorit.index==8] = 'PROB1'
wide_all_tissues_lv_copy_priorit['symbol'][wide_all_tissues_lv_copy_priorit.index==11] = 'HEY2'

wide_all_tissues_lv_copy_priorit.to_csv("fold_change_lv_medians_used_priorit_bulk_genes_for_top_53_vars.csv", index=False)


# >>> wide_all_tissues_lv_copy['fold_change_lv'].quantile(0.95)
# 2.67707081
# >>> wide_all_tissues_lv_copy['fold_change_lv'].describe()
# count    2.825800e+04
# mean              inf
# std               NaN
# min      0.000000e+00
# 25%      1.990022e-01
# 50%      3.284491e-01
# 75%      5.321955e-01
# max               inf
# Name: fold_change_lv, dtype: float64

# wide_all_tissues_lv_copy_f =wide_all_tissues_lv_copy[wide_all_tissues_lv_copy['Heart_LV_median_TPM']<=17.33]
# wide_all_tissues_lv_copy_f.to_csv("fold_change_lv_medians_used_only_within_95_quantiles.csv", index=False)

# wide_all_tissues_lv_copy_f2 =wide_all_tissues_lv_copy[wide_all_tissues_lv_copy['fold_change_lv']<=2.67707081]
# wide_all_tissues_lv_copy_f2.to_csv("fold_change_lv_medians_used_only_within_95_quantiles_fc.csv", index=False)


#

.libPaths('/well/PROCARDIS/domwest/R')

library(tidyverse)
#library(showtext)
#library(ggtext)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)

# -rw-r--r-- 1 wiq135 watkins  71286 Aug  9 10:06 DNASE_endothelial_cell_of_umbilical_vein_newborn_1_huvec.csv
# -rw-r--r-- 1 wiq135 watkins  71888 Aug  9 10:06 DNASE_endothelial_cell_of_umbilical_vein_newborn_2_huvec.csv

long_prediction_df <- read.csv('fold_change_lv_medians_used.csv')

p <- ggplot(long_prediction_df, aes(x=tissue, y=fold_change_lv)) +
geom_boxplot(outlier.colour="black", outlier.shape=1, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
scale_fill_brewer(palette="Dark2")+
scale_y_continuous(n.breaks = 10) +
theme(text=element_text(size=18), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = 'none') +
labs(x = "Fold change distribution across genes in Heart Left Ventricle") 

ggsave(filename="heart_lv_fold_change_distribution_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

p <- ggplot(long_prediction_df, aes(x=tissue, y=Heart_LV_median_TPM)) +
geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
scale_fill_brewer(palette="Dark2")+
scale_y_continuous(n.breaks = 10) +
theme(text=element_text(size=18), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = 'none') +
labs(x = "Median TPM distribution across genes in Heart Left Ventricle") 

ggsave(filename="heart_lv_median_tpm_distribution_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

long_prediction_df <- read.csv('fold_change_lv_medians_used_nonzero_median_tpm.csv')

p <- ggplot(long_prediction_df, aes(x=tissue, y=Heart_LV_median_TPM)) +
geom_boxplot(outlier.colour="black", outlier.shape=1, outlier.size=2) +
stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
scale_fill_brewer(palette="Dark2")+
scale_y_continuous(n.breaks = 10) +
theme(text=element_text(size=18), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = 'none') +
labs(x = "Median TPM distribution across genes in Heart Left Ventricle") 

ggsave(filename="heart_lv_median_tpm_distribution_boxplot_used_nonzero_median_tpm.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

long_prediction_df <- read.csv('fold_change_lv_medians_used_priorit_bulk_genes_for_top_53_vars.csv')

ggplot(long_prediction_df, aes(x=Heart_LV_median_TPM, y=fold_change_lv)) +
# geom_hline(yintercept = 1.3, color="#99CCCC", size =0.5) +
# geom_vline(xintercept = 1.0395, color="#99CCCC", size =0.5) +
geom_point(colour = "#336666") + 
geom_text_repel(aes(label = symbol)) +
labs(title = "Fold change vs median TPM for prioritised bulk genes identified across case studies",
   subtitle = NULL,
   tag = NULL, 
   x = "Median TPM",
   y= "Fold change values compared to other tissues",
   color = NULL) +
theme(text=element_text(size=18))
ggsave(filename="scatter_fold_change_lv_medians_used_priorit_bulk_genes_for_top_53_vars.tiff",width=18, height=15, dpi = 300)

#

# long_prediction_df <- read.csv('fold_change_lv_medians_used_only_within_95_quantiles.csv')

# p <- ggplot(long_prediction_df, aes(x=tissue, y=Heart_LV_median_TPM)) +
# geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
# stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
# scale_fill_brewer(palette="Dark2")+
# scale_y_continuous(n.breaks = 10) +
# theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = 'none') +
# labs(x = "Median TPM distribution across genes in Heart Left Ventricle") 

# ggsave(filename="heart_lv_median_tpm_distribution_boxplot_only_within_95_quantiles.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

#

# long_prediction_df <- read.csv('fold_change_aa_medians_used.csv')

# p <- ggplot(long_prediction_df, aes(x=tissue, y=fold_change_aa)) +
# geom_boxplot(outlier.colour="black", outlier.shape=1, outlier.size=2) +
# stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
# scale_fill_brewer(palette="Dark2")+
# scale_y_continuous(n.breaks = 10) +
# theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = 'none') +
# labs(x = "Fold change distribution across genes in Heart Atrial Appendage") 

# ggsave(filename="heart_aa_fold_change_distribution_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/

# p <- ggplot(long_prediction_df, aes(x=tissue, y=Heart_AA_median_TPM)) +
# geom_boxplot(outlier.colour="black", outlier.shape=2, outlier.size=2) +
# stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
# scale_fill_brewer(palette="Dark2")+
# scale_y_continuous(n.breaks = 10) +
# theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = 'none') +
# labs(x = "Median TPM distribution across genes in Heart Atrial Appendage") 

# ggsave(filename="heart_aa_median_tpm_distribution_boxplot.tiff", width=18, height=15, dpi = 300) #was: scatter_plots/all/


############################################################################################################################################################################################################

##IDENTIFYING TOP GENES

#AA:

##for identifying top expressed genes
#we want to look at those above the 95% percentile
#but would also be good to include the genes showing a fold change of 'inf' since these genes were positive in the heart tissue and showing a median of 0 across all other tissues
aa_filtered = wide_all_tissues_aa[(wide_all_tissues_aa['fold_change_aa'].notnull()) & (wide_all_tissues_aa['fold_change_aa']!=0)]
aa_filtered = aa_filtered[~np.isinf(aa_filtered['fold_change_aa'])]
Q95 = aa_filtered['fold_change_aa'].quantile(0.95) #2.490151739120483
above_Q95 = aa_filtered[(aa_filtered['fold_change_aa'] > Q95)]
#to just get the unique gene symbol names from this df:
genes_above_Q95 = above_Q95['gene_id'].unique() #1295 genes
#now let's add the inf genes too
#won't work if already applied aa_filtered = aa_filtered[~np.isinf(aa_filtered['fold_change_aa'])] so need to run this again first: 
aa_filtered = wide_all_tissues_aa[(wide_all_tissues_aa['fold_change_aa'].notnull()) & (wide_all_tissues_aa['fold_change_aa']!=0)]
inf_genes = aa_filtered['gene_id'][np.isinf(aa_filtered['fold_change_aa'])].unique() #911 genes
#add the 2 lists
genes_above_Q95_plus_inf = list(genes_above_Q95) + list(inf_genes) #gives 2206 genes

#LV:

##for identifying top expressed genes
#we want to look at those above the 95% percentile
#but would also be good to include the genes showing a fold change of 'inf' since these genes were positive in the heart tissue and showing a median of 0 across all other tissues
lv_filtered = wide_all_tissues_lv[(wide_all_tissues_lv['fold_change_lv'].notnull()) & (wide_all_tissues_lv['fold_change_lv']!=0)]
lv_filtered = lv_filtered[~np.isinf(lv_filtered['fold_change_lv'])]
Q95 = lv_filtered['fold_change_lv'].quantile(0.95) #1.3
above_Q95 = lv_filtered[(lv_filtered['fold_change_lv'] > Q95)]
#to just get the unique gene symbol names from this df:
genes_above_Q95 = above_Q95['gene_id'].unique() #1232 genes
#now let's add the inf genes too
#won't work if already applied aa_filtered = aa_filtered[~np.isinf(aa_filtered['fold_change_aa'])] so must first run: 
lv_filtered = wide_all_tissues_lv[(wide_all_tissues_lv['fold_change_lv'].notnull()) & (wide_all_tissues_lv['fold_change_lv']!=0)]
inf_genes = lv_filtered['gene_id'][np.isinf(lv_filtered['fold_change_lv'])].unique() #571 genes
#add the 2 lists
genes_above_Q95_plus_inf = list(genes_above_Q95) + list(inf_genes) #gives 1803 genes

##############################################################################################################################################################################################################

##FETCHING THE LOCI GENE SET SO THAT WE CAN CHECK IF THESE ARE COMMON TO THE TOP EXPR GENE LIST 

#using all comers loci text file from Anuj to determine regions for gene investigation
#looking 500kbp (500000bp) both upstream and downstream of the variant as the search regions
#determine what the search regions would be if one were to enter these in the genome viewer GRCh37h19 (https://www.ncbi.nlm.nih.gov/genome/gdv/browser/genome/?id=GCF_000001405.13)
#one variant in the open target slide Anuj sent me was missing from the original gws_allcomer_loci.txt file so needed to add this in:
#'sno', 'chr', 'pos_hg19', 'rsid', 'gene_openTarget', 'Tadros', 'Harper'
#22  11   44,901,082       MYBPC3               1          NaN    NaN
#this directory must be where the output files (after running the ensembl api perl code) are located
directory = 'ensembl_api_output'
gene_ID_list = []
for filename in os.listdir(directory):
	if filename.endswith(".txt"): 
        	# print(os.path.join(directory, filename))
		#print(filename)
		df = pd.read_csv('ensembl_api_output/' + filename, names=["gene_id"], sep='\t')
		df = df[(df['gene_id'].str.contains("Gene stable id:")) | (df['gene_id'].str.contains("Gene external name:")) | (df['gene_id'].str.contains("Position:"))]
		df = pd.DataFrame({'gene_id':df['gene_id'].loc[df['gene_id'].str.contains("Gene stable id:")].values, 'gene_symbol':df['gene_id'].loc[df['gene_id'].str.contains("Gene external name:")].values, 'gene_position':df['gene_id'].loc[df['gene_id'].str.contains("Position:")].values})
		df['gene_id']=df['gene_id'].str.replace('Gene stable id: ','')
		df['gene_symbol']=df['gene_symbol'].str.replace('Gene external name: ','')
		df['gene_position']=df['gene_position'].str.replace('Position: ','')
		df['variant'] = filename.split('.')[0]
		df['variant']=df['variant'].str.replace('variant_','')
		gene_ID_list.append(df)
		#position already provided in  ann_df 'start' 'end'
genes = pd.concat(gene_ID_list) 
genes = genes.sort_values('variant')
genes['gene_start'] = ''
genes['gene_end'] = ''
genes['gene_position'] = genes['gene_position'].str.lstrip('-')
genes[['gene_start', 'gene_end']] = genes['gene_position'].str.split('-', expand=True)
genes.rename({'variant':'evidence_gene'}, axis=1, inplace=True)
genes.drop(['gene_position'], axis=1, inplace=True)

genes[['evidence_gene', 'gene_symbol']].groupby(['evidence_gene']).size()

search_regions = pd.read_csv('mtag_loci_search_regions.csv')
search_regions.rename(columns={"rsid": "evidence_gene"}, inplace=True)

genes = genes.merge(search_regions[['chrom', 'position', 'evidence_gene']], on=['evidence_gene'])
genes =  genes[['chrom', 'evidence_gene', 'position', 'gene_id', 'gene_symbol', 'gene_start', 'gene_end']]

############################################################################################################################################################################################################

#ACTUALLY CHECKING IF ANY OF THE GENE SETS IN THE LOCI ARE CONTAINED WITHIN THE TOP EXPR GENES

#AA:

loci_gene_set = genes['gene_id'].unique() #1977 genes
loci_gene_set_as_set = set(loci_gene_set) #1977
inter = loci_gene_set_as_set.intersection(genes_above_Q95_plus_inf) #100
inter_as_list = list(inter) 
aa_top_expr_genes = pd.DataFrame()
aa_top_expr_genes['top_expr_genes'] = inter_as_list
# ['ENSG00000233117', 'ENSG00000140403', 'ENSG00000158286', 'ENSG00000104881', 'ENSG00000173641', 'ENSG00000224163', 'ENSG00000259910', 'ENSG00000124772', 'ENSG00000184471', 'ENSG00000185739', 'ENSG00000177791', 'ENSG00000136383', 'ENSG00000115602', 'ENSG00000243711', 'ENSG00000184908', 'ENSG00000228672', 'ENSG00000228485', 'ENSG00000239964', 'ENSG00000140986', 'ENSG00000170417', 'ENSG00000172572', 'ENSG00000266946', 'ENSG00000267257', 'ENSG00000175946', 'ENSG00000136378', 'ENSG00000267979', 'ENSG00000137251', 'ENSG00000198626', 'ENSG00000271579', 'ENSG00000238166', 'ENSG00000256287', 'ENSG00000111981', 'ENSG00000235050', 'ENSG00000250479', 'ENSG00000156222', 'ENSG00000146147', 'ENSG00000198523', 'ENSG00000188488', 'ENSG00000255394', 'ENSG00000137198', 'ENSG00000198796', 'ENSG00000232456', 'ENSG00000249816', 'ENSG00000136574', 'ENSG00000217330', 'ENSG00000230408', 'ENSG00000104879', 'ENSG00000163492', 'ENSG00000088726', 'ENSG00000155657', 'ENSG00000267784', 'ENSG00000242853', 'ENSG00000162009', 'ENSG00000128591', 'ENSG00000185519', 'ENSG00000260403', 'ENSG00000249526', 'ENSG00000010310', 'ENSG00000125740', 'ENSG00000114854', 'ENSG00000205221', 'ENSG00000215533', 'ENSG00000258930', 'ENSG00000252396', 'ENSG00000244020', 'ENSG00000166317', 'ENSG00000153531', 'ENSG00000220773', 'ENSG00000184608', 'ENSG00000197321', 'ENSG00000077522', 'ENSG00000198729', 'ENSG00000128285', 'ENSG00000160539', 'ENSG00000234123', 'ENSG00000242902', 'ENSG00000242849', 'ENSG00000167971', 'ENSG00000134775', 'ENSG00000251576', 'ENSG00000261713', 'ENSG00000196358', 'ENSG00000117707', 'ENSG00000226005', 'ENSG00000137441', 'ENSG00000196136', 'ENSG00000236426', 'ENSG00000082175', 'ENSG00000124743', 'ENSG00000162614', 'ENSG00000271848', 'ENSG00000237452', 'ENSG00000164530', 'ENSG00000126218', 'ENSG00000144712', 'ENSG00000272167', 'ENSG00000126882', 'ENSG00000236056', 'ENSG00000133477', 'ENSG00000251777']

#LV:

loci_gene_set = genes['gene_id'].unique() #1977 genes
loci_gene_set_as_set = set(loci_gene_set)
inter = loci_gene_set_as_set.intersection(genes_above_Q95_plus_inf) #93
inter_as_list = list(inter) 
lv_top_expr_genes = pd.DataFrame()
lv_top_expr_genes['top_expr_genes'] = inter_as_list
# ['ENSG00000134769', 'ENSG00000183888', 'ENSG00000140403', 'ENSG00000158286', 'ENSG00000237742', 'ENSG00000104881', 'ENSG00000173641', 'ENSG00000224163', 'ENSG00000259910', 'ENSG00000124772', 'ENSG00000185739', 'ENSG00000136383', 'ENSG00000184908', 'ENSG00000200378', 'ENSG00000228672', 'ENSG00000140986', 'ENSG00000170417', 'ENSG00000065054', 'ENSG00000172572', 'ENSG00000143499', 'ENSG00000267257', 'ENSG00000175946', 'ENSG00000136378', 'ENSG00000137251', 'ENSG00000023330', 'ENSG00000198626', 'ENSG00000271579', 'ENSG00000147684', 'ENSG00000235050', 'ENSG00000079215', 'ENSG00000250479', 'ENSG00000214814', 'ENSG00000122971', 'ENSG00000156222', 'ENSG00000146147', 'ENSG00000198523', 'ENSG00000255394', 'ENSG00000137198', 'ENSG00000198796', 'ENSG00000249816', 'ENSG00000136574', 'ENSG00000217330', 'ENSG00000104879', 'ENSG00000163492', 'ENSG00000088726', 'ENSG00000255146', 'ENSG00000144596', 'ENSG00000155657', 'ENSG00000267784', 'ENSG00000263096', 'ENSG00000162009', 'ENSG00000128591', 'ENSG00000185519', 'ENSG00000260403', 'ENSG00000114854', 'ENSG00000205221', 'ENSG00000215533', 'ENSG00000150401', 'ENSG00000252396', 'ENSG00000256849', 'ENSG00000104936', 'ENSG00000213560', 'ENSG00000244020', 'ENSG00000166317', 'ENSG00000153531', 'ENSG00000254294', 'ENSG00000184608', 'ENSG00000255303', 'ENSG00000197321', 'ENSG00000077522', 'ENSG00000198729', 'ENSG00000140990', 'ENSG00000160539', 'ENSG00000242902', 'ENSG00000242849', 'ENSG00000134775', 'ENSG00000251576', 'ENSG00000261713', 'ENSG00000117707', 'ENSG00000226005', 'ENSG00000137441', 'ENSG00000124743', 'ENSG00000162614', 'ENSG00000271848', 'ENSG00000144712', 'ENSG00000233954', 'ENSG00000272167', 'ENSG00000135547', 'ENSG00000223387', 'ENSG00000126882', 'ENSG00000236056', 'ENSG00000133477', 'ENSG00000251777']

#############################################################################################################################################################################################################

##SORTING AND RANKING ACCORDING TO EXPRESSION VALUES

#AA:

genes.drop(['gene_symbol'], axis=1, inplace=True)
combined_aa = genes.merge(long_initial_merge_aa, on='gene_id')

combined_aa.sort_values(['evidence_gene', 'tissue', 'median_expr'], ascending=[True, True, False], inplace=True)
combined_aa['expr_rank'] = combined_aa.groupby(['evidence_gene', 'tissue'])['median_expr'].rank(ascending=False)
new_df = pd.DataFrame(combined_aa.groupby(['evidence_gene', 'tissue'])['median_expr'].count()).reset_index()[['evidence_gene','median_expr']]
new_df.rename({'median_expr':'expr_rank_denom'}, axis=1, inplace=True)
combined_aa = combined_aa.merge(new_df, how='left', on='evidence_gene')
combined_aa['expr_rank'] = combined_aa['expr_rank'].astype(int)
combined_aa['expr_rank'] = combined_aa['expr_rank'].astype(str) + '/' + combined_aa['expr_rank_denom'].astype(str)
combined_aa.drop(['expr_rank_denom'], axis=1, inplace=True)

#LV:

#genes.drop(['gene_symbol'], axis=1, inplace=True)
combined_lv = genes.merge(long_initial_merge_lv, on='gene_id')

#sorting and ranking 
combined_lv.sort_values(['evidence_gene', 'tissue', 'median_expr'], ascending=[True, True, False], inplace=True)
combined_lv['expr_rank'] = combined_lv.groupby(['evidence_gene', 'tissue'])['median_expr'].rank(ascending=False)
new_df = pd.DataFrame(combined_lv.groupby(['evidence_gene', 'tissue'])['median_expr'].count()).reset_index()[['evidence_gene','median_expr']]
new_df.rename({'median_expr':'expr_rank_denom'}, axis=1, inplace=True)
combined_lv = combined_lv.merge(new_df, how='left', on='evidence_gene')
combined_lv['expr_rank'] = combined_lv['expr_rank'].astype(int)
combined_lv['expr_rank'] = combined_lv['expr_rank'].astype(str) + '/' + combined_lv['expr_rank_denom'].astype(str)
combined_lv.drop(['expr_rank_denom'], axis=1, inplace=True)

#############################################################################################################################################################################################################

##GIVING A PROXIMITY OF GENE TO EVIDENCE GENE RANKING

#AA:

fin = combined_aa

#calculate absolute distance from beginning of gene to variant position as well as end of gene to variant position. Because if the gene precedes the variant then the end distance will be the closest distance to the variant, whereas if it procedes the variant then the start distance will be the closest. So get the smallest distance between the start and end and then compare this absolute value to all the rest of the genes to rank closest to furthest from variant 
fin['gene_start'] = fin['gene_start'].astype(float)
fin['gene_end'] = fin['gene_end'].astype(float)
fin['position'] = fin['position'].astype(float)
fin['start_dist'] = abs(fin['gene_start'] - fin['position'])
fin['end_dist'] = abs(fin['gene_end'] - fin['position'])
fin['gene_to_evgene_dist'] = fin[['start_dist','end_dist']].min(axis=1)
fin.drop(['start_dist', 'end_dist'], axis=1, inplace=True)
fin['proximity_rank'] = fin.groupby(['evidence_gene', 'tissue'])['gene_to_evgene_dist'].rank(ascending=True)

final_merge_aa = fin.merge(wide_all_tissues_aa[['fold_change_aa', 'median_of_medians', 'gene_id']], on='gene_id')

#LV:

fin = combined_lv

#calculate absolute distance from beginning of gene to variant position as well as end of gene to variant position. Because if the gene precedes the variant then the end distance will be the closest distance to the variant, whereas if it procedes the variant then the start distance will be the closest. So get the smallest distance between the start and end and then compare this absolute value to all the rest of the genes to rank closest to furthest from variant 
fin['gene_start'] = fin['gene_start'].astype(float)
fin['gene_end'] = fin['gene_end'].astype(float)
fin['position'] = fin['position'].astype(float)
fin['start_dist'] = abs(fin['gene_start'] - fin['position'])
fin['end_dist'] = abs(fin['gene_end'] - fin['position'])
fin['gene_to_evgene_dist'] = fin[['start_dist','end_dist']].min(axis=1)
fin.drop(['start_dist', 'end_dist'], axis=1, inplace=True)
fin['proximity_rank'] = fin.groupby(['evidence_gene', 'tissue'])['gene_to_evgene_dist'].rank(ascending=True)
new_df = pd.DataFrame(fin.groupby(['evidence_gene', 'tissue'])['gene_to_evgene_dist'].count()).reset_index()[['evidence_gene','gene_to_evgene_dist']]
new_df.rename({'gene_to_evgene_dist':'proximity_rank_denom'}, axis=1, inplace=True)
fin = fin.merge(new_df, how='left', on='evidence_gene')
fin['proximity_rank'] = fin['proximity_rank'].astype(int)
fin['proximity_rank'] = fin['proximity_rank'].astype(str) + '/' + fin['proximity_rank_denom'].astype(str)
fin.drop(['proximity_rank_denom'], axis=1, inplace=True)

final_merge_lv = fin.merge(wide_all_tissues_lv[['fold_change_lv', 'median_of_medians', 'gene_id']], on='gene_id')

##############################################################################################################################################################################################################

##ENSURING THAT ONLY THE TOP EXPRESSED GENES (IDENTIFIED PREVIOUSLY IN generate_fold_changes_using_median.py) ARE INCLUDED IN THIS FINAL OUTPUT

#AA:

aa_top_expr_gene_ids = aa_top_expr_genes['top_expr_genes'].unique()
final_merge_aa_incl = final_merge_aa[final_merge_aa['gene_id'].isin(aa_top_expr_gene_ids)].sort_values(by=['evidence_gene', 'fold_change_aa'], ascending=[True, False])

#LV:

##now get only top expr genes in this df:

lv_top_expr_gene_ids = lv_top_expr_genes['top_expr_genes'].unique()
final_merge_lv_incl = final_merge_lv[final_merge_lv['gene_id'].isin(lv_top_expr_gene_ids)].sort_values(by=['evidence_gene', 'fold_change_lv'], ascending=[True, False])

#############################################################################################################################################################################################################

##GOING FORWARDS WITH ONLY LV NOT AA

##BRINGING ALL RELEVANT INFO TOGETHER TO HAVE DESIRED OUTPUT DF REPRESENTING THE TOP EXPR GENES ONLY

#we want evidence_gene | gene | expr_in_heart | median_of_medians | nr_tissues_in_calc | median_of_medians range | fold_change
#Getting nr other tissues involved in the calculation for each gene
other_tisues = ['Adipose - Subcutaneous', 'Adipose - Visceral (Omentum)', 'Adrenal Gland', 'Artery - Aorta', 'Artery - Coronary', 'Artery - Tibial', 'Bladder', 'Brain - Amygdala', 'Brain - Anterior cingulate cortex (BA24)', 'Brain - Caudate (basal ganglia)', 'Brain - Cerebellar Hemisphere', 'Brain - Cerebellum', 'Brain - Cortex', 'Brain - Frontal Cortex (BA9)', 'Brain - Hippocampus', 'Brain - Hypothalamus', 'Brain - Nucleus accumbens (basal ganglia)', 'Brain - Putamen (basal ganglia)', 'Brain - Spinal cord (cervical c-1)', 'Brain - Substantia nigra', 'Breast - Mammary Tissue', 'Cells - EBV-transformed lymphocytes', 'Cells - Transformed fibroblasts', 'Cervix - Ectocervix', 'Cervix - Endocervix', 'Colon - Sigmoid', 'Colon - Transverse', 'Esophagus - Gastroesophageal Junction', 'Esophagus - Mucosa', 'Esophagus - Muscularis', 'Fallopian Tube', 'Kidney - Cortex', 'Liver', 'Lung', 'Minor Salivary Gland', 'Muscle - Skeletal', 'Nerve - Tibial', 'Ovary', 'Pancreas', 'Pituitary', 'Prostate', 'Skin - Not Sun Exposed (Suprapubic)', 'Skin - Sun Exposed (Lower leg)', 'Small Intestine - Terminal Ileum', 'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina', 'Whole Blood']
wide_all_tissues_lv[other_tisues].count(axis=1) #all rows are 51, so no null values in all the tissue entries 
#Getting range of median_of_medians, as well as mean
#wide_all_tissues_lv['range_median_of_medians'] = "Min: " + str(wide_all_tissues_lv[other_tisues].min(axis = 1)) + ", Max: " + str(wide_all_tissues_lv[other_tisues].max(axis = 1)) + ", Mean: " +  str(wide_all_tissues_lv[other_tisues].mean(axis = 1))
wide_all_tissues_lv['min_other_tissue_medians'] = wide_all_tissues_lv[other_tisues].min(axis = 1)
wide_all_tissues_lv['max_other_tissue_medians'] = wide_all_tissues_lv[other_tisues].max(axis = 1)
wide_all_tissues_lv['mean_other_tissue_medians'] = wide_all_tissues_lv[other_tisues].mean(axis = 1)
#Getting expr in heart and median_of_medians in one df:
filt = wide_all_tissues_lv[['gene_id', 'gene_symbol', 'Heart_LV_median_TPM', 'median_of_medians', 'fold_change_lv', 'min_other_tissue_medians', 'max_other_tissue_medians', 'mean_other_tissue_medians']]
filt.rename(columns={"Heart_LV_median_TPM": "median_expr_lv", "Heart_AA_median_TPM":"Heart - Atrial Appendage"}, inplace=True)
filt.rename(columns={"median_of_medians": "median_other_tissue_medians"}, inplace=True)
#merge on the evidence_gene info
final = final_merge_lv_incl[['evidence_gene', 'gene_id', 'expr_rank', 'gene_to_evgene_dist', 'proximity_rank']].merge(filt, on=['gene_id'], how='left')
final = final[['evidence_gene', 'gene_id', 'gene_symbol', 'median_expr_lv', 'median_other_tissue_medians', 'min_other_tissue_medians', 'mean_other_tissue_medians', 'max_other_tissue_medians', 'fold_change_lv', 'expr_rank', 'gene_to_evgene_dist', 'proximity_rank']]
final.to_csv('top_expr_lv_genes_with_expr_fold_change_and_proximity.csv', index=False)

#############################################################################################################################################################################################################

##BRINGING ALL THE GENES BACK INTO THE DF NOT ONLY THE TOP EXPR GENES, BUT MAKING SURE TO MARK THESE GENES AS TOP EXPR 

##Creating a different output df where all genes within the loci are depicted regardless of whether or not they fall within the top 49 lv genes. One table showing them sorted by median_expr (in heart), another table showing them sorted by fold_change (comparison to other tissues), and possibly another table showing them sorted by a combo of both integrated into a score of some sort...
#First run generate_fold_changes_using_medians.py
#then run from fold_changes_to_top_expr_genes.py up until (and including) final_merge_lv[['gene_openTarget', 'gene_id', 'gene_symbol', 'median_expr']]
#now wrangle this df to get the info we want:
#first we need to implement a binary col which tells you whether or not the genes is included in the top 49 expr genes determined via top fold change values:

final_merge_lv['top_expr_gene_id'] = ''
final_merge_lv['top_expr_gene_id'][final_merge_lv['gene_id'].isin(lv_top_expr_gene_ids)] = 'Y'

#############################################################################################################################################################################################################

##SORTING DF BY 1). CUSTOM QUADRANTS IE 2,1,4,3 AND 2). DIAGONAL DISTANCE WITHIN THESE QUADRANTS

##New plan in light of inspecting the above images etc
#Basically, the whole point of doing those statistics like adj r, slope, std error, and p val was simply just for the sake of getting the slope -- the other values do not tell us anything regarding statistical significance because we already know that fold change and median expr are not independent variables (ie expression is included in the formula to reach fold change).. so p value in this context does not give significance. The only value of interest here is the slope since this tells us the higher the slope, the higher the fold change. 
#The way slope is calculated in linear regression --> b (slope) = r.Sy/Sx (where r is pearsons correlation coefficient and Sx and Sy are standard dev of x and y respectively). So this is a very different way to calcualting the slope ie using the Ordinary Least Square method versus the standard way in y = mx + c formula where it is change in y/change in x. 
#

loci_expr_and_foldchange = final_merge_lv[['evidence_gene', 'gene_id', 'gene_symbol', 'median_expr', 'fold_change_lv' , 'top_expr_gene_id']]
loci_expr_and_foldchange = loci_expr_and_foldchange[~(loci_expr_and_foldchange['median_expr']==0)]

#Quadrant for 3 loci:
#first determine the median of median_expr in the heart so that we know what separates the different quadrants eg we know that fold change >= 1.3 is considered high and that below would be considered low
#getting the df one with infinity fold change values included:
#we want to give the inf values a more accurate comparison and so wherever there is an inf value in the fold change column, take the highest fold change value in the df (besides inf) and give it that value instead
#loci_expr_and_foldchange['fold_change_lv'][~np.isinf(loci_expr_and_foldchange['fold_change_lv'])].max() #this gives 5873.9785
#loci_expr_and_foldchange['fold_change_lv'] = loci_expr_and_foldchange['fold_change_lv'].replace([np.inf], 5873.9785)
#trying a different approach here where we take the maximum fold change value within that particular locus rather than across all the loci:
testing = []
for un in loci_expr_and_foldchange['evidence_gene'].unique():
	temp = loci_expr_and_foldchange[loci_expr_and_foldchange['evidence_gene']==un]
	temp2 = temp.loc[temp['fold_change_lv'] != np.inf, 'fold_change_lv'].max()
	temp['fold_change_lv'].replace(np.inf,temp2,inplace=True)
	testing.append(temp)
final_testing = pd.concat(testing)
loci_expr_and_foldchange = final_testing

#should maybe consider doing this without including genes that had expression = 0 in the mix. Since we want to get an idea of the typical range of expression among those genes that are actually important in the heart in order to determine an expression cutoff in the quadrants. We also want to include all genes in the heart, so must go back to older df long_initial_merge_lv
long_initial_merge_lv = long_initial_merge_lv[~(long_initial_merge_lv['median_expr']==0)]
median_of_median_expr = long_initial_merge_lv["median_expr"].median() #1.0395

loci_expr_and_foldchange = loci_expr_and_foldchange.sort_values(['fold_change_lv', 'median_expr'], ascending = [False, False])

#using top gene metric to prioritise loci
#Sort genes within loci by BOTH fold change as well as median_expr
#the way it works is to sort by fold change and then expression
loci_expr_and_foldchange['gradient_per_gene'] = loci_expr_and_foldchange['fold_change_lv']/loci_expr_and_foldchange['median_expr']

#now determine which quadrant it is in:
loci_expr_and_foldchange['quadrant'] = ''
loci_expr_and_foldchange['quadrant'][(loci_expr_and_foldchange['fold_change_lv']>=1.3) & (loci_expr_and_foldchange['median_expr']<1.0395)] = 1
loci_expr_and_foldchange['quadrant'][(loci_expr_and_foldchange['fold_change_lv']>=1.3) & (loci_expr_and_foldchange['median_expr']>=1.0395)] = 2
loci_expr_and_foldchange['quadrant'][(loci_expr_and_foldchange['fold_change_lv']<1.3) & (loci_expr_and_foldchange['median_expr']<1.0395)] = 3
loci_expr_and_foldchange['quadrant'][(loci_expr_and_foldchange['fold_change_lv']<1.3) & (loci_expr_and_foldchange['median_expr']>=1.0395)] = 4
#adding diagonal distance (ie hypotenuse to fold change and median expr) to help sort in addition (and after) sorting by custom quadrant
loci_expr_and_foldchange['diagonal_distance'] = np.sqrt((loci_expr_and_foldchange['fold_change_lv'] ** 2) + (loci_expr_and_foldchange['median_expr'] ** 2))
loci_expr_and_foldchange.sort_values(by=['evidence_gene', 'quadrant', 'diagonal_distance'], inplace=True, ascending=[True,True,False])

#for ranking loci:
loci_expr_and_foldchange['locus_prioritised'] = ''
locus_prioritisation = pd.DataFrame()
for targ in loci_expr_and_foldchange['evidence_gene'].unique():
	quad2 = loci_expr_and_foldchange[(loci_expr_and_foldchange['evidence_gene']==targ) & (loci_expr_and_foldchange['quadrant']==2)]
	quad1 = loci_expr_and_foldchange[(loci_expr_and_foldchange['evidence_gene']==targ) & (loci_expr_and_foldchange['quadrant']==1)]
	quad3 = loci_expr_and_foldchange[(loci_expr_and_foldchange['evidence_gene']==targ) & (loci_expr_and_foldchange['quadrant']==3)]
	quad4 = loci_expr_and_foldchange[(loci_expr_and_foldchange['evidence_gene']==targ) & (loci_expr_and_foldchange['quadrant']==4)]
	if quad2.shape[0] > 0:
		#get top one ie already sorted by diagonal_distance, so this gets the longest distance
		
		loci_expr_and_foldchange['locus_prioritised'] = np.where( ( (loci_expr_and_foldchange['evidence_gene']==targ) & (loci_expr_and_foldchange['quadrant']==2) ), 'Y', 'N')
		prior_loci_2 = loci_expr_and_foldchange[loci_expr_and_foldchange['locus_prioritised']=='Y'].head(1)
		print(prior_loci_2)
		locus_prioritisation = pd.concat([locus_prioritisation, prior_loci_2])
		print(locus_prioritisation)
		#loci_expr_and_foldchange['locus_prioritised'][(loci_expr_and_foldchange['gene_openTarget']==targ) & (loci_expr_and_foldchange['quadrant']==2)].head(1) = 'Y'
		#loci_expr_and_foldchange['locus_prioritised'][~(loci_expr_and_foldchange['gene_openTarget']==targ) & (loci_expr_and_foldchange['quadrant']==2)] = 'N'
	elif quad1.shape[0] > 0:
		loci_expr_and_foldchange['locus_prioritised'] = np.where( ( (loci_expr_and_foldchange['evidence_gene']==targ) & (loci_expr_and_foldchange['quadrant']==1) ), 'Y', 'N')
		prior_loci_1 = loci_expr_and_foldchange[loci_expr_and_foldchange['locus_prioritised']=='Y'].head(1)	
		print(prior_loci_1)
		locus_prioritisation = pd.concat([locus_prioritisation, prior_loci_1])
		print(locus_prioritisation)
	elif quad4.shape[0] > 0:
		loci_expr_and_foldchange['locus_prioritised'] = np.where( ( (loci_expr_and_foldchange['evidence_gene']==targ) & (loci_expr_and_foldchange['quadrant']==4) ), 'Y', 'N')
		prior_loci_4 = loci_expr_and_foldchange[loci_expr_and_foldchange['locus_prioritised']=='Y'].head(1)
		print(prior_loci_4)
		locus_prioritisation = pd.concat([locus_prioritisation, prior_loci_4])
		print(locus_prioritisation)
	else:
		loci_expr_and_foldchange['locus_prioritised'][loci_expr_and_foldchange['evidence_gene']==targ] = 'Cannot prioritise locus due no gradients in quadrant 2, 1, or 4'
		cannot_prior = loci_expr_and_foldchange[loci_expr_and_foldchange['locus_prioritised']=='Cannot prioritise locus due no gradients in quadrant 2, 1, or 4']
		print(cannot_prior)
		locus_prioritisation = pd.concat([locus_prioritisation, cannot_prior])
		print(locus_prioritisation)
locus_prioritisation.drop(['locus_prioritised'], axis=1, inplace=True)

locus_prioritisation['quadrant'] = pd.Categorical(locus_prioritisation['quadrant'], [2, 1, 4, 3])
locus_prioritisation.sort_values("quadrant", inplace=True)

locus_prioritisation['median_expr'] = locus_prioritisation['median_expr'].round(decimals = 3)
locus_prioritisation['fold_change_lv'] = locus_prioritisation['fold_change_lv'].round(decimals = 3)
locus_prioritisation['gradient_per_gene'] = locus_prioritisation['gradient_per_gene'].round(decimals = 3)
locus_prioritisation['diagonal_distance'] = locus_prioritisation['diagonal_distance'].round(decimals = 3)

sorter = [2, 1, 4, 3]
x = pd.DataFrame({'quadrant': sorter})
x.index = x.index.set_names('number')
x = x.reset_index()
locus_prioritisation = pd.merge(locus_prioritisation, x, how='left', on='quadrant')
locus_prioritisation.sort_values(['number', 'diagonal_distance'], ascending = [True, False], inplace=True)
locus_prioritisation.drop(['number'], axis=1, inplace=True)

locus_prioritisation.to_csv('loci_prioritisation_sorted_by_quadrant.csv', index=False)

#perform lines 323 to 334 before running next command:
gene_prioritisation = loci_expr_and_foldchange
gene_prioritisation.sort_values(["evidence_gene", "quadrant", "diagonal_distance"], inplace=True, ascending=[True, True, False])
gene_prioritisation.drop(['top_expr_gene_id'], axis=1, inplace=True)
gene_prioritisation['median_expr'] = gene_prioritisation['median_expr'].round(decimals = 3)
gene_prioritisation['fold_change_lv'] = gene_prioritisation['fold_change_lv'].round(decimals = 3)
gene_prioritisation['gradient_per_gene'] = gene_prioritisation['gradient_per_gene'].round(decimals = 3)
gene_prioritisation['diagonal_distance'] = gene_prioritisation['diagonal_distance'].round(decimals = 3)
gene_prioritisation.to_csv("gene_prioritisation.csv", index=False)


