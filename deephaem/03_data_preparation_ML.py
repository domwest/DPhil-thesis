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
import numpy as np
import anndata
import scanpy as sc
import h5py, os, pybedtools, shutil
from sklearn import preprocessing
from tqdm import tqdm

if not os.path.exists('tmp'):
    os.makedirs('tmp')
pybedtools.set_tempdir('tmp')


# In[56]:


obj = sc.read_h5ad("/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/cm.h5ad")
#obj
#obj.var_names is cell_type
#obj.obs_names is chr:start-end
#obj.X is whether or not the specific peak is present in teh celltypes (binary so 1 is yes & 0 is no)

labels = pd.DataFrame(obj.X, columns=obj.var.index) 
labels.shape
# (321280, 65)

chrom = obj.obs.index.str.split("-").str[0].astype(str)
start = obj.obs.index.str.split("-").str[1].astype(float)
end   = obj.obs.index.str.split("-").str[2].astype(float)
#gives an array of all those different cols

df = pd.DataFrame([chrom, start, end], index=['chrom', 'start', 'end']).T
df['center'] = np.mean((df['start'], df['end']), axis=0)
df['start'] = df['center']-500
df['end'] = df['center']+500
df['diff'] = df['end']-df['start']
df
#321280

genome_fasta = "/well/PROCARDIS/domwest/upstr_processing/UpStreamPipeline_old/reference_genomes/upstream_pipelines_output/Download_Index_Genomes/hg38_bowtie2/ref_genome_analysis/results/hg38.fa"

regions = []
dfT = df[['chrom', 'start', 'end']]
dfT['start'] = dfT['start'].astype(int)
dfT['end'] = dfT['end'].astype(int)
dfT = dfT.T
for i in dfT:
    regions.append(dfT[i].tolist())
len(regions)
#321280

region_pybedtools = pybedtools.BedTool(regions)
region_pybedtools = region_pybedtools.sequence(fi=genome_fasta)
list_fasta = open(region_pybedtools.seqfn).read().split('\n')[:-1]
sequences = [list_fasta[i].upper() for i in range(1, df.shape[0]*2, 2)]
df.shape[0]*2
#642560
len(sequences)
#321280

##previous troubleshooting:
#len(list_fasta)
#488664
#list_fasta_df = pd.DataFrame()
#list_fasta_df[0] = list_fasta
#new_df = pd.DataFrame({'peak':list_fasta_df[0].iloc[::2].values, 'sequence':list_fasta_df[0].iloc[1::2].values})      
#new_df['chrom'] = new_df['peak'].str.split(':', expand=True)[0]
#new_df['chrom'] = new_df['chrom'].str.replace('>','')
#new_df['just_peak'] = new_df['peak'].str.split(':', expand=True)[1]
#new_df['start'] = new_df['just_peak'].str.split('-', expand=True)[0]
#new_df['end'] = new_df['just_peak'].str.split('-', expand=True)[1]
#new_df['start'] = new_df['start'].astype(float)
#new_df['end'] = new_df['end'].astype(float)
#new_df['seq_present'] = 'Y'
#244332 rows
#new_df2 = df.merge(new_df, how='left', on=['chrom', 'start', 'end'])
#new_df2[new_df2['seq_present']!='Y']
#new_df2[new_df2['seq_present']=='Y'] #122231 #Need to keep this nr rather than df.shape[0]*2 because there are some regions for which fasta sequences couldn't be fetched, so can't use the length of them all as a reference for how many seqs to expect (it is actually less)
#new_df2_with_seqs = new_df2[new_df2['seq_present']=='Y']
#new_df2_without_seqs = new_df2[new_df2['seq_present']!='Y']
#sequences = [list_fasta[i].upper() for i in range(1, 122231*2, 2)]
#sequences

label_encoder = preprocessing.LabelEncoder()
integer_encoded = label_encoder.fit(['A', 'C', 'G', 'T', 'N'])
integer_encoded = label_encoder.transform(['A', 'C', 'G', 'T', 'N'])

onehot_encoder  = preprocessing.OneHotEncoder(sparse=False)
integer_encoded = integer_encoded.reshape(len(integer_encoded), 1)
onehot_encoded = onehot_encoder.fit(integer_encoded)
onehot_encoded = onehot_encoder.transform(integer_encoded)

data = []
for x in tqdm(sequences):
    integer_encoded = label_encoder.transform([*x])
    integer_encoded = integer_encoded.reshape(len(integer_encoded), 1)
    onehot_encoded = onehot_encoder.transform(integer_encoded)
    data.append(np.delete(onehot_encoded, 3, 1)) #array, object, axis - so: onehot_encoded (array), object (3), axis (1 ie columns)
#integer_encoded
#onehot_encoded
#onehot_encoded.type()
#len(onehot_encoded)
#onehot_encoded[3] #array([0., 0., 0., 0., 1.]) #so this is the N entries?

data = np.array(data)
data.shape
# (321280, 1000, 4)

shutil.rmtree('tmp')

#data.shape #(38882, 1000, 4)
#data.type() #'numpy.ndarray' object
#for a 3D numpy nd array:
#   |       -- axis-2 ->
#   |    |
#   |  axis-1 [0, 1]
#   |    |    [2, 3]
#   |    V
#axis-0
#   |      -- axis-2 ->
#   |    |
#   |  axis-1 [4, 5]
#   |    |    [6, 7]
#   V    V

#len(data) #38882
#len(data[0]) #1000
#len(data[0][0]) #4

file = h5py.File('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/ML_training_sets/sequences.h5','w')
file.create_dataset('data', data = data)
file.close()
#so this is the actual hot encode format of the sequences (1000bp seqs) per peak

file = h5py.File('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/ML_training_sets/labels.h5','w')
file.create_dataset('label', data = labels.values)
file.close()
#so this is whether or not a specific peak is present in the different cell_types

labels.values
data
