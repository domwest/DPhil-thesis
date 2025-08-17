##IN TERMINAL!!!

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

gpu_device = 'GPU:0'

with tf.device(gpu_device):
    #
    labels = h5py.File('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/ML_training_sets/labels.h5', 'r')
    labels = labels['label'][:]
    labels.shape

with tf.device(gpu_device):
    #
    data = h5py.File('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/ML_training_sets/sequences.h5', 'r')
    data = data['data'][:]
    data.shape

with tf.device(gpu_device):
    #
    X_train, X_test, y_train, y_test = train_test_split(data,   labels, test_size=0.4, random_state=42, shuffle=True)
    X_val,   X_test, y_val,   y_test = train_test_split(X_test, y_test, test_size=0.5, random_state=42, shuffle=True)
    #
    print("Train:     ", X_train.shape, y_train.shape)
    print("Validation:", X_val.shape,   y_val.shape)
    print("Test:      ", X_test.shape,  y_test.shape)

with tf.device(gpu_device):
    #
    obj = sc.read_h5ad("/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/cm.h5ad")
    df_train = pd.DataFrame([np.sum(y_train, axis=0)], index=['count'], columns=obj.var.index.tolist()).T
    df_train['set'] = 'train_label'
    df_val   = pd.DataFrame([np.sum(y_val, axis=0)],   index=['count'], columns=obj.var.index.tolist()).T
    df_val['set'] = 'val_label'
    df_test  = pd.DataFrame([np.sum(y_test, axis=0)],  index=['count'], columns=obj.var.index.tolist()).T
    df_test['set'] = 'test_label'
    df = pd.concat([df_train, df_val, df_test])
    df['cell_type'] = df.index.tolist()

with tf.device(gpu_device):
    batch_size = 100
    epochs     = 100

with tf.device(gpu_device):
    input_shape  = (X_train.shape[-2], X_train.shape[-1])
    output_shape = labels.shape[-1]

del X_train, y_train, X_val, y_val

model = Sequential()
model.add(Input(shape=input_shape))
model.add(Conv1D(filters=4, kernel_size=(1), strides=(1), activation='relu', kernel_initializer=GlorotUniform(seed=42), padding="same"))
model.add(Dropout(0.2))
model.add(Conv1D(filters=300, kernel_size=(20), strides=(1), activation='relu', kernel_initializer=GlorotUniform(seed=42), padding="same"))
model.add(MaxPooling1D(pool_size=(3), strides=(3), padding="same"))
model.add(Dropout(0.2))
model.add(Conv1D(filters=600, kernel_size=(10), strides=1, activation='relu', kernel_initializer=GlorotUniform(seed=42), padding="same"))
model.add(MaxPooling1D(pool_size=(4), strides=(4), padding="same"))
model.add(Dropout(0.2))
model.add(Conv1D(filters=900, kernel_size=(8),  strides=1, activation='relu', kernel_initializer=GlorotUniform(seed=42), padding="same"))
model.add(MaxPooling1D(pool_size=(5), strides=(5), padding="same"))
model.add(Dropout(0.2))
model.add(Conv1D(filters=900, kernel_size=(4),  strides=1, activation='relu', kernel_initializer=GlorotUniform(seed=42), padding="same"))
model.add(MaxPooling1D(pool_size=(10), strides=(10), padding="same"))
model.add(Dropout(0.2))
model.add(Conv1D(filters=900, kernel_size=(8),  strides=1, activation='relu', kernel_initializer=GlorotUniform(seed=42), padding="same"))
model.add(MaxPooling1D(pool_size=(4), strides=(4), padding="same"))
model.add(Dropout(0.2))
model.add(Flatten())
model.add(Dense(output_shape, activation='sigmoid', kernel_regularizer=l2(0.001)))
model.summary()

with tf.device(gpu_device):
    #
    early_stop = EarlyStopping(monitor='val_loss', patience=10, verbose=1)
    #
    wBestModel = '/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/weights_new/model.h5'
    model.load_weights(wBestModel)
    #
    best_model = ModelCheckpoint(wBestModel, verbose=1, save_best_only=True)
    #
    hBestModel     = '/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/histories_new/model.csv'
    history_logger = CSVLogger(hBestModel, separator=",", append=True)
    #
    model.compile(Adam(learning_rate=0.0001),
                      loss    = "binary_crossentropy", 
                      metrics = [tf.keras.metrics.AUC()])

import os, pybedtools, shutil

if not os.path.isdir("/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/tmp"):
    os.mkdir("/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/tmp")
pybedtools.set_tempdir('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/tmp')

genome_fasta = "/well/PROCARDIS/domwest/upstr_processing/UpStreamPipeline_old/reference_genomes/upstream_pipelines_output/Download_Index_Genomes/hg38_bowtie2/ref_genome_analysis/results/hg38.fa"

##Center - Hidden layers

#Choosing regions in which to predict seqs for the classifiers...

#Region 1:
#Went back to 01_summary_dataset_creation notebook and looked at the df created (after adding pct col):
#Checked for a peak where there were quite a few classifiers in the name col ie pct was fairly high (in this case 0.79)
#So:
#df['name'][df['pct']==0.79][115]
#returned:
#'atrial_cm_ATAC,cardiac_muscle_cell_DNASE_paired,endothelial_ATAC,bulk_heart_LV_tissue_ATAC_paired,macrophage_ATAC,ventricular_cm_ATAC,cardiac_muscle_cell_DNASE_single,lymphocyte_ATAC,bulk_heart_LV_tissue_CHIP_HISTONE_H3K36me3_single,cardiac_muscle_cell_CHIP_HISTONE_H3K4me3_single,bulk_heart_LV_tissue_CHIP_TF_single,adipocyte_ATAC,bulk_heart_LV_tissue_CHIP_HISTONE_H3K4me3_single,bulk_heart_LV_tissue_CHIP_HISTONE_H3K4me1_single,bulk_heart_LV_tissue_CHIP_TF_paired,smooth_muscle_ATAC,bulk_heart_LV_tissue_DNASE_paired,bulk_heart_LV_tissue_CHIP_HISTONE_H3K27ac_single,fibroblast_ATAC'
#this is 19 out of the 24 classifers sharing this peak overlap...
#At this index (115) the peak is: chr1 	1396946 	1401081 #4135bp
#So let's extend this region but include that peak too:
#chr1 1379046 1418981 #39935bp
#Region 2:
#df['name'][df['pct']==0.12][227070]
#'cardiac_muscle_cell_DNASE_single,cardiac_muscle_cell_CHIP_TF_single,bulk_heart_LV_tissue_CHIP_TF_single'
#peak is: 227070 	chr9 	175817 	176695 #878
#chr9 163756 188756 #25000
#Region 3:
#df['name'][df['pct']==0.04][130605]
#'ventricular_cm_ATAC'
#peak is: 130605 	chr2 	242080657 	242081103 #446
#chr2 242067380 242094380 #27000

# # 01_ - Simone's region 1 was 40000bp
# regionOI = 'chr1:1379046-1418981' #39935bp

# 02_ - Simone's region 2 was 25000bp
#regionOI = 'chr9:163756-188756'  #25000bp

# 03_ - Simone's region 3 was 27000bp
regionOI = 'chr2:242067380-242094380' #27000bp

chrom = regionOI.split(":")[0]
start = regionOI.split(":")[-1].split("-")[0]
end   = regionOI.split(":")[-1].split("-")[-1]
print(chrom, start, end)

starts = np.arange(int(start), int(end)-1000)
ends   = np.arange(int(start)+1000, int(end))

regions = []
for s, e in zip(starts, ends):
    regions.append('%s:%s-%s'%(chrom, s, e))

for region in regions:
    print(region)
    d = {'chrom': [region.split(':')[0]], 'start': [region.split(':')[-1].split("-")[0]], 'end': [region.split(':')[-1].split("-")[-1]]}
    print(d)
    df = pd.DataFrame(data=d)
    print(df)
    break

#chr2:242067380-242094380
#first region: chr2:242067380-242068380 #1000
#last region: chr2:242093379-242094379 #1000

#GETTING PEAK FROM REGION 1 AND INSERTING THAT TO THE BEGINNING OF REGION3 FOR ALL CLASSIFIERS

###########################################################################################################################################################################

#all sequence

regionOI = 'chr2:242067380-242094380'

d = {'chrom': [regionOI.split(':')[0]], 'start': [regionOI.split(':')[-1].split("-")[0]], 'end': [regionOI.split(':')[-1].split("-")[-1]]}
df = pd.DataFrame(data=d)
df['start'] = df['start'].astype(int)
df['end'] = df['end'].astype(int)
df = df[['chrom', 'start', 'end']]

regions = []
dfT = df.T
for i in dfT:
    regions.append(dfT[i].tolist())

region_pybedtools = pybedtools.BedTool(regions)
region_pybedtools = region_pybedtools.sequence(fi=genome_fasta)
list_fasta = open(region_pybedtools.seqfn).read().split('\n')[:-1]
sequence = [list_fasta[i].upper() for i in range(1, df.shape[0]*2, 2)]

start_seq = 0
end_seq = int(regionOI.split(':')[-1].split('-')[-1])-int(regionOI.split(':')[-1].split('-')[0]) #39935

#region to edit ie peak

region_to_edit = 'chr2:242067680-242071815'

d = {'chrom': [region_to_edit.split(':')[0]], 'start': [region_to_edit.split(':')[-1].split("-")[0]], 'end': [region_to_edit.split(':')[-1].split("-")[-1]]}
df = pd.DataFrame(data=d)
df['start'] = df['start'].astype(int)
df['end'] = df['end'].astype(int)
df = df[['chrom', 'start', 'end']]

regions = []
dfT = df.T
for i in dfT:
    regions.append(dfT[i].tolist())

region_pybedtools = pybedtools.BedTool(regions)
region_pybedtools = region_pybedtools.sequence(fi=genome_fasta)
list_fasta = open(region_pybedtools.seqfn).read().split('\n')[:-1]
sequences_to_edit = [list_fasta[i].upper() for i in range(1, df.shape[0]*2, 2)]

start_seq_to_edit = int(region_to_edit.split(':')[-1].split('-')[0])-int(regionOI.split(':')[-1].split('-')[0]) #17900
end_seq_to_edit = start_seq_to_edit+(int(region_to_edit.split(':')[-1].split('-')[-1])-int(region_to_edit.split(':')[-1].split('-')[0])) #22035

#merge sequences

left_seq = sequence[0][:start_seq_to_edit]
len(left_seq)
#300

middle_seq = sequence[0][start_seq_to_edit:end_seq_to_edit]
len(middle_seq)
#4135

###########

# import random

# sequence_edited = []
# for i in range(len(middle_seq)):
#     sequence_edited.append(random.choice(list(set(['A','C','G','T']))))
# sequence_edited = ''.join(sequence_edited)
#4135

#GETTING PEAK FROM REGION 1 AND INSERTING THAT TO THE BEGINNING OF REGION3 FOR ALL CLASSIFIERS

d2 = {'chrom': 'chr1', 'start': 1396946, 'end': 1401081} #peak
df2 = pd.DataFrame(data=d, index=[0])

regions2 = [['chr1', 1396946, 1401081]] #peak
print(1401081 - 1396946)
lengthh = 1401081 - 1396946 #length of peak
region_pybedtools2 = pybedtools.BedTool(regions2)
region_pybedtools2 = region_pybedtools2.sequence(fi=genome_fasta)
list_fasta2 = open(region_pybedtools2.seqfn).read().split('\n')[:-1]
print(list_fasta2)
sequences2 = [list_fasta2[i] for i in range(1, df2.shape[0]*2, 2)] #sequences2 = [list_fasta2[i].upper() for i in range(1, df2.shape[0]*2, 2)]
print(sequences2)

print(len(sequences2))
for seq2 in sequences2:
   print(len(seq2)) ###should be 1?? also, has it assigned the variable properly?

sequence_edited = seq2 
#4135

###########

right_seq = sequence[0][end_seq_to_edit:]
len(right_seq)
#22565

final_sequence = left_seq+sequence_edited+right_seq
len(final_sequence)
#27000

#build data

chrom = regionOI.split(":")[0]
start = start_seq
end   = end_seq
print(chrom, start, end)
#chr1 0 39935

starts = np.arange(int(start), int(end)-1000)
ends   = np.arange(int(start)+1000, int(end))

regions_to_loop = []
for s, e in zip(starts, ends):
    regions_to_loop.append('%s:%s-%s'%(chrom, s, e))

from tqdm import tqdm

datas = []
for region in tqdm(regions_to_loop):
    d = {'chrom': [region.split(':')[0]], 'start': [region.split(':')[-1].split("-")[0]], 'end': [region.split(':')[-1].split("-")[-1]]}
    df = pd.DataFrame(data=d)
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    df['center'] = np.mean((df['start'], df['end']), axis=0)
    df['start'] = df['center']-500
    df['end'] = df['center']+500
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    df['diff'] = df['end']-df['start']
    df = df[['chrom', 'start', 'end']]
    #
    regions = []
    dfT = df.T
    for i in dfT:
        regions.append(dfT[i].tolist())
    #region_pybedtools = pybedtools.BedTool(regions)
    #region_pybedtools = region_pybedtools.sequence(fi=genome_fasta)
    #list_fasta = open(region_pybedtools.seqfn).read().split('\n')[:-1]
    #sequences = [list_fasta[i].upper() for i in range(1, df.shape[0]*2, 2)]
    sequences = [final_sequence[regions[0][1]:regions[0][2]].upper()]
    #
    #sequences = reverse_compl(sequences)
    #
    label_encoder = preprocessing.LabelEncoder()
    integer_encoded = label_encoder.fit(['A', 'C', 'G', 'T', 'N'])
    integer_encoded = label_encoder.transform(['A', 'C', 'G', 'T', 'N'])
    #
    onehot_encoder  = preprocessing.OneHotEncoder(sparse=False)
    integer_encoded = integer_encoded.reshape(len(integer_encoded), 1)
    onehot_encoded = onehot_encoder.fit(integer_encoded)
    onehot_encoded = onehot_encoder.transform(integer_encoded)
    #
    data = []
    for x in sequences:
        integer_encoded = label_encoder.transform([*x])
        integer_encoded = integer_encoded.reshape(len(integer_encoded), 1)
        onehot_encoded = onehot_encoder.transform(integer_encoded)
        data.append(np.delete(onehot_encoded, 3, 1))
    data = np.array(data)
    #
    datas.append(data)

datas_f = np.flip(datas)

predictions = []    
for data in tqdm(datas):
    prediction = model.predict(data, verbose=0)
    predictions.append(prediction)

predictions_f = []    
for data in tqdm(datas_f):
    prediction = model.predict(data, verbose=0)
    predictions_f.append(prediction)

target_names = obj.var.index.tolist()
n_classes = len(obj.var.index.tolist())

pred_df = []
for pred in predictions:
    p = pred[0]
    pred_df.append(p.tolist())

pred_df = pd.DataFrame(pred_df, columns=target_names)
pred_df.head()

pred_f_df = []
for pred in predictions_f:
    p = pred[0]
    pred_f_df.append(p.tolist())

pred_f_df = pd.DataFrame(pred_f_df, columns=target_names)[::-1]
pred_f_df.reset_index(inplace=True, drop=True)
pred_f_df.head()

from scipy.ndimage import gaussian_filter1d

pred_df_min = pd.DataFrame()
for col in pred_df.columns:
    tmp_value = np.min([pred_df[col], pred_f_df[col]], axis=0)
    value = gaussian_filter1d(tmp_value, sigma=15)
    pred_df_min[col] = value

from IPython.display import Image
from IPython.core.display import HTML 

region03 = 'chr2:242067380-242094380'

for col in pred_df.columns:
    fig, ax = plt.subplots(3,1,figsize=(150,30))
    ax[0] = pred_df[col].plot.line(title=col, fontsize=50, ax=ax[0])
    ax[0].get_xaxis().set_visible(False)
    ax[1] = pred_f_df[col].plot.line(title=col, fontsize=50, ax=ax[1])
    ax[0].set_title(col, fontsize = 60)
    ax[1].set_title("")
    ax[2] = pred_df_min[col].plot.line(title=col, fontsize=50, ax=ax[2])
    #ax[2].axhline(0.65, c='red', linewidth=0.3)
    ax[2].set_title("")
    plt.show()
    fig.savefig('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_modified_peak_rolling_window/validation_rolling_window_peak_creation_region3_prediction_' + col + '.png', bbox_inches='tight')
    print('@'*100)
    print('@'*100)

for col in pred_df.columns:
    fig, ax = plt.subplots(1,1,figsize=(150,30))
    ax = pred_df_min[col].plot.line(title=col, fontsize=50, ax=ax, linewidth=2, color='blue')
    ax.get_xaxis().set_visible(False)
    ax.set_title(col, fontsize = 60)
    plt.show()
    fig.savefig('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_modified_peak_rolling_window/validation_rolling_window_peak_creation_only_overlay_region3_prediction_' + col + '.png', bbox_inches='tight')
    print('@'*100)

#kill

# shutil.rmtree('tmp')