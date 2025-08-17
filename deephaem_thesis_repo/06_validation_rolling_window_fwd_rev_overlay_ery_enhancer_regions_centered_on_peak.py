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
from scipy.ndimage import gaussian_filter1d

from tensorflow.keras.layers import Dense, Flatten, Dropout, Input, Conv1D, MaxPooling1D, BatchNormalization
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint, CSVLogger
from tensorflow.keras.initializers import GlorotUniform
from tensorflow.keras.models import Sequential
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.regularizers import l2
from tensorflow.keras import backend

import tensorflow as tf
import keras

with tf.device("CPU"):
    labels = h5py.File('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/ML_training_sets/labels.h5', 'r')
    labels = labels['label'][:]
    labels.shape

with tf.device("CPU"):
    data = h5py.File('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/ML_training_sets/sequences.h5', 'r')
    data = data['data'][:]
    data.shape

with tf.device("CPU"):
    X_train, X_test, y_train, y_test = train_test_split(data,   labels, test_size=0.4, random_state=42, shuffle=True)
    X_val,   X_test, y_val,   y_test = train_test_split(X_test, y_test, test_size=0.5, random_state=42, shuffle=True)
    print("Train:     ", X_train.shape, y_train.shape)
    print("Validation:", X_val.shape,   y_val.shape)
    print("Test:      ", X_test.shape,  y_test.shape)

with tf.device("CPU"):
    obj = sc.read_h5ad("/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/cm.h5ad")
    df_train = pd.DataFrame([np.sum(y_train, axis=0)], index=['count'], columns=obj.var.index.tolist()).T
    df_train['set'] = 'train_label'
    df_val   = pd.DataFrame([np.sum(y_val, axis=0)],   index=['count'], columns=obj.var.index.tolist()).T
    df_val['set'] = 'val_label'
    df_test  = pd.DataFrame([np.sum(y_test, axis=0)],  index=['count'], columns=obj.var.index.tolist()).T
    df_test['set'] = 'test_label'
    df = pd.concat([df_train, df_val, df_test])
    df['cell_type'] = df.index.tolist()

with tf.device("CPU"):
    batch_size = 100
    epochs     = 100

with tf.device("CPU"):
    input_shape  = (X_train.shape[-2], X_train.shape[-1])
    output_shape = labels.shape[-1]

# del X_train, y_train, X_val, y_val

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

with tf.device("CPU"):
    early_stop = EarlyStopping(monitor='val_loss', patience=10, verbose=1)
    wBestModel = '/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/weights_new/model.h5'
    model.load_weights(wBestModel)
    best_model = ModelCheckpoint(wBestModel, verbose=1, save_best_only=True)
    hBestModel     = '/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/histories_new/model.csv'
    history_logger = CSVLogger(hBestModel, separator=",", append=True)
    model.compile(Adam(learning_rate=0.0001),
                      loss    = "binary_crossentropy", 
                      metrics = [tf.keras.metrics.AUC()])

import os, pybedtools, shutil

if not os.path.isdir("/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/tmp"):
    os.mkdir("/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/tmp")
pybedtools.set_tempdir('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/tmp')

genome_fasta = "/well/PROCARDIS/domwest/upstr_processing/UpStreamPipeline_old/reference_genomes/upstream_pipelines_output/Download_Index_Genomes/hg38_bowtie2/ref_genome_analysis/results/hg38.fa"

############################################################################################################################################################

#PART A).

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

#            rsid   chr  grch38_POS REF ALT  downstream   upstream                  regionOI
# 0      rs198415  chr1    11841535   T   C    11840535   11842535    chr1:11840535-11842535
# 400  rs12469283  chr2   209324187   C   T   209323187  209325187  chr2:209323187-209325187
# 600  rs73064343  chr5    31470831   C   A    31469831   31471831    chr5:31469831-31471831
# 700   rs7798211  chr7    25781981   T   A    25780981   25782981    chr7:25780981-25782981

# regionOI = 'chr1:11830921-11852113' #21192bp # actual peak: regionOI = 'chr1:11841171-11841863' #692

regionOI = 'chr12:69648877-69671712' #22835bp # actual peak: regionOI = 'chr12:69658877-69661712' #2835

snp_dataset = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/bi_allelic_Enhancer_erythroid_pos_control_snp_dataset.csv')
snp_dataset['chr'] = 'chr' + snp_dataset['chr'].astype(str)
snp_dataset['downstream'] = snp_dataset['grch38_POS'] - 1000
snp_dataset['upstream'] = snp_dataset['grch38_POS'] + 1000
snp_dataset['regionOI'] = regionOI
snp_dataset = snp_dataset.iloc[[176]] #rs11177764

#####

def multi_allelic_visualisations(rsid, regionOI, ref, alt):
    #
    chrom = regionOI.split(":")[0]
    start = regionOI.split(":")[-1].split("-")[0]
    end   = regionOI.split(":")[-1].split("-")[-1]
    print(chrom, start, end)
    #
    starts = np.arange(int(start), int(end)-1000)
    ends   = np.arange(int(start)+1000, int(end))
    #
    regions = []
    for s, e in zip(starts, ends):
        regions.append('%s:%s-%s'%(chrom, s, e))
    #
    for region in regions:
        print(region)
        d = {'chrom': [region.split(':')[0]], 'start': [region.split(':')[-1].split("-")[0]], 'end': [region.split(':')[-1].split("-")[-1]]}
        print(d)
        df = pd.DataFrame(data=d)
        print(df)
        break
    #
    from tqdm import tqdm
    #
    datas = []
    #
    counter = 0
    #
    for region in tqdm(regions):
        #
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
        #
        region_pybedtools = pybedtools.BedTool(regions)
        region_pybedtools = region_pybedtools.sequence(fi=genome_fasta)
        list_fasta = open(region_pybedtools.seqfn).read().split('\n')[:-1]
        sequences = [list_fasta[i].upper() for i in range(1, df.shape[0]*2, 2)]
        #
        # string = sequences[0]
        # index = 999 - counter
        # new_char = alt
        # string_list = list(string)
        # string_list[index] = new_char
        # new_string = "".join(string_list)
        # sequences = [new_string]
        # counter +=1 
        #
        label_encoder = preprocessing.LabelEncoder()
        integer_encoded = label_encoder.fit(['A', 'C', 'G', 'T', 'N'])
        integer_encoded = label_encoder.transform(['A', 'C', 'G', 'T', 'N'])
        #
        onehot_encoder  = preprocessing.OneHotEncoder(sparse_output=False)
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
        #
    datas_rc = np.flip(datas)
    #
    predictions = []
    for data in tqdm(datas):
        prediction = model.predict(data, verbose=0)
        predictions.append(prediction)
    #
    predictions_rc = []
    for data in tqdm(datas_rc):
        prediction = model.predict(data, verbose=0)
        predictions_rc.append(prediction)
    #
    ##my addition:
    #
    def flatten_list(nested_list):
        flat_list = []
        for sublist in nested_list:
            for item in sublist:
                flat_list.append(item)
        return flat_list
    #
    predictions_reformat = [l.tolist() for l in predictions]
    predictions_reformat = flatten_list(predictions_reformat)
    #
    predictions_rc_reformat = [l.tolist() for l in predictions_rc]
    predictions_rc_reformat = flatten_list(predictions_rc_reformat)
    #
    #
    target_names = obj.var.index.tolist()
    n_classes = len(obj.var.index.tolist())
    pred_df = pd.DataFrame(predictions_reformat, columns=target_names)
    pred_df.head()
    pred_df_rc = pd.DataFrame(predictions_rc_reformat, columns=target_names)[::-1]
    pred_df_rc.head()
    #pred_df.to_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/01_prediction40kb_forward.cvs', sep='\t')
    # pred_df.to_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/02_prediction25kb_forward.cvs', sep='\t')
    # pred_df.to_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/03_prediction27kb_forward.cvs', sep='\t')
    #pred_df_rc.to_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/01_prediction40kb_reverse.cvs', sep='\t')
    # pred_df_rc.to_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/02_prediction25kb_reverse.cvs', sep='\t')
    # pred_df_rc.to_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/03_prediction27kb_reverse.cvs', sep='\t')
    # pred_df = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/01_prediction40kb_forward.cvs', sep='\t', index_col=0)
    # pred_df = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/02_prediction25kb_forward.cvs', sep='\t', index_col=0)
    # pred_df = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/03_prediction27kb_forward.cvs', sep='\t', index_col=0)
    # pred_df_rc = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/01_prediction40kb_reverse.cvs', sep='\t', index_col=0)
    # pred_df_rc = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/02_prediction25kb_reverse.cvs', sep='\t', index_col=0)
    # pred_df_rc = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/03_prediction27kb_reverse.cvs', sep='\t', index_col=0)
    #pred_df = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/01_prediction40kb_forward.cvs', sep='\t')
    #pred_df.drop(['Unnamed: 0'], axis=1, inplace=True)
    #
    #uncomment for 1000bp:
    # pred_df = pred_df.iloc[399:600,:]
    #
    #pred_df_rc = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/01_prediction40kb_reverse.cvs', sep='\t')
    #pred_df_rc.drop(['Unnamed: 0'], axis=1, inplace=True)
    #
    #uncomment for 1000bp:
    # pred_df_rc = pred_df_rc.iloc[399:600,:]
    #
    #
    #smooth the signal:
    #
    #smooth the signal:
    #
    pred_df_min = pd.DataFrame()
    for col in pred_df.columns:
        tmp_value = np.min([pred_df[col], pred_df_rc[col]], axis=0)
        value = gaussian_filter1d(tmp_value, sigma=15)
        pred_df_min[col] = value
    #
    from IPython.display import Image
    from IPython.core.display import HTML 
    #
    # region01 = 'chr1:1379046-1418981'
    #
    for col in pred_df.columns:
        fig, ax = plt.subplots(3,1,figsize=(150,30)) #try making it longer, width is good
        ax[0] = pred_df[col].plot.line(title=col, fontsize=50, ax=ax[0])
        ax[0].get_xaxis().set_visible(False)
        ax[1] = pred_df_rc[col].plot.line(title=col, fontsize=50, ax=ax[1]) #try taking col out here
        ax[0].set_title(col, fontsize = 60)
        ax[1].set_title("")
        ax[2] = pred_df_min[col].plot.line(title=col, fontsize=50, ax=ax[2]) #try taking col out here
        #ax[2].axhline(0.65, c='red', linewidth=0.3)
        ax[2].set_title("")
        plt.show()
        fig.savefig('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/sigmoid_65_validation_rolling_window_fwd_rev_prediction_' + col + '_' + rsid + '_' + ref + '_' + alt + '_Enhancer_erythroid_pos_control_without_vars_peak_centered.png', bbox_inches='tight')
        print('@'*100)
    #
    for col in pred_df.columns:
        fig, ax = plt.subplots(1,1,figsize=(150,30))
        ax = pred_df_min[col].plot.line(title=col, fontsize=50, ax=ax, linewidth=2, color='blue')
        ax.get_xaxis().set_visible(False)
        ax.set_title(col, fontsize = 60)
        plt.show()
        fig.savefig('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/sigmoid_65_validation_rolling_window_only_overlay_prediction_' + col + '_' + rsid + '_' + ref + '_' + alt + '_Enhancer_erythroid_pos_control_without_vars_peak_centered.png', bbox_inches='tight')
        print('@'*100)


for index, row in snp_dataset.iterrows():
    #print(index)
    rsid = row['rsid']
    regionOI = row['regionOI']
    alt = row['ALT']
    ref = row['REF']
    multi_allelic_visualisations(rsid, regionOI, ref, alt)

######################################################################################################################################################

############################################################################################################################################################

#PART F).

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
#At this index (115) the peak is: chr1  1396946     1401081 #4135bp
#So let's extend this region but include that peak too:
#chr1 1379046 1418981 #39935bp
#Region 2:
#df['name'][df['pct']==0.12][227070]
#'cardiac_muscle_cell_DNASE_single,cardiac_muscle_cell_CHIP_TF_single,bulk_heart_LV_tissue_CHIP_TF_single'
#peak is: 227070    chr9    175817  176695 #878
#chr9 163756 188756 #25000
#Region 3:
#df['name'][df['pct']==0.04][130605]
#'ventricular_cm_ATAC'
#peak is: 130605    chr2    242080657   242081103 #446
#chr2 242067380 242094380 #27000

#            rsid   chr  grch38_POS REF ALT  downstream   upstream                  regionOI
# 0      rs198415  chr1    11841535   T   C    11840535   11842535    chr1:11840535-11842535
# 400  rs12469283  chr2   209324187   C   T   209323187  209325187  chr2:209323187-209325187
# 600  rs73064343  chr5    31470831   C   A    31469831   31471831    chr5:31469831-31471831
# 700   rs7798211  chr7    25781981   T   A    25780981   25782981    chr7:25780981-25782981

# regionOI = 'chr1:11830921-11852113' #21192bp # actual peak: regionOI = 'chr1:11841171-11841863' #692

regionOI = 'chr12:69648877-69671712' #22835bp # actual peak: regionOI = 'chr12:69658877-69661712' #2835

snp_dataset = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/bi_allelic_Enhancer_erythroid_pos_control_snp_dataset.csv')
snp_dataset['chr'] = 'chr' + snp_dataset['chr'].astype(str)
snp_dataset['downstream'] = snp_dataset['grch38_POS'] - 1000
snp_dataset['upstream'] = snp_dataset['grch38_POS'] + 1000
snp_dataset['regionOI'] = regionOI
snp_dataset = snp_dataset.iloc[[176]] #rs11177764

# genome_fasta_with_vars = "/well/PROCARDIS/domwest/upstr_processing/UpStreamPipeline_old/reference_genomes/upstream_pipelines_output/Download_Index_Genomes/hg38_bowtie2/ref_genome_analysis/results/hg38_edited_for_ery_enhancer_pos_control_rs198415.fa" #hg38_edited_for_ten_pos_control_vars.fa

genome_fasta_with_vars = "/well/PROCARDIS/domwest/upstr_processing/UpStreamPipeline_old/reference_genomes/upstream_pipelines_output/Download_Index_Genomes/hg38_bowtie2/ref_genome_analysis/results/hg38_edited_for_ery_enhancer_pos_control_rs11177764.fa" #hg38_edited_for_ten_pos_control_vars.fa

#####

def multi_allelic_visualisations(rsid, regionOI, ref, alt):
    #
    chrom = regionOI.split(":")[0]
    start = regionOI.split(":")[-1].split("-")[0]
    end   = regionOI.split(":")[-1].split("-")[-1]
    print(chrom, start, end)
    #
    starts = np.arange(int(start), int(end)-1000)
    ends   = np.arange(int(start)+1000, int(end))
    #
    regions = []
    for s, e in zip(starts, ends):
        regions.append('%s:%s-%s'%(chrom, s, e))
    #
    for region in regions:
        print(region)
        d = {'chrom': [region.split(':')[0]], 'start': [region.split(':')[-1].split("-")[0]], 'end': [region.split(':')[-1].split("-")[-1]]}
        print(d)
        df = pd.DataFrame(data=d)
        print(df)
        break
    #
    from tqdm import tqdm
    #
    datas = []
    #
    # counter = 0
    #
    for region in tqdm(regions):
        #
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
        #
        region_pybedtools = pybedtools.BedTool(regions)
        region_pybedtools = region_pybedtools.sequence(fi=genome_fasta_with_vars)
        list_fasta = open(region_pybedtools.seqfn).read().split('\n')[:-1]
        sequences = [list_fasta[i].upper() for i in range(1, df.shape[0]*2, 2)]
        # #
        # string = sequences[0]
        # index = 999 - counter
        # new_char = alt
        # string_list = list(string)
        # string_list[index] = new_char
        # new_string = "".join(string_list)
        # sequences = [new_string]
        # counter +=1 
        # #
        label_encoder = preprocessing.LabelEncoder()
        integer_encoded = label_encoder.fit(['A', 'C', 'G', 'T', 'N'])
        integer_encoded = label_encoder.transform(['A', 'C', 'G', 'T', 'N'])
        #
        onehot_encoder  = preprocessing.OneHotEncoder(sparse_output=False)
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
        #
    datas_rc = np.flip(datas)
    #
    predictions = []
    for data in tqdm(datas):
        prediction = model.predict(data, verbose=0)
        predictions.append(prediction)
    #
    predictions_rc = []
    for data in tqdm(datas_rc):
        prediction = model.predict(data, verbose=0)
        predictions_rc.append(prediction)
    #
    ##my addition:
    #
    def flatten_list(nested_list):
        flat_list = []
        for sublist in nested_list:
            for item in sublist:
                flat_list.append(item)
        return flat_list
    #
    predictions_reformat = [l.tolist() for l in predictions]
    predictions_reformat = flatten_list(predictions_reformat)
    #
    predictions_rc_reformat = [l.tolist() for l in predictions_rc]
    predictions_rc_reformat = flatten_list(predictions_rc_reformat)
    #
    #
    target_names = obj.var.index.tolist()
    n_classes = len(obj.var.index.tolist())
    pred_df = pd.DataFrame(predictions_reformat, columns=target_names)
    pred_df.head()
    pred_df_rc = pd.DataFrame(predictions_rc_reformat, columns=target_names)[::-1]
    pred_df_rc.head()
    #pred_df.to_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/01_prediction40kb_forward.cvs', sep='\t')
    # pred_df.to_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/02_prediction25kb_forward.cvs', sep='\t')
    # pred_df.to_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/03_prediction27kb_forward.cvs', sep='\t')
    #pred_df_rc.to_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/01_prediction40kb_reverse.cvs', sep='\t')
    # pred_df_rc.to_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/02_prediction25kb_reverse.cvs', sep='\t')
    # pred_df_rc.to_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/03_prediction27kb_reverse.cvs', sep='\t')
    # pred_df = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/01_prediction40kb_forward.cvs', sep='\t', index_col=0)
    # pred_df = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/02_prediction25kb_forward.cvs', sep='\t', index_col=0)
    # pred_df = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/03_prediction27kb_forward.cvs', sep='\t', index_col=0)
    # pred_df_rc = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/01_prediction40kb_reverse.cvs', sep='\t', index_col=0)
    # pred_df_rc = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/02_prediction25kb_reverse.cvs', sep='\t', index_col=0)
    # pred_df_rc = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/03_prediction27kb_reverse.cvs', sep='\t', index_col=0)
    #pred_df = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/01_prediction40kb_forward.cvs', sep='\t')
    #pred_df.drop(['Unnamed: 0'], axis=1, inplace=True)
    #
    #uncomment for 1000bp:
    # pred_df = pred_df.iloc[399:600,:]
    #
    #pred_df_rc = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/01_prediction40kb_reverse.cvs', sep='\t')
    #pred_df_rc.drop(['Unnamed: 0'], axis=1, inplace=True)
    #
    #uncomment for 1000bp:
    # pred_df_rc = pred_df_rc.iloc[399:600,:]
    #
    #
    #smooth the signal:
    #
    #smooth the signal:
    #
    pred_df_min = pd.DataFrame()
    for col in pred_df.columns:
        tmp_value = np.min([pred_df[col], pred_df_rc[col]], axis=0)
        value = gaussian_filter1d(tmp_value, sigma=15)
        pred_df_min[col] = value
    #
    from IPython.display import Image
    from IPython.core.display import HTML 
    #
    # region01 = 'chr1:1379046-1418981'
    #
    for col in pred_df.columns:
        fig, ax = plt.subplots(3,1,figsize=(150,30)) #try making it longer, width is good
        ax[0] = pred_df[col].plot.line(title=col, fontsize=50, ax=ax[0])
        ax[0].get_xaxis().set_visible(False)
        ax[1] = pred_df_rc[col].plot.line(title=col, fontsize=50, ax=ax[1]) #try taking col out here
        ax[0].set_title(col, fontsize = 60)
        ax[1].set_title("")
        ax[2] = pred_df_min[col].plot.line(title=col, fontsize=50, ax=ax[2]) #try taking col out here
        #ax[2].axhline(0.65, c='red', linewidth=0.3)
        ax[2].set_title("")
        plt.show()
        fig.savefig('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/sigmoid_65_validation_rolling_window_fwd_rev_prediction_' + col + '_' + rsid + '_' + ref + '_' + alt + '_Enhancer_erythroid_pos_control_with_vars_peak_centered.png', bbox_inches='tight')
        print('@'*100)
    #
    for col in pred_df.columns:
        fig, ax = plt.subplots(1,1,figsize=(150,30))
        ax = pred_df_min[col].plot.line(title=col, fontsize=50, ax=ax, linewidth=2, color='blue')
        ax.get_xaxis().set_visible(False)
        ax.set_title(col, fontsize = 60)
        plt.show()
        fig.savefig('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_fwd_rev_rolling_window_ery_enhancer/sigmoid_65_validation_rolling_window_only_overlay_prediction_' + col + '_' + rsid + '_' + ref + '_' + alt + '_Enhancer_erythroid_pos_control_with_vars_peak_centered.png', bbox_inches='tight')
        print('@'*100)


for index, row in snp_dataset.iterrows():
    #print(index)
    rsid = row['rsid']
    regionOI = row['regionOI']
    alt = row['ALT']
    ref = row['REF']
    multi_allelic_visualisations(rsid, regionOI, ref, alt)
