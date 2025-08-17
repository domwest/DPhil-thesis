#this paper: https://www.sciencedirect.com/science/article/pii/S0092867416314635?via%3Dihub#app2

#msgs form matt:
# So - I don't really work on the RBC stuff - but I think people who do use GWAS from a combination of 3 sources - but the primary one is the Astle 2016 paper:
# https://pubmed.ncbi.nlm.nih.gov/27863252/
# The other studies are van der Harst 2012 and Vuckovic 2020, but I think you should be fine just to use the Astle dataset
# I'm sure Simone has a vcf of the snps saved somewhere but I'm not sure where to look - If you are under time pressure you could look on GWAS catalogue for the Astle dataset

# I downloaded Table S4 (Summary of Associated Variants and Their Consequences) - there seem to be about ~6700 variants
# What are the erythrocyte indices?
# The indices include: Average red blood cell size (MCV) Hemoglobin amount per red blood cell (MCH) The amount of hemoglobin relative to the size of the cell (hemoglobin concentration) per red blood cell (MCHC)
#so filtered for these and end up with: 

# #first in python:
# import pandas as pd
# df = pd.read_excel('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/erythroid_pos_control_vars_metadata.xlsx')
# df_filt = df[df['Associated Blood Index'].isin(['MLV', 'MCH', 'MCHC'])]
# df_filt2.rename({'Chr (GRCh37)':'CHROM', 'BP (GRCh37)':'POS'}, axis=1, inplace=True)
# df_filt2.to_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/erythroid_pos_control_vars_for_liftover_to_38.csv', index=False)

# #first in python but redo:
# import pandas as pd
# df = pd.read_excel('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/erythroid_pos_control_vars_metadata.xlsx')
# # df_filt = df[df['Associated Blood Index'].isin(['MLV', 'MCH', 'MCHC'])]
# df.rename({'Chr (GRCh37)':'CHROM', 'BP (GRCh37)':'POS'}, axis=1, inplace=True)
# df[['Unique Variant ID']].drop_duplicates() #3755 unique variants
# df.to_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/erythroid_pos_control_vars_for_liftover_to_38.csv', index=False) #6736

# #first get build37 to build38:

# #CONVERT FROM B37 TO B38

# # [wiq135@rescomp2 input_files_for_deephaem]$ module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2
# # [wiq135@rescomp2 input_files_for_deephaem]$ module load R/4.2.2-foss-2022b
# # R

# # install and load rtracklayer
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("rtracklayer")

# library(rtracklayer)
# library(dplyr)

# # import the chain file

# # download.file("https://hgdownload.soe.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz", "hg19ToHg38.over.chain.gz")
# # system("gzip -d hg19ToHg38.over.chain.gz")
# # #import.chain("hg19ToHg38.over.chain.gz")
# chainObject <- import.chain("hg19ToHg38.over.chain")

# tmp = read.csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/erythroid_pos_control_vars_for_liftover_to_38.csv')
# data = list()
# for(i in 1:nrow(tmp)) {
#     row <- tmp[i,]
#     chrom <- paste(c("chr", row$CHROM), collapse = "")
#     x <- row$POS
#     y <- row$POS
#     grObject <- GRanges(seqnames=chrom, ranges=IRanges(start=x, end=y))
#     #print(grObject)
#     # run liftOver
#     results <- as.data.frame(liftOver(grObject, chainObject))
#     if (nrow(results) > 0) {
#         names(results)[names(results) == "start"] <- "grch38_POS"
#         to_merge <- select(results, grch38_POS)
#         #
#         row$grch38_POS <- to_merge$grch38_POS
#         data <- append(data, row)
#         }
#     else {
#         names(results)[names(results) == "start"] <- "grch38_POS"
#         grch38_POS <- c('')
#         results[1:1,"grch38_POS"] = ''
#         to_merge <- select(results, grch38_POS)
#         #
#         row$grch38_POS <- to_merge$grch38_POS
#         data <- append(data, row)
#         }
#     }
# fin2 <- data.frame(matrix(unlist(data), nrow=nrow(tmp), byrow=TRUE))
# columns <- colnames(tmp)
# columnss <- c(columns, 'grch38_POS')
# colnames(fin2) <- columnss   
# write.csv(fin2, '/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/erythroid_pos_control_vars_post_liftover.csv', row.names = FALSE)

# #now again in python:
# import pandas as pd
# df = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/erythroid_pos_control_vars_post_liftover.csv')
# df.rename({'rsID..where.available.':'rsid', 'CHROM':'chr'}, axis=1, inplace=True)
# df[['REF', 'ALT']] = df['REF.ALT..GRCh37.'].str.split('/', expand=True)
# df = df[['rsid', 'chr', 'grch38_POS', 'REF', 'ALT']].drop_duplicates()
# df[df['grch38_POS'].isnull()]
# #           rsid  chr  grch38_POS REF ALT
# # 3654  rs687289    9         NaN   G   A
# df = df[df['grch38_POS'].notnull()]
# df['grch38_POS'] = df['grch38_POS'].astype(int)
# df.to_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/bi_allelic_3754_erythroid_pos_control_snp_dataset.csv', index=False)

#they all appear to be bi-allelic

#new actual erythroid positive control files: CTCF, CTCF, CTCF
#we want:
import pandas as pd 
import os 
import glob 
import numpy as np  

banger_files = glob.glob('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/*bangers*')
for file in banger_files:
    typ = file.split('/')[-1].split('.bangers.tsv')[0]
    df = pd.read_csv(file, sep='\t')
    df = df[['TEST.SNP.ID', 'TEST.SNP.CHROM', 'TEST.SNP.POS', 'TEST.SNP.REF.ALLELE', 'TEST.SNP.ALT.ALLELE']]
    df['TEST.SNP.CHROM'] = df['TEST.SNP.CHROM'].str.replace('chr','')
    df.rename({'TEST.SNP.ID':'rsid', 'TEST.SNP.CHROM':'chr', 'TEST.SNP.POS':'grch38_POS', 'TEST.SNP.REF.ALLELE':'REF', 'TEST.SNP.ALT.ALLELE':'ALT'}, axis=1, inplace=True)
    df.to_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/bi_allelic_' + typ + '_erythroid_pos_control_snp_dataset.csv', index=False)

##BI-ALLELIC ONLY

import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
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
    os.mkdir("tmp")
pybedtools.set_tempdir('tmp')

genome_fasta_normal = "/well/PROCARDIS/domwest/upstr_processing/UpStreamPipeline_old/reference_genomes/upstream_pipelines_output/Download_Index_Genomes/hg38_bowtie2/ref_genome_analysis/results/hg38.fa"

##Center - Hidden layers

#Let's do regions of 1000bp each:

snp_dataset = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/bi_allelic_CTCF_erythroid_pos_control_snp_dataset.csv')
snp_dataset['chr'] = 'chr' + snp_dataset['chr'].astype(str)
snp_dataset['downstream'] = snp_dataset['grch38_POS'] - 500
snp_dataset['upstream'] = snp_dataset['grch38_POS'] + 500
snp_dataset['regionOI'] = snp_dataset['chr'].astype(str) + ':' + snp_dataset['downstream'].astype(str) + '-' + snp_dataset['upstream'].astype(str)

regionOIs_df = []

for index, row in snp_dataset.iterrows():
    #print(index)
    rsid = row['rsid']
    regionOI = row['regionOI']
    alt = row['ALT']
    ref = row['REF']
    multi_allelic_predictions(rsid, regionOI, ref, alt)

def multi_allelic_predictions(rsid, regionOI, ref, alt):
    #
    chrom = regionOI.split(":")[0]
    start = regionOI.split(":")[-1].split("-")[0]
    end   = regionOI.split(":")[-1].split("-")[-1]
    print(chrom, start, end)
    #
    starts = np.array([int(regionOI.split(':')[1].split('-')[0])])
    ends = np.array([int(regionOI.split(':')[1].split('-')[1])])
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
    predictions = []
    #
    for region in tqdm(regions):
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
        regions = []
        dfT = df.T
        for i in dfT:
            regions.append(dfT[i].tolist())
        region_pybedtools = pybedtools.BedTool(regions)
        region_pybedtools = region_pybedtools.sequence(fi=genome_fasta_normal)
        list_fasta = open(region_pybedtools.seqfn).read().split('\n')[:-1]
        sequences = [list_fasta[i].upper() for i in range(1, df.shape[0]*2, 2)]
        #comment below out for generating without_vars
        # string = sequences[0]
        # index = 499
        # new_char = alt
        # string_list = list(string)
        # string_list[index] = new_char
        # new_string = "".join(string_list)
        # sequences = [new_string]
        #comment above out for generating without_vars
        label_encoder = preprocessing.LabelEncoder()
        integer_encoded = label_encoder.fit(['A', 'C', 'G', 'T', 'N'])
        integer_encoded = label_encoder.transform(['A', 'C', 'G', 'T', 'N'])
        onehot_encoder  = preprocessing.OneHotEncoder(sparse=False)
        integer_encoded = integer_encoded.reshape(len(integer_encoded), 1)
        onehot_encoded = onehot_encoder.fit(integer_encoded)
        onehot_encoded = onehot_encoder.transform(integer_encoded)
        data = []
        for x in sequences:
            integer_encoded = label_encoder.transform([*x])
            integer_encoded = integer_encoded.reshape(len(integer_encoded), 1)
            onehot_encoded = onehot_encoder.transform(integer_encoded)
            data.append(np.delete(onehot_encoded, 3, 1))
        data = np.array(data)
        prediction = model.predict(data, verbose=0)
        predictions.append(prediction)
    #
    pred_df = []
    for pred in predictions:
        p = pred[0]
        pred_df.append(p.tolist())
    #
    target_names = obj.var.index.tolist()
    n_classes = len(obj.var.index.tolist())
    #
    pred_df = pd.DataFrame(pred_df, columns=target_names)
    pred_df['variant'] = rsid
    #
    pred_df['var_type'] = 'bi'
    pred_df['REF'] = ref 
    pred_df['ALT'] = alt
    #
    regionOIs_df.append(pred_df)

regionOIs_df = pd.concat(regionOIs_df)

# regionOIs_df.to_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_shortlist_vars_prediction/with_vars_present_CTCF_erythroid_pos_control_snp_prediction.csv', sep='\t')
regionOIs_df.to_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/redo_shortlist_vars_prediction/without_vars_present_CTCF_erythroid_pos_control_snp_prediction.csv', sep='\t')


################################################################################################################################################################################
