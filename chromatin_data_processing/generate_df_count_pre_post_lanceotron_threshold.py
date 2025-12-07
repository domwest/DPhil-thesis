import pandas as pd
import os
import numpy as np
import pickle
import glob

bed_files_path = '/well/PROCARDIS/domwest/upstr_processing/UpStreamPipeline_success_single_merged_SC/genetics/CATCH-UP/analysis/results/09_peak_calling/*.bed'
indiv_beds = glob.glob(bed_files_path)
counts_df = pd.DataFrame()
filenames = []
original_lengths = []
new_lengths = []
for file in indiv_beds:
	filename = file.split('/')[11]
	filename_edited = filename.split('.')[0]
	print(filename_edited)
	filenames.append(filename)
	df = pd.read_csv(file, sep='\t')
	original_length = df.shape[0]
	original_lengths.append(original_length)
	filt_df = df[df['overall_peak_score']>=0.5]
	filt_df = filt_df[['chrom', 'start', 'end']]
	filt_df.to_csv('/well/PROCARDIS/domwest/upstr_processing/UpStreamPipeline_success_single_merged_SC/genetics/CATCH-UP/analysis/results/09_peak_calling/' + filename_edited + '_post_filtering_peaks.bed', sep='\t', header=None)
	new_length = filt_df.shape[0]
	new_lengths.append(new_length)

counts_df['filename'] = filenames
counts_df['data'] = counts_df['filename'].str.split('_Don', expand=True)[0]
counts_df['pre_0.5_threshold'] = original_lengths
counts_df['post_0.5_threshold'] = new_lengths
counts_df['%_leftover'] = (counts_df['post_0.5_threshold']/counts_df['pre_0.5_threshold'])*100
counts_df = counts_df[['data', 'filename', 'pre_0.5_threshold', 'post_0.5_threshold', '%_leftover']]
counts_df = counts_df.sort_values(by=['data', 'filename', '%_leftover'], ascending=[True, True, False])
counts_df['cell_type_or_tissue'] = 'cardiac muscle cell'

counts_df.to_csv('/well/PROCARDIS/domwest/upstr_processing/UpStreamPipeline_success_single_merged_SC/genetics/CATCH-UP/analysis/results/09_peak_calling/peak_counts_single_SC_merged_data.csv', index=False)
