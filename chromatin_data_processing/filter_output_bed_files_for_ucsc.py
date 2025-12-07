import pandas as pd
import os
import numpy as np
import pickle
import glob

bed_files_path = '/well/PROCARDIS/domwest/upstr_processing/UpStreamPipeline_success_single/genetics/CATCH-UP/analysis/results/09_peak_calling/*.bed'
indiv_beds = glob.glob(bed_files_path)
for file in indiv_beds:
	filename = file.split('/')[11]
	df = pd.read_csv(file, sep='\t')
	filt_df = df[df['overall_peak_score']>=0.5]
	filt_df = filt_df.sort_values(by=['chrom', 'start'], ascending=[True, True])
	filt_df['filt_score'] = filt_df['overall_peak_score']*1000
	filt_df.rename({'chrom':'col1', 'start':'col2', 'end':'col3', 'filt_score':'col5', 'enrichment_score':'col8'}, axis=1, inplace=True)
	filt_df['col4'] = '.'
	filt_df['col6'] = '.'
	filt_df['col7'] = '1.0'
	filt_df['col9'] = '-1'
	filt_df['col10'] = '-1'
	filt_df = filt_df[['col1', 'col2', 'col3', 'col4', 'col5', 'col6', 'col7', 'col8', 'col9', 'col10']]
	filt_df.to_csv('filt_' + filename, index=False, sep='\t', header=False)
	

#add to the top of each file: 
#track type=narrowPeak visibility=3 db=hg38 name="nPk" description="filtered for >=0.5 LanceOtron score peaks" useScore=1

cd /well/PROCARDIS/domwest/upstr_processing/UpStreamPipeline_success_single/genetics/CATCH-UP/analysis/results/09_peak_calling/
for i in filt_*.bed; do
	sed -i '1i track type=narrowPeak visibility=3 db=hg38 name="nPk" description="filtered for >=0.5 LanceOtron score peaks" useScore=1' $i
done




