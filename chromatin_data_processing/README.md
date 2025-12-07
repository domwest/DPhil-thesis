**How to:**

**Important links:**
- [https://github.com/Genome-Function-Initiative-Oxford/UpStreamPipeline]
- [https://github.com/Genome-Function-Initiative-Oxford/UpStreamPipeline/tree/main/genetics/CATCH-UP]

- 1). The first link explains how to set up the upstream pipeline conda environment.
- 2). The second link details how to use the CATCHUP pipeline (for chromatin data QC and standardisation) within the upstream environment.

**Methodology as described in thesis:**

> This pipeline was run as a precursor to DeepHaem with the aim to standardise and QC the data beforehand.

> Open chromatin data from bulk tissue heart left ventricle as well as single-cell heart datasets (including cardiac muscle cell, fibroblast, lymphocyte, macrophage, ventricular cardiomyocyte, atrial cardiomyocyte, adipocyte, smooth muscle cell, and endothelial cell) were downloaded from various databank sources. Most of the data was downloaded from ENCODE (Dunham et al., 2012), some from the cardiac atlas of regulatory elements (CARE) portal, and one dataset from the national center for biotechnology information sequence read archive (NCBI SRA). 

> The following methodologies were captured across the various datasets: ATAC-Seq, DNase-Seq, TF ChIP-seq, and Histone ChIP-seq (H3K27ac, H3K27me3, H3K36me3, H3K4me1, H3K4me3, H3K9me3). All TF ChIP-seq datasets were CTCF transcription factor site data. See the tracking sheet at [github.com/domwest/DPhil-thesis/blob/main/data_tracking_sheet_for_thesis.xlsx] for dataset related details.

> Both the command-line and web version of CATCH-UP were employed depending on the dataset file formats. Input files need to be in the form of fastq files for the command-line version to run. A combination of single and paired-end fastq files were downloaded wherever possible and used as input to the command-line version of this pipeline. In the event that fastq files were not available but bed and bigwig files were, LanceOtron web for peak-calling was run on the multi-locus view (MLV) interface (Sergeant et al., 2021) instead. The command line version was carried out on the BMRC cluster by accessing the Upstream pipeline conda environment where all relevant packages were pre-installed.

> Before running CATCH-UP, the reference genome (build 38) was downloaded directly from UCSC and indexed by Bowtie2 (integrated in CATCH-UP). When executing the CATCH-UP pipeline, certain parameters needed to be introduced to inform the pipeline as to whether concatenation of lanes, merging of samples, and/or adapter trimming was required. Lane concatenation was not implemented as there were no sample replicates in the input and adapter trimming was not performed due to the use of the reliable aligner, bowtie. The files were merged on single versus paired end sequencing, bulk tissue versus individual cell type data, and according to the methodology employed (ATAC versus DNase etc.). In addition, the files were merged across donors. The pipeline output bed, bigwig, and bam files per merged dataset.

**Results as described in thesis:**

> Having selected the heart tissue and cell type datasets for inclusion as input to machine learning (post-cell type enrichment), these datasets needed to undergo quality control and standardisation by an upstream pipeline called CATCH-UP (with LanceOtron integrated) before being fed to machine learning for functional variant prioritisation. This was essential for noise filtering as well as peak inspection – ensuring accurate downstream DeepHaem predictions.

> The quality control and standardisation of ATAC-Seq, ChIP-seq and DNase-Seq data found LanceOtron to surpass the established gold standard tools through enhancing selectivity and obtaining flawless sensitivity (Hentges et al., 2022). For this reason, LanceOtron was implemented for peak calling across the datasets.

> Running CATCH-UP led to more refined and therefore accurate chromatin datasets as shown in the example below. In figure 31, the cardiac muscle cell DNase peaks showed a cleaner track with more stringently assigned peaks after running CATCH-UP with LanceOtron (overall peak scores ≥ 0.5) compared to before. Overall peak scores ≥ 0.5 were chosen as this is the widely accepted threshold (as seen by Gao et al. (2023) and Crump et al. (2023)).

<img width="1920" height="864" alt="catchup" src="https://github.com/user-attachments/assets/fafd1564-fb54-4d78-bc83-9d8693b19190" />


**A few additional notes:**

If copying an UpstreamPipeline folder to re-do some kind of analysis (eg by tweaking fastq files to run on) then make sure to delete the 'analysis' folder containing all the previous results. Otherwise this won't work

Commands I used to execute the analysis:
- If running the job for the first time:
```
conda activate upstream
snakemake --configfile=config/analysis.yaml all --cores 8 --rerun-incomplete
```
- If running the job after it was interrupted ie through connection being lost:
```
conda activate upstream
snakemake --configfile=config/analysis.yaml all --cores 8 --rerun-incomplete --unlock
snakemake --configfile=config/analysis.yaml all --cores 8 --rerun-incomplete
```
