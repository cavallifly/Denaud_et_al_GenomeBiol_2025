# Dependencies #
We suggest to install [conda](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html) and create an [environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html). 

Ensure that you have a working version of R. The scripts in this repository have been tested for version 4.0.5 (2021-03-31).
Install the following R packages: R.utils, Rsubread, xlsx, tidyverse, data.table, DESeq2, EnhancedVolcano.
The packages subread, samtools, and deeptools.

# Input data #
Next, you should download the .fastq files from the GEO entry GSE274159.

Now, you are ready to run the scripts in ./scripts. To do so, you can access the directory RNAseq_analysis of this repository using
```
cd RNAseq_analysis
```
and run one after the other the following commands following the messages:
```
scriptsDir=./scripts/

bash    ${scriptsDir}/01_fastaqc_analysis.sh                     &>> 01_fastaqc_analysis.out
bash    ${scriptsDir}/02_multiqc_analysis.sh                     &>> 02_multiqc_analysis.out
bash    ${scriptsDir}/03_alignReads_using_subread.sh             &>> 03_alignReads_using_subread.out
bash    ${scriptsDir}/04_indexBam_using_samtools.sh              &>> 04_indexBam_using_samtools.out
bash    ${scriptsDir}/05_generationCountTable_using_subreads.sh  &>> 05_generationCountTable_using_subreads.out
Rscript ${scriptsDir}/06_DEanalysis_using_DEseq2.R               &>> 06_DEanalysis_using_DEseq2.out
bash    ${scriptsDir}/07_merge_and_generate_bigWigs.sh           &>> 07_merge_and_generate_bigWigs.out
```

Once the scripts are finished, you will obtain the panels and the data points used to obtain them for all the Figures in the paper. 
To obtain the final version of the figures the panels have been assembled using the PowerPoint program.

## Contributions ##
The code in this repository has been developed at the [Cavalli Lab](https://www.igh.cnrs.fr/en/research/departments/genome-dynamics/chromatin-and-cell-biology) with the contributions of Marco Di Stefano, Gonzalo Sabaris, and Giorgio L. Papadopoulos.
