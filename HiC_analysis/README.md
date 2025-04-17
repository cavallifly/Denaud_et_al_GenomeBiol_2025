# Dependencies #
We suggest to install [conda](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html) and create an [environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)

Ensure that you have a working version of R. The scripts in this repository have been tested for version 4.0.5 (2021-03-31).
Install the following R packages: misha, shaman, ggpubr, dplyr, doParallel.
Install the following bash commands: convert (part of imagemagick) and identify.
For example, you can use the commands:
```
conda install -y bioconda::r-misha
conda install -y bioconda::r-shaman
conda install -y conda-forge::r-ggpubr
conda install -y conda-forge::r-dplyr
conda install -y conda-forge::r-doparallel
conda install -y conda-forge::imagemagick
conda install -y conda-forge::identify
```

Download the script [hicpro2huicebox.sh](https://github.com/nservant/HiC-Pro/blob/master/bin/utils/hicpro2juicebox.sh) from the hicpro repository and the juicer_tools_1.22.01.jar using 
```
wget https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar
```
and place them in the ./scripts/ folder.

# Input data #
Next, you should download the misha tracks from GEO for the observed counts in the directory HiC_analysis/mishaDB/trackdb/dm6/tracks/hic/ of this repository.
- GSM7888066	Larv WT HiC Repli1
- GSM7888067	Larv WT HiC Repli2
- GSM7888068	Larv WT HiC Repli3
- GSM7888077	Larv gypsy2 insertion HiC Repli1
- GSM7888078	Larv gypsy2 insertion HiC Repli2
- GSM7888079	Larv gypsy2 insertion HiC Repli3
- GSM8444487     Larv Double HiC Repli1
- GSM8444488     Larv Double HiC Repli2
- GSM8444489     Larv Double HiC Repli3
- GSM8444490     Larv deltaPRE1 HiC Repli1
- GSM8444491     Larv delatPRE1 HiC Repli2
- GSM8444492     Larv PRE1_UP HiC Repli1
- GSM8444493     Larv PRE1_UP HiC Repli2

These tracks have been obtained using the "scHiC2" pipeline in T. Nagano et al., Cell-cycle dynamics of chromosomal organization at single-cell resolution. Nature 547, 61–67 (2017). For the
samples of the PRE1_UP condition, you should run the following command to generate the modified dm6 assembly:
```
cd HiC_analysis
bash     ./scripts/01_editFasta.sh &> 01_editFasta.out
```

Now, you are ready to run the scripts in HiC_analysis/scripts. To do so, you can access the directory HiC_analysis of this repository using
```
cd HiC_analysis
```
and run one after the other the following commands:
```
scriptsDir=./scripts/

# Compute the shuffled (expected) tracks and the shaman Hi-C scores
Rscript ${scriptsDir}/02_generateShuffleTrack.R &>> 02_generateShuffleTrack.out
Rscript ${scriptsDir}/03_computeScoreTrack.R    &>> 03_computeScoreTrack.out

# Analysis of the Hi-C maps
Rscript ${scriptsDir}/04_computeInsulationTracks.R        &>> 04_computeInsulationTracks.out
Rscript ${scriptsDir}/05_computeDownSampledTracks.R       &>> 05_computeDownSampledTracks.out
Rscript ${scriptsDir}/06_computeDiffScoreTrack_parallel.R &>> 06_computeDiffScoreTrack_parallel.out

# Generate data and Figures for each figure panel
Rscript ${scriptsDir}/07_ScoreMapsPlot_Fig1C.R                 &>> 07_ScoreMapsPlot_Fig1C.out
Rscript ${scriptsDir}/08_INSplot_Fig1D.R                       &>> 08_INSplot_Fig1D.out
Rscript ${scriptsDir}/09_INSquantification_Fig1E.R             &>> 09_INSquantification_Fig1E.out
Rscript ${scriptsDir}/10_DiffScoreMapsPlot_Fig1F.R             &>> 10_DiffScoreMapsPlot_Fig1F.out
Rscript ${scriptsDir}/11_ScoreMapsQuantification_Fig1G.R       &>> 11_ScoreMapsQuantification_Fig1G.out
Rscript ${scriptsDir}/12_ScoreMapsPlot_Fig2C.R                 &>> 12_ScoreMapsPlot_Fig2C.out
Rscript ${scriptsDir}/13_INSplot_Fig2D.R                       &>> 13_INSplot_Fig2D.out
Rscript ${scriptsDir}/14_INSquantification_Fig2E.R             &>> 14_INSquantification_Fig2E.out
Rscript ${scriptsDir}/15_DiffScoreMapsPlot_Fig2G.R             &>> 15_DiffScoreMapsPlot_Fig2G.out
Rscript ${scriptsDir}/16a_ScoreMapsQuantification_Fig2Fa.R     &>> 16a_ScoreMapsQuantification_Fig2Fa.out
Rscript ${scriptsDir}/16b_ScoreMapsQuantification_Fig2Fb.R     &>> 16b_ScoreMapsQuantification_Fig2Fb.out
bash    ${scriptsDir}/17_retrieve_points_from_mishaDB.sh       &>> 17_retrieve_points_from_mishaDB.out
bash    ${scriptsDir}/18_convert_ints2hic.sh                   &>> 18_convert_ints2hic.out 
bash    ${scriptsDir}/19_dump_scaleNormalizedMaps.sh           &>> 19_dump_scaleNormalizedMaps.out
bash    ${scriptsDir}/20_plot_scaleNormalizedMaps.sh           &>> 20_plot_scaleNormalizedMaps.out

# Obtain the triangular maps
bash ${scriptsDir}/21_generate_triangular_pngs.sh &> 21_generate_triangular_pngs.out

# Move results in folders
mkdir -p Data_for_figures
mkdir -p Figure_panels
mkdir -p intsFiles
mkdir -p hicFiles

mv *.png *.pdf Figure_panels
mv *.tab *.tsv Data_for_figures/
mv *.ints intsFiles
mv *.hic  hicFiles
```

Once the scripts are finished, you will obtain the panels and the data points used to obtain them for all the Figures in Denaud_at_al_2024.
To obtain the final version of the figures the panels have been assembled using the PowerPoint program.

## Contributions ##
The code in this repository has been developed at the [Cavalli Lab](https://www.igh.cnrs.fr/en/research/departments/genome-dynamics/chromatin-and-cell-biology) with the contributions of Marco Di Stefano, Gonzalo Sabaris, and Giorgio L. Papadopoulos.
