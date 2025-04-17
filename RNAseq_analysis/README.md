# Dependencies #
We suggest to install [conda](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html) and create an [environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)

# Input data #
Next, you should download the misha tracks from GEO for the .fastq files.
- GSM8444494	WT_Repli1
- GSM8444495	WT_Repli2
- GSM8444496	WT_Repli3
- GSM8444497  Double Repli1
- GSM8444498  Double Repli2
- GSM8444499  Double Repli3

Now, you are ready to run the scripts in ./scripts. To do so, you can access the directory RNAseq_analysis of this repository using
```
cd RNAseq_analysis
```
and run one after the other the following commands:
```
scriptsDir=./scripts/

# 
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
