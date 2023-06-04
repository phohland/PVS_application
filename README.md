# hRELSA_application

This GitHub repository enables the user to apply the hRELSA code with the given data as shown in the initial publication.

## Table of contents

- [Get started](#Get started)
- [Output and figures](#Output and figures)

## Get started

The raw data used for calculation can be found in output/raw.csv. If you are interested in how the raw data table was created, you can take a look into R/data_setup.R. For further calculation and analysis, you can start the R script R/main.R.

In the main R script you can see step by step, how the results as shown in the initial publication were calculated. The script calls different functions from the hRELSA package. The hRELSA package can be found in the regarding GitHub repository: https://github.com/phohland/hRELSA.

The main R script itself generates both the figures and the tables containing the information about the calculations.

## Outputs and figures

The output folder contains both dat.csv and pre.csv. The dat file contains a formatted version of the raw data, the pre file contains the normalized version of the data. To calculate the final severity weights and severity values, calculate it yourself with the main R script. A file called final.csv will be created inside of the output folder showing all severity calculations based on the previous data sets. It may be a bit bigger.
