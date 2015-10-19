# wtc3_flux_share
Whole tree flux dataset and processing scripts. Hawkesbury Institute for the Environment, University of Western Sydney.

This repository contains the data and R code to replicate the figures and analyses presented in the following manuscript. 

"Does thermal acclimation to climate warming stabilize the ratio of canopy respiration to photosynthesis?"
Drake JE, Tjoelker MG, Aspinwall MJ, Reich PB, Barton CVM, Medlyn BE, Duursma RA 

The data files are moderately large (~25 MB total), so cloning this repository may take longer than is typical. I recommend that you clone this repository into an Rstudio project by copying the SSH link on the right, and pasting the link into a new Rstudio project URL (File > New Project > Version Control > Git). This process is described in useful detail here:  http://www.molecularecologist.com/2013/11/using-github-with-r-and-rstudio/ . The data will be downloaded as csv files in the folder "data/" and three R-scripts will be downloaded into the folder "R/". The "main_script.R" contains code to recreate all of the manuscript figures; set the "export" variable to "T" to create pdfs in folder "output/". The "statistical_analysis.R" script contains statistical analyses to recreate Table 1; note that these take considerable time to run.

This code was developed in R v3.2.2 "Fire Safety".
