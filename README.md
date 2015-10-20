# wtc3_flux_share
Whole tree flux dataset and processing scripts. This repository contains the data and R code to replicate the figures and analyses presented in the following manuscript:

"Does thermal acclimation to climate warming stabilize the ratio of canopy respiration to photosynthesis?"
Drake JE, Tjoelker MG, Aspinwall MJ, Reich PB, Barton CVM, Medlyn BE, Duursma RA 

There are two ways to download and use this code.
(1) I recommend that you clone this repository into an Rstudio project by copying the SSH link on the right, and pasting the link into a new Rstudio project URL (File > New Project > Version Control > Git). You will also need git. This process is described in useful detail here:  http://www.molecularecologist.com/2013/11/using-github-with-r-and-rstudio/ . The data files are moderately large (~32 MB total), so cloning this repository may take longer than is typical.

(2) You could choose to download a zip file of this repository using the link on the right. Upzip it in a location of your choosing and update the setwd() command on line 20 of “main_script.R” accordingly.

The data can be found as csv files in the folder "data/", and three R-scripts will be in the folder "R/". The "main_script.R" contains code to recreate all of the manuscript figures; set the "export" variable to "T" to create pdfs in folder "output/". The "statistical_analysis.R" script contains statistical analyses to recreate Table 1; note that these take considerable time to run. The "functions.R" script contains all of the custom analysis and plotting functions that do all of the actual work; these functions are called by "main_script.R".

This code was developed in R v3.2.2 "Fire Safety" with Rstudio version 0.99.669 on Windows 7. It has been successfully tested on several machines. However the code would need modification to run on other platforms (e.g. OSX users would need to replace windows() calls with X11() ).
