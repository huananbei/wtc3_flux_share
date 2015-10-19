#- required libraries

if (!require("mgcv")) install.packages("mgcv")
if (!require("scales")) install.packages("scales")
if (!require("gplots")) install.packages("gplots")
if (!require("magicaxis")) install.packages("magicaxis")
if (!require("lubridate")) install.packages("lubridate")
if (!require("doBy")) install.packages("doBy")
if (!require("Hmisc")) install.packages("Hmisc")
if (!require("zoo")) install.packages("zoo")
if (!require("hexbin")) install.packages("hexbin")
if (!require("nlme")) install.packages("nlme")
if (!require("lsmeans")) install.packages("lsmeans")
if (!require("car")) install.packages("car")
if (!require("data.table")) install.packages("data.table")
if (!require("calibrate")) install.packages("calibrate")


library(mgcv)
library(scales)
library(gplots)
library(magicaxis)
library(lubridate)
library(doBy)
library(Hmisc)
library(zoo)
library(hexbin)
library(nlme)
library(lsmeans)
library(car)
library(data.table)
library(calibrate)

#- the following libraries aren't on CRAN, but can be installed from github or bitbucket with devtools
if (!require("devtools")) install.packages("devtools")
if (!require("plantecophys")) install_bitbucket("remkoduursma/plantecophys")
if (!require("HIEv")) install_bitbucket("remkoduursma/HIEv")
if (!require("plotBy")) install_bitbucket("remkoduursma/plotBy")
if (!require("piecewiseSEM")) install_github("jslefche/piecewiseSEM")
library(devtools)
library(piecewiseSEM) # for estimating r2 value in mixed-effects models
library(plantecophys) # for modeling leaf-level gas exchange
library(HIEv)  
library(plotBy)