#- required libraries

if (require("mgcv")==F) install.packages("mgcv")
if (require("scales")==F) install.packages("scales")
if (require("gplots")==F) install.packages("gplots")
if (require("magicaxis")==F) install.packages("magicaxis")
if (require("lubridate")==F) install.packages("lubridate")
if (require("doBy")==F) install.packages("doBy")
if (require("Hmisc")==F) install.packages("Hmisc")
if (require("zoo")==F) install.packages("zoo")
if (require("hexbin")==F) install.packages("hexbin")
if (require("nlme")==F) install.packages("nlme")
if (require("lsmeans")==F) install.packages("lsmeans")
if (require("car")==F) install.packages("car")
if (require("data.table")==F) install.packages("data.table")
if (require("calibrate")==F) install.packages("calibrate")


library(mgcv,quietly=T)
library(scales,quietly=T)
library(gplots,quietly=T)
library(magicaxis,quietly=T)
library(lubridate,quietly=T)
library(doBy,quietly=T)
library(Hmisc,quietly=T)
library(zoo,quietly=T)
library(hexbin,quietly=T)
library(nlme,quietly=T)
library(lsmeans,quietly=T)
library(car,quietly=T)
library(data.table,quietly=T)
library(calibrate,quietly=T)

#- the following libraries aren't on CRAN, but can be installed from github or bitbucket with devtools
if (require("devtools")==F) install.packages("devtools")
library(devtools,quietly=T)

if (require("plantecophys")==F) install_bitbucket("remkoduursma/plantecophys")
if (require("HIEv")==F) install_bitbucket("remkoduursma/HIEv")
if (require("plotBy")==F){ install_bitbucket("remkoduursma/plotBy")}
if (require("piecewiseSEM")==F) install_github("jslefche/piecewiseSEM")

library(piecewiseSEM,quietly=T) # for estimating r2 value in mixed-effects models
library(plantecophys,quietly=T) # for modeling leaf-level gas exchange
library(HIEv,quietly=T)  
library(plotBy,quietly=T)
