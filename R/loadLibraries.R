#- required libraries

if (require("mgcv",quietly=T)==F) install.packages("mgcv")
if (require("scales",quietly=T)==F) install.packages("scales")
if (require("gplots",quietly=T)==F) install.packages("gplots")
if (require("magicaxis",quietly=T)==F) install.packages("magicaxis")
if (require("lubridate",quietly=T)==F) install.packages("lubridate")
if (require("doBy",quietly=T)==F) install.packages("doBy")
if (require("Hmisc",quietly=T)==F) install.packages("Hmisc")
if (require("zoo",quietly=T)==F) install.packages("zoo")
if (require("hexbin",quietly=T)==F) install.packages("hexbin")
if (require("nlme",quietly=T)==F) install.packages("nlme")
if (require("lsmeans",quietly=T)==F) install.packages("lsmeans")
if (require("car",quietly=T)==F) install.packages("car")
if (require("data.table",quietly=T)==F) install.packages("data.table")
if (require("calibrate",quietly=T)==F) install.packages("calibrate")

# 
# library(mgcv,quietly=T)
# library(scales,quietly=T)
# library(gplots,quietly=T)
# library(magicaxis,quietly=T)
# library(lubridate,quietly=T)
# library(doBy,quietly=T)
# library(Hmisc,quietly=T)
# library(zoo,quietly=T)
# library(hexbin,quietly=T)
# library(nlme,quietly=T)
# library(lsmeans,quietly=T)
# library(car,quietly=T)
# library(data.table,quietly=T)
# library(calibrate,quietly=T)

#- the following libraries aren't on CRAN, but can be installed from github or bitbucket with devtools
if (require("devtools")==F) install.packages("devtools")
library(devtools,quietly=T)

if (require("plantecophys")==F) install_bitbucket("remkoduursma/plantecophys")
if (require("HIEv")==F) install_bitbucket("remkoduursma/HIEv")
if (require("plotBy")==F){ install_bitbucket("remkoduursma/plotBy")}
if (require("piecewiseSEM")==F) install_github("jslefche/piecewiseSEM")

# library(piecewiseSEM,quietly=T) # for estimating r2 value in mixed-effects models
# library(plantecophys,quietly=T) # for modeling leaf-level gas exchange
# library(HIEv,quietly=T)  
# library(plotBy,quietly=T)
