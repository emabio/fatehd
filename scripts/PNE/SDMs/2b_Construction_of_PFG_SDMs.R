### HEADER #####################################################################
##' @title Build PFG SDMs from individual species SDMs
##'
##' @author Damien G. & Maya G.
##' @contact damien.georges2 at gmail.com
##' 
##' @date 22/01/2016
##' 
##' @description Here we will build PFGs SDMs fom individual determinant species 
##'   SDMs (~ 270 determinant species) that have been built via the script
##'   "2a_Application_biomod2_species_by_species.R".
##'   We will test several way to construct the meta SDMs:
##'     - max of individual sdms
##'     - quantiles of individual sdms
##'     - based on weigted mean or on comittee averaging
##'   
##' @log 
##' 
##' @licencing GPL
##'     Copyright (C) 2015  Damien G.
##' 
##'     This program is free software: you can redistribute it and/or modify
##'     it under the terms of the GNU General Public License as published by
##'     the Free Software Foundation, either version 3 of the License, or
##'     (at your option) any later version.
##' 
##'     This program is distributed in the hope that it will be useful,
##'     but WITHOUT ANY WARRANTY; without even the implied warranty of
##'     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##'     GNU General Public License for more details.
##' 
##'     You should have received a copy of the GNU General Public License
##'     along with this program.  If not, see <http://www.gnu.org/licenses/>.
## END OF HEADER ###############################################################

## Script initialisation -------------------------------------------------------
rm(list=ls())

## Retireve input args ---------------------------------------------------------
args <- commandArgs(trailingOnly=TRUE)
pfg.name <- as.character(args[1]) ## give the name of the PFG
# pfg.name <- "P6"

## Constants definition --------------------------------------------------------
user = "luke" ## the id of user/machine which do the analyses
sce = "AUST" ## the vlimatic environmental variable source ("NICK" or "AUST")
proj.name = "ParcEcrins_current"
mod.pattern = "EMwmeanByTSS" ## the selected model ("EMcaByTSS", or "EMwmeanByTSS")

## Paths to data definition ----------------------------------------------------
if(user == "maya"){
  path_input <- "~/Documents/_BIOMOVE/EX_BIOMOD2/_NEW_VERSION/_INPUT_DATA/"
  path_output <- "~/Documents/_BIOMOVE/EX_BIOMOD2/_NEW_VERSION/_OUTPUT_DATA/"
} else if (user == "damien"){
  path_input <- "~/Work/FATEHD/data/scripts_fatehd_january_2016_by_Maya/SDMs_pourDamien/_INPUT_DATA/"
  path_output <- "~/Work/FATEHD/data/scripts_fatehd_january_2016_by_Maya/SDMs_pourDamien/_OUTPUT_DATA_BIS/"
} else if (user == "ftp"){
  path_input <- "/media/ftp-public/GUEGUEN_Maya/_SP_VERSION/_INPUT_DATA/"
  path_output <- "/media/ftp-public/GUEGUEN_Maya/_SP_VERSION/_OUTPUT_DATA/"
} else if (user == "luke"){
  path_input <- "/nfs_scratch2/emabio/FATEHD/_SP_VERSION/_INPUT_DATA/"
  path_output <- "/nfs_scratch2/emabio/FATEHD/_SP_VERSION/_OUTPUT_DATA_NEW_ENV/"
  .libPaths('/nfs_scratch2/emabio/R_PKG_LUKE') ## here are the shared library installed on luke
} else stop("Unsupported 'user' value")

## Load needed packages --------------------------------------------------------
library(raster)
library(rasterVis)

## Load the list of species belonging to PFGs ----------------------------------
determ <- get(load(file.path(path_input, "determinantes")))

## define the specific output directory ----------------------------------------
path_output_sdm <- file.path(path_output, paste0("DATA_", sce), pfg.name, proj.name, mod.pattern)
dir.create(path_output_sdm, showWarnings = FALSE, recursive = TRUE)

## Build PFG species SDM stacks ------------------------------------------------
PFG.sp.sdm <- lapply(determ[pfg.name], 
  function(sps_){
    cat("\n> creating SDM stack of", sps_)
    sps_files <- list.files(paste0(path_output, "DATA_",sce, "/", sps_, "/proj_", 
                                   proj.name, "/individual_projections"), 
                            pattern = paste0(mod.pattern, ".*img$"), full.names = TRUE)
    return(raster::stack(sps_files))
    })

##' @todo Implement the parallel version of this function or make the full scrip
##'   runable PFG by PFG because it is far too long to compute this way
PFG.sdm <- lapply(names(PFG.sp.sdm), 
  function(pfg.name_){
    cat("\nDeal with", pfg.name_)
    sdm.stk_ <- PFG.sp.sdm[[pfg.name_]]
    cat("\n> calc quantile 70 of sdms...")
    sdm.q70 <- calc(sdm.stk_, fun = function(x){ quantile(x, probs = 0.7, na.rm=TRUE) }, filename = file.path(path_output_sdm, paste0(pfg.name_, "_q70_SDM.img")), overwrite = TRUE)
    cat("\n> calc quantile 80 of sdms...")
    sdm.q80 <- calc(sdm.stk_, fun = function(x){ quantile(x, probs = 0.8, na.rm=TRUE) }, filename = file.path(path_output_sdm, paste0(pfg.name_, "_q80_SDM.img")), overwrite = TRUE)
    cat("\n> calc quantile 90 of sdms...")
    sdm.q90 <- calc(sdm.stk_, fun = function(x){ quantile(x, probs = 0.9, na.rm=TRUE) }, filename = file.path(path_output_sdm, paste0(pfg.name_, "_q90_SDM.img")), overwrite = TRUE)
    cat("\n> calc cv of sdms...")
    sdm.cv <- calc(sdm.stk_, fun = raster::cv, filename = file.path(path_output_sdm, paste0(pfg.name_, "_cv_SDM.img")), overwrite = TRUE)
    return(raster::stack(sdm.q70, sdm.q80, sdm.q80, sdm.cv))
  })
names(PFG.sdm) <- names(PFG.sp.sdm)

# ## Produce gaphs of the PFG SDM maps -------------------------------------------
# ras.theme <- rasterTheme(region = brewer.pal(9, 'Greens'))
# for(pfg.name_ in names(PFG.sdm)){
#   cat("\nProducing SDM for", pfg.name_)
#   png(file.path(path_output_sdm, paste0(pfg.name_, "_SDM.png")))
#   print(rasterVis::levelplot(PFG.sdm[[pfg.name_]], par.settings = ras.theme))
#   dev.off()
# }

# ## Save the chosen HS in a "old fashion" suitable directory to be make it compatible 
# ## with all scripts developped previously
# 
# SDMhab.dir = paste0(path_output,"DATA_",sce,"/HS_Current/cartes") # directory of PFG SDMs
# dir.create(SDMhab.dir, showWarnings = FALSE, recursive = TRUE)
# 
# for(pfg.name_ in names(PFG.sdm)){
#   writeRaster(subset(PFG.sdm[[pfg.name_]], 2), file.path(SDMhab.dir, paste0("HS_f0_", pfg.name_, ".asc")), overwrite = TRUE)
# }

q("no")