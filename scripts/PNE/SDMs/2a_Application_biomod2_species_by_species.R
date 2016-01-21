### HEADER #####################################################################
##' @title Build PFG determinant species SDM
##'
##' @author Damien G. & Maya G.
##' @contact damien.georges2 at gmail.com
##' 
##' @date 16/01/2016
##' 
##' @description Here we will build SDMs of our ~ 270 determinant species. We 
##'   are interested in the final ensemble model as weighted mean and comitee 
##'   averaging. PFGDMs will then be constructed taking the max probs accros
##'   species of the PFG. 
##'   
##' @note 
##'   - Here the philosophy is quite different from Isa's Paper where the PFGDMs
##'     were constructed using merged occurences of all the species that bellong it.
##'   - A peace of code to generate parameters file to run script on luke is avalable
##'     at the end of the script.
##' 
##' @log 
##'   - 17/01/2016 (Damien)
##'       - do projection with rasters instead of tables
##'       - make it gridded friendly way
##'       - add the possiblity to check if models/ensemble models/ensemble models 
##'         projections have been produced not to start from scratch
##'       - give the species name as input parameter (not the id of the species)
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
sp.name <- as.character(args[1]) ## give the name of the species
# sp.name <- "X10352"


## Constants definition --------------------------------------------------------
user = "luke" ## the id of user/machine which do the analyses
sce = "AUST" ## the vlimatic environmental variable source ("NICK" or "AUST")
nrep = 5 ## number of repetitions in Biomod
env.var.names <- c("bio_3", "bio_4", "bio_7", "bio_11", "bio_12", "carbon", "slope") ## variables we are interested in
nb.extr.abs.glacier <- 200 ## if > 0, the number of absences that will be added from glacier area
check.computed <- TRUE ## check if models/ensemble models/ensemble models projections have been computed and load them directly if it is the case

## Print brief summary of the script args --------------------------------------
cat("\nBuild PFG determinant species SDM -------------------------------------")
cat("\nstart at:", format(Sys.time(), "%a %d %b %Y %X"))
cat("\n species:", sp.name)
cat("\n user:", user)
cat("\n nrep:", nrep)
cat("\n env.var.names:", env.var.names)
cat("\n nb.extr.abs.glacier:", nb.extr.abs.glacier)
cat("check.computed:", check.computed)
cat("\n-----------------------------------------------------------------------")

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
  path_output <- "/nfs_scratch2/emabio/FATEHD/_SP_VERSION/_OUTPUT_DATA/"
  .libPaths('/nfs_scratch2/emabio/R_PKG_LUKE') ## here are the shared library installed on luke
} else stop("Unsupported 'user' value")

## load needed packages --------------------------------------------------------
library(biomod2)

## create the output directory
path_sce <- paste(path_output,"DATA_",sce,"/",sep="")
dir.create(path_sce, showWarnings = FALSE, recursive = TRUE)
setwd(path_sce)


## get the data we need --------------------------------------------------------
load(paste(path_output,"PLOTS_Ecrins",sep=""))
coord_XY <- input[,c("PlotID_KEY","X_ETRS89","Y_ETRS89")]

if(sce == "NICK"){
  ## releves environment for models creation
  data.env <- get(load(paste(path_output,"data.env.NICK",sep="")))
  ## full environment of PNE for projections
  PNE.env.files <- file.path(path_input,"ENV_ALL_ECRINS", "CURRENT", paste0(env.var.names, ".img"))
  PNE.env.stk <- raster::stack(PNE.env.files)
} else if(sce == "AUST"){
  ## releves environement for models creation
  data.env <- get(load(paste(path_output,"data.env.AUST",sep="")))
  ## full environment of PNE for projections
  PNE.env.files <- file.path(path_input,"ENV_ALL_ECRINS", "EOBS_1970_2005", paste0(env.var.names, ".img"))
  PNE.env.stk <- raster::stack(PNE.env.files)
} else stop("Unsupported 'sce' value")

###############################################################################################
### Implementation of the models in BIOMOD and projections onto Ecrins Park
###############################################################################################

load(paste(path_sce,"PFT_occ/OCC_",sp.name,sep="")) # responses vector
load(paste(path_sce,"PFT_env/ENV_",sp.name,sep="")) # explanatory variables

## Add Glacier absences
if(nb.extr.abs.glacier > 0){
  site_gla <- coord_XY$PlotID_KEY[grep("GLACIER-",coord_XY$PlotID_KEY)]
  abs_add <- rep(0, nb.extr.abs.glacier)
  names(abs_add) <- sample(site_gla, nb.extr.abs.glacier)
  pft.occ <- c(pft.occ,abs_add)
  pft.env <- rbind(pft.env,data.env[which(rownames(data.env) %in% names(abs_add)),])
}

xy <- coord_XY[which(names(pft.occ) %in% coord_XY$PlotID_KEY),c("X_ETRS89","Y_ETRS89")]

## formating data in a biomod2 friendly way ------------------------------------
bm.form <- BIOMOD_FormatingData(resp.var=pft.occ, expl.var=pft.env[,c(1:7)], resp.xy=xy, 
                                  resp.name=sp.name)

## define models options -------------------------------------------------------
bm.opt <- BIOMOD_ModelingOptions(GLM=list(type="polynomial",test="AIC"), GBM=list(n.trees=3000), GAM=list(k=4))


## run single models -----------------------------------------------------------

if(check.computed){
  bm.mod.file <- list.files(paste0(path_output,"DATA_",sce,"/", sp.name), pattern = "mod1.models.out$", full.names = TRUE)
  bm.em.file <- list.files(paste0(path_output,"DATA_",sce,"/", sp.name), pattern = "ensemble.models.out$", full.names = TRUE)
  bm.ef.file <- list.files(paste0(path_output,"DATA_",sce,"/", sp.name, "/ParcEcrins_current"), pattern = "ensemble.projection.out$", full.names = TRUE)
} 

if(check.computed & length(bm.mod.file)){
  cat("\n loading previous version of bm.mod..")
  bm.mod <- get(load(bm.mod.file))
} else {
  bm.mod <- BIOMOD_Modeling(data = bm.form, 
                            models = c('RF', 'MARS', 'GLM', 'GAM', 'GBM'), 
                            models.options = bm.opt,
                            NbRunEval = nrep, DataSplit = 70,
                            Prevalence = 0.5, VarImport=0, 
                            models.eval.meth = c('TSS','ROC'), 
                            do.full.models = TRUE, modeling.id = 'mod1')
}

## run ensemble models ---------------------------------------------------------
if(check.computed & length(bm.em.file)){
  cat("\n loading previous version of bm.em..")
  bm.em <- get(load(bm.em.file))
} else {
  bm.em <- BIOMOD_EnsembleModeling(modeling.output = bm.mod, 
                                   chosen.models = "all",
                                   em.by = "all",
                                   eval.metric = c('TSS'),
                                   eval.metric.quality.threshold = 0.4, 
                                   models.eval.meth = c('TSS', 'ROC'),
                                   prob.mean = FALSE,
                                   prob.mean.weight = TRUE, 
                                   prob.mean.weight.decay = 'proportional',
                                   committee.averaging = TRUE)
}

## project ensemble models -----------------------------------------------------
if(check.computed & length(bm.ef.file)){
  cat("\n loading previous version of bm.ef..")
  bm.ef <- get(load(bm.ef.file))
} else{
  bm.ef <- BIOMOD_EnsembleForecasting(EM.output = bm.em, 
                                      new.env = PNE.env.stk, 
                                      output.format = ".img",
                                      proj.name = "ParcEcrins_current",
                                      selected.models="all", 
                                      binary.meth=c('TSS'))
}

cat("\n\nCompleted!")
cat("\nended at:", format(Sys.time(), "%a %d %b %Y %X"))

## end of script ---------------------------------------------------------------

# ## generate params for grid computing ------------------------------------------
# sp.all <- gtools::mixedsort(list.files("/media/ftp-public/GUEGUEN_Maya/_SP_VERSION/_OUTPUT_DATA/DATA_AUST/", "^X"))
# cat(sp.all, sep = "\n", file = "~/Work/FATEHD/fatehdhub/scripts/PNE/SDMs/2a_Application_biomod2_species_by_species.params")
# ## end generate params for grid computing --------------------------------------
