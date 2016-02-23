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
##'   - 28/01/2016 (Damien)
##'       - remove MARS models for 23 species because it randomly failed
##'         [1] "X457"  "X479"  "X561"  "X786"  "X796"  "X846"  "X902"  "X976"  "X1080" "X1085" "X1087" "X1455"
##'         [13] "X1544" "X4826" "X4837" "X4969" "X4999" "X5025" "X5029" "X5059" "X5185" "X5327" "X5343"
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
# sp.name <- "P8"


## Constants definition --------------------------------------------------------
user = "luke" ## the id of user/machine which do the analyses
sce = "NICK_VAR_ISA_ALPS" ## the vlimatic environmental variable source ("NICK" or "AUST")
version.name = "_PFG_VERSION" ## the type of model unit ("_SP_VERSION" or "_PFG_VERSION")
nrep = 2 ## number of repetitions in Biomod
# env.var.names <- c("bio_6", "bio_9", "bio_12", "bio_15", "carbon", "slope") ## new variables we are interested in
env.var.names <- c("bio_3", "bio_4", "bio_7", "bio_11", "bio_12", "carbon", "slope") ## isa's variables
nb.extr.abs.glacier <- 200 ## if > 0, the number of absences that will be added from glacier area
check.computed <- FALSE ## check if models/ensemble models/ensemble models projections have been computed and load them directly if it is the case

## Print brief summary of the script args --------------------------------------
cat("\nBuild PFG determinant species SDM -------------------------------------")
cat("\nstart at:", format(Sys.time(), "%a %d %b %Y %X"))
cat("\n species:", sp.name)
cat("\n user:", user)
cat("\n version.name:", version.name)
cat("\n nrep:", nrep)
cat("\n env.var.names:", env.var.names)
cat("\n nb.extr.abs.glacier:", nb.extr.abs.glacier)
cat("check.computed:", check.computed)
cat("\n-----------------------------------------------------------------------")

## Paths to data definition ----------------------------------------------------
if (user == "ftp"){
  path_input <- paste0("/media/ftp-public/GUEGUEN_Maya/", version.name , "/_INPUT_DATA/")
  path_output <- paste0("/media/ftp-public/GUEGUEN_Maya/", version.name , "/_OUTPUT_DATA/")
} else if (user == "luke"){
  path_input <- "/nfs_scratch2/emabio/FATEHD/"
  path_output <- paste0("/nfs_scratch2/emabio/FATEHD/", version.name , "/_OUTPUT_DATA_NEW_ENV/")
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

## load calibration environment
(load("/nfs_scratch2/emabio/FATEHD/_ISA_VERSION/bioclimatic_data"))

## to get slope
slope <- raster("/nfs_scratch2/emabio/FATEHD/ENV_ALL_ALPS/EOBS_1970_2005/slope.img")
input$slope <- extract(slope, coord_XY[input$PlotID_KEY , c("X_ETRS89","Y_ETRS89")])
rownames(input) <- input$PlotID_KEY

## load projection environment
(load("/nfs_scratch2/emabio/FATEHD/_ISA_VERSION/ENV_ALL_ECRINS/CURRENT/ParcEcrinsSG_env.table"))


###############################################################################################
### Implementation of the models in BIOMOD and projections onto Ecrins Park
###############################################################################################

# load(paste(path_sce,"PFT_occ/OCC_",sp.name,sep="")) # responses vector
# load(paste("/nfs_scratch2/emabio/FATEHD/_PFG_VERSION/_OUTPUT_DATA_NEW_ENV/DATA_AUST/PFT_occ/OCC_",sp.name,sep="")) # responses vector
load(paste("/nfs_scratch2/emabio/FATEHD/PFT_OCC_FR_Alps/OCC_",sp.name,sep=""))

# ## Add Glacier absences
# if(nb.extr.abs.glacier > 0){
#   site_gla <- coord_XY$PlotID_KEY[grep("GLACIER-",coord_XY$PlotID_KEY)]
#   abs_add <- rep(0, nb.extr.abs.glacier)
#   names(abs_add) <- sample(site_gla, nb.extr.abs.glacier)
#   pft.occ <- c(pft.occ,abs_add)
# }


xy <- input[names(pft.occ), c("X_ALBERS",  "Y_ALBERS")]

# # Weight vector initialization to balance presences/absences
# poids <- c()
# k = length(pft.occ[pft.occ==1])/length(pft.occ[pft.occ==0])
# poids[pft.occ==1] <- 1
# poids[pft.occ==0] <- k

## formating data in a biomod2 friendly way ------------------------------------
bm.form <- BIOMOD_FormatingData(resp.var = pft.occ, 
                                expl.var = input[names(pft.occ), sub("_", "", env.var.names)], 
                                resp.xy = xy, 
                                resp.name = sp.name)

## define models options -------------------------------------------------------
bm.opt <- BIOMOD_ModelingOptions(GLM=list(type="polynomial",test="AIC"), GBM=list(n.trees=3000), GAM=list(k=4))


## run single models -----------------------------------------------------------

if(check.computed){
  bm.mod.file <- list.files(paste0(path_output,"DATA_",sce,"/", sp.name), pattern = "mod1.models.out$", full.names = TRUE)
  bm.em.file <- list.files(paste0(path_output,"DATA_",sce,"/", sp.name), pattern = "ensemble.models.out$", full.names = TRUE)
  bm.ef.file <- list.files(paste0(path_output,"DATA_",sce,"/", sp.name, "/proj_ParcEcrins_current"), pattern = "ensemble.projection.out$", full.names = TRUE)
} else {
  bm.mod.file <- bm.em.file <- bm.ef.file <- NULL
}

if(check.computed & length(bm.mod.file)){
  cat("\n loading previous version of bm.mod..")
  bm.mod <- get(load(bm.mod.file))
} else {
  bm.mod <- BIOMOD_Modeling(data = bm.form, 
                            models = c('RF', 'GLM', 'GBM'), #c('RF', 'MARS', 'GLM', 'GAM', 'GBM'), 
                            models.options = bm.opt,
                            NbRunEval = nrep, DataSplit = 70,
#                             Yweights = poids,
                            Prevalence = 0.5, 
                            VarImport=0, 
                            models.eval.meth = c('TSS','ROC'), 
                            do.full.models = FALSE, modeling.id = 'mod1')
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
                                   eval.metric.quality.threshold = 0.35, 
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
                                      new.env = ParcEcrinsSG_env,
                                      proj.name = "ParcEcrins_current",
                                      selected.models="all", 
                                      binary.meth=c('TSS'))
}

cat("\n\nCompleted!")
cat("\nended at:", format(Sys.time(), "%a %d %b %Y %X"))
q("no")

## end of script ---------------------------------------------------------------

# ## generate params for grid computing ------------------------------------------
# ## for "_SP_VERSION" 
# sp.all <- sub("OCC_", "", gtools::mixedsort(list.files("/nfs_scratch2/emabio/FATEHD/_SP_VERSION/_OUTPUT_DATA_NEW_ENV/DATA_AUST/PFT_occ/", "^OCC_X[0-9]{1,5}$")))
# cat(sp.all, sep = "\n", file = "/nfs_scratch2/emabio/FATEHD/fatehdhub/scripts/PNE/SDMs/2a_Application_biomod2_species_by_species.params")
# ## for "_PFG_VERSION"
# pfg.all <- sub("OCC_", "", gtools::mixedsort(list.files("/nfs_scratch2/emabio/FATEHD/_PFG_VERSION/_OUTPUT_DATA_NEW_ENV/DATA_AUST/PFT_occ", "^OCC_[C,H,P][0-9]{1,2}$")))
# cat(pfg.all, sep = "\n", file = "/nfs_scratch2/emabio/FATEHD/fatehdhub/scripts/PNE/SDMs/2a_Application_biomod2_species_by_species.params2")
# ## end generate params for grid computing --------------------------------------

# ## check campaain status -------------------------------------------------------
# ## define path to param and output files
# param.file <- "/nfs_scratch2/emabio/FATEHD/fatehdhub/scripts/PNE/SDMs/2a_Application_biomod2_species_by_species.params"
# output.file <- "/nfs_scratch2/emabio/FATEHD/_SP_VERSION/_OUTPUT_DATA_NEW_ENV/DATA_AUST"
# # param.file <- "/nfs_scratch2/emabio/FATEHD/fatehdhub/scripts/PNE/SDMs/2a_Application_biomod2_species_by_species.params2"
# # output.file <- "/nfs_scratch2/emabio/FATEHD/_PFG_VERSION/_OUTPUT_DATA_NEW_ENV/DATA_AUST"
# 
# ## get the names of species we want to model
# sp.all <- as.character(read.table(param.file, stringsAsFactors = FALSE)[, 1])
# 
# ## check if models/ensemble models/ensemble models projections have been produced
# sp.mod.test <- sapply(sp.all, function(sp_){file.exists(file.path(output.file, sp_, paste0(sp_, ".mod1.models.out")))})
# sp.em.test <- sapply(sp.all, function(sp_){file.exists(file.path(output.file, sp_, paste0(sp_, ".mod1ensemble.models.out")))})
# sp.ef.test <- sapply(sp.all, function(sp_){file.exists(file.path(output.file, sp_, "proj_ParcEcrins_current", paste0(sp_, ".ParcEcrins_current.ensemble.projection.out")))})
# 
# ## print the campain status summary
# cat("\ncampain stats:",
#     "\n\tsingle models:", sum(sp.mod.test), "/", length(sp.mod.test),
#     "\n\tensemble models:", sum(sp.em.test), "/", length(sp.em.test),
#     "\n\tensemble models projection:", sum(sp.ef.test), "/", length(sp.ef.test))
# ## end check campaain status ---------------------------------------------------
