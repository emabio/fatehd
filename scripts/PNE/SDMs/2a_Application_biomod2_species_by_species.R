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
##' @note Here the philosophy is quite different from Isa's Paper where the PFGDMs
##'   were constructed using merged occurences of all the species that bellong it.
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
library(biomod2)


### CAREFUL !! IF USING GBM, parallelization sometimes goes wrong...
args <- commandArgs(trailingOnly=TRUE)
i = as.numeric(args[1])

## Constants definition --------------------------------------------------------
sce = "AUST"
user = "ftp"
nrep = 5 # number of repetitions in Biomod


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
}


## create the output directory
path_sce <- paste(path_output,"DATA_",sce,"/",sep="")
dir.create(path_sce, showWarnings = FALSE, recursive = TRUE)
setwd(path_sce)


## get the data we need --------------------------------------------------------
sp = gtools::mixedsort(sub("ENV_","",list.files(paste(path_output,"DATA_",sce,"/PFT_env/",sep=""))))
nbgr = length(sp)  # number of PFT
load(paste(path_output,"PLOTS_Ecrins",sep=""))
coord_XY <- input[,c("PlotID_KEY","X_ETRS89","Y_ETRS89")]

if(sce=="NICK"){
  load(paste(path_input,"ENV_ALL_ECRINS/CURRENT/ParcEcrinsSG_env.table",sep="")) # environment for projections
  data.env <- get(load(paste(path_output,"data.env.NICK",sep="")))
} else {
  load(paste(path_input,"ENV_ALL_ECRINS/EOBS_1970_2005/ParcEcrinsSG_env.table",sep="")) # environment for projections
  data.env <- get(load(paste(path_output,"data.env.AUST",sep="")))
}

mask <- get(load(paste(path_output,"maskParcEcrinsSansGlacierNiLac.raster",sep="")))
xy_proj <- xyFromCell(mask,which(!is.na(mask[])))


###############################################################################################
### Implementation of the models in BIOMOD and projections onto Ecrins Park
###############################################################################################

print(sp[i])
load(paste(path_sce,"PFT_occ/OCC_",sp[i],sep="")) # responses vector
load(paste(path_sce,"PFT_env/ENV_",sp[i],sep="")) # explanatory variables

## Add Glacier absences
site_gla <- coord_XY$PlotID_KEY[grep("GLACIER-",coord_XY$PlotID_KEY)]
abs_add <- rep(0,200)
names(abs_add) <- sample(site_gla,200)
pft.occ <- c(pft.occ,abs_add)
pft.env <- rbind(pft.env,data.env[which(rownames(data.env) %in% names(abs_add)),])

xy <- coord_XY[which(names(pft.occ) %in% coord_XY$PlotID_KEY),c("X_ETRS89","Y_ETRS89")]

## formating data in a biomod2 friendly way ------------------------------------
bm.form <- BIOMOD_FormatingData(resp.var=pft.occ, expl.var=pft.env[,c(1:7)], resp.xy=xy, 
                                  resp.name=sp[i])

## define models options -------------------------------------------------------
bm.opt <- BIOMOD_ModelingOptions(GLM=list(type="polynomial",test="AIC"), GBM=list(n.trees=3000), GAM=list(k=4))

## run single models -----------------------------------------------------------
bm.mod <- BIOMOD_Modeling(data=bm.form, models=c('RF', 'MARS', 'GLM', 'GAM', 'GBM'), models.options=bm.opt,
                        NbRunEval=nrep, DataSplit=70, 
                        Prevalence = 0.5,
                        VarImport=0, models.eval.meth=c('TSS','ROC'), do.full.models = TRUE, modeling.id = 'mod1')

## run ensemble models ---------------------------------------------------------
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

## project ensemble models -----------------------------------------------------
bm.ef <- BIOMOD_EnsembleForecasting(EM.output = bm.em, 
                                    new.env = as.data.frame(ParcEcrinsSG_env[,1:7]), 
                                    proj.name = "ParcEcrins_current",
                                    xy.new.env = xy_proj, 
                                    selected.models="all", 
                                    binary.meth=c('TSS'))

## end of script ---------------------------------------------------------------


