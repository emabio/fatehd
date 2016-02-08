### HEADER #####################################################################
##' @title FATE-HD Simulations post-treatment : PFG maps
##'
##' @author Damien G. & Maya G., Isabelle B., Marta C., Ceres B.
##' @contact damien.georges2 at gmail.com, maya.gueguen at gmail.com
##' 
##' @date 21/01/2016
##' 
##' @description This script is designed to plot for each PFG both maps from SDM
##'              and FATE-HD simulations.
##'              
##' @note Parameters required :
##'              - Simulation name
##'              - ParamSimul file name
##'              
##'              - SIMUL FOLDER is hard written at the beginning of the script...
##'                maybe should be given with the simulation name...
##'   
## END OF HEADER ###############################################################


## Script initialisation -------------------------------------------------------
rm(list=ls())
gc(reset=TRUE)

library(raster)
library(sp)
library(rgdal)
library(foreign)
library(Hmisc)

#################################################################################################
### PARAMETERS OF THE SIMULATION (FATE)
#################################################################################################

## Retrieve input args ---------------------------------------------------------
args <- commandArgs(trailingOnly=TRUE)
simul_name = as.character(args[1]) ## give the simulation name
file_name = as.character(args[2]) ## give the paramSimul file name

## Set directories -------------------------------------------------------------
simul_folder = paste("/nfs_scratch/mayagueguen/FATE_newHS/", simul_name,"/",sep="")
setwd(sub(simul_name,"",simul_folder))
output_folder = "resultsVALID_graphs/" ## directory to plot graphs
dir.create(file.path(simul_folder,output_folder))

## Retrieve name of the simul replication --------------------------------------
cat("\n PARAM SIMUL FILE : ",file_name)
i = sub(".txt","",sub(paste(simul_name,"/PARAM_SIMUL/paramSimul_",sep=""),"",file_name))
cat("\n Rep : ", i)

## Get parameters from GLOBAL PARAMETERS FILE ----------------------------------
cat("\n GLOBAL PARAM FILE : Global_parameters_",i,".txt\n")
if(file.exists(paste(simul_name,"/DATA/Global_parameters_",i,".txt",sep=""))){
  param_file <- file(paste(simul_name,"/DATA/Global_parameters_",i,".txt",sep=""),"r")
  param_lines <- readLines(param_file)
  do_soil = as.numeric(strsplit(param_lines[which(sapply(param_lines,function(x) strsplit(x,split=" ")[[1]][1])=="DO_SOIL_COMPETITION")],split=" ")[[1]][2])
  nb_PFGS = as.numeric(strsplit(param_lines[which(sapply(param_lines,function(x) strsplit(x,split=" ")[[1]][1])=="NB_FG")],split=" ")[[1]][2])
  if(do_soil){
    nb_soil = as.numeric(strsplit(param_lines[which(sapply(param_lines,function(x) strsplit(x,split=" ")[[1]][1])=="NB_SOIL_CATEGORIES")],split=" ")[[1]][2])
    soil_thres = as.numeric(strsplit(param_lines[which(sapply(param_lines,function(x) strsplit(x,split=" ")[[1]][1])=="SOIL_CATEGORIES_TRESHOLDS")],split=" ")[[1]][2:(2+nb_soil)])
  }
} else { cat("\n GLOBAL PARAMETERS FILE of replication ",i," does not exist !!\n") }

## Get parameters from PARAM SIMUL FILE ----------------------------------
cat("\n PARAM SIMUL FILE : ",file_name,"\n")
if(file.exists(file_name)){
  paramSim_file <- file(file_name,"r")
  paramSim_lines <- readLines(paramSim_file)
  
  ## Get simulation mask name
  mask_name = paramSim_lines[which(paramSim_lines=="--MASK--")+1]
  cat("\n MASK NAME : ",mask_name,"\n")
  mask <- raster(mask_name)
  mask[which(mask[] == 0) ] = NA
 
  ## Get PFG succession files
  PFG = paramSim_lines[(which(paramSim_lines=="--PFG_LIFE_HISTORY_PARAMS--")+1):(which(paramSim_lines=="--PFG_LIFE_HISTORY_PARAMS--")+nb_PFGS)]
  
  ## Get for each PFG the soil contribution
  if(do_soil){
    soil_contrib = vector()
    for(fg in PFG){
      pfg_file = file(fg,"r")
      pfg_lines = readLines(pfg_file)
      soil_contrib = c(soil_contrib,as.numeric(strsplit(pfg_lines[which(sapply(pfg_lines,function(x) strsplit(x,split=" ")[[1]][1])=="SOIL_CONTRIB")],split=" ")[[1]][2]))
    }
  }
  
  ## Get PFG names
  PFG_files = paramSim_lines[(which(paramSim_lines=="--PFG_ENVSUIT--")+1):(which(paramSim_lines=="--PFG_ENVSUIT--")+nb_PFGS)]
  pattern = paste(simul_name,"/DATA/PFGS/SUCC/SUCC_",sep="")
  PFG = sub(".txt","",sub(pattern,"",PFG))
} else { cat("\n PARAM SIMUL FILE ",file_name," does not exist !!\n") }

#################################################################################################
### CREATE PFG MAPS
#################################################################################################

## Load FATE-HD results --------------------------------------------------------
arrayPFG <- get(load(paste(simul_name,"/outputsTables/arrayPFG_year840_rep",i,sep="")))
cat("\n",colnames(arrayPFG),"\n")
cat("\n",colSums(arrayPFG),"\n")

## PLOT MAPS SDM vs FATE-HD RESULTS --------------------------------------------
pdf(paste(simul_name,"/",output_folder,"MapsPFGdistri_SDM_vs_",simul_name,"_rep",i,".pdf",sep=""), width=12, height=12)
for(pfg in 1:length(PFG)){
  cat("\n PFG : ",PFG[pfg],"\n")
# pdf(paste(simul_name,"/",output_folder,"MapsPFGdistri_SDM_vs_",simul_name,"_rep",i,"_",PFG[pfg],".pdf",sep=""), width=12, height=12) 
  ## Get PFG Habitat Suitability map
  habitat = raster(PFG_files[pfg])
  cat("\n FILE loaded :",PFG_files[pfg],"\n")

  ## Create FATE-HD predicted abundances maps
  simul = mask
  simul[] = arrayPFG[,PFG[pfg]]


  op = par(mfrow=c(1,2))
  plot(habitat, main="Habitat Suitability")
  plot(simul, main=paste("FATE-HD simulation - PFG ",PFG[pfg]),
       legend.width=2.5,
       legend.mar=12,
       col=rev(terrain.colors(length(seq(1, max(c(100,simul[]), na.rm=TRUE), 1))-1)),
       breaks=seq(1, max(c(100,simul[]), na.rm=TRUE), 1),
       axis.args=list(at= seq(1, max(c(100,simul[]), na.rm=TRUE), 200)))
  if(do_soil){legend("topleft", legend=paste("Soil contribution = ",soil_contrib), bty="n")
    par(op)
    
  }
}
dev.off()

