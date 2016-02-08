### HEADER #####################################################################
##' @title FATE-HD Simulations post-treatment : outputs tables
##'
##' @author Damien G. & Maya G.
##' @contact damien.georges2 at gmail.com, maya.gueguen at gmail.com
##' 
##' @date 21/01/2016
##' 
##' @description This script is designed to construct from FATE-HD fullArrays outputs
##'              yearly results of PFG abundances in each pixel summarised by PFG or
##'              by stratum. This new tables have the advantage to be easier to use.
##'              In the second part of the script we will show an example to
##'              produce the evolution of % of coverage of each PFG.
##'              
##' @note Parameters required :
##'              - Simulation name
##'              - ParamSimul file name
##'              - Number of resources to be used
##'              - Do you want PDF maps to be printed ? (much longer)
##'              
##'              - SIMUL FOLDER is hard written at the beginning of the script...
##'                maybe should be given with the simulation name...
##'   
## END OF HEADER ###############################################################


## Script initialisation -------------------------------------------------------
rm(list=ls())
library(parallel)
library(raster)

#################################################################################################
### PARAMETERS OF THE SIMULATION (FATE)
#################################################################################################

## Retrieve input args ---------------------------------------------------------
args <- commandArgs(trailingOnly=TRUE)
simul_name = as.character(args[1]) ## give the simulation name
file_name = as.character(args[2]) ## give the paramSimul file name
nbCores = as.character(args[3]) ## give the number of resources to create outputs tables
maps = as.character(args[4]) ## do you want to print pdf maps ?

## Set directories -------------------------------------------------------------
simul_folder = paste("/nfs_scratch/mayagueguen/FATE_newHS/", simul_name,"/",sep="")
setwd(sub(simul_name,"",simul_folder))
output_folder = "outputsTables/"
dir.create(file.path(simul_folder,output_folder))

#################################################################################################
### FUNCTION OF YEARS
#################################################################################################

createTables <- function(y){
  print(paste("Create tables year",y))
  ## unzip and read the arrays saved
  fileTmp = paste(simul_folder, "RESULTS/", i, "/FullArray_", y, ".txt", sep="")
  system(paste("gunzip ",paste(fileTmp, ".gz", sep=""), sep=""))
  mytable = read.table(fileTmp, header=TRUE, sep= " ", row.names=1)
  
  ## table strata : sites x strata
  arrayStrata = data.frame(matrix(NA, nrow=nrow(mytable), ncol=nb_strat))
  colnames(arrayStrata) = 1:nb_strat
  for(st in 1:nb_strat){
    selectCol = mytable[,grep(paste("str", st, sep=""), colnames(mytable))]
    arrayStrata[,st] = apply(selectCol, 1, sum)
  }
  
  ## table PFGs : sites x pfgs
  arrayPFG = data.frame(matrix(NA, nrow=nrow(mytable), ncol=length(PFG)))
  colnames(arrayPFG) = PFG
  for (gp in PFG){
    selectCol = mytable[,grep(gp, colnames(mytable))]
    arrayPFG[,gp] = apply(selectCol, 1, sum) # SUM or MAX, that is the question
  }
  
  ## save and rezip the arrays
  save(arrayStrata, file=paste(simul_folder, output_folder, "arrayStrata_year",y, "_rep",i, sep=""))
  save(arrayPFG, file=paste(simul_folder, output_folder, "arrayPFG_year",y, "_rep",i, sep=""))
  system(paste("gzip -9 ",fileTmp, sep=""))
  
  ## table soil : sites x years
  if(do_soil){
    fileTmpSoil = paste(simul_folder, "RESULTS/", i, "/FullSoilArray_", y, ".txt", sep="")
    system(paste("gunzip ",paste(fileTmpSoil, ".gz", sep=""), sep=""))
    mytableSoil = read.table(fileTmpSoil, header=TRUE, sep= " ", row.names=1)
    system(paste("gzip -9 ",fileTmpSoil, sep=""))
  }
  
  ### IF ASKED : PDF MAPS ###############################################################
  if (maps){
    pdf(paste(simul_folder, output_folder, "strata_year",y, "_rep",i, ".pdf", sep=""))
    for (st in 1:nb_strat){
      mask[] = arrayStrata[,st]
      mask[which(mask[]>10000)] = 10000
      plot(mask, main=paste('strata', st), breaks=c(1,3000,7000,10000), col=c("lightgreen", "orange", "darkgreen"))
    }
    dev.off()
    
    pdf(paste(simul_folder, output_folder, "PFG_year",y, "_rep",i, ".pdf", sep=""))
    for (gp in 1:length(PFG)){
      mask[which(!is.na(mask[]))] = arrayPFG[,gp]
      mask[which(mask[]>10000)] = 10000
      plot(mask, main=paste('PFG', colnames(arrayPFG)[gp]), breaks=c(0,1,3000,7000,10000), col=c("grey","lightgreen", "orange", "darkgreen"))
    }
    dev.off()
  }
  #######################################################################################
  
  if(do_soil) return(mytableSoil)
  else return(y)
}


#################################################################################################
### SIMULATION PARAMETERS
#################################################################################################

## Retrieve name of the simul replication --------------------------------------
cat("\n PARAM SIMUL FILE : ",file_name)
i = sub(".txt","",sub(paste(simul_name,"/PARAM_SIMUL/paramSimul_",sep=""),"",file_name))
cat("\n Rep : ", i)

## Get parameters from GLOBAL PARAMETERS FILE ----------------------------------
cat("\n GLOBAL PARAM FILE : Global_parameters_",i,".txt\n")
if(file.exists(paste(simul_name,"/DATA/Global_parameters_",i,".txt",sep=""))){
  param_file <- file(paste(simul_name,"/DATA/Global_parameters_",i,".txt",sep=""),"r")
  param_lines <- readLines(param_file)
  nb_strat = as.numeric(strsplit(param_lines[which(sapply(param_lines,function(x) strsplit(x,split=" ")[[1]][1])=="NB_STRATUM")],split=" ")[[1]][2])
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
  if(maps){
    mask <- raster(mask_name)
    mask[which(mask[] == 0) ] = NA
  }
  
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
  pattern = paste(simul_name,"/DATA/PFGS/SUCC/SUCC_",sep="")
  PFG = sub(".txt","",sub(paste(simul_name,"/DATA/PFGS/SUCC/SUCC_",sep=""),"",PFG))
} else { cat("\n PARAM SIMUL FILE ",file_name," does not exist !!\n") }


#################################################################################################
### EXTRACT FATE-HD OUTPUTS
#################################################################################################

## Get list of arrays and extract years of simulation --------------------------
fullArrays <- grep("FullArray_", list.files(paste(simul_name,"/RESULTS/",i,sep="")), value=TRUE)
years <- sort(as.numeric(sub("FullArray_", "", sub(".txt","", sub(".gz", "", fullArrays)))))

## Zip files that were unzipped ------------------------------------------------
if(length(fullArrays[-grep('.txt.gz', fullArrays)])>0){
  warnings("mixed zipped and unzipped files")
  toCorrect = fullArrays[-grep('.gz', fullArrays)]
  for(cc in toCorrect){ system(paste("gzip -9 ",paste(simul_folder,"RESULTS/",i,"/",cc,sep=""),sep="")) }
}

## Create outputs tables -------------------------------------------------------
cat("\n BEGINNING OF OUTPUTS CREATION : ",date(),"\n")
res <- mclapply(years, createTables, mc.cores=nbCores)
cat("\n END OF OUTPUTS CREATION : ",date(),"\n")

## If soil, create soil array --------------------------------------------------
if(do_soil){
  cat("\n BEGINNING OF SOIL OUTPUT CREATION : ",date(),"\n")
  arraySoil = matrix(unlist(res),nrow=length(years),byrow=TRUE) ## one row for each year, with pixels in columns
  rownames(arraySoil) = years
  save(arraySoil, file=paste(simul_folder, output_folder, "arraySoil_rep",i, sep=""))
  cat("\n END OF SOIL OUTPUT CREATION : ",date(),"\n")
}

#################################################################################################
### CREATE COVERAGE PLOT
#################################################################################################

distri = array(NA, dim=c(length(years),length(PFG)), dimnames=list(years,PFG) )
distriAbund = array(NA, dim=c(length(years),length(PFG)), dimnames=list(years,PFG) )
for (y in years){
  print(y)
  ## load the post treated FATE-HD output for considered year
  load(paste(simul_folder,"outputsTables/arrayPFG_year",y,"_rep",i, sep=""))
  ## calculate the % of cover of each PPFG
  distri[as.character(y),] = apply(arrayPFG, 2, function(x){ sum (x>0, na.rm=T) }) / nrow(arrayPFG)
  distriAbund[as.character(y),] = colSums(arrayPFG)#/sum(colSums(arrayPFG))
}

## produce the plot
pdf(paste(simul_name,"/PFG_coverage_evolution_", i, ".pdf",sep=""))
par(mfrow = c(ceiling(sqrt(nb_PFGS)),ceiling(sqrt(nb_PFGS))),mar=c(2,4,2,1))
for (pfg in dimnames(distri)[[2]]){
  plot(years, distri[,pfg][which(!is.na(distri[,pfg]))], 
       type="l", 
       main=strsplit(pfg,"_")[[1]][1] , 
       xlab = "Time (year)",
       ylab = "% occupation",
       ylim = c(0,1))
}
dev.off()

