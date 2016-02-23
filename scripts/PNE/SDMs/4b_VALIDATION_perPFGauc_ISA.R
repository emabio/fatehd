###########################################################################################
# Validation de la structure de la vegetation
# From Isabelle Boulangeat files - 03/12/12
# Rearranged by Maya Gueguen and Marta Carboni - 13/07/15
# Cleaned/finalised Ceres - 11/08/2017
###########################################################################################

gc(reset=TRUE)
rm(list=ls())

library(raster)
library(sp)
library(rgdal)
library(foreign)
library(Hmisc)

###########################################################################################

#===============================================
# 0. Define directories, simulation to validate and map
#===============================================

## Constants definition --------------------------------------------------------
user = "luke" ## the id of user/machine which do the analyses
sce = "AUST" ## the vlimatic environmental variable source ("NICK" or "AUST")

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

projETRS89 <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"

# maskSimul <- get(load(paste(path_output,"maskParcEcrinsSansGlacierNiLacCONSENSUS.raster",sep="")))
maskSimul <- get(load(file=paste(path_output,"maskParcEcrins.raster",sep="")))

# for(sce in c("NICK","AUST")){
  
  validObj.dir = paste(path_input,"Objects/",sep="") # directory of validation objects
  SDMhab.dir = paste(path_output,"DATA_",sce,"/HS_Current/cartes/",sep="") # directory of PFG SDMs
  results.dir = paste(path_output,"DATA_",sce,"/SDM_validation/",sep="") # directory to save results
  # dir.create(results.dir)

  #===============================================
  # 1. expectedMaps / Observations
  #===============================================
  
  # MAPS ---------
  
  # polygon ID per pixel
  polygonsID = raster(file.path(validObj.dir,"DELPHINE/polygonsID.img")) # polygon ID per pixel; 18915 different polygons for 484058 cells
  polygonsID <- projectRaster(from=polygonsID,res=100,crs=CRS(projETRS89))
  
  # Correspondance table for PFGs present in representativity (86 habitats, for further details see tableau PFG_reluCD.xls)
  # 86 DELPHINE CODE x 24 PFGs (presence/absence)
  # 0 = absent; 1 = potentially present; 2 = present (only 0s and 2s are used for validation - conservative approach)
  gpementsTable = read.table(file.path(validObj.dir,"groupementsETpfg.txt"),sep="\t",h=T)
  rownames(gpementsTable) = toupper(gpementsTable$code.delphine)
  head(gpementsTable)
  rownames(gpementsTable)
  
  load(paste(validObj.dir,"SAVED_analyzedPix",sep=""))
  load(paste(validObj.dir,"SAVED_PAobs_list",sep=""))
  
  
  #===============================================
  # 3. Comparisons 
  #===============================================
  
  #------------------------------------------------
  # CALCULATING PREDICTED PRESENCES AND ABSENCES HABITAT P/As TO POLYGON P/As
  # + accuracy statistics
  
  if(!file.exists(paste(results.dir,"res_STAT_SDMS",sep="")) || !file.exists(paste(results.dir,"res_STAT_POLYGONS",sep=""))){
    # Tables to store false positive (FP) and false negative (FN)
    # Polygons x PFG
    FP_Poly = array(NA, dim=c(length(unique(polygonsID[analyzedPix])),24), dimnames=list(unique(polygonsID[analyzedPix]),colnames(gpementsTable)[-1]))
    FN_Poly = array(NA, dim=c(length(unique(polygonsID[analyzedPix])),24), dimnames=list(unique(polygonsID[analyzedPix]),colnames(gpementsTable)[-1]))
    
    #------------------------------------------------
    
    # Create vectors to store accuracy stats per PFG for polygons in FATE and SDMs 
    STAT_AUC_POLY = STAT_SPECIF_POLY = STAT_SENSIB_POLY = ErRate_POLY =  nOBS_0_POLY = nOBS_2_POLY =  rep(NA, 24)
    names(STAT_AUC_POLY) = names(STAT_SENSIB_POLY) = names(STAT_SPECIF_POLY) = names(ErRate_POLY) = names(nOBS_0_POLY) = names(nOBS_2_POLY) = colnames(gpementsTable)[-1]
    
    STAT_AUC_SDM = STAT_SPECIF_SDM = STAT_SENSIB_SDM = ErRate_SDM =  rep(NA, 24)
    names(STAT_AUC_SDM) = names(STAT_SENSIB_SDM) = names(STAT_SPECIF_SDM) = names(ErRate_SDM)  = colnames(gpementsTable)[-1]
    
    
    # Loop over PFGs
    for(pfg in sort(colnames(gpementsTable)[-1])){
      print(pfg)
      
      if(file.exists(paste(SDMhab.dir,"HS_f0_",pfg,".asc",sep=""))){
      
      # pfg = strsplit(pfg, "_")[[1]][1]   # Get PFG name
      PAobs = PAobs_list[[pfg]]          # Get observations for this PFG (based on 84 delphi codes)    
      
      # MERGE for the selected pixels : df(polygonsID) + df(sites of PAobs, PAobs) by POLYGONS --> PA data
      tt1 = merge(data.frame(polygons = as.character(polygonsID[analyzedPix])), data.frame(polygons = names(PAobs), PAobs), by="polygons")
      
      # FOR SDMs ----------------------------------------------------------------------------------------------------------------------------
      habitat = raster(paste(SDMhab.dir,"HS_f0_",pfg,".asc",sep=""))
      ABpred_SDM = habitat[analyzedPix]
      
      # PRESENCE-ABSENCE data by PIXELS
      # Transform abundance data in presence-absence + remove presence when abundance < % observed absence
      PApred_SDM = ifelse(ABpred_SDM > quantile(ABpred_SDM,sum(tt1$PAobs==0)/nrow(tt1),na.rm=T),1,0)
      print(quantile(ABpred_SDM, sum(tt1$PAobs==0)/nrow(tt1),na.rm=T))
      
      # PRESENCE-ABSENCE data by POLYGONS
      PApred_SDM = unlist(lapply(split(PApred_SDM, as.character(polygonsID[analyzedPix])), function(x){ifelse(sum(x)==0,0,1)}))
      # ABUNDANCE (mean) data by POLYGONS
      ABpred_SDM = unlist(lapply(split(ABpred_SDM, as.character(polygonsID[analyzedPix])), mean))
      
      # MERGE for the selected polygons : df(polygons of PApred_SDM, PApred_SDM) + df(sites of PAobs, PAobs) by POLYGONS --> PA data
      DAT_PA_SDM = merge(data.frame(polygons=names(PApred_SDM), PApred=PApred_SDM),data.frame(polygons=names(PAobs),PAobs),by="polygons")
      # MERGE for the selected polygons : df(polygons of ABpred_SDM, ABpred_SDM) + df(sites of PAobs, PAobs) by POLYGONS --> ABUND data
      DAT_ABUND_SDM = merge(data.frame(polygons=names(ABpred_SDM), ABpred=ABpred_SDM),data.frame(polygons=names(PAobs),PAobs),by="polygons")
      
      # AUC CALCULATION FOR SDMs ------------------------------------------------------------------------------------------------------------
      
      # ["C"] : means that we want the ROC area
      STAT_AUC_SDM[pfg] = somers2(DAT_ABUND_SDM$ABpred[which(DAT_ABUND_SDM$PAobs!=1)], DAT_ABUND_SDM$PAobs[which(DAT_ABUND_SDM$PAobs!=1)]/2)["C"]
      
      CN = sum(DAT_PA_SDM$PApred==0 & DAT_PA_SDM$PAobs==0,na.rm=T) # TRUE NEGATIVE
      CP = sum(DAT_PA_SDM$PApred==1 & DAT_PA_SDM$PAobs==2,na.rm=T) # TRUE POSITIVE
      FN = sum(DAT_PA_SDM$PApred==0 & DAT_PA_SDM$PAobs==2,na.rm=T) # FALSE NEGATIVE
      FP = sum(DAT_PA_SDM$PApred==1 & DAT_PA_SDM$PAobs==0,na.rm=T) # FALSE POSITIVE
      STAT_SPECIF_SDM[pfg] = CN/(CN+FP) # SPECIFICITY : True negative rate, proportion of negatives which are correctly identified as such
      STAT_SENSIB_SDM[pfg] = CP/(CP+FN) # SENSITIVITY : True positive rate, proportion of positives which are correctly identified as such
      
      ErRate_SDM[pfg] = (FN + FP) / (CN + FP + FN + CP) # Error rate
      
    }#end pfg
    
    # saving results ------------------------------------------------
    res_sdm = data.frame(STAT_AUC_SDM, STAT_SPECIF_SDM, STAT_SENSIB_SDM, ErRate_SDM)
    save(res_sdm,file=paste(results.dir,"res_STAT_SDMS",sep=""))
    }
  }
  
  #========================================================================
  # 4. Comparisons STATS - HABITAT SUITABILITY vs FATE for all PFGS
  #       GRAPHS
  #========================================================================
  
  res_sdm = get(load(paste(results.dir,"res_STAT_SDMS",sep="")))
  pdf(paste(results.dir,"GRAPH_PFG_AUC_SENS_SPEC_HS.pdf",sep=""),width=12,height=6)
  
  for(i in 1:ncol(res_sdm)){
    print(colnames(res_sdm)[i])
    barplot(res_sdm[7:17,i],col="grey",ylab=toupper(colnames(res_sdm[i])),main="herbaceous PFGs",ylim=c(0,1),names.arg=rownames(res_sdm)[7:17])
    abline(h=seq(0,1,0.2),lty=2,col="grey")
    barplot(res_sdm[c(1:6,18:24),i],col="grey",ylab=toupper(colnames(res_sdm[i])),main="woody PFGs",ylim=c(0,1),names.arg=rownames(res_sdm)[c(1:6,18:24)])
    abline(h=seq(0,1,0.2),lty=2,col="grey")
  }
  
#   # AUC ------------------------------------------------------------------------------------------
#   barplot(res_sdm$STAT_AUC_SDM[1:10],col="grey",ylab="AUC",main="herbaceous PFGs",ylim=c(0,1),names.arg=rownames(res_sdm)[1:10])
#   abline(h=seq(0,1,0.2),lty=2,col="grey")
#   barplot(res_sdm$STAT_AUC_SDM[11:24],col="grey",ylab="AUC",main="woody PFGs",ylim=c(0,1),names.arg=rownames(res_sdm)[11:24])
#   abline(h=seq(0,1,0.2),lty=2,col="grey")
#   
#   # SENSITIVITY ---------------------------------------------------------------------------------------
#   barplot(res_sdm$STAT_SENSIB_SDM[1:10],col="grey",ylab="Sensitivity",main="herbaceous PFGs",ylim=c(0,1),names.arg=rownames(res_sdm)[1:10])
#   abline(h=seq(0,1,0.2),lty=2,col="grey")
#   barplot(res_sdm$STAT_SENSIB_SDM[11:24],col="grey",ylab="Sensitivity",main="woody PFGs",ylim=c(0,1),names.arg=rownames(res_sdm)[11:24])
#   abline(h=seq(0,1,0.2),lty=2,col="grey")
#   
#   # SPECIFICITY ---------------------------------------------------------------------------------------
#   barplot(res_sdm$STAT_SPECIF_SDM[1:10],col="grey",ylab="Specificity",main="herbaceous PFGs",ylim=c(0,1),names.arg=rownames(res_sdm)[1:10])
#   abline(h=seq(0,1,0.2),lty=2,col="grey")
#   barplot(res_sdm$STAT_SPECIF_SDM[11:24],col="grey",ylab="Specificity",main="woody PFGs",ylim=c(0,1),names.arg=rownames(res_sdm)[11:24])
#   abline(h=seq(0,1,0.2),lty=2,col="grey")
  
  dev.off()

# }


