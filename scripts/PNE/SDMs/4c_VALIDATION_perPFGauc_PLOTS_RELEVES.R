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
library(PresenceAbsence)

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
projLCC <- "+proj=lcc +lat_1=45.89891888888889 +lat_2=47.69601444444444 +lat_0=46.8 +lon_0=2.337229166666667 +x_0=600000 +y_0=2200000 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

# maskSimul <- get(load(paste(path_output,"maskParcEcrinsSansGlacierNiLacCONSENSUS.raster",sep="")))
maskSimul <- get(load(file=paste(path_output,"maskParcEcrins.raster",sep="")))

setwd(paste(path_output,"DATA_",sce,sep=""))

validObj.dir = paste(path_input,"Objects/ABUND_DATA/",sep="") # directory of validation objects
SDMhab.dir = paste(path_output,"DATA_",sce,"/HS_Current/cartes/",sep="") # directory of PFG SDMs
results.dir = paste(path_output,"DATA_",sce,"/SDM_validation/",sep="") # directory to save results
dir.create(results.dir)

#===============================================
# 1. expectedMaps / Observations
#===============================================

# MAPS ---------

#Recoded delphine habitats from 26-code categories (check Isa's script)
#recoded into 7 large habitat categories (see HabCorrespTab)
#only includes habitats in polygons selected for validation
ValidHabitats = raster(paste(validObj.dir,"delphCODE_validationHabs.img",sep=""))
HabCorrespTab <- read.table(paste(validObj.dir,"HabCorrespTab.txt",sep=""),head=T)     
habitats <- raster(paste(validObj.dir,"RECODED_categories_physio.asc",sep=""))
plot(habitats)
#  0 = glaciers et neiges éternelles (= 0 in RECODED)
# 11 = lacs et mares (= 0 in RECODED)
# 14 = ravines, eaux vives (= 0 in RECODED)
# 20 = marais, eaux stagnantes (= 0 in RECODED)
# 31 = roches non colonisées (grouped with 36 in RECODED)
# 36 = roches en voie de colonisation (grouped with 31 in RECODED)
# 40 = pelouses et prairies
# 50 = landes basses
# 60 = milieux ouverts, broussailles
# 70 = milieux semi-fermés
# 81 = milieux fermés
# 83 = forêts
# 90 = milieux fortement artificialisés (= 0 in RECODED)


#===============================================
# Observed releves data
#===============================================

# Observed PFGs Releves ---------
plots = read.dbf(paste(validObj.dir,'plots.dbf',sep=""))
rownames(plots) <- plots$numchrono
# PFG_AB_XY = get(load(paste(validObj.dir,"PFG_AB_XY_uni",sep=""))) #One species name per CBNA_code (but PFG might be repeated)
PFGxPlots = get(load(paste(validObj.dir,"PFGxPlots",sep="")))  
plots <- plots[rownames(PFGxPlots),]

# Transform coordinates in good projection system
xy_in <- SpatialPoints(as.data.frame(plots[,c("X_LII","Y_LII")]), proj4string=CRS(projLCC) )
xy_out <- spTransform(xy_in, CRS(projETRS89))

# Plot the data to check ---------
plot(maskSimul)
plot(xy_out,add=T)
sp_plots = as.data.frame(xy_out)
sp_plots$delphCODE <- extract(habitats, xy_in)
sp_plots$ValidHab <- extract(ValidHabitats, xy_in)
sp_plots$ValidHab[sp_plots$ValidHab==0] <- NA
rownames(sp_plots) <- plots$numchrono

# Habitat codes ---------
hab.labs <- c("Excluded","Rock", "Grasslands", "Lowlands","Open_habitats", "Semi-closed_habitats", "Closed_habitats","Forests")
sp_plots$delphCODE <- as.factor(sp_plots$delphCODE)
hab.tab <- data.frame(code=levels(sp_plots$delphCODE),labs=hab.labs)

#===============================================
# 3. Comparisons 
#===============================================

#===============================================
# SIMULATED DATA: Calculate Abundance based goodness of prediction
#===============================================

pfgs = c(paste("C",1:6,sep=""),paste("H",1:10,sep=""),paste("P",1:8,sep=""))
#create summary table to store stats
sumTab<-matrix(NA,nrow=length(pfgs),ncol=2+2*length(hab.labs),dimnames=list(pfgs,paste(rep(c("Overall",hab.labs), each=2),c(".r",".p"), sep="")))

for(fg in pfgs){
  print(fg)
  
  if(file.exists(paste(SDMhab.dir,"HS_f0_",fg,".asc",sep=""))){
    # load simulated PFG distribution in releves
    sim_PFGxPlots <- raster(paste(SDMhab.dir,"HS_f0_",fg,".asc",sep=""))
    sim_PFGxPlots <- extract(sim_PFGxPlots,xy_out)*100
    
    #Overall performance
    tst<-cor.test(PFGxPlots[,fg],sim_PFGxPlots)
    sumTab[fg,1] <-tst$estimate
    sumTab[fg,2] <-tst$p.value
    
    #Per habitat
    for(hab in levels(plots$delphCODE)){
      tst<-cor.test(PFGxPlots[plots$delphCODE==hab,fg],sim_PFGxPlots[plots$delphCODE==hab,pfg])
      sumTab[pfg,paste(hab.tab[hab.tab$code==hab,"labs"],".r", sep="")] <-tst$estimate
      sumTab[pfg,paste(hab.tab[hab.tab$code==hab,"labs"],".p", sep="")] <-tst$p.value
    }
  }
}
write.table(sumTab,paste(results.dir,"PFG_STATS_Releves_Abund.txt",sep=""))
write.table(sumTab[,grep(".r",colnames(sumTab),fixed=T)], paste(results.dir,"PFG_STATS_Releves_Abund_r.txt",sep=""))

#===============================================
# SIMULATED DATA: Calculate P/A based goodness of prediction
#===============================================      

#create summary tables to store stats
sumTab <- matrix(NA,nrow=length(pfgs),ncol=5,dimnames=list(pfgs,c("sensitivity","specificity","tss","AUC","AUC_ab")))
AUC_ab.hab <- matrix(NA,nrow=length(pfgs),ncol=length(hab.labs)+1,dimnames=list(pfgs,c("overall",hab.labs)))

# PFG predicted distribution in observed releves
for(fg in pfgs){
  print(fg)
  
  if(file.exists(paste(SDMhab.dir,"HS_f0_",fg,".asc",sep=""))){
    # load simulated PFG distribution in releves
    sim_PFGxPlots <- raster(paste(SDMhab.dir,"HS_f0_",fg,".asc",sep=""))
    sim_PFGxPlots <- extract(sim_PFGxPlots,xy_out)*100
    
    # PRESENCE-ABSENCE data in simulated data
    # keeping prevalence as in observed data
    sim_PA = ifelse(sim_PFGxPlots >= quantile(sim_PFGxPlots,sum(PFGxPlots[,fg]==0)/nrow(PFGxPlots), na.rm=T), 1, 0)
    print(quantile(sim_PFGxPlots,sum(PFGxPlots[,fg]==0)/nrow(PFGxPlots),na.rm=T)) ## percentage of observed absences and the value of abundance at this percentage
    
    # bind with observed
    dat <- data.frame(id=rownames(PFGxPlots),PAobs=as.numeric(PFGxPlots[,fg]>0),PApred=sim_PA,ABobs=PFGxPlots[,fg],ABpred=sim_PFGxPlots)
    
    #Overall performance
    test <- data.frame(presence.absence.accuracy(dat[,1:3],st.dev=F,na.rm=T))
    test$tss <- (test$sensitivity+test$specificity)-1
    AUC_ab.hab[fg,"overall"] <- test$AUC_ab <- somers2(dat$ABpred,dat$PAobs)["C"] ## Predicted abundances
    sumTab[fg,] <- as.numeric(test[,colnames(sumTab)])
    
    #Per habitat
    for(hab in levels(plots$delphCODE)){
      AUC_ab.hab[fg,as.character(hab.tab[hab.tab$code==hab,"labs"])] <- somers2(dat$ABpred[plots$delphCODE==hab], dat$PAobs[plots$delphCODE==hab])["C"]
    }
  }
}
write.table(sumTab,paste(results.dir,"PFG_STATS_Releves_PA.txt",sep="")) 
write.table(AUC_ab.hab,paste(results.dir,"PFG_AUCab_Releves_Habitats.txt",sep=""))  



#========================================================================
# 4. Comparisons STATS - HABITAT SUITABILITY vs FATE for all PFGS
#       GRAPHS
#========================================================================

res_sdm = read.table(paste(results.dir,"PFG_STATS_Releves_PA.txt",sep=""))
pdf(paste(results.dir,"GRAPH_PFG_AUC_SENS_SPEC_HS_releves.pdf",sep=""),width=12,height=6)

for(i in 1:ncol(res_sdm)){
  barplot(res_sdm[7:17,i],col="grey",ylab=toupper(colnames(res_sdm[i])),main="herbaceous PFGs",ylim=c(0,1),names.arg=rownames(res_sdm)[7:17])
  abline(h=seq(0,1,0.2),lty=2,col="grey")
  barplot(res_sdm[c(1:6,18:24),i],col="grey",ylab=toupper(colnames(res_sdm[i])),main="woody PFGs",ylim=c(0,1),names.arg=rownames(res_sdm)[c(1:6,18:24)])
  abline(h=seq(0,1,0.2),lty=2,col="grey")
}

dev.off()
