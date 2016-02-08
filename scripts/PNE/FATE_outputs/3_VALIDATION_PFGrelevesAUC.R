### HEADER #####################################################################
##' @title FATE-HD Simulations post-treatment : Validation with releves data
##'
##' @author Marta C. & Maya G.
##' @contact maya.gueguen at gmail.com
##' 
##' @date 21/01/2016
##' 
##' @description This script is designed to 
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
library(foreign)
library(Hmisc)
library(PresenceAbsence)

#################################################################################################
### PARAMETERS OF THE SIMULATION (FATE)
#################################################################################################

projETRS89 <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
projLCC <- "+proj=lcc +lat_1=45.89891888888889 +lat_2=47.69601444444444 +lat_0=46.8 +lon_0=2.337229166666667 +x_0=600000 +y_0=2200000 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

## Retrieve input args ---------------------------------------------------------
args <- commandArgs(trailingOnly=TRUE)
simul_name = as.character(args[1]) ## give the simulation name
file_name = as.character(args[2]) ## give the paramSimul file name

## Set directories -------------------------------------------------------------
simul_folder = paste("/nfs_scratch/mayagueguen/FATE_newHS/", simul_name,"/",sep="")
setwd(sub(simul_name,"",simul_folder))

validObj.dir = "Objects/" ## directory of objects to be used for validation
results.dir = "/resultsVALID_objects/" ## directories where validation results are saved
resultsgraph.dir = "/resultsVALID_graphs/" ## directory to plot graphs
dir.create(paste(simul_name, results.dir, sep = ""))
dir.create(paste(simul_name, resultsgraph.dir, sep=""))


#################################################################################################
### PARAMETERS OF THE SIMULATION (FATE)
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
  
  ## Get PFG names
  PFG = paramSim_lines[(which(paramSim_lines=="--PFG_LIFE_HISTORY_PARAMS--")+1):(which(paramSim_lines=="--PFG_LIFE_HISTORY_PARAMS--")+nb_PFGS)]
  pattern = paste(simul_name,"/DATA/PFGS/SUCC/SUCC_",sep="")
  PFG = sub(".txt","",sub(paste(simul_name,"/DATA/PFGS/SUCC/SUCC_",sep=""),"",PFG))
  PFG_files = paramSim_lines[(which(paramSim_lines=="--PFG_ENVSUIT--")+1):(which(paramSim_lines=="--PFG_ENVSUIT--")+nb_PFGS)]
} else { cat("\n PARAM SIMUL FILE ",file_name," does not exist !!\n") }


#################################################################################################
### VALIDATION OBJECTS
#################################################################################################

elevation = raster(paste(validObj.dir,"ABUND_DATA/elevation_ecrins.asc",sep=""))

## Delphine habitat from 26-code categories recoded into 7 large habitat categories (see HabCorrespTab)
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

## Only includes habitats in polygons selected for validation ??
ValidHabitats <- raster(paste(validObj.dir,"ABUND_DATA/delphCODE_validationHabs.img",sep=""))
HabCorrespTab <- read.table(paste(validObj.dir,"ABUND_DATA/HabCorrespTab.txt",sep=""),head=T)     
habitats <- raster(paste(validObj.dir,"ABUND_DATA/RECODED_categories_physio.asc",sep=""))

## Get observed PFG releves ----------------------------------------------------
plots <- read.dbf(paste(validObj.dir,'ABUND_DATA/plots.dbf',sep="")) ## coordinates
rownames(plots) <- plots$numchrono
PFGxPlots <- get(load(paste(validObj.dir,"ABUND_DATA/PFGxPlots",sep=""))) ## PFGxSITES matrix
plots <- plots[rownames(PFGxPlots),]

## Transform coordinates into the right projection system
xy_in <- SpatialPoints(as.data.frame(plots[,c("X_LII","Y_LII")]), proj4string=CRS(projLCC) )
xy_out <- spTransform(xy_in, CRS(projETRS89))
XY <- xy_out

# Plot the data to check ---------
# plot(mask)
# plot(XY,add=T)

## Gather data : coordinates, habitats types -----------------------------------
sp_plots = as.data.frame(XY)
sp_plots$delphCODE <- extract(habitats, xy_in)
sp_plots$ValidHab <- extract(ValidHabitats, xy_in)
sp_plots$ValidHab[sp_plots$ValidHab==0] <- NA
rownames(sp_plots) <- plots$numchrono

## Get Habitat codes -----------------------------------------------------------
hab.labs <- c("Excluded","Rock", "Grasslands", "Lowlands","Open_habitats", "Semi-closed_habitats", "Closed_habitats","Forests")
sp_plots$delphCODE <- as.factor(sp_plots$delphCODE)
hab.tab <- data.frame(code=levels(sp_plots$delphCODE),labs=hab.labs)

#################################################################################################
### PREDICTED DISTRIBUTION FOR SDM and FATE-HD SIMULATION
#################################################################################################

## Load FATE-HD results --------------------------------------------------------
arrayPFG <- get(load(paste(simul_name,"/outputsTables/arrayPFG_year840_rep",i,sep="")))
cat("\n",colnames(arrayPFG),"\n")

## Rearrange colnames (PFG names) of output of FATE-HD simulation --------------
if(length(strsplit(colnames(arrayPFG), "_")[[1]])>2) {
  cut <- length(strsplit(colnames(arrayPFG), "_")[[1]]) - 2
  pattern2 <- paste(strsplit(colnames(arrayPFG), "_")[[1]][1:cut], collapse="_")
  colnames(arrayPFG) = sub(paste(pattern2,"_",sep=""),"",colnames(arrayPFG))
}
cat("\n",colnames(arrayPFG),"\n")

# Create empty spatial points objects to record predicted PFGs -----------------
PFG_mat <- matrix(NA,nrow=nrow(sp_plots),ncol=nb_PFGS,dimnames=list(rownames(sp_plots),PFG))
sp_PFGs_SDM <- sp_PFGs_simu <- data.frame(PFG_mat)


# PFG predicted distribution in observed releves -------------------------------
for(pfg in 1:length(PFG)){
  cat("\n PFG : ",PFG[pfg],"\n")
  
  ## Load SDM predicted PFG abundances:
  SDM = raster(PFG_files[pfg])
  
  ## Create FATE predicted abundances maps
  simul = mask
  simul[]= arrayPFG[,PFG[pfg]]
  
  # Record simulated/predicted abundances in releves
  sp_PFGs_simu[,PFG[pfg]]  <- extract(simul, XY)
  sp_PFGs_SDM[,PFG[pfg]]  <- extract(SDM, XY)
}

#Save
write.table(sp_PFGs_simu,file=paste(simul_name,results.dir,"PFG_distri_Releves_",i,".txt",sep=""))
write.table(sp_PFGs_SDM, file=paste(simul_name,results.dir,"PFG_distri_Releves_SDM.txt",sep=""))	

#################################################################################################
### SIMULATED DATA: Calculate Abundance based goodness of prediction
#################################################################################################

## Load simulated PFG distribution for releves ---------------------------------
sim_PFGxPlots <- read.table( file=paste(simul_name,results.dir,"PFG_distri_Releves_",i,".txt",sep=""),head=T)

## Create summary table to store stats -----------------------------------------
sumTab <- matrix(NA,nrow=ncol(sim_PFGxPlots),ncol=2+2*length(hab.labs),dimnames=list(colnames(sim_PFGxPlots),paste(rep(c("Overall",hab.labs), each=2),c(".r",".p"), sep="")))

## PFG predicted distribution performance --------------------------------------
for(pfg in colnames(sim_PFGxPlots)){
  cat("\n PFG : ",pfg,"\n")
  subpfg = strsplit(pfg,"_")[[1]][1]
  
  ## Overall performance
  tst <- cor.test(PFGxPlots[,subpfg],sim_PFGxPlots[,pfg])
  sumTab[pfg,1] <- tst$estimate
  sumTab[pfg,2] <- tst$p.value
  
  ## Per habitat
  for(hab in levels(plots$delphCODE)){
    tst <- cor.test(PFGxPlots[plots$delphCODE==hab,subpfg],sim_PFGxPlots[plots$delphCODE==hab,pfg])
    sumTab[pfg,paste(hab.tab[hab.tab$code==hab,"labs"],".r", sep="")] <- tst$estimate
    sumTab[pfg,paste(hab.tab[hab.tab$code==hab,"labs"],".p", sep="")] <- tst$p.value
  }
}
write.table(sumTab,paste(simul_name,results.dir,"PFG_STATS_Releves_Abund_",i,".txt",sep="")) 
write.table(sumTab[,grep(".r",colnames(sumTab),fixed=T)],paste(simul_name,results.dir,"PFG_STATS_Releves_Abund_ESTIMATES_",i,".txt",sep=""))

#################################################################################################
### SIMULATED DATA: Calculate P/A based goodness of prediction
#################################################################################################

## Load simulated PFG distribution for releves ---------------------------------
sim_PFGxPlots <- read.table( file=paste(simul_name,results.dir,"PFG_distri_Releves_",i,".txt",sep=""),head=T)

## Create summary table to store stats -----------------------------------------
sumTab <- matrix(NA,nrow=ncol(sim_PFGxPlots),ncol=5,dimnames=list(colnames(sim_PFGxPlots),c("sensitivity","specificity","tss","AUC","AUC_ab")))
AUC_ab.hab <- matrix(NA,nrow=ncol(sim_PFGxPlots),ncol=length(hab.labs)+1,dimnames=list(colnames(sim_PFGxPlots),c("overall",hab.labs)))

## PFG predicted distribution performance --------------------------------------
for(pfg in colnames(sim_PFGxPlots)){
  cat("\n PFG : ",pfg,"\n")
  subpfg = strsplit(pfg,"_")[[1]][1]
  
  ## PRESENCE-ABSENCE data in simulated data : keeping prevalence as in observed data
  sim_PA = ifelse(sim_PFGxPlots[,pfg] >= quantile(sim_PFGxPlots[,pfg], sum(PFGxPlots[,subpfg]==0)/nrow(PFGxPlots), na.rm=T), 1, 0)
  cat("\n QUANTILE : ",quantile(sim_PFGxPlots[,pfg], sum(PFGxPlots[,subpfg]==0)/nrow(PFGxPlots),na.rm=T),"\n") ## percentage of observed absences and the value of abundance at this percentage
  
  ## Bind with observed
  dat <- data.frame(id=rownames(PFGxPlots), PAobs=as.numeric(PFGxPlots[,subpfg]>0), PApred=sim_PA,ABobs=PFGxPlots[,subpfg], ABpred=sim_PFGxPlots[,pfg])
  
  ## Overall performance
  test <- data.frame(presence.absence.accuracy(dat[,1:3], st.dev=F))
  test$tss <- ((test$sensitivity+test$specificity)-1)
  AUC_ab.hab[pfg,"overall"] <- test$AUC_ab <- somers2(dat$ABpred, dat$PAobs)["C"] ## Predicted abundances
  sumTab[pfg,] <- as.numeric(test[,colnames(sumTab)])
  
  ## Per habitat
  for(hab in levels(plots$delphCODE)){
    AUC_ab.hab[pfg,as.character(hab.tab[hab.tab$code==hab,"labs"])] <- somers2(dat$ABpred[plots$delphCODE==hab], dat$PAobs[plots$delphCODE==hab])["C"]
  }
}
write.table(sumTab, paste(simul_name, results.dir, "PFG_STATS_Releves_PA_",i,".txt",sep="")) 
write.table(AUC_ab.hab, paste(simul_name, results.dir, "PFG_AUCab_Releves_Habitats_",i,".txt",sep=""))  

#################################################################################################
### PLOT : Comparisons STATS - HABITAT SUITABILITY vs FATE for all PFGS
#################################################################################################

res_sdm = read.table(paste(simul_name, results.dir,"PFG_STATS_Releves_PA_",i,".txt",sep=""))
pdf(paste(simul_name, results.dir,"GRAPH_PFG_AUC_SENS_SPEC_HS_releves.pdf",sep=""),width=12,height=6)
for(stati in 1:ncol(res_sdm)){
  barplot(res_sdm[grep("H",rownames(res_sdm)),stati],col="grey",ylab=toupper(colnames(res_sdm[stati])),main="herbaceous PFGs",ylim=c(0,1),names.arg=rownames(res_sdm)[grep("H",rownames(res_sdm))])
  abline(h=seq(0,1,0.2),lty=2,col="grey")
  barplot(res_sdm[-grep("H",rownames(res_sdm)),stati],col="grey",ylab=toupper(colnames(res_sdm[stati])),main="woody PFGs",ylim=c(0,1),names.arg=rownames(res_sdm)[-grep("H",rownames(res_sdm))])
  abline(h=seq(0,1,0.2),lty=2,col="grey")
}
dev.off()


#################################################################################################
### PLOT : SIMULATED ABUNDANCE PER HABITAT
#################################################################################################

## Load simulated PFG distribution for releves ---------------------------------
sim_PFGxPlots <- read.table( file=paste(simul_name,results.dir,"PFG_distri_Releves_",i,".txt",sep=""),head=T)

## Boxplots in relation to habitat (Ceres Habitats)
pdf(paste(simul_name,resultsgraph.dir,"PFG_AbundxHabitat_RELEVESvsSimu",i,".pdf",sep=""),height=14)
for(pfg in colnames(sim_PFGxPlots)){
  subpfg = strsplit(pfg,"_")[[1]][1]
  
  if(sum(sim_PFGxPlots[,pfg]>0)>0){
    op = par(mfrow=c(2,1),mar=c(6,4,4,2))
    boxplot(PFGxPlots[PFGxPlots[,subpfg]>0,subpfg] ~ sp_plots$delphCODE[PFGxPlots[,subpfg]>0],
            col="grey", ylab=pfg, main= "Observed abundance",
            names=c("Excluded","Rock", "Grasslands", "Lowlands","Open\nhabitats", "Semi-closed\nhabitats", "Closed\nhabitats","Forests"), las=2)
    boxplot(sim_PFGxPlots[sim_PFGxPlots[,pfg]>0,pfg]~ sp_plots$delphCODE[sim_PFGxPlots[,pfg]>0],
            col="grey", ylab=pfg, main= "Simulated abundance",
            names=c("Excluded","Rock", "Grasslands", "Lowlands","Open\nhabitats", "Semi-closed\nhabitats", "Closed\nhabitats","Forests"), las=2)
    par(op)
  }
}
dev.off()  

## Boxplots in relation to habitat (Ceres Habitats) WITH THE SAME COVER SCALE --
sim_PFGxPlots_cov <- (sim_PFGxPlots/rowSums(sim_PFGxPlots))*100            
PFGxPlots_cov <- (PFGxPlots/rowSums(PFGxPlots))*100 ## standardize cover of observed

pdf(paste(simul_name,resultsgraph.dir,"PFG_CoverxHabitat_RelevesvsSimu",i,".pdf",sep=""),height=14)
for(pfg in colnames(sim_PFGxPlots)){
  subpfg = strsplit(pfg,"_")[[1]][1]
  
  if(sum(sim_PFGxPlots_cov[,pfg]>0, na.rm=T)>0){
    op=par(mfrow=c(2,1), mar=c(6,4,4,2))
    #observed
    boxplot(PFGxPlots_cov[PFGxPlots_cov[,subpfg]>0,subpfg] ~ sp_plots$delphCODE[PFGxPlots_cov[,subpfg]>0],
            col="grey", ylab=pfg, main= "Observed abundance",
            names=c("Excluded","Rock", "Grasslands", "Lowlands","Open\nhabitats", "Semi-closed\nhabitats", "Closed\nhabitats","Forests"),
            las=2)
    #simulated
    boxplot(sim_PFGxPlots_cov[sim_PFGxPlots_cov[,pfg]>0,pfg]~ sp_plots$delphCODE[sim_PFGxPlots_cov[,pfg]>0], 
            col="grey", ylab=pfg, main= "Simulated abundance", 
            names=c("Excluded","Rock", "Grasslands", "Lowlands","Open\nhabitats", "Semi-closed\nhabitats", "Closed\nhabitats","Forests"), 
            las=2)
    par(op)
  }
}
dev.off()  



# In relation to habitat - proportion of cover across habitats
# (in which habitat is the PFG most abundant? or in other words,
# of all the abundance of the PFG in the releves what proportion is in each habitat?)
sim_PFGxHab_prop <- obs_PFGxHab_prop <- matrix(NA,nrow=length(levels(sp_plots$delphCODE)),ncol=ncol(sim_PFGxPlots),dimnames=list(levels(sp_plots$delphCODE),colnames(sim_PFGxPlots)))

## Per PFG ---------------------------------------------------------------------
pdf(paste(simul_name, resultsgraph.dir, "PFGxHabitat_Prop_Releves_",i,".pdf",sep=""))
for(pfg in colnames(sim_PFGxPlots)){   
  subpfg = strsplit(pfg,"_")[[1]][1]
  
  sim_PFGxHab_prop[,pfg] <- tapply(sim_PFGxPlots[,pfg],sp_plots$delphCODE,function(x){(sum(x)/sum(sim_PFGxPlots[,pfg])*100)})
  obs_PFGxHab_prop[,pfg] <- tapply(PFGxPlots[,subpfg],sp_plots$delphCODE,function(x){(sum(x)/sum(PFGxPlots[,subpfg])*100)})
  
  #Plot
  op = par(mar=c(6,5,4,2))
  tmp <- cbind(sim=sim_PFGxHab_prop[,pfg],obs=obs_PFGxHab_prop[,pfg])
  barplot(t(tmp), beside=T, las=2, legend.text = T, main=pfg,
          names.arg=c("Excluded","Rock", "Grasslands", "Lowlands","Open\nhabitats", "Semi-closed\nhabitats", "Closed\nhabitats","Forests"), 
          ylab="Proportion of Abundance",
          args.legend=list(x="topleft", bty="n"), 
  )
  par(op)
}
dev.off()  

## Overall trend ---------------------------------------------------------------
diff <- sim_PFGxHab_prop - obs_PFGxHab_prop 
pdf(paste(simul_name, resultsgraph.dir, "PFGxHabitat_Error_Releves_",i,".pdf",sep=""), width=14)
palette(rainbow(8))	
op=par(mar=c(8,5,4,2))  
barplot(diff, beside=T, ylim=c(-40,40), las=2, col=1:8, ylab="Underprediction          Overprediction")
legend("topleft", fill=1:8, legend=c("Excluded","Rock", "Grasslands", "Lowlands","Open\nhabitats", "Semi-closed\nhabitats", "Closed\nhabitats","Forests"), bty="n",ncol=4, cex=0.9)
par(op)
dev.off() 

write.table(diff ,file=paste(simul_name, results.dir, "PFGxHabitat_Error_Releves_",i,".txt",sep=""))


#################################################################################################
### PLOT : SOIL VALUES
#################################################################################################

if(do_soil){
  
  ## PFG soil parameters ---------------------------------------------------------
  params <- read.csv2("editParamsC24.csv",sep=",")
  rownames(params) <- unlist(lapply(strsplit(as.character(params$group),"_"), function(x){x[[1]]}))
  params$Soil_tol_ave <- rowMeans(cbind(as.numeric(as.character(params[,"Soil_tol_min"])), as.numeric(as.character(params[,"Soil_tol_max"]))))
  params$Soil_contrib <- as.numeric(as.character(params$Soil_contrib)) 
  
  ## CALCULATE CWM OF OBSERVED SOIL CONTRIBUTION ---------------------------------
  .libPaths("/home/mayagueguen/luke_libs")
  library(FD)
  cwms <- functcomp(as.matrix(params[,c("Soil_contrib","Soil_tol_ave")]), as.matrix(PFGxPlots[,-25]))
  # cor.test(cwms[,1],cwms[,2])
  sp_plots$Soil_contrib <- cwms$Soil_contrib
  sp_plots$Soil_tol_ave <- cwms$Soil_tol_ave
  
  ## Potential soil values compared to Delphine codes, BASED ON OBSERVED ABUNDANCES
  pdf(paste(simul_name,resultsgraph.dir,"Potential_Soil_per_Habitat_DElph.pdf",sep=""),width=14)
  op=par(mfrow=c(1,2), mar=c(6,4,4,2))
  boxplot(sp_plots$Soil_contrib ~ sp_plots$delphCODE, col="grey", ylab="Soil Value", main= "Soil Contribution (LNC)", names=hab.tab$labs, las=2)
  boxplot(sp_plots$Soil_tol_ave ~ sp_plots$delphCODE, col="grey", ylab="Soil Value", main= "Soil Tolerance (Landolt)", names=hab.tab$labs, las=2)
  par(op)
  dev.off()
  
  
  ## Load the full soil values across years and extract 1 year
  Soil_val <- get(load(paste(simul_name,"/outputsTables/arraySoil_rep",i,sep="")))
  Soil_val <- Soil_val["840",] #Extract 1 year
  
  simul = mask ## create maps for soil
  simul[mask[]==1]= Soil_val ## create FATE predicted soil maps 
  sp_plots$SoilVal_simu <- extract(simul, XY) ## store in releves
  
  ## SIMULATED SOIL VALUES -------------------------------------------------------
  maskSoil = mask
  maskSoil[] = sp_plots$SoilVal_simu
  pdf(paste(simul_name,resultsgraph.dir,"Simulated_Soil_values_Map.pdf",sep=""),width=14)
  op = par(mfrow=c(1,1))
  plot(maskSoil)
  contour(elevation, add=T)
  par(op)
  dev.off()
  
  
  ## SIMULATED SOIL VALUES PER HABITAT -------------------------------------------
  pdf(paste(simul_name,resultsgraph.dir,"/Simulated_Soil_per_Habitat_DElph_Releves_",i,".pdf",sep=""),width=14)
  op = par(mar=c(9,4,4,2))
  boxplot(sp_plots$SoilVal_simu ~ sp_plots$delphCODE,col="grey",ylab="Soil Value",main=paste("Simu SoilVvalue",simul_name),names=hab.tab$labs,las=2)
  par(op)
  dev.off()
  
  ## OBSERVED SOIL CONTRIB & TOLERANCE vs SIMULATED SOIL VALUES ------------------
  pdf(paste(simul_name,resultsgraph.dir,"Potential_Soil_values_Map_Observed.pdf",sep=""),width=14)
  op=par(mfrow=c(1,3))
  for(sol in c("Soil_contrib","Soil_tol_ave","SoilVal_simu")){
    palette(heat.colors(5))
    image(mask, col="grey", asp=1, xlab="", ylab="", xaxt="n", yaxt="n", bty="n", bg=1, main=paste("Observed", sol))
    # plot(elevation, asp=1, xlab="", ylab="", xaxt="n", yaxt="n", bty="n", bg=1, main=paste("Observed", sol))
    points(sp_plots$X_LII, sp_plots$Y_LII, pch=19, col=cut(sp_plots[,sol],5))
    contour(elevation, add=T)
    tst <- cor.test(sp_plots[,sol], extract(elevation, XY))
    legend("topleft",paste("R =", round(tst$estimate, 2), ifelse(tst$p.value<0.05,"*"," ")), bty="n", title="Correlation with altitude")
  }
  par(op)
  dev.off()
  
}
