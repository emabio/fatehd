##
# Data creation for BIOMOD
# run april 2012 - Damien G.
# update 6/9/12 - Isa B.
# update with biomod2 7/12/15 - Maya G.
##

rm(list=ls())
library(raster)
library(rgdal)
library(sqldf)

path_script <- "~/Documents/_BIOMOVE/EX_BIOMOD2/_SP_VERSION/_SCRIPTS/"
path_input <- "~/Documents/_BIOMOVE/EX_BIOMOD2/_SP_VERSION/_INPUT_DATA/"
path_output <- "~/Documents/_BIOMOVE/EX_BIOMOD2/_SP_VERSION/_OUTPUT_DATA/"

projETRS89 <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
projETRS89_bis <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
projAlbers <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs"
projLCC <- "+proj=lcc +lat_1=45.89891888888889 +lat_2=47.69601444444444 +lat_0=46.8 +lon_0=2.337229166666667 +x_0=600000 +y_0=2200000 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

###################################################################################################
# All data used is located in ../_INPUT_DATA
###################################################################################################

###############################
## 1. raster of study area

# DELPHINE map over the entire Ecrins park
r <- raster(paste(path_input,"Delphine_maps/Categories_physio/",sep=""))
r <- projectRaster(from=r,res=100,crs=CRS(projETRS89))

# Ecrins park : whole area
maskParcEcrins <- reclassify(r, c(-Inf, Inf, 1))
save(maskParcEcrins,file=paste(path_output,"maskParcEcrins.raster",sep=""))

# Ecrins park : without glaciers or lakes
maskParcEcrinsSG <- reclassify(r,c(-Inf,12,NA, 12,Inf,1))
save(maskParcEcrinsSG,file=paste(path_output,"maskParcEcrinsSansGlacierNiLac.raster",sep=""))

# Ecrins park : glaciers or lakes
maskGlacier <- reclassify(r,c(-Inf,12,1, 12,Inf,NA))
save(maskGlacier,file=paste(path_output,"maskGlacier.raster",sep=""))

maskSimul <- maskParcEcrinsSG


###################################################################################################
# 0. SELECT OCCURENCES DATA and CORRESPONDING BIOCLIM VARIABLES

# occurence data for all species
load(paste(path_input,"matrice",sep=""))

# bioclimatic data
load(paste(path_input,"PLOTS_FR_Alps",sep=""))
input <- input[rownames(matrice),]

# get coordinates and project into the adequate projection system (here from Albers to ETRS89)
xy_in <- SpatialPoints(as.data.frame(input[,c("X_ALBERS","Y_ALBERS")]), proj4string = CRS(projAlbers) )
xy_out <- spTransform(xy_in, CRS(projETRS89))

input$X_ETRS89 <- xy_out@coords[,1]
input$Y_ETRS89 <- xy_out@coords[,2]

# Add glaciers and lakes coordinates
xy_gla <- xyFromCell(maskGlacier,which(maskGlacier[]==1))
xy1 <- SpatialPoints(xy_gla,proj4string=CRS(projETRS89) )
xy2 <- spTransform(xy1,CRS(projAlbers))
xy_gla <- data.frame(PlotID_KEY=paste("GLACIER-",1:nrow(xy_gla),sep=""),
                     X_ALBERS=xy2@coords[,1],Y_ALBERS=xy2@coords[,2],X_ETRS89=xy1@coords[,1],Y_ETRS89=xy1@coords[,2])
rownames(xy_gla) <- xy_gla$PlotID_KEY
input = rbind(input,xy_gla)

# # Add glaciers and lakes absences
# toAdd = matrix("NA",ncol=ncol(matrice),nrow=length(grep("GLACIER-",input$PlotID_KEY)),
#                dimnames=list(rownames(input)[grep("GLACIER-",input$PlotID_KEY)],colnames(matrice)))
# listeGrp = get(load(paste(path_input,"determinantes",sep=""))) # list species / group
# sp <- sapply(listeGrp,function(x) sub('X','',unlist(x))[1]) # determinantes of group considered
# toAdd[,sp] = apply(toAdd[,sp],2,function(x){
#   tmp = x
#   tmp[sample(1:length(x),50)] = "0"
#   return(tmp)
# })
# matrice = rbind(matrice,toAdd)

# reduce extent
input = sqldf(paste('select * from input where X_ETRS89<=',EXT[2],' AND X_ETRS89>=',EXT[1],'
                     AND Y_ETRS89<=',EXT[4],' AND Y_ETRS89>=',EXT[3],sep=''))
rownames(input) <- input$PlotID_KEY
matrice = matrice[which(rownames(matrice) %in% input$PlotID_KEY),]

# save new input data
save(matrice,file=paste(path_output,"matrice",sep=""))
save(input,file=paste(path_output,"PLOTS_Ecrins",sep=""))

###################################################################################################
# 1. REPROJECT and EXTRACT VALUES FOR OCCURENCES DATA

load(paste(path_output,"PLOTS_Ecrins",sep=""))

# convert coordinates of sites in LCC projection
xyAlb <- SpatialPoints(as.data.frame(input[,c("X_ALBERS","Y_ALBERS")]), proj4string=CRS(projAlbers) )
xyLCC <- spTransform(xyAlb,CRS(projLCC))

varEnv <- c("bio_3","bio_4","bio_7","bio_11","bio_12","carbon","spring","winter","slope","twi","rock")

## load shared data
RAScarbon <- raster(paste(path_input,"carbon/",sep=""))
RAScarbon <- projectRaster(RAScarbon,crs=CRS(projLCC))
RASslope <- raster(paste(path_input,"slope/",sep=""))
RAStwi <- raster(paste(path_input,"twi/",sep=""))
RASrock <- raster(paste(path_input,"rock.img",sep=""))
RASrock[which(is.na(RASrock[]))] <- 0

###---------------------------------------------------------------------------------------------###
### NICK's DATA

path_input_bio <- paste(path_input,"ENV_FR_ALPS/CURRENT/",sep="")
path_input_spr_win <- paste(path_input,"ENV_FR_ALPS/CURRENT/",sep="")

# slope, twi smoothed, bare rock, spring & winter temperature (LCC projection)
for(var in c("bio_3","bio_4","bio_7","bio_11","bio_12")){
  eval(parse(text=paste("RAS",var," <- raster(paste(path_input_bio,var,sep=''))",sep="")))
}
RASspring <- raster(paste(path_input_spr_win,"tave_456.img",sep=""))
RASwinter <- raster(paste(path_input_spr_win,"tave_120102.img",sep=""))

# extract slope, twi and spring temp for points of interest
for(var in varEnv){
  print(var)
  print(projection(get(paste("RAS",var,sep=""))))
  if(projection(get(paste("RAS",var,sep="")))==projETRS89 || projection(get(paste("RAS",var,sep="")))==projETRS89_bis){
    eval(parse(text=paste(var," = extract(RAS",var,",input[,c('X_ETRS89','Y_ETRS89')])",sep="")))
  } else if(projection(get(paste("RAS",var,sep="")))==projLCC){
    eval(parse(text=paste(var," = extract(RAS",var,",xyLCC)",sep="")))
  } else if(projection(get(paste("RAS",var,sep="")))==projAlbers){
    eval(parse(text=paste(var," = extract(RAS",var,",xyAlb)",sep="")))
  }
}

# sharing of data
data.env <- data.frame(bio_3,bio_4,bio_7,bio_11,bio_12,carbon,slope,twi,spring,winter)
data.env$rock <- as.factor(rock)
rownames(data.env) <- rownames(input)
save(data.env,file=paste(path_output,"data.env.NICK",sep=""))

###---------------------------------------------------------------------------------------------###
### AUSTRIA's DATA

path_input_bio <- "/media/gueguen/equipes/emabio/GIS_DATA/Alpes/PHYSIQUE/CLIMATE/Alp_dwnscl_100m_15yr_newprec/EOBS_1970_2005/bio/"
path_input_spr_win <- "/home/gueguen/Documents/_DATA/DATA_FATE_Ecrins/RASTERS_clim/Spring_winter_temp/_AUTR_Scenario_EOBS_rcp26_45_85/"

# slope, twi smoothed, bare rock, spring & winter temperature (LCC projection)
# for(var in c("bio_3","bio_4","bio_7","bio_11","bio_12")){
varEnv <- c("bio_6", "bio_9", "bio_12", "bio_15", "carbon", "slope")
for(var in grep("bio_", varEnv, value = TRUE)){ 
  eval(parse(text=paste("RAS",var," <- raster(paste(path_input_bio,var,'.asc',sep=''))",sep="")))
  eval(parse(text=paste("projection(RAS",var,") <- projETRS89",sep="")))
}
RASspring <- raster(paste(path_input_spr_win,"EOBS_1970_2005_tave_456.img",sep=""))
RASwinter <- raster(paste(path_input_spr_win,"EOBS_1970_2005_tave_120102.img",sep=""))
projection(RASspring) = projection(RASwinter) = projETRS89

# extract slope, twi and spring temp for points of interest
for(var in varEnv){
  print(var)
  print(projection(get(paste("RAS",var,sep=""))))
  if(projection(get(paste("RAS",var,sep="")))==projETRS89 || projection(get(paste("RAS",var,sep="")))==projETRS89_bis){
    eval(parse(text=paste(var," = extract(RAS",var,",input[,c('X_ETRS89','Y_ETRS89')])",sep="")))
  } else if(projection(get(paste("RAS",var,sep="")))==projLCC){
    eval(parse(text=paste(var," = extract(RAS",var,",xyLCC)",sep="")))
  } else if(projection(get(paste("RAS",var,sep="")))==projAlbers){
    eval(parse(text=paste(var," = extract(RAS",var,",xyAlb)",sep="")))
  }
}

# sharing of data
data.env <- data.frame(bio_3,bio_4,bio_7,bio_11,bio_12,carbon,slope,twi,spring,winter)
data.env$rock <- as.factor(rock)
rownames(data.env) <- rownames(input)
save(data.env,file=paste(path_output,"data.env.AUST",sep=""))


###################################################################################################
# 2. FOR EACH PFG, EXTRACT CORRESPONDING OCCURENCES and ENVIRONMENTAL VARIABLES

load(paste(path_output,"PLOTS_Ecrins",sep=""))
# load(paste(path_output,"matrice",sep=""))

# list species / group
listeGrp = get(load(paste(path_input,"determinantes",sep="")))
(nbgr <- length(listeGrp))
listeSP = unlist(listeGrp)

for(sce in c("AUST")){
  # selection
  dir.create(paste(path_output,"DATA_",sce,"/PFT_occ/",sep=""),recursive=T)
  dir.create(paste(path_output,"DATA_",sce,"/PFT_env/",sep=""),recursive=T)
  
  for(sp in listeSP){
    cat("\n", sp, "\n")
    
    ## occurences
    sp <- sub('X','',sp) # determinantes of group considered
    tab.occ <- matrice[,which(colnames(matrice) %in% sp)]
    
    if( length(sp) > 1){
      pft.occ <- apply(tab.occ, 1, function(x){
        if(sum(!is.na(as.numeric(x)))>0){ # at least one value defined
          if(sum(as.numeric(x)>0, na.rm=T) > 0){ # at least one presence
            return(1)
          } else if( sum(as.numeric(x)==0, na.rm=T) == length(x)){ # only 0 values
            return(0)
          } else{ return (NA) }
        } else { return (NA) }
      })
    } else{
      pft.occ <- tab.occ
      n <- names(pft.occ)
      pft.occ <- as.numeric(pft.occ)
      names(pft.occ) <- n
    }
    
    pft.occ <- na.omit(pft.occ) ## occurences
    load(paste(path_output,"data.env.",sce,sep=""))
    pft.env <- na.omit(data.env[which(rownames(data.env) %in% names(pft.occ)),]) ## environment
    
    ## we only keep the intersection
    pft.occ <- pft.occ[which(names(pft.occ) %in% intersect(rownames(pft.env), names(pft.occ)))]
    pft.env <- pft.env[which(rownames(pft.env) %in% intersect(rownames(pft.env), names(pft.occ))),]
    
    save(pft.occ, file=paste(path_output,"DATA_",sce,"/PFT_occ/OCC_X",sp,sep=""))
    save(pft.env, file=paste(path_output,"DATA_",sce,"/PFT_env/ENV_X",sp,sep=""))
  }
}


###################################################################################################
# Creation of tables of environmental data for a defined area for BIOMOD projections
###################################################################################################

###############################
## 2. Variables transformation : project into another coordinate system and on another area

path_ecrins <- paste(path_input,"ENV_ALL_ECRINS/",sep="")
dir.create(path_ecrins)
varEnv <- c("bio_3","bio_4","bio_7","bio_11","bio_12","carbon","slope","twi","spring","winter","rock")
maskSimul <- maskParcEcrinsSG

for(sce in c("AUST")){
  if(sce=="NICK"){
    path_fr_alps <- paste(path_input,"ENV_FR_ALPS/",sep="")
    path_fr_alps_spr_win <- paste(path_input,"ENV_FR_ALPS/",sep="")
    time_frame <- c("CURRENT/","ssmhi_rca30_ccsm_ar4_a1b_clim2150/","ssmhi_rca30_ccsm_ar4_a1b_clim5180/","ssmhi_rca30_ccsm_ar4_a1b_clim9120/")
  } else {
    path_fr_alps <- "/media/gueguen/equipes/emabio/GIS_DATA/Alpes/PHYSIQUE/CLIMATE/Alp_dwnscl_100m_15yr_newprec/"
    path_fr_alps_spr_win <- "/home/gueguen/Documents/_DATA/DATA_FATE_Ecrins/RASTERS_clim/Spring_winter_temp/_AUTR_Scenario_EOBS_rcp26_45_85/"
    # years = seq(2020,2090,10)
    years = 2090
    time_frame <- c("rcp26_ICHEC-EC-EARTH_SMHI-RCA4/","rcp45_CNRM-CERFACS-CNRM-CM5_SMHI-RCA4/","rcp85_ICHEC-EC-EARTH_DMI-HIRHAM5/")
    # time_frame <- sapply(time_frame,function(x) paste(x,"bio",years,"_",sep=""))
    time_frame <- c("EOBS_1970_2005/",time_frame)
  }
  
  for(ti in time_frame){
    dir.create(paste(path_ecrins,ti,sep=""),recursive=T)
    for(v in varEnv){
      cat("\n", ti, v, "\n")
      if(v=="spring"){
        if(sce=="NICK"){ rasName = paste(path_fr_alps_spr_win,ti,"tave_456.img",sep="")
        } else { rasName = paste(path_fr_alps_spr_win,sub("/","_",ti),"tave_456.img",sep="") }
      } else if(v=="winter"){
        if(sce=="NICK"){ rasName = paste(path_fr_alps_spr_win,ti,"tave_120102.img",sep="")
        } else { rasName = paste(path_fr_alps_spr_win,sub("/","_",ti),"tave_120102.img",sep="") }
      } else {
        if(sce=="NICK"){ rasName = paste(path_fr_alps,ti,v,sep="")
        } else {
          if(ti=="EOBS_1970_2005/"){ rasName = paste(path_fr_alps,ti,"bio/",v,".asc",sep="")
          } else { rasName = paste(path_fr_alps,ti,"bio/bio",years,"_",strsplit(v,"_")[[1]][2],sep="") }
        }
      }
      if(v %in% c("carbon","slope","twi")){ rasName = paste(path_input,v,sep="") }
      if(v=="rock"){ rasName = paste(path_input,v,".img",sep="") }
      ras = raster(rasName)
      if(sce=="AUST" && !(v %in% c("carbon","slope","twi"))) projection(ras) <- projETRS89
      if(v=="carbon") ras[ras[]>100] <- NA
      if(v=="rock"){
        ras[which(is.na(ras[]))] <- 0
        ras <- crop(projectRaster(from=ras,to=maskSimul,crs=projETRS89,res=100,methode='ngb'),extent(maskSimul))*maskSimul
        ras[which(ras[]<0.5)] <- 0
        ras[intersect(which(ras[]>0.5),which(ras[]<1.5))] <- 1
        ras[which(ras[]>1.5)] <- 2
      } else { ras <- crop(projectRaster(from=ras,to=maskSimul,crs=projETRS89,res=100,methode=bilinear),extent(maskSimul))*maskSimul }
      writeRaster(ras,file=paste(path_ecrins,ti,v,".img",sep=""))
    }
  }
}

###############################
# 3. Table creation

# standardization of all masks
path_ecrins <- paste(path_input,"ENV_ALL_ECRINS/",sep="")
load(paste(path_output,"maskParcEcrinsSansGlacierNiLac.raster",sep=""))
# time_frame <- c("CURRENT/","ssmhi_rca30_ccsm_ar4_a1b_clim2150/","ssmhi_rca30_ccsm_ar4_a1b_clim5180/",
#                 "ssmhi_rca30_ccsm_ar4_a1b_clim9120/","EOBS_1970_2005/","rcp26_ICHEC-EC-EARTH_SMHI-RCA4/",
#                 "rcp45_CNRM-CERFACS-CNRM-CM5_SMHI-RCA4/","rcp85_ICHEC-EC-EARTH_DMI-HIRHAM5/")
time_frame <- c("EOBS_1970_2005/")
varEnv <- c("bio_3","bio_4","bio_7","bio_11","bio_12","carbon","slope","twi","spring","winter","rock")

maskSimul <- maskParcEcrinsSG
for(v in varEnv){ eval(parse(text=paste(v," = raster(paste(path_ecrins,time_frame[1],v,'.img',sep=''))",sep=""))) }
eval(parse(text=paste("maskSimul <- reclassify(",paste(varEnv,collapse=' * ')," * maskSimul,c(-Inf,Inf,1))",sep="")))
save(maskSimul ,file=paste(path_output,"maskSimulCONSENSUS.raster",sep=""))

for(ti in time_frame){
  for(v in varEnv){
    cat("\n", ti, v, "\n")
    ras = raster(paste(path_ecrins,ti,v,".img",sep=""))
    ras <- ras * maskSimul
    writeRaster(ras,file=paste(path_ecrins,ti,v,'.img',sep=''),overwrite=T)
  }
}







