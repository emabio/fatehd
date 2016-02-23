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

validObj.dir = paste(path_input,"Objects/",sep="") # directory of validation objects
results.dir = paste(path_output,"DATA_",sce,"/SDM_validation/",sep="") # directory to save results
dir.create(results.dir)


#===============================================
# 1. expectedMaps / Observations
#===============================================

# MAPS ---------

# polygon ID per pixel
polygonsID = raster(file.path(validObj.dir,"DELPHINE/polygonsID.img")) # polygon ID per pixel; 18915 different polygons for 484058 cells
polygonsID <- projectRaster(from=polygonsID,res=100,crs=CRS(projETRS89))

# more general DELPHINE codes (26, after NAing 0, 128; for further details see ambiancesPFGs.xls)
delphCODE = raster(file.path(validObj.dir,"DELPHINE/delphCODEreclass.img"))
delphCODE [delphCODE [] %in% c(0,128)] = NA
delphCODE <- projectRaster(from=delphCODE,res=100,crs=CRS(projETRS89),method="ngb")

#any(!is.na(polygonsID[]) == is.na(maskSimul[]))   # some polygons have extra NAs

# Creating new mask, keeping cells within mask AND polygons map
maskCons = maskSimul
maskCons[is.na(polygonsID[])] = NA             
# creating a list of polygons in pixels
polygonPix = polygonsID[which(maskCons[] == 1)]


# HABITAT AND PFG CORRESPONDENCE --------

# Correspondence table for PFGs present in delphCODE, "ambiance" level classification
# 26 habitats x 24 PFGs -  Habitat code, description, and presence/absence of each PFG
compatibTable = read.table(file.path(validObj.dir, "ambiances_v2.txt"), sep="\t", h=T)
# This is not really used in this script. It is just the reference for the values in the delphCODE map which is later used to select the pixels in the analysis and exclude certain habitats
compatibTable$description = c("Snow and ice",
                              "Snow-pack",
                              "Alpine crests on calcareous",
                              "Alpine crests on silica",
                              "Other alpine grasslands",
                              "Rich grazed subalpine grasslands",
                              "Poor grazed subalpine grasslands",
                              "Subalpine natural ubacs",
                              "Subalpine natural adrets",
                              "Subalpine disturbed",
                              "Mowned grasslands",
                              "Mountainous mature ubac forests",
                              "Moutainous young ubac forests",
                              "Moutainous grasslands & heathlands",
                              "Mountainous adret forests",
                              "Grazed alpine grasslands",
                              "Collinean shrublands",
                              "Mediterranean heathlands & forests",
                              "Dry low grasslands",
                              "Seasonal wetlands",
                              "Riverside vegetation",
                              "Formerly mowned grasslands",
                              "Wetlands",
                              "Water vegetation",
                              "Rocks and screes",
                              "Anthropogenic disturbed habitat")


# Correspondance table for PFGs present in representativity (86 habitats, for further details see tableau PFG_reluCD.xls)
# 86 DELPHINE CODE x 24 PFGs (presence/absence)
# 0 = absent; 1 = potentially present; 2 = present (only 0s and 2s are used for validation - conservative approach)
gpementsTable = read.table(file.path(validObj.dir,"groupementsETpfg.txt"),sep="\t",h=T)
rownames(gpementsTable) = toupper(gpementsTable$code.delphine)
head(gpementsTable)
rownames(gpementsTable)


# POLYGON INFORMATION --------
# loading table of polygon IDs (19959) and their characteristics (55)
# both SSP and FICHE are forms of polygon IDs, "fiche" meaning the actual hard document used to classify the polygons. 
tab = read.dbf(file.path(validObj.dir, "DELPHINE/delphTABLE.dbf"))
rownames(tab) = tab$SSP # only SSP (polygon ID) and FICHE (polygon info/record file) will be used in this script 
length(unique(tab$SSP)); length(unique(polygonPix)); length(unique(polygonsID[])) # more polygons in tab.
head(tab)


# DELPHINE HABITATS WITHIN POLYGONS----------
# (for gpementsTable and representativity map) -

# 1 - Merging tables
# getting the fiche (polygonID2) information - CodeDelphine (percentage of habitats) 
gpements = read.table(file.path(validObj.dir,"DELPHINE/GroupementVegetaux.txt"),h=T) # Fiche, CodeDelphine, Surface
length(unique(gpements$Fiche)); length(unique(tab$FICHE)) # less FICHE have associated CodeDelphine than the total no. of fiche in tab

# Building a table with the CodeDelphine and the SSP polygonID (which corresponds to our polygonID rasters)
# for the'fiche' ID common to both tables there are several SSP polygons
lapply(split(tab$SSP[which(tab$FICHE %in% gpements$Fiche)], tab$FICHE[which(tab$FICHE %in% gpements$Fiche)]), function(x) length(x))

# merging tables with reference to fiche (to which (possibly) corresponds the polygon characterisation)
gpements2 = merge(gpements,tab[,c("SSP","FICHE")],by.x="Fiche",by.y="FICHE",all.x=T,all.y=F) # Fiche, CodeDelphine, Surface, SSP. 
length(unique(gpements2$SSP)); length(unique(tab$SSP))  

# 2 - Selecting final habitats to keep
# 2a - identifying the habitats present in each polygon
# 2b - renaming habitat codes (keeping only 2 first characters)
# 2c - building a list on unwated habitat codes ("foireux"? Further details in tableau PFG_reluCD.xls: habitat to be erased in grey/white)
# 2d - only keeping the habitats which do not belong to the unwanted habitats, when their surface is >20% 
# AND
# 2e - only keeping the habitats which have >20% surface in total (avoids having many very small habitats in a polygon)

polyList = split(gpements2[,c("CodeDelphine","Surface")],as.character(gpements2$SSP)) # 2a

delphCode = lapply(polyList, function(x){
  newCode = unlist(lapply(x$CodeDelphine, substr, start=1, stop=2)) # 2b
  habCode = c("A","A1","A2","A3",
              "O4","N4","Q2","Q5","T1",
              "U1","U2","U3","U4","U5",
              "W","W1","W2","W3","W4",
              "X","X1","X2",
              "Z1","Z2","Z3","Z4","Z5","Z6","Z7","Z9") # 2c
  if(sum(x[newCode %in% habCode,"Surface"]>20)) { # 2d
    return(NA)
  } else { return(newCode[which(x$Surface>20)]) } # 2e
})

delphCode = lapply(delphCode, unique)
delphCode = delphCode[-which(is.na(delphCode))] 

#===============================================
# 3. Comparisons 
#===============================================

# SELECTED PIXELS FOR DISTRIBUTION validation (based on 26 delphine codes (compatibTable))

# Keep pixels in mask : maskSimul[]==1
# Keep pixels in mask which have polygons : maskCons[]==1
# Keep pixels not in humid or anthropic areas : !delphCODE[]%in%c(0,1,24,20,21,23,26)

# Analyzed pixels:
analyzedPix = which(maskCons[] == 1 & (!delphCODE[]%in%c(0,1,24,20,21,23,26))) 

#Analyzed Polygons:
sum(!unique(tab$SSP) %in% unique(polygonsID[])) # 1045 SSP which are not in polygonsID
analyzedSSP = intersect(tab$SSP, as.character(polygonsID[])) # 19959 over 18915 polygonsID 
analyzedSSP = intersect(analyzedSSP, names(delphCode)) # 8835 SSP comon to 9244 delphCode

length(unique(polygonsID[analyzedPix])) # 8760 different polygons

#Analyzed pixels for which we have polygon and habitat (delphcode) info
analyzedPix = intersect(which(as.character(polygonsID[]) %in% analyzedSSP),analyzedPix) # 124045 for soil; 124064 drought init base
analyzedRow = which(which(maskSimul[]==1) %in% analyzedPix) # 124045 for soil; 124064 drought init base (all the analyzedPix are in the simulation mask)

save(analyzedPix,file=paste(validObj.dir,"SAVED_analyzedPix",sep=""))

# ------------------------------------------------
# SELECTED PIXELS FOR STRUCTURE and DISTRIBUTION (based on 84 delphine codes (gpementsTable)) 
selPix = maskCons
sum(!unique(gpements2$SSP) %in% unique(polygonsID[]))  # 527 SSP which are not in polygonsID
selPix[which(as.character(polygonsID[]) %in% unique(gpements2$SSP))] = 100 # Yellow, pixels for which we have information both about (Fiche,CodeDelphine,Surface) and SSP (info based on 84codes, and thus on vegetation groups - analyzed for structure)
selPix[analyzedPix] = 200 # Green, pixels above (info based on 26codes - analyzed for distribution)
table(selPix[])

pdf(paste(results.dir,"SAMPLING_analyzedPixels_Str_Dist.pdf",sep=""))
image(selPix,col=c("beige","grey50","black"),asp=1,ylab="",xlab="",yaxt="n",xaxt="n",bty="n")
legend("bottomleft", legend=c("Simulation pixels", "Analyzed pixels for structure", "Analyzed pixels for distribution"), bty="n", pch=22, pt.bg=c("beige","grey50","black"))
axis(1, line=-2,at=c(905000, 915000), labels=c("0", "10km"))
dev.off()


#--------------------------------------------------
# CALCULATING OBSERVED PRESENCES AND ABSENCES

# Converting PFG habitat presence/absence info. to polygon obs. presences absences:
# For each PFG, habitats present in the polygon are firstly identified and then compared to the habitat x PFG
# correspondance table (gpementsTable), using a conservative approach
# A polygon will be attibuted:
# 0 (absence) if the PFG is absent in all it's habitats
# 1 (ptentially present) if the PFG is 1 for at least one of the habitats in the polygon
# 2 (surely present) if the PFG is 2 for all habitats in the polygon

# NotF: for validation purposes, only 0s and 2s (absences and sure presences) are used

# Based on 84 delphi codes (STRUCTURE)
PAobs_list = list()
for(pfg in colnames(gpementsTable)[-1]){      # For each PFG
  
  compatib = gpementsTable[,pfg]              # P/A of the PFG per habitat
  names(compatib) = rownames(gpementsTable)   # attributing habitat codes as names
  
  PAobs_list[[pfg]] = unlist(lapply(delphCode[as.character(analyzedSSP)], function(x){ # For each SSP, which delphCode are corresponding
    xx = compatib[x]                                        
    resu = 1                                                
    if(sum(xx==0, na.rm=T) == length(na.omit(xx))) resu = 0 
    if(sum(xx==2, na.rm=T) == length(na.omit(xx))) resu = 2 
    return(resu)
  }))
}

save(PAobs_list,file=paste(validObj.dir,"SAVED_PAobs_list",sep=""))



