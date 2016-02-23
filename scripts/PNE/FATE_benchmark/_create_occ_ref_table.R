## reevaluate fatehd outputs
# .libPaths("/nfs_scratch2/emabio/R_PKG_LUKE/")
library(gridExtra)
library(raster)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rasterVis)


rm(list = ls())

sim.dir <- "/nfs_scratch2/dgeorges/FATE_newHS_devel/"
setwd(file.path(sim.dir, "workdir"))


path_input_data <- file.path(sim.dir, "data")
path_to_mask <- file.path(sim.dir, "SIMUL_6STRATA/DATA/MASK/maskEcrins.asc")
  
projETRS89 <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
projLCC <- "+proj=lcc +lat_1=45.89891888888889 +lat_2=47.69601444444444 +lat_0=46.8 +lon_0=2.337229166666667 +x_0=600000 +y_0=2200000 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

listeGrp <- get(load(file.path(path_input_data, "determinantes")))
pfg.list <- names(listeGrp) 
  
maskSimul = raster(path_to_mask, crs = CRS(projETRS89))
source("/nfs_scratch2/dgeorges/PolygonFromExtent.R")
poly.simul <- PolygonFromExtent(extent(maskSimul), crs = CRS(projETRS89), asSpatial = TRUE) 
  
  
## laod PFGS observations over the ecrins
occ.obj.list <- occ.ras.list <- vector(mode = "list", length = length(pfg.list))
names(occ.obj.list) <- names(occ.ras.list) <- pfg.list
  
cat("\n> getting pfg occurences..")
for(pfg_ in pfg.list){
    cat("\t", pfg_)
    (load(file.path(path_input_data, paste0("data.", pfg_, ".aust.newSetOfVar.austExt.RData"))))
    pts.in.simul <- pft.xyETRS89 %over% poly.simul
    pts.in.simul <- !is.na(pts.in.simul)
    occ.obj.list[[pfg_]] <- SpatialPointsDataFrame(pft.xyETRS89[pts.in.simul], data.frame(pft.occ = pft.occ[pts.in.simul]))
    occ.ras.list[[pfg_]] <- maskSimul
    occ.ras.list[[pfg_]][] <- NA ## set all to NA 
    occ.ras.list[[pfg_]][cellFromXY(occ.ras.list[[pfg_]], occ.obj.list[[pfg_]][occ.obj.list[[pfg_]]$pft.occ == 1, ])] <- 1
    occ.ras.list[[pfg_]][cellFromXY(occ.ras.list[[pfg_]], occ.obj.list[[pfg_]][occ.obj.list[[pfg_]]$pft.occ == 0, ])] <- 0
}

## create the observation matrix in a arrayPFG shape
occ.ras.tab <- sapply(occ.ras.list, function(r) r[])
head(occ.ras.tab)
save(occ.ras.tab, file = "occ.ras.tab.RData")

# occ.ras.stk <- stack(occ.ras.list)
# levelplot(occ.ras.stk)
