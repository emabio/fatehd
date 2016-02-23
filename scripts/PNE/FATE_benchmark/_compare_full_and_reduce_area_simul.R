## compare full and reduce version of fatehd simulaitons ##
rm(list = ls())
setwd("/nfs_scratch2/dgeorges/FATE_newHS_devel/workdir/")

library(raster)
library(rasterVis)


projETRS89 <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"

rep.id <- 4
year.id <- 850

## load full simul attributes
path.to.full.sim <- "/nfs_scratch2/dgeorges/FATE_newHS/SIMUL_6STRATA/outputsTables/"
path.to.mask.full.sim <- "/nfs_scratch2/dgeorges/FATE_newHS/SIMUL_6STRATA/DATA/MASK/maskEcrins.asc"
mask.full.sim <- raster(path.to.mask.full.sim, crs = CRS(projETRS89))

## load reduced simule attributes
path.to.red.sim <- "/nfs_scratch2/dgeorges/FATE_newHS_devel/SIMUL_6STRATA/outputsTables/" 
path.to.mask.red.sim <- "/nfs_scratch2/dgeorges/FATE_newHS_devel/SIMUL_6STRATA/DATA/MASK/maskEcrins.asc"
mask.red.sim <- raster(path.to.mask.red.sim, crs = CRS(projETRS89))
red.ext <- extent(mask.red.sim)

## load simul outputs
sim.full.df <- get(load(paste0(path.to.full.sim, "arrayPFG_year", year.id, "_rep", rep.id )))
sim.red.df <- get(load(paste0(path.to.red.sim, "arrayPFG_year", year.id, "_rep", rep.id )))

pfg_ <- "H1"

pdf("compare_full_and_reduced_area_simulations.pdf")

for(pfg_ in colnames(sim.full.df)){
  cat("\n>", pfg_)
  
  ## deal with full area simul
  sim.full.ras <- mask.full.sim
  sim.full.ras[] <- sim.full.df[, pfg_]
  sim.full.ras <- crop(sim.full.ras, red.ext)
  
  ## deal with reduced area simul
  sim.red.ras <- mask.red.sim
  sim.red.ras[] <- sim.red.df[, pfg_]
  
  ## produce plots
  sim.stk <- stack(sim.full.ras, sim.red.ras) 
  names(sim.stk) <- c("full_area", "reduced_area")
  print(levelplot(sim.stk, main = pfg_))
}
dev.off()



