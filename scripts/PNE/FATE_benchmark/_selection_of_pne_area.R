## define an subarea of PNE to try some parameter settings

setwd("/nfs_scratch2/dgeorges/FATE_newHS/workdir/")
rm(list = ls())

library(raster)
library(rasterVis)
library(rgdal)
library(GSIF)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

projETRS89 <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"

ref.obj.path <- "/nfs_scratch2/dgeorges/FATE_newHS/SIMUL_6STRATA/outputsTables/arrayPFG_year850_rep4"
load(ref.obj.path)

ref.mask <- raster("/nfs_scratch2/dgeorges/FATE_newHS/SIMUL_6STRATA/DATA/MASK/maskEcrins.asc", crs = CRS(projETRS89))


dim(arrayPFG)
ncell(ref.mask)

arrayPFG.bin <- arrayPFG > 0
arrayPFG.sum <- apply(arrayPFG.bin, 1, sum)

arrayPFG.sum.ras <- ref.mask
arrayPFG.sum.ras[] <- arrayPFG.sum
arrayPFG.sum.ras <- arrayPFG.sum.ras * ref.mask

## load couple of masks we want to check
dist.ras <- raster("/nfs_scratch2/dgeorges/FATE_newHS/SIMUL_6STRATA/DATA/MASK/deforestation.asc", crs = CRS(projETRS89))
alti.ras <- raster("/nfs_scratch2/dgeorges/FATE_newHS/SIMUL_6STRATA/DATA/MASK/slope.asc", crs = CRS(projETRS89))


###############################
PolygonFromExtent <-
  function(ext, asSpatial=T, crs=CRS(NA), id=1)
  {
    x1 <- ext@xmin
    x2 <- ext@xmax
    y1 <- ext@ymin
    y2<-ext@ymax
    
    coords <- matrix(c(x1, y1, 
                       x1, y2,
                       x2, y2,
                       x2, y1,
                       x1, y1), ncol=2, byrow=T)
    
    poly <- Polygon(coords)
    if(asSpatial)
    {
      spPoly <- SpatialPolygons(list(Polygons(list(poly), ID=id)), proj4string=crs)
      return(spPoly)
      
    }
    return(poly)
    
  }
##############################

## define our extent of 20000 cells
obj <- GDALinfo("/nfs_scratch2/dgeorges/FATE_newHS/SIMUL_6STRATA/DATA/MASK/maskEcrins.asc")
ras.lst <- getSpatialTiles(obj, block.x=10000)

ext.list <- poly.list <-  bp.alti.list <- bp.nb.pfg.list <- nb.pix.list <- vector(mode = "list", nrow(ras.lst))

for(id in 1:nrow(ras.lst)){
  cat("\t", id)
  ext.tmp <- extent(ras.lst[id,"xl"],ras.lst[id,"xu"],ras.lst[id,"yl"],ras.lst[id,"yu"])
  ext.list[[id]] <- ext.tmp
  poly.list[[id]] <- PolygonFromExtent(ext.tmp, crs = CRS(projETRS89), asSpatial = FALSE)
  bp.alti.list[[id]] <- boxplot(crop(alti.ras, ext.tmp), plot = FALSE)$stat
  bp.nb.pfg.list[[id]] <- boxplot(crop(arrayPFG.sum.ras, ext.tmp), plot = FALSE)$stat
  nb.pix.list[[id]] <- sum(!is.na(crop(arrayPFG.sum.ras, ext.tmp)[]))
}
 
bp.alti.df <- do.call(cbind, bp.alti.list) %>% t %>% as.data.frame 
bp.alti.df$ext.id <- 1:nrow(bp.alti.df)
colnames(bp.alti.df) <- c("q005", "q025", "q050", "q075", "q095", "ext.id")
bp.alti.df <- bp.alti.df %>% gather(quant.name, quant.val, c(q005, q025, q050, q075, q095))

bp.nb.pfg.df <- do.call(cbind, bp.nb.pfg.list) %>% t %>% as.data.frame 
bp.nb.pfg.df$ext.id <- 1:nrow(bp.nb.pfg.df)
colnames(bp.nb.pfg.df) <- c("q005", "q025", "q050", "q075", "q095", "ext.id")
bp.nb.pfg.df <- bp.nb.pfg.df %>% gather(quant.name, quant.val, c(q005, q025, q050, q075, q095))

nb.pix.df <- do.call(cbind, nb.pix.list) %>% t %>% as.data.frame 
nb.pix.df$ext.id <- 1:nrow(nb.pix.df)
colnames(nb.pix.df) <- c("nb.pix", "ext.id")

grid.layer <- layer(sp.polygons(SpatialPolygons(list(Polygons(poly.list, ID = 1)), proj4string = CRS(projETRS89))))

lp.alti <- levelplot(alti.ras, margin = FALSE, main = "alti", par.settings = rasterVis::BuRdTheme) + grid.layer
lp.nb.pfg <- levelplot(arrayPFG.sum.ras, margin = FALSE, main = "nb PFG", par.settings = rasterVis::BuRdTheme) + grid.layer
lp.nb.pix <- levelplot(ref.mask, margin = FALSE, main = "nb pix", par.settings = rasterVis::BuRdTheme) + grid.layer
  
gg.alti <- ggplot(data = bp.alti.df, aes(x = factor(ext.id), y = quant.val)) + geom_boxplot()
gg.nb.pfg <- ggplot(data = bp.nb.pfg.df, aes(x = factor(ext.id), y = quant.val)) + geom_boxplot()
gg.nb.pix <- ggplot(data = nb.pix.df, aes(x = factor(ext.id), y = nb.pix)) + geom_point()

png('PNE_sub_area_selection.png', width = 650, height = 720)
print(grid.arrange(lp.nb.pfg, gg.nb.pfg,
                   lp.alti, gg.alti,
                   lp.nb.pix, gg.nb.pix, ncol = 2))
dev.off()

##' The area 20/27/36 seems to be quite optimal to conduct our tests

selected.area <- c(20,27,43)
ras.files <- list.files("/nfs_scratch2/dgeorges/FATE_newHS/SIMUL_6STRATA/DATA", ".asc$", full.names = TRUE, recursive = TRUE)
## remove useless files
ras.files <- grep("_archiveMaya|_FORMAL_HS|_FORMAL_HS_BIN|_SCALED_HS_0.75|maskDemo.asc", ras.files, value = TRUE, invert = TRUE)

for(i in 1:length(selected.area)){
  cat('\n\n> Zone', i)
  ext.pne.reduce <- ext.list[[selected.area[i]]]
  save(ext.pne.reduce, file = paste0("ext.pne.reduce.zone", i, ".RData"))
  
  ##' create the new fate subarea
  for(ras.files_ in ras.files){
    cat("\n", ras.files_)
    r <- raster(ras.files_, crs = CRS(projETRS89))
    r <- try(crop(r, ext.pne.reduce))
    if(!inherits(r, 'try-error')){
      fn_ <- file.path(dirname(sub("/SIMUL_6STRATA/", "/SIMUL_6STRATA_TEST/", ras.files_)),
                       paste0("zone", i), basename(ras.files_))
      dir.create(dirname(fn_), showWarnings = FALSE, recursive = TRUE)
      writeRaster(r, filename = fn_, overwrite = TRUE)
    }
  }
}
