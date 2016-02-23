##' ----------------------------------------------------------------------------
##' @title check aliens invasion simulation parameters
##' @description the aim of this script is to understand why simul of introduction
##' of aliens species are returning stange results.
##' @date 2016-02-23
##' @author damien g
##' @licence GPL-2
##' ----------------------------------------------------------------------------

rm(list = ls())
work.dir <- "~/Work/FATEHD/workdir/TEST_ALIEN/INTRODUCE_AFG/workdir"
# work.dir <- "~/Work/FATEHD/workdir/TEST_ALIEN/INTRODUCE_NATs/workdir"
dir.create(work.dir, showWarnings = FALSE, recursive = TRUE)
setwd(work.dir)

##' 1. check that all masks are aligned and defines on the right coordinate system

library(raster)
library(dplyr)

## define the system of projections
projLCC <- "+proj=lcc +lat_1=45.89891888888889 +lat_2=47.69601444444444 +lat_0=46.8 +lon_0=2.337229166666667 +x_0=600000 +y_0=2200000 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

map.files <- list.files("../DATA_all_ecrins/", ".asc$", recursive = TRUE, full.names = TRUE)
map.ras.list <- lapply(map.files, function(m_) raster(m_, crs = CRS(projLCC)))

## get all maps attributes
map.ras.attr <- lapply(map.ras.list, function(r_){
  ext_ <- extent(r_)
  data.frame(file = filename(r_),
    xmin = ext_@xmin, xmax = ext_@xmax, ymin = ext_@ymin, ymax = ext_@ymax,
    ncell = ncell(r_), nNA = sum(is.na(r_[])), stringsAsFactors = FALSE)
})

map.ras.attr <- do.call(rbind, map.ras.attr)
## check if they are all equals
summary(map.ras.attr)

## remove useless masks
map.ras.attr <- map.ras.attr %>% dplyr::filter(!grepl("maskDemo.asc$|maskEcrins0.asc$|/demo/", file))
for(i in 2:ncol(map.ras.attr)){
  cat("\n>", colnames(map.ras.attr)[i], ":", paste(unique(map.ras.attr[,i]), collapse = ", "))
}

##' @note: potential issue: The alien HS maps don't have the same number of NA's than other pfgs.


##' 2. define a reduce area to run simulations

hs.ref.file <- grep("H2a_alien.asc$", map.ras.attr$file, value = TRUE)
hs.ref.ras <- raster(hs.ref.file)
plot(hs.ref.ras)

##' first extent
# ext.H2a_alien <- drawExtent()
# > ext.H2a_alien
# class       : Extent 
# xmin        : 880016.5 
# xmax        : 896589.4 
# ymin        : 2004261 
# ymax        : 2020636 
# save(ext.H2a_alien, file = "ext.H2a_alien.RData")

##' second extent
# ext.H2a_alien <- drawExtent()
# > ext.H2a_alien
# class       : Extent 
# xmin        : 886330 
# xmax        : 891459.7 
# ymin        : 2011758 
# ymax        : 2015507 
# save(ext.H2a_alien, file = "ext.H2a_alien_2.RData")

# load("~/Work/FATEHD/workdir/TEST_ALIEN/INTRODUCE_AFG/workdir/ext.H2a_alien.RData")

##' full extent
ext.H2a_alien <- raster::extent(hs.ref.ras)
plot(ext.H2a_alien, add = TRUE)

## construct the min common area mask
ref.mask <- raster::stack(as.character(map.ras.attr$file))
ref.mask <- crop(ref.mask, ext.H2a_alien)
ref.mask <- calc(ref.mask, fun = sum, na.rm = FALSE)
ref.mask <- reclassify(ref.mask, c(-Inf, Inf, 1))
plot(ref.mask)

## clip all rasters
pdf("all_masks.pdf")
fn_ <- map.ras.attr$file[1]
new.ras.files.list <- lapply(map.ras.attr$file, function(fn_){
  cat("\n>", fn_)
  r_ <- raster(fn_)
  r_ <- crop(r_ , ref.mask) * ref.mask
  plot(r_, main = basename(fn_))
  fn__ <- sub("DATA_all_ecrins", "DATA", fn_)
  writeRaster(r_, filename = fn__, overwrite = TRUE, NAflag = -3.4e+38)
  ## remove the ".000000000000000" pattern if exists
  f_ <- readLines(fn__)
  f_ <- sub(".000000000000000", "", f_)
  writeLines(text = f_, con = fn__)
  return(fn__)
})
dev.off()

## check that all raster attributes are matching
map.ras.list <- lapply(unlist(new.ras.files.list), function(m_) raster(m_))
map.ras.attr <- do.call(rbind, map.ras.attr)
## check if they are all equals
summary(map.ras.attr)

##' @note Should be ok with this version!!


