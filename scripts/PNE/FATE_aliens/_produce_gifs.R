##' ----------------------------------------------------------------------------
##' @title create PNE alien invasion gifs
##' @description we will produce here some animated images of the invation of 
##'   PNE by some alien species
##' @date 2016-02-23
##' @author damien g
##' @licence GPL-2
##' ----------------------------------------------------------------------------

rm(list = ls())

## load libraries --------------------------------------------------------------
library(raster)
library(rasterVis)
library(animation)
library(dichromat)

## laod simul  mask ------------------------------------------------------------
mask.ras <- raster("INTRODUCE_AFG/DATA/MASK/maskEcrins.asc")  

## define species and simul params and load abund arrays -----------------------
species.name <- "H2a_alien"
array.output.pattern <- "repH2a_alien_PNE"
array.output.dir <- "INTRODUCE_AFG/outputsTables/"
output.dir <- "figs"
dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)

array.output.files <- list.files(array.output.dir, paste0("arrayPFG_.*", array.output.pattern), full.names = TRUE)
array.output.files <- gtools::mixedsort(array.output.files)

sp.ras.list <- lapply(array.output.files, function(ar_){
  cat("\n>", ar_)
  load(ar_)
  ras.tmp <- mask.ras
  ras.tmp[] <- arrayPFG[, species.name]
  ras.tmp <- ras.tmp * mask.ras 
  return(ras.tmp)
}) 
sp.stk <- do.call(raster::stack, sp.ras.list)
names(sp.ras.list) <- names(sp.stk) <- sub("_.*$", "", sub("^.*_year", "t", array.output.files))

## define graphical custom params ----------------------------------------------

## transform the 0-10000 scale into a 0-100 one
my.at <- seq(0,10000, 2000)
myColorkey <- list(at = seq(0,10000, 100), ## where the colors change
                   labels = list(
                     at = my.at, ## where to print labels
                     labels = my.at/100
                   ))
## define a easy to read color theme
myTheme <- rasterTheme(region=rev(dichromat(terrain.colors(15))))

## produce the levelplot -------------------------------------------------------
png(file.path(output.dir, paste0(species.name, "_PNE.png")), width = 960)
print(levelplot(sp.stk, at = my.at, colorkey = myColorkey, par.settings = myTheme, main = species.name))
dev.off()

## produce the gif -------------------------------------------------------------
saveGIF({
  for(t_ in names(sp.ras.list)){
    print(levelplot(sp.ras.list[[t_]], at = my.at, colorkey = myColorkey, main = list(t_), par.settings = myTheme))
  }
}, movie.name = file.path(output.dir, paste0(species.name, "_PNE.gif")), interval = 0.8, nmax = 50, ani.width = 600,
ani.height = 600)

