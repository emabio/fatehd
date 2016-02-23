### HEADER #####################################################################
##' @title Check PFG SDM quality
##'
##' @author Damien G. & Maya G.
##' @contact damien.georges2 at gmail.com
##' 
##' @date 29/01/2016
##' 
##' @description Here we will compare SDMs of our 23 PFG SDMs.
##'  - We will compare Isa's Paper SDM's to the new version where climatic 
##'    variables have been changed 
##'   
##' @note 
##'   - We made some tests considering quantile of determinante species to define
##'     PFG SDMs, but because it not seems to be better in therme of performances
##'     (TSS comparation and FATEHD simul outputs) and because it is much more
##'     computive intensive to produce SDM maps, we decided not to keep working
##'     with this alternative.
##' 
##' @log 
##' 
##' @licencing GPL
##'     Copyright (C) 2015  Damien G.
##' 
##'     This program is free software: you can redistribute it and/or modify
##'     it under the terms of the GNU General Public License as published by
##'     the Free Software Foundation, either version 3 of the License, or
##'     (at your option) any later version.
##' 
##'     This program is distributed in the hope that it will be useful,
##'     but WITHOUT ANY WARRANTY; without even the implied warranty of
##'     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##'     GNU General Public License for more details.
##' 
##'     You should have received a copy of the GNU General Public License
##'     along with this program.  If not, see <http://www.gnu.org/licenses/>.
## END OF HEADER ###############################################################

rm(list = ls())

library(raster)
library(rasterVis)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

## evaluate PFG models
dat.path <- "/nfs_scratch2/emabio/FATEHD/_PFG_VERSION/_OUTPUT_DATA_NEW_ENV/DATA_AUST/PFT_occ"
# mod.path.pfg <- "/nfs_scratch2/emabio/FATEHD/_PFG_VERSION/_OUTPUT_DATA_NEW_ENV/DATA_AUST"
mod.path.pfg <- "/nfs_scratch2/emabio/FATEHD/_PFG_VERSION/_OUTPUT_DATA_NEW_ENV/DATA_AUST_VAR_ISA_ALPS/"
mod.path.isa <- "/nfs_scratch2/emabio/FATEHD/_ISA_VERSION"

pfg.list <- list.files(mod.path.pfg, "[C,H,P][0-9]{1,2}")

# pfg_ <- "H5"
for(pfg_ in pfg.list){
  cat("\n\n>", pfg_)
  
  cat("\n\tloading pfg based models")
  ## the SDM based on PFG occurence
  pfg.mod.pred.ras <- raster(list.files(file.path(mod.path.pfg, pfg_, "proj_ParcEcrins_current", "individual_projections"), "_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData.img$", full.names = TRUE))
    
  cat("\n\tloading isa's models")
  ## the formal SDM
  pfg.mod.pred.formal <- raster(list.files(mod.path.isa, paste0("_", pfg_, "_"), full.names = TRUE))
  projection(pfg.mod.pred.formal) <- CRS("+proj=lcc +lat_1=45.89891888888889 +lat_2=47.69601444444444 +lat_0=46.8 +lon_0=2.337229166666667 +x_0=600000 +y_0=2200000 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
  pfg.mod.pred.formal <- projectRaster(pfg.mod.pred.formal, pfg.mod.pred.ras)
  
  ## create the stack of predictions
  pred.stk <- stack(pfg.mod.pred.formal, pfg.mod.pred.ras / 1000)
  names(pred.stk) <- c("isa sdm", "new sdm")
  
  cat("\n\tcreate level plot")
  lp <- levelplot(pred.stk, layout = c(2,1), main = pfg_, colorkey = list(space = "bottom")) 
  
  pfg.dat <- get(load(list.files(dat.path, paste0("OCC_", pfg_, "$"), full.names = TRUE)))
  
  load(file.path(dirname(mod.path.pfg),"PLOTS_Ecrins"))
  coord_XY <- input[,c("PlotID_KEY","X_ETRS89","Y_ETRS89")]
  
  pfg.mod.pred.df <- raster::extract(pred.stk, coord_XY[, c("X_ETRS89","Y_ETRS89")])
  
  to.keep <- unique(which(!is.na(pfg.mod.pred.df), arr.ind = TRUE)[, 1])
  
  cat("\n\tevaluating models")
  eval.test.df <- NULL
  for(mod.name in colnames(pfg.mod.pred.df)){
    eval.test.list <- lapply(seq(0, 1, 0.02), function(x){
      biomod2::Find.Optim.Stat(Stat='TSS',
                               pfg.mod.pred.df[to.keep, mod.name],
                               pfg.dat[to.keep],
                               Fixed.thresh = x)
    })
    eval.test.df <- bind_rows(eval.test.df, data.frame(do.call(rbind, eval.test.list), mod.name = mod.name, row.names = NULL, stringsAsFactors = FALSE))
  }
  
  
  ## scale variables
  eval.test.df$cutoff <- eval.test.df$cutoff
  eval.test.df$sensitivity <- eval.test.df$sensitivity / 100
  eval.test.df$specificity <- eval.test.df$specificity / 100
  
  eval.test.df <- eval.test.df %>% gather(val.name, val, c(best.stat, sensitivity, specificity))
  
  cat("\n\tproduce evaluation graph")
  gg <- ggplot(data = eval.test.df, aes(x = cutoff, y = val, colour = val.name)) + geom_line() + facet_grid(~mod.name) + coord_cartesian(ylim = c(0,1)) + theme(legend.position="bottom")
  
  fig.out.dir <- "/nfs_scratch2/emabio/FATEHD/workdir/pfg_envsuit_test_graphs"
  dir.create(fig.out.dir, showWarnings = FALSE, recursive = TRUE)
  
  cat("\n\tproduce final graph")
  png(file.path(fig.out.dir, paste0(pfg_, "_envsuit_test.png")), width = 800, height = 600)
  print(grid.arrange(lp, gg, ncol = 1))
  dev.off()
}



####
plot(pfg.mod.pred.formal)
points(coord_XY[intersect(coord_XY$PlotID_KEY, names(pfg.dat)[pfg.dat == 1]), c("X_ETRS89","Y_ETRS89")], col = "red", pch = 1)
points(coord_XY[intersect(coord_XY$PlotID_KEY, names(pfg.dat)[pfg.dat == 0]), c("X_ETRS89","Y_ETRS89")], col = "blue", pch = 1)


plot(pfg.mod.pred.ras)
points(coord_XY[intersect(coord_XY$PlotID_KEY, names(pfg.dat)[pfg.dat == 1]), c("X_ETRS89","Y_ETRS89")], col = "red", pch = 1)
points(coord_XY[intersect(coord_XY$PlotID_KEY, names(pfg.dat)[pfg.dat == 0]), c("X_ETRS89","Y_ETRS89")], col = "blue", pch = 1)


points(coord_XY[intersect(coord_XY$PlotID_KEY, names(pfg.dat)[pfg.dat == 0]), c("X_ETRS89","Y_ETRS89")], col = "black", pch = 1)
