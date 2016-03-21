### compute the graph of sens, spe change in function of time
rm(list = ls())

## Retrieve input args ---------------------------------------------------------
args <- commandArgs(trailingOnly=TRUE)
path.to.dat <- as.character(args[1])
file_name <- as.character(args[2]) ## give the paramSimul file name
nbCores <- as.numeric(args[3]) ## give the number of resources to create outputs tables
maps <- as.logical(args[4]) ## do you want to print pdf maps ?
lib.dir <- as.character(args[5])

## tests ------
# setwd("/var/tmp/dgeorges/fhdpfgpt/20")
# path.to.dat <- "/home/dgeorges/fhdpfgpt" 
# file_name <- "ParamsSimul_531.txt"
# nbCores <- 1 
# maps <- FALSE 
# # lib.dir <- "/home/dgeorges/R_PKG_LUKE"
## end tests -----

## Load libraries --------------------------------------------------------------
if(length(lib.dir)){
  .libPaths(lib.dir)
}

library(raster)
library(dplyr)
library(PresenceAbsence)
library(foreign)


## print simul parameters ------------------------------------------------------
cat("\nInput parameters --------")
cat("\n- path.to.dat <-", path.to.dat)
cat("\n- file_name <-", file_name)
cat("\n- nbCores <-", nbCores)
cat("\n- maps <-", maps)
cat("\n- lib.dir <-", lib.dir)
cat("\n-------------------------")

## Retrieve name of the simul replication --------------------------------------
cat("\n PARAM SIMUL FILE : ",file_name)

## get the simul results directory name
source(file.path(path.to.dat, "Functions", "_get_parameters_utils.R"))

simul.param.list <- readLines(file_name, warn = FALSE)
simul.res.dir <- get.char.param(pl = simul.param.list, flag = "--SAVE_DIR--")
global.param.file <- get.char.param(pl = simul.param.list, flag = "--GLOBAL_PARAMS--")
ref.mask.path <- get.char.param(simul.param.list, "--MASK--")

## define the id of the simul as the saving directory name
simul.id = basename(simul.res.dir)

path.input.data <- file.path(path.to.dat, "Data")
pattern.input.data <- ".aust.newSetOfVar.austExt.RData"

path.output <- file.path("simul_comparaisons", simul.id)
dir.create(path.output, showWarnings = FALSE, recursive = TRUE)
# path.output.tab <- file.path(path.output, "table")
# path.output.graph <- file.path(path.output, "graph")
# dir.create(path.output.tab, showWarnings = FALSE, recursive = TRUE)
# dir.create(path.output.graph, showWarnings = FALSE, recursive = TRUE)

## define some projection system
projETRS89 <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
projLCC <- "+proj=lcc +lat_1=45.89891888888889 +lat_2=47.69601444444444 +lat_0=46.8 +lon_0=2.337229166666667 +x_0=600000 +y_0=2200000 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

abund.thresh <- 0.0001 ## abundance use to convert FATEHD abundances into PA

test.simul <- file.path("outputsTables", simul.id)
test.array.PFG.files <- list.files(test.simul, paste0("arrayPFG"), full.names = TRUE)
test.times <- sort(as.numeric(sub("^.*_year", "", test.array.PFG.files)))

inter.times <- test.times

## test that objects dimentions are compatible
test.array.test <- get(load(test.array.PFG.files[1]))
pfg.list <- colnames(test.array.test)


test.array.PFG <- array(NA, dim = c(dim(test.array.test), length(inter.times)),
                        dimnames = c(dimnames(test.array.test), list(paste0("year", inter.times))))  

cat("\n> loading array PFG:")
for(yr_ in dimnames(test.array.PFG)[[3]]){
  cat("\t", yr_)
  test.array.PFG[, , yr_] <- as.matrix(get(load(grep(paste0("_", yr_, "$"), test.array.PFG.files, value = TRUE))))
}

## transform abundances into presences absences
test.array.PFG.bin <- test.array.PFG > abund.thresh

## load the releve observation data
source(file.path(path.to.dat, "Functions", "PolygonFromExtent.R"))
mask.ref <- raster(ref.mask.path, crs = CRS(projETRS89))
poly.simul <- PolygonFromExtent(extent(mask.ref), crs = CRS(projETRS89), asSpatial = TRUE) 


## laod PFGS observations over the ecrins
occ.obj.list <- occ.ras.list <- vector(mode = "list", length = length(pfg.list))
names(occ.obj.list) <- names(occ.ras.list) <- pfg.list

cat("\n> getting pfg occurences..")
for(pfg_ in pfg.list){
  cat("\t", pfg_)
  (load(file.path(path.input.data, paste0("data.", sub("_.*$", "", pfg_), pattern.input.data))))
  pts.in.simul <- pft.xyETRS89 %over% poly.simul
  pts.in.simul <- !is.na(pts.in.simul)
  occ.obj.list[[pfg_]] <- SpatialPointsDataFrame(pft.xyETRS89[pts.in.simul], data.frame(pft.occ = pft.occ[pts.in.simul]))
  occ.ras.list[[pfg_]] <- mask.ref
  occ.ras.list[[pfg_]][] <- NA ## set all to NA 
  occ.ras.list[[pfg_]][cellFromXY(occ.ras.list[[pfg_]], occ.obj.list[[pfg_]][occ.obj.list[[pfg_]]$pft.occ == 0, ])] <- 0
  occ.ras.list[[pfg_]][cellFromXY(occ.ras.list[[pfg_]], occ.obj.list[[pfg_]][occ.obj.list[[pfg_]]$pft.occ == 1, ])] <- 1
}

## create the observation matrix in a arrayPFG shape
occ.array <- sapply(occ.ras.list, function(r) r[])

## get habitats specificity
cat("\n> getting habitat specifications..")
habitats <- raster(file.path(path.input.data, "Objects", "ABUND_DATA", "RECODED_categories_physio.asc"), crs = CRS(projLCC))
habitats.vect <- raster::extract(habitats, xyFromCell(mask.ref, 1:ncell(mask.ref), spatial = TRUE))


## calculate some quality indicies

## test
# yr_ <- "year100"
# pfg_ <- pfg.list[1]
## end test

cat("\n> evaluating test simulation...")
eval.occ.test <- do.call(rbind, lapply(pfg.list, function(pfg_){
  do.call(rbind, lapply(dimnames(test.array.PFG)[[3]], function(yr_){
    ## deal with the full area
    dat_ <- data.frame(id = 1:nrow(occ.array), 
                       obs = occ.array[, pfg_], 
                       pred = as.numeric(test.array.PFG[, pfg_, yr_]),
                       hab = habitats.vect)
    cont.tab <- cmx(na.omit(dat_), threshold = abund.thresh, na.rm = TRUE)
    nb.occ <- sum(occ.array[, pfg_], na.rm = TRUE)
    nb.abs <- sum(occ.array[, pfg_] == 0, na.rm = TRUE)
    nb.pred.occ <- sum(test.array.PFG.bin[, pfg_, yr_], na.rm = TRUE)
    nb.pred.abs <- sum(!test.array.PFG.bin[, pfg_, yr_], na.rm = TRUE)
    
    df.out <- data.frame(
      simul.id = simul.id,
      dat.source = "releves",
      habitat = "all",
      pfg = as.character(pfg_),
      year = as.character(yr_),
      nb.occ = nb.occ, nb.abs = nb.abs, 
      nb.pred.occ = nb.pred.occ, nb.pred.abs = nb.pred.abs,
      sens = sensitivity(cont.tab, st.dev=FALSE), 
      spec = specificity(cont.tab, st.dev=FALSE), 
      pcc = pcc(cont.tab, st.dev = FALSE),
      kappa = Kappa(cont.tab, st.dev=FALSE),
      auc = auc(na.omit(dat_), st.dev=FALSE),
      stringsAsFactors = FALSE)
    
    ## deal with each habitat
    hab.list <- na.omit(unique(dat_$hab))
    for(hab_ in hab.list){
      sub_ <- dat_$hab == hab_
      if(length(sub_)){
        dat__ <- na.omit(dat_[sub_, , drop = FALSE])
        if(nrow(dat__)){
          cont.tab <- cmx(dat__, threshold = abund.thresh, na.rm = TRUE)
          nb.occ <- sum(occ.array[sub_, pfg_], na.rm = TRUE)
          nb.abs <- sum(occ.array[sub_, pfg_] == 0, na.rm = TRUE)
          nb.pred.occ <- sum(test.array.PFG.bin[sub_, pfg_, yr_], na.rm = TRUE)
          nb.pred.abs <- sum(!test.array.PFG.bin[sub_, pfg_, yr_], na.rm = TRUE)
          
          df.out <- bind_rows(df.out, data.frame(
            simul.id = simul.id,
            dat.source = "releves",
            habitat = as.character(hab_),
            pfg = as.character(pfg_),
            year = as.character(yr_),
            nb.occ = nb.occ, nb.abs = nb.abs, 
            nb.pred.occ = nb.pred.occ, nb.pred.abs = nb.pred.abs,
            sens = sensitivity(cont.tab, st.dev=FALSE), 
            spec = specificity(cont.tab, st.dev=FALSE), 
            pcc = pcc(cont.tab, st.dev = FALSE),
            kappa = Kappa(cont.tab, st.dev=FALSE),
            auc = auc(dat__, st.dev=FALSE, na.rm = TRUE),
            stringsAsFactors = FALSE))        
        }
      }
    }
    return(df.out)
  }))
}))

## deal with abundances metrics

## Get observed PFG releves ----------------------------------------------------
plots <- read.dbf(file.path(path.input.data, "Objects", "ABUND_DATA", "plots.dbf")) ## coordinates
rownames(plots) <- plots$numchrono
PFGxPlots <- get(load(file.path(path.input.data, "Objects", "ABUND_DATA", "PFGxPlots"))) ## PFGxSITES matrix
plots <- plots[rownames(PFGxPlots),]

## create the relative PFG abund relative to other PFG
PFGxPlots.rel <- PFGxPlots %>% select(-Others) 
PFGxPlots.rowsums <- rowSums(PFGxPlots.rel)
PFGxPlots.rel <- PFGxPlots.rel / PFGxPlots.rowsums


## Transform coordinates into the right projection system
plots.xy.lcc <- SpatialPoints(as.data.frame(plots[,c("X_LII","Y_LII")]), proj4string=CRS(projLCC))
plots.xy.etrs89 <- spTransform(plots.xy.lcc, CRS(projETRS89))
plots.rows <- cellFromXY(mask.ref, plots.xy.etrs89)

test.array.PFG.plots <- test.array.PFG[plots.rows, , , drop = FALSE]
test.array.PFG.bin.plots <- test.array.PFG.bin[plots.rows, , , drop = FALSE]
rownames(test.array.PFG.plots) <- rownames(test.array.PFG.bin.plots) <- rownames(PFGxPlots)

## create relative abundance table
test.array.PFG.plots.rel <- test.array.PFG.plots
for(yr_ in dimnames(test.array.PFG)[[3]]){
  test.array.PFG.plots.rowsums <- rowSums(test.array.PFG.plots[, , yr_])
  test.array.PFG.plots.rel[, , yr_] <- test.array.PFG.plots[, , yr_] / test.array.PFG.plots.rowsums  
}



habitats.vect.plots <- habitats.vect[plots.rows] 



eval.abund.test <- do.call(rbind, lapply(pfg.list, function(pfg_){
  do.call(rbind, lapply(dimnames(test.array.PFG)[[3]], function(yr_){
    pfg.short_ <- sub("_.*$", "", pfg_)
    ## deal with the full area
    dat_ <- data.frame(id = rownames(PFGxPlots), 
                       obs = as.numeric(PFGxPlots[, pfg.short_] >= abund.thresh), 
                       pred = as.numeric(test.array.PFG.plots[, pfg_, yr_]),
                       hab = habitats.vect.plots)
    cont.tab <- cmx(na.omit(dat_), threshold = abund.thresh, na.rm = TRUE)
    nb.occ <- sum(PFGxPlots[, pfg.short_] >= abund.thresh, na.rm = TRUE)
    nb.abs <- sum(PFGxPlots[, pfg.short_] < abund.thresh, na.rm = TRUE)
    nb.pred.occ <- sum(test.array.PFG.bin.plots[, pfg_, yr_], na.rm = TRUE)
    nb.pred.abs <- sum(!test.array.PFG.bin.plots[, pfg_, yr_], na.rm = TRUE)
    
    df.out <- data.frame(
      simul.id = simul.id,
      dat.source = "phyto",
      habitat = "all",
      pfg = as.character(pfg_),
      year = as.character(yr_),
      nb.occ = nb.occ, nb.abs = nb.abs, 
      nb.pred.occ = nb.pred.occ, nb.pred.abs = nb.pred.abs,
      sens = sensitivity(cont.tab, st.dev=FALSE), 
      spec = specificity(cont.tab, st.dev=FALSE), 
      pcc = pcc(cont.tab, st.dev = FALSE),
      kappa = Kappa(cont.tab, st.dev=FALSE),
      auc = auc(na.omit(dat_), st.dev=FALSE),
      stringsAsFactors = FALSE)
    
    ## deal with each habitat
    hab.list <- na.omit(unique(dat_$hab))
    for(hab_ in hab.list){
      sub_ <- dat_$hab == hab_
      if(length(sub_)){
        dat__ <- na.omit(dat_[sub_, , drop = FALSE])
        if(nrow(dat__)){
          cont.tab <- cmx(dat__, threshold = abund.thresh, na.rm = TRUE)
          nb.occ <- sum(PFGxPlots[sub_, pfg.short_] >= abund.thresh, na.rm = TRUE)
          nb.abs <- sum(PFGxPlots[sub_, pfg.short_] < abund.thresh, na.rm = TRUE)
          nb.pred.occ <- sum(test.array.PFG.bin.plots[sub_, pfg_, yr_], na.rm = TRUE)
          nb.pred.abs <- sum(!test.array.PFG.bin.plots[sub_, pfg_, yr_], na.rm = TRUE)
          
          df.out <- bind_rows(df.out, data.frame(
            simul.id = simul.id,
            dat.source = "phyto",
            habitat = as.character(hab_),
            pfg = as.character(pfg_),
            year = as.character(yr_),
            nb.occ = nb.occ, nb.abs = nb.abs, 
            nb.pred.occ = nb.pred.occ, nb.pred.abs = nb.pred.abs,
            sens = sensitivity(cont.tab, st.dev=FALSE), 
            spec = specificity(cont.tab, st.dev=FALSE), 
            pcc = pcc(cont.tab, st.dev = FALSE),
            kappa = Kappa(cont.tab, st.dev=FALSE),
            auc = auc(dat__, st.dev=FALSE, na.rm = TRUE),
            stringsAsFactors = FALSE))        
        }
      }
    }
    return(df.out)
  }))
}))


eval.abund.test.bis <- do.call(rbind, lapply(pfg.list, function(pfg_){
  do.call(rbind, lapply(dimnames(test.array.PFG)[[3]], function(yr_){
    pfg.short_ <- sub("_.*$", "", pfg_)
    ## deal with the full area
    dat_ <- na.omit(data.frame(id = rownames(PFGxPlots), 
                       obs = PFGxPlots[, pfg.short_], 
                       pred = test.array.PFG.plots[, pfg_, yr_],
                       hab = habitats.vect.plots,
                       stringsAsFactors = FALSE)) 
    
    dat.rel_ <- na.omit(data.frame(id = rownames(PFGxPlots.rel), 
                                  obs = PFGxPlots.rel[, pfg.short_], 
                                  pred = test.array.PFG.plots.rel[, pfg_, yr_],
                                  hab = habitats.vect.plots,
                                  stringsAsFactors = FALSE))
    df.out <- data.frame(
      simul.id = simul.id,
      dat.source = "phyto",
      habitat = "all",
      pfg = as.character(pfg_),
      year = as.character(yr_),
      cor.spear.abund = cor(dat_$obs, dat_$pred, method = "spearman"), ## spearman correlation
      mean.abund.obs = mean(dat_$obs, na.rm = TRUE), ## mean abundance in obs data
      mean.abund.fhd = mean(dat_$pred, na.rm = TRUE),
      median.abund.obs = median(dat_$obs, na.rm = TRUE),
      median.abund.fhd = median(dat_$pred, na.rm = TRUE),
      mean.cov.obs = mean(dat.rel_$obs, na.rm = TRUE),
      median.cov.obs = median(dat.rel_$obs, na.rm = TRUE),
      mean.cov.fhd = mean(dat.rel_$pred, na.rm = TRUE),
      median.cov.fhd = median(dat.rel_$pred, na.rm = TRUE),
      hab.prop.obs = 100,
      hab.prop.fhd = 100,
      stringsAsFactors = FALSE)
    
    ## deal with each habitat
    hab.list <- na.omit(unique(dat_$hab))
    for(hab_ in hab.list){
      sub_ <- dat_$hab == hab_
      if(length(sub_)){
        dat__ <- na.omit(dat_[sub_, , drop = FALSE])
        dat.rel__ <- na.omit(dat.rel_[sub_, , drop = FALSE])
        if(nrow(dat__)){
          df.out <- bind_rows(df.out, data.frame(
            simul.id = simul.id,
            dat.source = "phyto",
            habitat = as.character(hab_),
            pfg = as.character(pfg_),
            year = as.character(yr_),
            cor.spear.abund = cor(dat__$obs, dat__$pred, method = "spearman"), ## spearman correlation
            mean.abund.obs = mean(dat__$obs, na.rm = TRUE), ## mean abundance in obs data
            mean.abund.fhd = mean(dat__$pred, na.rm = TRUE),
            median.abund.obs = median(dat__$obs, na.rm = TRUE),
            median.abund.fhd = median(dat__$pred, na.rm = TRUE),
            mean.cov.obs = mean(dat.rel__$obs, na.rm = TRUE),
            median.cov.obs = median(dat.rel__$obs, na.rm = TRUE),
            mean.cov.fhd = mean(dat.rel__$pred, na.rm = TRUE),
            median.cov.fhd = median(dat.rel__$pred, na.rm = TRUE),
            hab.prop.obs = sum(dat__$obs) / sum(dat_$obs) * 100,
            hab.prop.fhd = sum(dat__$pred) / sum(dat_$pred) * 100,
            stringsAsFactors = FALSE))       
        }
      }
    }
    return(df.out)
  }))
}))

eval.abund.test.ter <- full_join(eval.abund.test, eval.abund.test.bis)


## save the output tables
cat("\n> Saving tables...")

write.table(eval.occ.test, file = file.path(path.output, paste0("eval.occ.test.", simul.id, ".txt")), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(eval.abund.test.ter, file = file.path(path.output, paste0("eval.abund.test.", simul.id, ".txt")), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

q("no")