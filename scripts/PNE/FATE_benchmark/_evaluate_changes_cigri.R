### compute the graph of sens, spe change in function of time
rm(list = ls())

## Retrieve input args ---------------------------------------------------------
args <- commandArgs(trailingOnly=TRUE)
path.to.dat <- as.character(args[1])
file_name <- as.character(args[2]) ## give the paramSimul file name
file_name_ref = as.character(args[3])
nbCores <- as.numeric(args[4]) ## give the number of resources to create outputs tables
maps <- as.logical(args[5]) ## do you want to print pdf maps ?
lib.dir <- as.character(args[6])

## Load libraries --------------------------------------------------------------
if(length(lib.dir)){
  .libPaths(lib.dir)
}

library(gridExtra)
library(raster)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rasterVis)


# path.to.dat <- "/home/dgeorges/fhdmpt" 
# file_name <- "/home/dgeorges/fhdmpt/SIMUL_6STRATA_TEST/PARAM_SIMUL/ParamsSimul_1.txt"
# file_name_ref <- "/home/dgeorges/fhdmpt/SIMUL_6STRATA_REF/PARAM_SIMUL/ParamsSimul_ref1.txt"
# nbCores <- 1 
# maps <- FALSE 
# lib.dir <- "/home/dgeorges/R_PKG_LUKE"

## print simul parameters ------------------------------------------------------
cat("\nInput parameters --------")
cat("\n- path.to.dat <-", path.to.dat)
cat("\n- file_name <-", file_name)
cat("\n- file_name_ref <-", file_name_ref)
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

## define the id of the simul as the saving directory name
simul_id = basename(simul.res.dir)

path.input.data <- file.path(path.to.dat, "Data")
pattern.input.data <- ".aust.newSetOfVar.austExt.RData"

path.output <- file.path("simul_comparaisons", simul_id)
path.output.tab <- file.path(path.output, "table")
path.output.graph <- file.path(path.output, "graph")
dir.create(path.output.tab, showWarnings = FALSE, recursive = TRUE)
dir.create(path.output.graph, showWarnings = FALSE, recursive = TRUE)

## define some projection system
projETRS89 <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
projLCC <- "+proj=lcc +lat_1=45.89891888888889 +lat_2=47.69601444444444 +lat_0=46.8 +lon_0=2.337229166666667 +x_0=600000 +y_0=2200000 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

abund.thresh <- 0 ## abundance use to convert FATEHD abundances into PA

ref.simul.dir <- dirname(dirname(file_name_ref))

ref.simul.file <- file(file_name_ref, "r")
ref.simul.file.pl <- readLines(ref.simul.file, warn = FALSE)
ref.rep <- basename(get.char.param(ref.simul.file.pl, "--SAVE_DIR--"))


ref.mask.path <- get.char.param(simul.param.list, "--MASK--")

ref.array.PFG.files <- list.files(file.path(ref.simul.dir, "outputsTables"), paste0("arrayPFG.*_rep", ref.rep), full.names = TRUE)
ref.times <- sort(as.numeric(sub("^.*_year", "", sub("_rep.*$", "", ref.array.PFG.files))))


test.simul <- file.path("outputsTables", simul_id)
test.array.PFG.files <- list.files(test.simul, paste0("arrayPFG"), full.names = TRUE)
test.times <- sort(as.numeric(sub("^.*_year", "", test.array.PFG.files)))

inter.times <- intersect(ref.times, test.times)

## test that objects dimentions are compatible
ref.array.test <- get(load(ref.array.PFG.files[1]))
test.array.test <- get(load(test.array.PFG.files[1]))
if(!all(dim(ref.array.test) == dim(test.array.test))) stop("ref and test arrayPFG have incompatible sizes")
pfg.list <- colnames(ref.array.test)


ref.array.PFG <- test.array.PFG <- array(NA, dim = c(dim(ref.array.test), length(inter.times)),
                                         dimnames = c(dimnames(ref.array.test), list(paste0("year", inter.times))))  

cat("\n> loading array PFG:")
for(yr_ in dimnames(ref.array.PFG)[[3]]){
  cat("\t", yr_)
  ref.array.PFG[, , yr_] <- as.matrix(get(load(grep(paste0("_", yr_, "_"), ref.array.PFG.files, value = TRUE))))
  test.array.PFG[, , yr_] <- as.matrix(get(load(grep(paste0("_", yr_, "$"), test.array.PFG.files, value = TRUE))))
}

## transform abundances into presences absences
ref.array.PFG.bin <- ref.array.PFG > abund.thresh
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

## calculate some quality indicies
cat("\n> evaluating reference simulation...")
eval.occ.ref <- do.call(rbind, lapply(pfg.list, function(pfg_){
  do.call(rbind, lapply(dimnames(ref.array.PFG.bin)[[3]], function(yr_){
    nb.occ <- sum(occ.array[, pfg_], na.rm = TRUE)
    nb.abs <- sum(occ.array[, pfg_] == 0, na.rm = TRUE)
    nb.pred.occ <- sum(ref.array.PFG.bin[, pfg_, yr_], na.rm = TRUE)
    nb.pred.abs <- sum(!ref.array.PFG.bin[, pfg_, yr_], na.rm = TRUE)
    cont.tab <- table(ref.array.PFG.bin[, pfg_, yr_], occ.array[, pfg_])
    sens = spec = ppv = npv = fpr = fnr = fdr = acc = f1 = NA
    if(all(dim(cont.tab) == c(2,2))){
      TP = cont.tab['TRUE', '1']
      TN = cont.tab['FALSE', '0']
      FP = cont.tab['TRUE', '0']
      FN = cont.tab['FALSE', '1']
      
      sens = TP / (TP + FN)
      spec = TN / (TN + FP)
      ppv =  TP / (TP + FP)
      npv = TN / (TN + FN)
      fpr = 1 - spec
      fnr = 1 - sens
      fdr = 1 - ppv
      acc = (TP + TN) / (TP + FP + FN + TN)
      f1 = 2*TP / (2*TP + FP + FN)    
    }
    return(data.frame(pfg = as.character(pfg_), year = as.character(yr_), nb.occ = nb.occ, nb.abs = nb.abs, nb.pred.occ = nb.pred.occ, nb.pred.abs = nb.pred.abs,
             sens = sens, spec = spec, ppv = ppv, npv = npv,
             fpr = fpr, fnr = fnr, fdr = fdr, acc = acc, f1 = f1))
  }))
}))

cat("\n> evaluating test simulation...")
eval.occ.test <- do.call(rbind, lapply(pfg.list, function(pfg_){
  do.call(rbind, lapply(dimnames(test.array.PFG.bin)[[3]], function(yr_){
    nb.occ <- sum(occ.array[, pfg_], na.rm = TRUE)
    nb.abs <- sum(occ.array[, pfg_] == 0, na.rm = TRUE)
    nb.pred.occ <- sum(test.array.PFG.bin[, pfg_, yr_], na.rm = TRUE)
    nb.pred.abs <- sum(!test.array.PFG.bin[, pfg_, yr_], na.rm = TRUE)
    cont.tab <- table(test.array.PFG.bin[, pfg_, yr_], occ.array[, pfg_])
    sens = spec = ppv = npv = fpr = fnr = fdr = acc = f1 = NA
    if(all(dim(cont.tab) == c(2,2))){
      TP = cont.tab['TRUE', '1']
      TN = cont.tab['FALSE', '0']
      FP = cont.tab['TRUE', '0']
      FN = cont.tab['FALSE', '1']
      
      sens = TP / (TP + FN)
      spec = TN / (TN + FP)
      ppv =  TP / (TP + FP)
      npv = TN / (TN + FN)
      fpr = 1 - spec
      fnr = 1 - sens
      fdr = 1 - ppv
      acc = (TP + TN) / (TP + FP + FN + TN)
      f1 = 2*TP / (2*TP + FP + FN)    
    }
    return(data.frame(pfg = as.character(pfg_), year = as.character(yr_), nb.occ = nb.occ, nb.abs = nb.abs, nb.pred.occ = nb.pred.occ, nb.pred.abs = nb.pred.abs,
             sens = sens, spec = spec, ppv = ppv, npv = npv,
             fpr = fpr, fnr = fnr, fdr = fdr, acc = acc, f1 = f1))
  }))
}))



##' comp.test.ref will contains:
##'   - the diff of each index
##'   - the correlation score btw both simulations

cat("\n> comparing ref and test simulations...")
comp.test.ref <- data.frame(eval.occ.ref[, 1:2], eval.occ.test[, -c(1:2)] - eval.occ.ref[, -c(1:2)])

cat("\n> summerizing differences btw simuls...")
comp.test.ref.summ <- comp.test.ref %>% group_by(year) %>% dplyr::summarize(nb.occ = mean(nb.occ, na.rm = TRUE),
                                                                     nb.abs = mean(nb.abs, na.rm = TRUE),
                                                                     nb.pred.occ = mean(nb.pred.occ, na.rm = TRUE),
                                                                     nb.pred.abs = mean(nb.pred.abs, na.rm = TRUE),
                                                                     sens = mean(sens, na.rm = TRUE),
                                                                     spec = mean(spec, na.rm = TRUE),
                                                                     ppv = mean(ppv, na.rm = TRUE),
                                                                     npv = mean(npv, na.rm = TRUE),
                                                                     fpr = mean(fpr, na.rm = TRUE),
                                                                     fnr = mean(fnr, na.rm = TRUE),
                                                                     fdr = mean(fdr, na.rm = TRUE),
                                                                     acc = mean(acc, na.rm = TRUE),
                                                                     f1 = mean(f1, na.rm = TRUE),
                                                                     pfg = 'all')

## produce summary graphs
comp.test.ref_ <- comp.test.ref %>% bind_rows(comp.test.ref.summ) %>% select(pfg, year, sens, spec, ppv, npv, acc) %>% 
  gather("metric.name", "metric.val", c(sens, spec, ppv, npv, acc)) %>%
  mutate(year = as.numeric(sub("year", "", year)))

# comp.test.ref_ <- comp.test.ref %>% bind_rows(comp.test.ref.summ) %>% select(pfg, year, sens) %>% 
#   gather("metric.name", "metric.val", c(sens)) %>%
#   mutate(year = as.numeric(sub("year", "", year)))

gg.diff <- ggplot(data = comp.test.ref_, aes(x = year, y = metric.val, colour = metric.name)) + 
  geom_line() + facet_wrap(~pfg) + coord_cartesian(ylim = c(-1,1)) + ggtitle(paste0("simul ", simul_id)) + 
  geom_vline(xintercept = c(300, 600, 800), lty = 2, colour = 'grey') + ylab("") + xlab("") + 
  theme(legend.title=element_blank())

## save graph
pdf(file.path(path.output, paste0("comp_ref_test_PA_", simul_id, ".pdf")), width = 9)
print(gg.diff)
dev.off()

## save tables
write.table(eval.occ.test, file = file.path(path.output, paste0("eval.occ.test.", simul_id, ".txt")), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(comp.test.ref, file = file.path(path.output, paste0("comp.test.ref.", simul_id, ".txt")), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(comp.test.ref.summ, file = file.path(path.output, paste0("comp.test.ref.summ.", simul_id, ".txt")), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

## copy the table on the home because irods is dead!!!
path.output.home <- "~/fhdmpt_simul_comparaisons/"
dir.create(path.output.home, showWarnings = FALSE, recursive = TRUE)
file.copy(path.output, path.output.home, recursive = TRUE, overwrite = TRUE)


q("no")