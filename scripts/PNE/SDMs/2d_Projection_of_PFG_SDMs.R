## HEADER #####################################################################
##' @title Project SDMs in the future
##'
##' @author Damien G. & Maya G.
##' @contact damien.georges2 at gmail.com
##' 
##' @date 16/01/2016
##' 
##' @description Here we will do the projection of our PFGs SDMs
##'   
##' @note 
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


## Script initialisation -------------------------------------------------------
rm(list=ls())

## Retireve input args ---------------------------------------------------------
args <- commandArgs(trailingOnly=TRUE)
sp.name <- as.character(args[1]) ## give the name of the species
# sp.name <- "P8"


## Constants definition --------------------------------------------------------
user = "ftp" ## the id of user/machine which do the analyses
sce = "AUST" ## the vlimatic environmental variable source ("NICK" or "AUST")
version.name = "_PFG_VERSION" ## the type of model unit ("_SP_VERSION" or "_PFG_VERSION")
env.var.names <- c("bio_6", "bio_9", "bio_12", "bio_15", "carbon", "slope") ## variables we are interested in

## Print brief summary of the script args --------------------------------------
cat("\nBuild PFG determinant species SDM -------------------------------------")
cat("\nstart at:", format(Sys.time(), "%a %d %b %Y %X"))
cat("\n species:", sp.name)
cat("\n user:", user)
cat("\n version.name:", version.name)
cat("\n env.var.names:", env.var.names)
cat("\n-----------------------------------------------------------------------")

## Paths to data definition ----------------------------------------------------
if (user == "ftp"){
  path_input <- path_output <- paste0("/media/ftp-public/GUEGUEN_Maya/", version.name , "/_OUTPUT_DATA/")
} else if (user == "luke"){
  path_input <- path_output <- paste0("/nfs_scratch2/emabio/FATEHD/", version.name , "/_OUTPUT_DATA_NEW_ENV/")
  .libPaths('/nfs_scratch2/emabio/R_PKG_LUKE') ## here are the shared library installed on luke
} else stop("Unsupported 'user' value")

## load needed packages --------------------------------------------------------
library(biomod2)

## create the output directory
path_sce <- paste(path_output,"DATA_",sce,"/",sep="")
dir.create(path_sce, showWarnings = FALSE, recursive = TRUE)
setwd(path_sce)

## load ensemble models object
cat("\n loading previous version of bm.em..")
bm.em.file <- list.files(paste0(path_input,"DATA_",sce,"/", sp.name), pattern = "ensemble.models.out$", full.names = TRUE)
bm.em <- get(load(bm.em.file))

## 


