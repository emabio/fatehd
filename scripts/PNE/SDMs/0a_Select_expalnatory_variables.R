### HEADER #####################################################################
##' @title Sslect new explanatory variables for PNE PFG modelling
##'
##' @author Damien G. & Maya G.
##' @contact damien.georges2 at gmail.com
##' 
##' @date 25/01/2016
##' 
##' @description Here we will conduct a little analysis to select a set of explanatory
##'   variables for FATEHD PFG SDMs
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

## enable caching 
knitr::opts_chunk$set(cache = TRUE, fig.width = 13, fig.height = 7)

## load needed packages --------------------------------------------------------
library(dplyr)
library(raster)
library(rgdal)
library(ggplot2)
library(ade4)



## Constants definition --------------------------------------------------------
user = "ftp" ## the id of user/machine which do the analyses
sce = "AUST" ## the vlimatic environmental variable source ("NICK" or "AUST")
# ext = "ALL_ECRINS" ## the extent of the study ("ALL_ECRINS" or "ALL_ALPS")
env.var.names <- gtools::mixedsort(c(paste0("bio_", 1:19), "spring", "winter", "carbon", "slope", "rock", "twi")) ## variables we are interested in

##+ ggplot, fig.width = 13, fig.height = 7
for(ext in c("ALL_ECRINS", "ALL_ALPS")){
  ## Print brief summary of the script args --------------------------------------
  cat("\nVariable selection procedure for SDMs -------------------------------")
  cat("\nstart at:", format(Sys.time(), "%a %d %b %Y %X"))
  cat("\n user:", user)
  cat("\n env.var.names:", env.var.names)
  cat("\n ext:", ext)
  cat("\n-----------------------------------------------------------------------")
  
  ## Paths to data definition ----------------------------------------------------
  if (user == "ftp"){
    path_input <- "/media/ftp-public/GUEGUEN_Maya/"
    path_output <- "/media/ftp-public/GUEGUEN_Maya/VAR_SELECTION"
  } else stop("Unsupported 'user' value")
  
  dir.create(path_output, showWarnings = FALSE, recursive = TRUE)
  
  ## load the environemetal variables stack --------------------------------------
  stk.cur <- stack(gtools::mixedsort(list.files(file.path(path_input, paste0("ENV_", ext), "EOBS_1970_2005"), ".img$", full.names = TRUE)))
  
  ##  compute correlation table btw variables
  stk.cur.extr <- extract(stk.cur, which(!is.na(subset(stk.cur, 'carbon')[])))
  expl.var.cor <- cor(na.omit(stk.cur.extr), method = "spearman")
  
  ## plot correlations
  dat_ <- data.frame(var1 = rownames(expl.var.cor), expl.var.cor, stringsAsFactors = FALSE) %>% tidyr::gather(var2, spear.cor, -1) %>% 
    mutate(col = spear.cor > 0, var1 = factor(var1, levels = env.var.names),  var2 = factor(var2, levels = env.var.names)) %>% 
    na.omit %>% rowwise() %>% mutate(tt = paste(gtools::mixedsort(c(var1, var2)), collapse = "-")) # %>% distinct(tt) 
  dat_$spear.cor[dat_$spear.cor > 0.7] <- NA
  
  print(ggplot(dat_, aes(x = var1, y = var2, size = 1 - abs(spear.cor), col = spear.cor, label = round(spear.cor, 2))) + geom_text() + scale_colour_continuous(limits=c(-1, 1), low = 'blue', high = 'red') +  guides(colour = FALSE, size = FALSE))

  
  
  ## use of pca to choose variables
  z <- dudi.pca(na.omit(stk.cur.extr), center = T, scale = T, scannf = F)
  
  inertie <- z$eig / sum(z$eig) * 100
  barplot( inertie ,ylab = "% d'inertie", names.arg = round(inertie,2))
  title("Eboulis des valeurs propres en %")
  
  print(inertia.dudi(z, col.inertia = T)$col.abs)
  
  s.label(z$li, clab = 0, xax = 1, yax = 2, sub = "axe1 vs axe2")
  # ajout des coordonn ́es des variables, m^me plan
  s.arrow(3 * z$co, xax = 1, yax = 2, clab = 0.8, add.plot = TRUE)
  s.corcircle(z$co)
}

##' ## Check future variations of our climatic variables in the future

## only done on PNE
ext = "ALL_ECRINS"

## load the environemetal variables stack --------------------------------------
stk.cur <- stack(gtools::mixedsort(list.files(file.path(path_input, paste0("ENV_", ext), "EOBS_1970_2005"), ".img$", full.names = TRUE)))


stat.out <- NULL

##+ , eval = FALSE
for(sc_ in list.files(file.path(path_input, paste0("ENV_", ext)), 'rcp[0-9]{2}_')){
  for(yr_ in seq(2020, 2090, 10)){
    cat("\n stk.name:", paste0("stk.", sc_, ".", yr_))
    stk.tmp <- stack(gtools::mixedsort(list.files(file.path(path_input, paste0("ENV_", ext), sc_), 
                                              paste0(yr_, ".img$"), full.names = TRUE)))
    stat.mean.tmp <- data.frame(area = ext, year = yr_, scenario = sc_, stat = 'mean', data.frame(stat_value = cellStats(stk.tmp, 'mean')), stringsAsFactors = FALSE)
    stat.sd.tmp <- data.frame(area = ext, year = yr_, scenario = sc_, stat = 'sd', data.frame(stat_value = cellStats(stk.tmp, 'sd')), stringsAsFactors = FALSE)
    stat.tmp <- bind_rows(stat.mean.tmp, stat.sd.tmp)
    stat.tmp$variable <- sub("_20[0-9]0$", "", rownames(stat.sd.tmp))
    stat.out <- bind_rows(stat.out, stat.tmp)
  }
}
write.table(stat.out, file = "~/Work/FATEHD/workdir/var_fut_trends.txt", col.names = TRUE, row.names = FALSE, sep = "\t")

##+ , eval = TRUE 
stat.out <- read.table("~/Work/FATEHD/workdir/var_fut_trends.txt", h = TRUE, stringsAsFactors = FALSE)
stat.out <- stat.out %>% tidyr::spread(stat, stat_value)
stat.out$variable <- factor(stat.out$variable, levels = gtools::mixedsort(unique(stat.out$variable)))


##+ , fig.width = 13, fig.height = 7
gg <- ggplot(data = stat.out, aes(x = year, y = mean, colour = scenario)) +
  geom_line() + facet_wrap(~variable, scales = 'free_y')
# gg <- ggplot(data = stat.out, aes(x = year, y = mean, ymin = mean - sd, ymax = mean + sd, colour = scenario)) +
#   geom_line() + geom_linerange(aes(alpha = .3)) + facet_wrap(~variable, scales = 'free_y')
print(gg)
