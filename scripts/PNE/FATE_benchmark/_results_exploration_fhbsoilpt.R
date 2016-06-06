
##' ---
##' title : FATEHD parameter exploration results explorration
##' author : damien g.
##' ---

##' # Descrtiption -------------------------------------------------------------
##' 
##' In this script we will setup some visual tools to explore what combination
##' of parameters should lead to the best model preformances.
##' A grid campain have been run on cigri to see on 3 subarea of PNE the influence
##' of :
##' 
##'   - envsuit.option 
##'   - seeding.params
##'   - dispers.mode
##'   - global.abund
##'   - global.resource.thresh
##'   - max.by.cohort
##'   

##' # Get simul results --------------------------------------------------------
##' 
##' Simul outputs are stored on each cigri cluters, we will merge all results on
##' luke in the directory "~/fhdmpt_simul_comparaisons_merged/"
##' 

#+ knitr option, echo = FALSE 
library(knitr)
knitr::opts_chunk$set(fig.width=12, fig.height=8, fig.path='Figs_fhdpfgpt/',
                      echo=TRUE, warning=FALSE, message=FALSE,
                      cache  = TRUE)

#+ start script, 
rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)

work.dir <- "~/Work/FATEHD/benchmarking/fatehd_params_test/grid/hs_params_test/fhdhspt/workdir/"
input.dir <- "~/Work/FATEHD/benchmarking/fatehd_params_test/grid/hs_params_test/fhdhspt/simul_comparaisons/"
output.dir <- "~/Work/FATEHD/benchmarking/fatehd_params_test/grid/hs_params_test/fhdhspt/figures/"

dir.create(work.dir, showWarnings = FALSE, recursive = TRUE)
dir.create(input.dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)

setwd(work.dir)

##' execute the scp command to retrieave results from all cigri clusters
files.to.load <- list.files(input.dir, full.names = TRUE, recursive = TRUE)
year_to_keep <- paste0("year", seq(100, 600, 50))

## create the summary table of jobs
res.occ.list <- lapply(grep("eval.occ", files.to.load, value = TRUE), function(fid_){
  res1_ <- read.table(fid_, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  ## keep only a subset of the data
  res_ <- res1_ %>% filter(is.element(year, year_to_keep))
  return(res_)
})

eval.occ.df <- bind_rows(res.occ.list)

save(eval.occ.df, file = file.path(work.dir, "eval.occ.df_fhdhspt.RData"))
# (load(file.path(work.dir, "eval.occ.df_fhdhspt.RData")))

## create the summary table of jobs
res.abu.list <- lapply(grep("eval.abund", files.to.load, value = TRUE), function(fid_){
  res1_ <- read.table(fid_, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  ## keep only a subset of the data
  res_ <- res1_ %>% filter(is.element(year, year_to_keep))
  return(res_)
})

eval.abu.df <- bind_rows(res.abu.list)
save(eval.abu.df, file = file.path(work.dir, "eval.abu.df_fhdhspt.RData"))
# (load(file.path(work.dir, "eval.abu.df_fhdhspt.RData")))

##' # Explore campain results

##' ## Which parfameters are important for global tss?



gg.dat <- eval.abu.df %>% mutate(tss = sens + spec -1, simul.id = sub("^.*hs_", "", simul.id)) %>%
  filter(year == "year600") 

gg <- ggplot(gg.dat, aes(y = tss, fill = simul.id, x = pfg)) + 
  geom_bar(stat = "identity", position = "dodge") + facet_grid(habitat ~ .) +
  geom_hline(data = gg.dat %>% group_by(simul.id, habitat) %>% summarise(mean.tss = mean(tss, na.rm = TRUE)),
             aes(yintercept = mean.tss, colour = simul.id), lty = 2)

gg

gg.dat <- eval.occ.df %>% mutate(tss = sens + spec -1, simul.id = sub("^.*hs_", "", simul.id)) %>%
  filter(year == "year600") 
gg.dat <- eval.abu.df %>% mutate(tss = sens + spec -1, simul.id = sub("^.*hs_", "", simul.id)) %>%
  filter(year == "year600") 

gg <- ggplot(gg.dat, aes(y = auc, fill = simul.id, x = pfg)) + 
  geom_bar(stat = "identity", position = "dodge") + facet_grid(habitat ~ .) +
  geom_hline(data = gg.dat %>% group_by(simul.id, habitat) %>% summarise(mean.auc = mean(auc, na.rm = TRUE)),
             aes(yintercept = mean.auc, colour = simul.id), lty = 2) + coord_cartesian(ylim = c(0.5,1))

gg

g



library(rasterVis)
stk.formal <- stack(list.files("~/Work/FATEHD/benchmarking/fatehd_params_test/grid/hs_params_test/fhdhspt/SIMUL_6STRATA_TEST/DATA/PFGS/ENVSUIT/_SCALED_SENS_95/formal/", ".asc", full.names = TRUE))
x11()
levelplot(stk.formal)

stk.rock.cont <- stack(list.files("~/Work/FATEHD/benchmarking/fatehd_params_test/grid/hs_params_test/fhdhspt/SIMUL_6STRATA_TEST/DATA/PFGS/ENVSUIT/_SCALED_SENS_95/mask_rock_cont/", ".asc", full.names = TRUE))
x11()
levelplot(stk.rock.cont)

stk.rock.cont.slope.bin <- stack(list.files("~/Work/FATEHD/benchmarking/fatehd_params_test/grid/hs_params_test/fhdhspt/SIMUL_6STRATA_TEST/DATA/PFGS/ENVSUIT/_SCALED_SENS_95/mask_rock_cont_slope_bin/", ".asc", full.names = TRUE))
x11()
levelplot(stk.rock.cont.slope.bin)

stk.slope.bin <- stack(list.files("~/Work/FATEHD/benchmarking/fatehd_params_test/grid/hs_params_test/fhdhspt/SIMUL_6STRATA_TEST/DATA/PFGS/ENVSUIT/_SCALED_SENS_95/mask_slope_bin/", ".asc", full.names = TRUE))
x11()
levelplot(stk.slope.bin)





gg.dat <- eval.abu.df %>% mutate(tss = sens + spec -1, simul.id = sub("^.*hs_", "", simul.id)) %>%
  filter(year == "year600") 

gg <- ggplot(gg.dat, aes(y = cor.spear.abund, fill = simul.id, x = pfg)) + 
  geom_bar(stat = "identity", position = "dodge") + facet_grid(habitat ~ .) +
  geom_hline(data = gg.dat %>% group_by(simul.id, habitat) %>% summarise(mean.cor.spear.abund = mean(cor.spear.abund, na.rm = TRUE)),
             aes(yintercept = mean.cor.spear.abund, colour = simul.id), lty = 2) + coord_cartesian(ylim = c(0,1))

gg.dat <- eval.abu.df %>% mutate(hab.prop.rel = hab.prop.fhd - hab.prop.obs , simul.id = sub("^.*hs_", "", simul.id)) %>%
  filter(year == "year600") 

gg <- ggplot(gg.dat, aes(y = hab.prop.rel, fill = simul.id, x = pfg)) + 
  geom_bar(stat = "identity", position = "dodge") + facet_grid(habitat ~ .) +
  geom_hline(data = gg.dat %>% group_by(simul.id, habitat) %>% summarise(mean.hab.prop.rel = mean(hab.prop.rel, na.rm = TRUE)),
             aes(yintercept = mean.hab.prop.rel, colour = simul.id), lty = 2)
gg
