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
##' Simul outputs are stored on each cigri cluters, we will merge All results on
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
library(parallel)

work.dir <- "~/Work/FATEHD/benchmarking/workdir/"
input.dir <- "~/Work/FATEHD/benchmarking/results/"
output.dir <- "~/Work/FATEHD/benchmarking/figures/"

dir.create(work.dir, showWarnings = FALSE, recursive = TRUE)
dir.create(input.dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)

setwd(work.dir)

##' execute the scp command to retrieave results from All cigri clusters

## get input parameters file
params <- read.csv(file.path(input.dir, "uncoupledLHS.csv"), row.names = 1, stringsAsFactors = FALSE)
params$simul.id <- rownames(params)
head(params)

## check jobs that have been completed
params$job.status <- file.exists(file.path(input.dir, "fhdbsoilpt", "5671", "9918", params$simul.id)) 
sum(params$job.status)
length(params$job.status)

## reshape params
params <- params %>% mutate(hs.type = factor(ceiling(hs.type)),
                            area = factor(ceiling(area)))

year_to_keep <- paste0("year", c(830,840,850))

# ## create the summary table of jobs
# res.occ.list <- mclapply(params$simul.id[params$job.status], function(jid_){
#   res1_ <- read.table(file.path(input.dir, "fhdbsoilpt", "5671", "9918", jid_, "simul_comparaisons", paste0("531_", jid_), paste0("eval.occ.test.531_", jid_,".txt")), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#   ## keep only a subset of the data
#   res_ <- res1_ %>% filter(is.element(year, year_to_keep)) %>%  mutate(simul.id = jid_)
#   return(res_)
# }, mc.cores = 16)
# 
# eval.occ.df <- bind_rows(res.occ.list)
# ## merge with parameters
# eval.occ.df <- eval.occ.df %>% left_join(params)
# ## define correctly the habitat 
# eval.occ.df <- eval.occ.df %>% mutate(habitat = factor(habitat, levels = c("0", "31", "40", "50", "60", "70", "81", "83", "All"),
#                                                        labels = c("Excluded","Rock", "Grasslands", "Lowlands","Open_habitats", "Semi-closed_habitats", "Closed_habitats","Forests", "All")))
# 
# save(eval.occ.df, file = file.path(work.dir, "eval.occ.df_fhdbsoilpt.RData"))
(load(file.path(work.dir, "eval.occ.df_fhdbsoilpt.RData")))

# ## create the summary table of jobs
# res.abu.list <- mclapply(params$simul.id[params$job.status], function(jid_){
#   res1_ <- read.table(file.path(input.dir, "fhdbsoilpt", "5671", "9918", jid_, "simul_comparaisons", paste0("531_", jid_), paste0("eval.abund.test.531_", jid_,".txt")), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#   ## keep only a subset of the data
#   res_ <- res1_ %>% filter(is.element(year, year_to_keep)) %>%  mutate(simul.id = jid_)
#   return(res_)
# }, mc.cores = 16)
# 
# eval.abu.df <- bind_rows(res.abu.list)
# ## merge with parameters
# eval.abu.df <- eval.abu.df %>% left_join(params)
# eval.abu.df <- eval.abu.df %>% mutate(habitat = factor(habitat, levels = c("0", "31", "40", "50", "60", "70", "81", "83", "All"),
#                                                        labels = c("Excluded","Rock", "Grasslands", "Lowlands","Open_habitats", "Semi-closed_habitats", "Closed_habitats","Forests", "All")))
# 
# save(eval.abu.df, file = file.path(work.dir, "eval.abu.df_fhdbsoilpt.RData"))
(load(file.path(work.dir, "eval.abu.df_fhdbsoilpt.RData")))

##' # Explore campain results

##' ## Which parfameters are important for global tss?

dat.name <- c('eval.occ.df', 'eval.abu.df')
# library(lme4)

## check the pfg cohexistance
eval.occ.df %>% filter(habitat == 'All', year == "year850") %>% group_by(simul.id, dat.source) %>% summarise(nb.pfg = sum(nb.pred.occ > 0))

##' only the simul 4407, 8841, 90181, 9421 maintains the 24 PFGs
best.cohexitence.simul <- c(4407, 8841, 9081, 9421)
eval.occ.df %>% filter(habitat == 'All', year == "year850", is.element(simul.id, best.cohexitence.simul)) %>% 
  group_by(simul.id, dat.source) %>% summarise(tss.mean = mean(spec + sens -1))

gg.dat <- eval.occ.df %>% filter(habitat == 'All', year == "year850", is.element(simul.id, best.cohexitence.simul)) %>% 
  group_by(simul.id, dat.source, pfg) %>% summarise(tss = spec + sens -1)

gg <- ggplot(data = gg.dat, aes(x = pfg, y = tss, fill = tss)) + geom_bar(stat = 'identity') + facet_grid(~simul.id) + coord_flip()
gg




## check the pfg cohexistance
eval.abu.df %>% filter(habitat == 'All', year == "year850", is.element(simul.id, best.cohexitence.simul)) %>% 
  group_by(simul.id, dat.source) %>% summarise(tss.mean = mean(spec + sens -1))

gg.dat <- eval.abu.df %>% filter(habitat == 'All', year == "year850", is.element(simul.id, best.cohexitence.simul)) %>% 
  group_by(simul.id, dat.source, pfg) %>% summarise(tss = spec + sens -1)

gg <- ggplot(data = gg.dat, aes(x = pfg, y = tss, fill = tss)) + geom_bar(stat = 'identity') + facet_grid(~simul.id) + coord_flip()
gg


eval.abu.df %>% filter(habitat == 'All', year == "year850", is.element(simul.id, best.cohexitence.simul)) %>% 
  group_by(simul.id, dat.source) %>% summarise(cor.mean = mean(cor.spear.abund))

gg.dat <- eval.abu.df %>% filter(habitat == 'All', year == "year850", is.element(simul.id, best.cohexitence.simul)) %>% 
  group_by(simul.id, dat.source, pfg) %>% summarise(cor = cor.spear.abund, tss = spec + sens -1)

gg <- ggplot(data = gg.dat, aes(x = pfg, y = tss, fill = cor)) + geom_bar(stat = 'identity') + facet_grid(~simul.id) + coord_flip()
gg

## lets keep only the simul 9421 an do the same analysis by environment
gg.dat <- eval.abu.df %>% filter(year == "year850", simul.id == 9421) %>% 
  group_by(simul.id, dat.source, pfg, habitat) %>% summarise(cor = cor.spear.abund, tss = spec + sens -1)

gg <- ggplot(data = gg.dat, aes(x = pfg, y = tss, fill = cor)) + geom_bar(stat = 'identity') + facet_grid(~habitat) + coord_flip()
gg


