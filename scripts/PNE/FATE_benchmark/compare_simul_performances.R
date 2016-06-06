##' ----------------------------------------------------------------------------
##' @title Compare Simulation Perforamnce
##' @description The aim of this script is to check that we made progress in
##'   PNE vegetation modelling by model parameters tuning
##' @date 06/04/2016
##' @author damien georges
##' ----------------------------------------------------------------------------

## script initialisation -------------------------------------------------------
rm(list = ls())

work.dir <- "~/Work/FATEHD/benchmarking/workdir/"
# setwd(work.dir)

library(dplyr)
library(tidyr)
library(ggplot2)

year_to_keep <- paste0("year", c(850))

## get the initial version of FATEHD run on PNE --------------------------------
eval.abund.ref <- read.table("/home/georgeda/Work/FATEHD/benchmarking/fatehd_params_test/FATE_newHS/SIMUL_6STRATA/simul_comparaisons/1/eval.abund.test.1.txt", h = TRUE)
eval.abund.ref <- eval.abund.ref %>% 
  mutate(habitat = factor(habitat, levels = c("0", "31", "40", "50", "60", "70", "81", "83", "All"),
                          labels = c("Excluded","Rock", "Grasslands", "Lowlands","Open_habitats", "Semi-closed_habitats", "Closed_habitats","Forests", "All"))) %>%
  filter(is.element(year, year_to_keep))


eval.occ.ref <- read.table("/home/georgeda/Work/FATEHD/benchmarking/fatehd_params_test/FATE_newHS/SIMUL_6STRATA/simul_comparaisons/1/eval.occ.test.1.txt", h = TRUE)
eval.occ.ref <- eval.occ.ref %>% 
  mutate(habitat = factor(habitat, levels = c("0", "31", "40", "50", "60", "70", "81", "83", "All"),
                          labels = c("Excluded","Rock", "Grasslands", "Lowlands","Open_habitats", "Semi-closed_habitats", "Closed_habitats","Forests", "All"))) %>%
  filter(is.element(year, year_to_keep))

## get the soil optimized version of FATEHD -----------------------------------
eval.abund.bs <- get(load(file.path(work.dir, "eval.abu.df_fhdbsoilpt.RData"))) %>% 
  filter(simul.id == 9421, year == "year850") 
eval.occ.bs <- get(load(file.path(work.dir, "eval.occ.df_fhdbsoilpt.RData"))) %>% 
  filter(simul.id == 9421, year == "year850") 
## simul comparaison ----------------------------------------------------------

### 

gg.dat.ref_ <- eval.abund.ref %>% 
  mutate(tss = sens + spec - 1) %>% 
  select(dat.source, habitat, pfg, year, nb.occ, nb.abs, nb.pred.occ, nb.pred.abs, sens, spec, pcc, kappa, auc, tss) %>%
  gather("metric.name", "metric.val.ref", - c(1:4))

gg.dat.bs_ <- eval.abund.bs %>% 
  mutate(tss = sens + spec - 1) %>% 
  select(dat.source, habitat, pfg, year, nb.occ, nb.abs, nb.pred.occ, nb.pred.abs, sens, spec, pcc, kappa, auc, tss) %>%
  gather("metric.name", "metric.val.bs", - c(1:4))

gg.dat <- inner_join(gg.dat.ref_, gg.dat.bs_)

gg.dat <- gg.dat %>% 
  mutate(metric.val.diff = metric.val.bs - metric.val.ref,
         pfg.short = sub("_.*$", "", pfg)) %>%
  filter(is.element(metric.name, c("sens", "spec", "pcc", "auc", "tss")))

head(gg.dat)


# gg <- ggplot(data = gg.dat, aes(y = metric.val.diff, x = metric.name, fill = metric.name)) +
#   geom_bar(stat = "identity") + facet_grid(pfg ~ habitat) + coord_flip()
# gg

gg <- ggplot(data = gg.dat, aes(metric.name, metric.val.diff, fill = metric.name)) +
  geom_bar(stat = "identity", position = "dodge") + facet_grid(pfg.short ~ habitat) + coord_flip() + 
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + 
  ylab("metric improvment in soil simul") + xlab("") + scale_fill_discrete(name="")

#+ , fig.height = 25, fig.width = 15
print(gg)


gg <- ggplot(data = gg.dat, aes(habitat, metric.val.diff, fill = habitat)) +
  geom_bar(stat = "identity", position = "dodge") + facet_grid(pfg.short ~ metric.name) + coord_flip() + 
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + 
  ylab("metric improvment in soil simul") + xlab("") + scale_fill_discrete(name="")

#+ , fig.height = 25, fig.width = 15
print(gg)





