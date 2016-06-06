##' ---
##' title : FATEHD get the best models parameters from full models parameters sensitivity analysis
##' author : damien g.
##' ---

##' # Descrtiption 
##' 
##' We made a run on 3 different subarea from PNE to determine the influence
##' of models parameters (constants) and to determine the optimal value of 
##' each of this parameters.
##' 
##' Tested parameters are:
##' 
##' | parameter name  | description  |
##' | ------------- | ------------- |
##' | "global.low.abund" | the number of units of light that should consume a low abundance pfg |
##' | "global.medium.abund" | the number of units of light that should consume a medium abundance pfg |
##' | "global.high.abund" | the number of units of light that should consume a high abundance pfg |
##' | "global.low.resources.tresh" | the number of units of light above wich we are in low light condition |
##' | "global.medium.resources.tresh" | the number of units of light above wich we are in low medium condition |
##' | "global.max.by.cohort" | th maximal units of light a single cohort should consume |
##' | "global.full.soil.coverage" | the number of light units that should be occupied to consider that we are in full coverage  |
##' | "envsuit.option" | the way the annual envsuit threshold is drawn (1: (default) a random number by cell by year; 2: a random mean and sd by year, then the random hs trhesh is drawn in a normal distribution) |
##' | "seeding.timestep" | the seeding timestep |
##' | "seeding.duration" | the seeding duration |
##' | "soil.default.value" | the starting value of soil |
##' | "soil.percent" | the percentage of decreasing of the soil value after a disturbance occurs |
##' | "soil.categories.threshold.low" | the soil value under which we are in low soil contition |
##' | "soil.categories.threshold.medium" | the soil value above which we are in high soil contition |
##' | "disturbance.changing.times.1" | the timing of the first deforestation event |
##' | "disturbance.changing.times.2" | the timing of the second deforestation event |
##' | "hs.type" | the type of habitat maps (1: default; 2: ; 3: ; 4: ) |
##' | "area" | the PNE sub-area considered |
##' 

##' # Get simul results 
##' 
##' Simul outputs are stored on each cigri cluters, results are stored in 
##' luke in the directory "~/fhdmfpt"
##' 

#+ knitr option, echo = FALSE 
library(knitr)
knitr::opts_chunk$set(fig.width=12, fig.height=8, fig.path='Figs_fhdmfpt/',
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

# setwd(work.dir)

##' execute the scp command to retrieave results from All cigri clusters

## get input parameters file
params <- read.csv(file.path(input.dir, "params.fhdmfpt.txt"), row.names = NULL, sep = "\t", stringsAsFactors = FALSE)
params$simul.id <- rownames(params)
head(params)

## check jobs that have been completed
params$job.status <- file.exists(file.path(input.dir, "fhdmfpt", "5864", params$simul.id)) 
sum(params$job.status)
length(params$job.status)

# year_to_keep <- paste0("year", c(830,840,850))

# ## create the summary table of jobs
# res.occ.list <- mclapply(params$simul.id[params$job.status], function(jid_){
#   res1_ <- read.table(file.path(input.dir, "fhdmfpt", "5864", jid_, paste0("eval.occ.test.531_", jid_,".txt")), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#   ## keep only a subset of the data
#   res_ <- res1_ %>% 
#     # filter(is.element(year, year_to_keep)) %>%  
#     mutate(simul.id = jid_)
#   return(res_)
# }, mc.cores = 16)
# 
# eval.occ.df <- bind_rows(res.occ.list)
# ## merge with parameters
# eval.occ.df <- eval.occ.df %>% left_join(params)
# eval.occ.df <- eval.occ.df %>% 
#   mutate(habitat = factor(habitat, levels = c("0", "31", "40", "50", "60", "70", "81", "83", "All"),
#                           labels = c("Excluded","Rock", "Grasslands", "Lowlands","Open_habitats", "Semi-closed_habitats", "Closed_habitats","Forests", "All")),
#          tss = spec + sens - 1)
# 
# 
# save(eval.occ.df, file = file.path(work.dir, "eval.occ.df_fhdmfpt.RData"))
# eval.occ.df <- eval.occ.df %>% filter(habitat == "All")
# save(eval.occ.df, file = file.path(work.dir, "eval.occ.df.hab.all_fhdmfpt.RData"))
# (load(file.path(work.dir, "eval.occ.df_fhdmfpt.RData")))
(load(file.path(work.dir, "eval.occ.df.hab.all_fhdmfpt.RData")))

# ## create the summary table of jobs
# res.abu.list <- mclapply(params$simul.id[params$job.status], function(jid_){
#   res1_ <- read.table(file.path(input.dir, "fhdmfpt", "5864", jid_, paste0("eval.abund.test.531_", jid_,".txt")), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#   ## keep only a subset of the data
#   res_ <- res1_ %>% 
#     # filter(is.element(year, year_to_keep)) %>%  
#     mutate(simul.id = jid_)
#   return(res_)
# }, mc.cores = 16)
# 
# eval.abu.df <- bind_rows(res.abu.list)
# ## merge with parameters
# eval.abu.df <- eval.abu.df %>% left_join(params) 
# eval.abu.df <- eval.abu.df %>% 
#   mutate(habitat = factor(habitat, levels = c("0", "31", "40", "50", "60", "70", "81", "83", "All"),
#                           labels = c("Excluded","Rock", "Grasslands", "Lowlands","Open_habitats", "Semi-closed_habitats", "Closed_habitats","Forests", "All")),
#          tss = spec + sens - 1)
# 
# save(eval.abu.df, file = file.path(work.dir, "eval.abu.df_fhdmfpt.RData"))
# eval.abu.df <- eval.abu.df %>% filter(habitat == "All")
# save(eval.abu.df, file = file.path(work.dir, "eval.abu.df.hab.all_fhdmfpt.RData"))
# (load(file.path(work.dir, "eval.abu.df_fhdmfpt.RData")))
(load(file.path(work.dir, "eval.abu.df.hab.all_fhdmfpt.RData")))

##' # Explore campain results

##' ## Which parfameters are important for global tss?

dat.name <- c('eval.occ.df', 'eval.abu.df')
# library(lme4)

# dat.name_ <- dat.name[1]
for(dat.name_ in dat.name){
  eval.df_ <- get(dat.name_) %>% filter(habitat == "All") %>% mutate(tss = sens + spec - 1) 
  eval.df_$tss[is.na(eval.df_$tss)] <- 0
  eval.df_ <- eval.df_ %>% group_by(simul.id, dat.source, habitat) %>% summarize(tss = mean(tss)) %>% left_join(params)
  
  head(eval.df_)
  # tss.mod <- lmer(as.formula(paste0("tss ~ ", paste0(colnames(params %>% select(- c(simul.id, job.status, area))), collapse = " + "), paste0(" + (1|area)"))), data = eval.df_)
  tss.lower.formula <- as.formula("tss ~ 1")
  tss.upper.formula <- as.formula(paste0("tss ~ ", paste0(colnames(params %>% select(- c(simul.id, job.status, area))), collapse = " + ")))
  tss.mod <- lm(tss.lower.formula, data = eval.df_)
  tss.mod <- step(tss.mod, scope = tss.upper.formula)
  print(summary(tss.mod))
  
  tss.vars <- attr(attr(tss.mod$model, "terms"), "term.labels")
  tss.coef <- coefficients(tss.mod)[-1]
  # names(tss.coef) <- tss.vars
  tss.vi <- biomod2::variables_importance(tss.mod, data = eval.df_, nb_rand = 3)$mat[tss.vars, ]
  
  gg.dat <- data.frame(importance = rowMeans(tss.vi), effect = c("decrease", "increase")[(tss.coef > 0) + 1], var = tss.vars, var.grp = sub("[[:punct:]].*$", "", tss.vars)) 
  
  gg <- ggplot(gg.dat, aes(x = factor(var), y = importance, fill = effect)) +
    ylim(0,1) + 
    geom_bar(aes(group = var.grp), stat = "identity") + 
    labs(title = paste0("variable importance on global tss (", dat.name_, ")"), x = "", y ="") +
    coord_flip()
  print(gg)
  
  ## distribution of parameters in the most performant simulations
  
  gg.dat <- eval.df_ %>% group_by(area) %>% filter(tss > quantile(tss, 0.95)) %>% mutate(subset = 'best.5.pc') %>% bind_rows(eval.df_ %>% mutate(subset = 'All')) %>% gather(param.name, param.value, grep("pfg|glob", names(.)))
  gg <- ggplot(gg.dat, aes(x = param.value, colour = area, lty = subset)) + geom_density() + facet_wrap(~param.name, scales = 'free') + 
    geom_vline(data = gg.dat %>% group_by(area, subset, param.name) %>% summarise(param.value.median = median(param.value, na.rm = TRUE)),
                aes(xintercept = param.value.median, colour = area), lty = 3)
  print(gg)
  
}

##' ## Which parfameters are important for global abundances correlations ?

dat.name <- c('eval.abu.df')
for(dat.name_ in dat.name){
  eval.df_ <- get(dat.name_) %>% filter(habitat == "All")
  eval.df_ <- eval.df_ %>% group_by(simul.id, dat.source) %>% summarize(cor.spear.abund = mean(cor.spear.abund, na.rm = TRUE)) %>% left_join(params)
  
  head(eval.df_)
  # csa.mod <- lmer(as.formula(paste0("cor.spear.abund ~ ", paste0(colnames(params %>% select(- c(simul.id, job.status, area))), collapse = " + "), paste0(" + (1|area)"))), data = eval.df_)
  csa.lower.formula <- as.formula("cor.spear.abund ~ 1")
  csa.upper.formula <- as.formula(paste0("cor.spear.abund ~ ", paste0(colnames(params %>% select(- c(simul.id, job.status, area))), collapse = " + ")))
  csa.mod <- lm(csa.lower.formula, data = eval.df_)
  csa.mod <- step(csa.mod, scope = csa.upper.formula)
  print(summary(csa.mod))
  
  csa.vars <- attr(attr(csa.mod$model, "terms"), "term.labels")
  csa.coef <- coefficients(csa.mod)[-1]
  # names(csa.coef) <- csa.vars
  csa.vi <- biomod2::variables_importance(csa.mod, data = eval.df_, nb_rand = 3)$mat[csa.vars, ]
  
  gg.dat <- data.frame(importance = rowMeans(csa.vi), effect = c("decrease", "increase")[(csa.coef > 0) + 1], var = csa.vars, var.grp = sub("[[:punct:]].*$", "", csa.vars)) 
  
  gg <- ggplot(gg.dat, aes(x = factor(var), y = importance, fill = effect)) +
    ylim(0,1) + 
    geom_bar(aes(group = var.grp), stat = "identity") + 
    labs(title = paste0("variable importance on global abundances correlation (", dat.name_, ")"), x = "", y ="") +
    coord_flip()
  print(gg)
  
  ## distribution of parameters in the most performant simulations
  
  gg.dat <- eval.df_ %>% group_by(area) %>% filter(cor.spear.abund > quantile(cor.spear.abund, 0.95, na.rm = TRUE)) %>% mutate(subset = 'best.5.pc') %>% bind_rows(eval.df_ %>% mutate(subset = 'All')) %>% gather(param.name, param.value, grep("pfg|glob", names(.)))
  gg <- ggplot(gg.dat, aes(x = param.value, colour = area, lty = subset)) + geom_density() + facet_wrap(~param.name, scales = 'free') + 
    geom_vline(data = gg.dat %>% group_by(area, subset, param.name) %>% summarise(param.value.median = median(param.value, na.rm = TRUE)),
               aes(xintercept = param.value.median, colour = area), lty = 3)
  print(gg)
  
  gg.dat <- eval.df_  %>% group_by(area) %>% filter(cor.spear.abund > quantile(cor.spear.abund, 0.95, na.rm = TRUE)) %>% mutate(subset = 'best.5.pc') %>% bind_rows(eval.df_ %>% mutate(subset = 'All')) %>% gather(param.name, param.value, grep("pfg|glob", names(.))) %>% filter(is.element(param.name, csa.vars))
  gg <- ggplot(gg.dat, aes(x = param.value, colour = area, lty = subset)) + geom_density() + facet_wrap(~param.name, scales = 'free') + 
    geom_vline(data = gg.dat %>% group_by(area, subset, param.name) %>% summarise(param.value.median = median(param.value, na.rm = TRUE)),
               aes(xintercept = param.value.median, colour = area, lty = subset))
  print(gg)
  
  gg.dat <- eval.df_  %>% group_by(area)  %>% filter(hs.type == 1) %>% filter(cor.spear.abund > quantile(cor.spear.abund, 0.95, na.rm = TRUE)) %>% mutate(subset = 'best.5.pc') %>% bind_rows(eval.df_ %>% mutate(subset = 'All')) %>% gather(param.name, param.value, grep("pfg|glob", names(.))) %>% filter(is.element(param.name, csa.vars))
  gg <- ggplot(gg.dat, aes(x = param.value, colour = area, lty = subset)) + geom_density() + facet_wrap(~param.name, scales = 'free') + 
    geom_vline(data = gg.dat %>% group_by(area, subset, param.name) %>% summarise(param.value.median = median(param.value, na.rm = TRUE)),
               aes(xintercept = param.value.median, colour = area, lty = subset))
  print(gg)
  
}

##' # get the best parameters combination
##' 
##' To define which are the best combination of parameters, we will (by PNE subarea) rank the models in term of best mean TSS
##' acording to both PA and abundance data, same vith the correlation of abundances and get the 5 best simulation  by area
##' according to the sum of ranks

## rank based on tss occ
occ.tss.ranked <- eval.occ.df %>% filter(habitat == "All") %>%
  group_by(area, simul.id) %>% summarize(tss.occ = mean(tss, na.rm = TRUE)) %>% 
  ungroup() %>% group_by(area) %>%
  mutate(tss.occ.rank = rank(tss.occ))

## rank based on tss abu
abu.tss.ranked <- eval.abu.df %>% filter(habitat == "All") %>%
  group_by(area, simul.id) %>% summarize(tss.abu = mean(tss, na.rm = TRUE)) %>% 
  ungroup() %>% group_by(area) %>%
  mutate(tss.abu.rank = rank(tss.abu))

## rank based on spearman cor abu
abu.cor.ranked <- eval.abu.df %>% filter(habitat == "All") %>%
  group_by(area, simul.id) %>% summarize(cor.spear.abun = mean(cor.spear.abund, na.rm = TRUE)) %>% 
  ungroup() %>% group_by(area) %>%
  mutate(cor.spear.abun.rank = rank(cor.spear.abun))

## calculate the global rank score
dat.rank_ <- occ.tss.ranked %>% full_join(abu.tss.ranked) %>% full_join(abu.cor.ranked) %>% 
  group_by(area, simul.id) %>% mutate(glob.rank.score = tss.occ.rank + tss.abu.rank + cor.spear.abun.rank) %>%
  ungroup() %>% group_by(area) %>%
  mutate(glob.rank = rank(glob.rank.score))

## get the 6 best simul by area
nb.simul.selected <- 30
dat.best.rank_ <- dat.rank_ %>% mutate(glob.inv.rank = max(glob.rank) + 1 - glob.rank) %>% filter(glob.inv.rank <= nb.simul.selected) %>% as.data.frame
dat.best.rank_

## represent the distribution of parameters for the 30 best simul by area
dat.best.rank_ <- dat.best.rank_ %>% left_join(params)

gg.dat_ <- dat.best.rank_ %>% gather( param.name, param.value, 
                                      c(global.low.abund, global.medium.abund, 
                                        global.high.abund, global.low.resources.tresh, 
                                        global.medium.resources.tresh, global.max.by.cohort, 
                                        global.full.soil.coverage, envsuit.option, seeding.timestep, 
                                        seeding.duration, soil.default.value, soil.percent, 
                                        soil.categories.threshold.low, soil.categories.threshold.medium, 
                                        disturbance.changing.times.1, disturbance.changing.times.2, hs.type)) %>%
  mutate(glob.inv.rank.bin = cut(glob.inv.rank, seq(0,nb.simul.selected, 5)))

gg <- ggplot(data = gg.dat_, aes(param.value,  fill = as.factor(glob.inv.rank.bin))) + #colour = as.factor(area),
  geom_histogram() + facet_wrap(~param.name, scales = "free") + scale_fill_manual(values=colorRampPalette(c("red", "white"))( nb.simul.selected %/% 5))
print(gg)

gg <- ggplot(data = gg.dat_, aes(param.value,  fill = as.factor(area))) + 
  geom_histogram() + facet_wrap(~param.name, scales = "free") 
print(gg)

best.simul.subset <- dat.best.rank_ %>% group_by(area) %>% filter(glob.inv.rank <= 3) %>% as.data.frame 
kable(best.simul.subset)

## create the best parameters file for grid computation
grid.params <- best.simul.subset$simul.id
writeLines(grid.params, con = "~/Work/FATEHD/benchmarking/fatehd_params_test/grid/model_full_params_test_PNE/params_fhdmfpne.txt")
