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

work.dir <- "~/Work/FATEHD/benchmarking/workdir/"
input.dir <- "~/Work/FATEHD/benchmarking/results/"
output.dir <- "~/Work/FATEHD/benchmarking/figures/"

dir.create(work.dir, showWarnings = FALSE, recursive = TRUE)
dir.create(input.dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)

setwd(work.dir)

##' execute the scp command to retrieave results from all cigri clusters

## get input parameters file
params <- read.csv(file.path(input.dir, "uncoupledLHS.csv"), row.names = 1, stringsAsFactors = FALSE)
params$simul.id <- rownames(params)
head(params)

## check jobs that have been completed
params$job.status <- file.exists(file.path(input.dir, "fhdsoilpt", "5590", params$simul.id)) 
sum(params$job.status)
length(params$job.status)

## reshape params
params <- params %>% mutate(hs.type = factor(ceiling(hs.type)),
                            area = factor(ceiling(area)))

year_to_keep <- paste0("year", c(830,840,850))

# ## create the summary table of jobs
# res.occ.list <- mclapply(params$simul.id[params$job.status], function(jid_){
#   res1_ <- read.table(file.path(input.dir, "fhdsoilpt", "5590", jid_, paste0("eval.occ.test.531_", jid_,".txt")), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#   ## keep only a subset of the data
#   res_ <- res1_ %>% filter(is.element(year, year_to_keep)) %>%  mutate(simul.id = jid_)
#   return(res_)
# }, mc.cores = 16)
# 
# eval.occ.df <- bind_rows(res.occ.list)
# ## merge with parameters
# eval.occ.df <- eval.occ.df %>% left_join(params)
# 
# save(eval.occ.df, file = file.path(work.dir, "eval.occ.df_fhdsoilpt.RData"))
(load(file.path(work.dir, "eval.occ.df_fhdsoilpt.RData")))

# ## create the summary table of jobs
# library(parallel)
# res.abu.list <- mclapply(params$simul.id[params$job.status], function(jid_){
#   res1_ <- read.table(file.path(input.dir, "fhdsoilpt", "5590", jid_, paste0("eval.abund.test.531_", jid_,".txt")), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#   ## keep only a subset of the data
#   res_ <- res1_ %>% filter(is.element(year, year_to_keep)) %>%  mutate(simul.id = jid_)
#   return(res_)
# }, mc.cores = 16)
# 
# eval.abu.df <- bind_rows(res.abu.list)
# ## merge with parameters
# eval.abu.df <- eval.abu.df %>% left_join(params)
# 
# save(eval.abu.df, file = file.path(work.dir, "eval.abu.df_fhdsoilpt.RData"))
(load(file.path(work.dir, "eval.abu.df_fhdsoilpt.RData")))

##' # Explore campain results

##' ## Which parfameters are important for global tss?

dat.name <- c('eval.occ.df', 'eval.abu.df')
# library(lme4)

# dat.name_ <- dat.name[2]
for(dat.name_ in dat.name){
  eval.df_ <- get(dat.name_) %>% filter(habitat == "all") %>% mutate(tss = sens + spec - 1) 
  eval.df_$tss[is.na(eval.df_$tss)] <- 0
  eval.df_ <- eval.df_ %>% group_by(simul.id, dat.source) %>% summarize(tss = mean(tss)) %>% left_join(params)
  
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
  
  gg.dat <- eval.df_ %>% group_by(area) %>% filter(tss > quantile(tss, 0.95)) %>% mutate(subset = 'best.5.pc') %>% bind_rows(eval.df_ %>% mutate(subset = 'all')) %>% gather(param.name, param.value, grep("pfg|glob", names(.)))
  gg <- ggplot(gg.dat, aes(x = param.value, colour = area, lty = subset)) + geom_density() + facet_wrap(~param.name, scales = 'free') + 
    geom_vline(data = gg.dat %>% group_by(area, subset, param.name) %>% summarise(param.value.median = median(param.value, na.rm = TRUE)),
                aes(xintercept = param.value.median, colour = area), lty = 3)
  print(gg)
  
}

##' ## Which parfameters are important for global abundances correlations ?
for(dat.name_ in dat.name){
  eval.df_ <- get(dat.name_) %>% filter(habitat == "all")
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
  
  gg.dat <- eval.df_ %>% group_by(area) %>% filter(cor.spear.abund > quantile(cor.spear.abund, 0.95, na.rm = TRUE)) %>% mutate(subset = 'best.5.pc') %>% bind_rows(eval.df_ %>% mutate(subset = 'all')) %>% gather(param.name, param.value, grep("pfg|glob", names(.)))
  gg <- ggplot(gg.dat, aes(x = param.value, colour = area, lty = subset)) + geom_density() + facet_wrap(~param.name, scales = 'free') + 
    geom_vline(data = gg.dat %>% group_by(area, subset, param.name) %>% summarise(param.value.median = median(param.value, na.rm = TRUE)),
               aes(xintercept = param.value.median, colour = area), lty = 3)
  print(gg)
  
  gg.dat <- eval.df_  %>% group_by(area) %>% filter(cor.spear.abund > quantile(cor.spear.abund, 0.95, na.rm = TRUE)) %>% mutate(subset = 'best.5.pc') %>% bind_rows(eval.df_ %>% mutate(subset = 'all')) %>% gather(param.name, param.value, grep("pfg|glob", names(.))) %>% filter(is.element(param.name, csa.vars))
  gg <- ggplot(gg.dat, aes(x = param.value, colour = area, lty = subset)) + geom_density() + facet_wrap(~param.name, scales = 'free') + 
    geom_vline(data = gg.dat %>% group_by(area, subset, param.name) %>% summarise(param.value.median = median(param.value, na.rm = TRUE)),
               aes(xintercept = param.value.median, colour = area, lty = subset))
  print(gg)
  
  gg.dat <- eval.df_  %>% group_by(area)  %>% filter(hs.type == 1) %>% filter(cor.spear.abund > quantile(cor.spear.abund, 0.95, na.rm = TRUE)) %>% mutate(subset = 'best.5.pc') %>% bind_rows(eval.df_ %>% mutate(subset = 'all')) %>% gather(param.name, param.value, grep("pfg|glob", names(.))) %>% filter(is.element(param.name, csa.vars))
  gg <- ggplot(gg.dat, aes(x = param.value, colour = area, lty = subset)) + geom_density() + facet_wrap(~param.name, scales = 'free') + 
    geom_vline(data = gg.dat %>% group_by(area, subset, param.name) %>% summarise(param.value.median = median(param.value, na.rm = TRUE)),
               aes(xintercept = param.value.median, colour = area, lty = subset))
  print(gg)
  
}

## get the best parameters for tss optim

dat_ <- eval.abu.df %>% filter(habitat == "all") %>% 
  mutate(area = ceiling(area), hs.type = ceiling(hs.type), tss = sens + spec - 1) %>% 
  group_by(area, simul.id) %>% summarize(tss = mean(tss, na.rm = TRUE), cor.spear.abund = mean(cor.spear.abund, na.rm = TRUE), sum.tss.cor = tss + cor.spear.abund) %>% ungroup() %>% group_by(area) 

dat_ %>% filter(tss == max(tss, na.rm = TRUE)) %>% left_join( dplyr::select(params, - area)) ## %>% as.data.frame
dat_ %>% filter(cor.spear.abund == max(cor.spear.abund, na.rm = TRUE)) %>% left_join( dplyr::select(params, - area)) ## %>% as.data.frame
dat_ %>% filter(sum.tss.cor == max(sum.tss.cor, na.rm = TRUE)) %>% left_join( dplyr::select(params, - area)) ## %>% as.data.frame
