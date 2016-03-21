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
knitr::opts_chunk$set(fig.width=12, fig.height=8, fig.path='Figs/',
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
params <- read.table(file.path(input.dir, "params_formal_fhdmpt.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(params)

## check jobs that have been completed
params$job.status <- file.exists(file.path(input.dir, "fhdmpt_simul_comparaisons", params$simul.id)) 
sum(params$job.status)
length(params$job.status)

## create the summary table of jobs
res.list <- lapply(params$simul.id[params$job.status], function(jid_){
  res1_ <- read.table(file.path(input.dir, "fhdmpt_simul_comparaisons", jid_, paste0("comp.test.ref.summ.", jid_,".txt")), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  res2_ <- read.table(file.path(input.dir, "fhdmpt_simul_comparaisons", jid_, paste0("comp.test.ref.", jid_,".txt")), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  res_ <- bind_rows(res1_, res2_)
  res_$simul.id <- jid_
  return(res_)
})

res.df <- do.call(rbind, res.list)

## merge with parameters
res.df <- res.df %>% left_join(params)

save(res.df, file = "res.df.RData")

## create the eval table of jobs
eval.list <- lapply(params$simul.id[params$job.status], function(jid_){
  res_ <- read.table(file.path(input.dir, "fhdmpt_simul_comparaisons", jid_, paste0("eval.occ.test.", jid_,".txt")), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  res_$simul.id <- jid_
  return(res_)
})

eval.df <- do.call(rbind, eval.list)

## merge with parameters
eval.df <- eval.df %>% left_join(params)

save(eval.df, file = "eval.df.RData")

##' # Explore campain results

##' ## Find best models parameters combination

eval.df_ <- eval.df %>% group_by(year, pfg) %>% summarise(best.job = simul.id[which.max(sens + spec)[1]], 
                                                          sens = sens[which.max(sens + spec)[1]],
                                                          spec = spec[which.max(sens + spec)[1]])

##' ## importance of each parameters

params.to.test <- c("envsuit.option", "seeding.params", "dispers.mode", "global.abund", "global.resource.thresh", "max.by.cohort", "area")
param.to.test <- "envsuit.option"

# pdf(file.path(output.dir, paste0("eval_fhd_densplot.pdf")), width = 7 , height = 14)

#+ eval_fhd_densplot, dev = c('png'), fig.width = 7, fig.height = 14
for(param.to.test in c(params.to.test)){
  cat("\n>", param.to.test)
  gg.dat <- res.df %>% filter(is.element(year, c("year840", "year850"))) %>% 
    select_("pfg", "year", "sens", "spec", param.to.test) %>% 
    mutate(sens_spec = sens + spec, pfg_short = sub("_.*$", "", pfg)) %>%
    mutate_(.dots = setNames(param.to.test, "param.val")) %>%
    gather(stat.name, stat.val, sens, spec, sens_spec)

  gg <- ggplot(gg.dat, aes(x = stat.val, colour = factor(param.val))) + geom_density() + facet_grid(pfg_short~stat.name, scales = "free") + 
    scale_color_discrete(name = param.to.test) + xlab("")
  
  print(gg)
}
# dev.off()

# pdf(file.path(output.dir, paste0("eval_fhd_boxplot.pdf")), width = 7 , height = 14)

#+ eval_fhd_boxplot, dev = c('png'), fig.width = 7, fig.height = 14
for(param.to.test in c(params.to.test)){
  cat("\n>", param.to.test)
  gg.dat <- res.df %>% filter(is.element(year, c("year840", "year850"))) %>% 
    select_("pfg", "year", "sens", "spec", param.to.test) %>% 
    mutate(sens_spec = sens + spec, pfg_short = sub("_.*$", "", pfg)) %>%
    mutate_(.dots = setNames(param.to.test, "param.val")) %>%
    gather(stat.name, stat.val, sens, spec, sens_spec)

  gg <- ggplot(gg.dat, aes(param.val, stat.val, fill = factor(param.val))) + geom_boxplot(varwidth = TRUE) + facet_grid(pfg_short~stat.name, scales = "free") + 
    scale_fill_discrete(name = param.to.test) + xlab("")
  
  print(gg)
}
# dev.off()

##' detect the distribution of parameters that lead to the 5 percent best simulations

##' criterria: max sens + spec
gg.dat <- res.df %>% filter(is.element(year, c("year830", "year840", "year850"))) %>% 
  filter(pfg == "all") %>% 
  select(pfg, year, sens, spec, one_of(params.to.test)) %>% 
  mutate(sens_spec = sens + spec, pfg_short = sub("_.*$", "", pfg)) %>%
  filter(sens_spec >= quantile(sens_spec, 0.95)) %>%
  gather(stat.name, stat.val, sens, spec, sens_spec) %>%
  gather_("param.name", "param.value", params.to.test) %>%
  distinct()


# pdf(file.path(output.dir, paste0("eval_fhd_barplot_top_5pc_sens_spec.pdf")), width = 7 , height = 14)

#+ eval_fhd_barplot_top_5pc_sens_spec, dev = c('png'), fig.width = 7, fig.height = 14
  gg <- ggplot(gg.dat %>% filter(stat.name == 'sens_spec'), aes(x = param.value, fill = factor(param.value))) + geom_bar() + facet_grid(param.name~stat.name) + 
    scale_fill_discrete(name = "parameter value") + xlab("")
  print(gg)
# dev.off()

# pdf(file.path(output.dir, paste0("eval_fhd_boxplot_top_5pc_sens_spec.pdf")), width = 7 , height = 14)

#+ eval_fhd_boxplot_top_5pc_sens_spec, dev = c('png'), fig.width = 7, fig.height = 21
  gg <- ggplot(gg.dat, aes(param.value, stat.val, fill = factor(param.value))) + geom_boxplot(varwidth = TRUE) + facet_grid(param.name~stat.name) + 
    scale_fill_discrete(name = "parameter value") + xlab("")
  print(gg)
# dev.off()

##' criterria: max sens 
gg.dat <- res.df %>% filter(is.element(year, c("year830", "year840", "year850"))) %>% 
  filter(pfg == "all") %>% 
  select(pfg, year, sens, spec, one_of(params.to.test)) %>% 
  mutate(sens_spec = sens + spec, pfg_short = sub("_.*$", "", pfg)) %>%
  filter(sens >= quantile(sens, 0.95)) %>%
  gather(stat.name, stat.val, sens, spec, sens_spec) %>%
  gather_("param.name", "param.value", params.to.test) %>%
  distinct()


# pdf(file.path(output.dir, paste0("eval_fhd_barplot_top_5pc_sens.pdf")), width = 7 , height = 14)

#+ eval_fhd_barplot_top_5pc_sens, dev = c('png'), fig.width = 7, fig.height = 14
gg <- ggplot(gg.dat %>% filter(stat.name == 'sens'), aes(x = param.value, fill = factor(param.value))) + geom_bar() + facet_grid(param.name~stat.name) + 
  scale_fill_discrete(name = "parameter value") + xlab("")
print(gg)
# dev.off()

# pdf(file.path(output.dir, paste0("eval_fhd_boxplot_top_5pc_max_coocc_sens.pdf")), width = 7 , height = 14)

#+ eval_fhd_boxplot_top_5pc_max_coocc_sens, dev = c('png'), fig.width = 7, fig.height = 14
gg <- ggplot(gg.dat, aes(param.value, stat.val, fill = factor(param.value))) + geom_boxplot(varwidth = TRUE) + facet_grid(param.name~stat.name) + 
  scale_fill_discrete(name = "parameter value") + xlab("")
print(gg)
# dev.off()

##' get the 20 best simulations by area
gg.dat <- res.df %>% group_by(area) %>% filter(is.element(year, c("year830", "year840", "year850"))) %>% 
  filter(pfg == "all") %>% 
  select(simul.id, pfg, year, sens, spec, one_of(params.to.test)) %>% 
  mutate(sens_spec = sens + spec, pfg_short = sub("_.*$", "", pfg)) %>%
  filter(is.element(sens_spec,  tail(sort(sens_spec), 25))) 

for(v_ in params.to.test){
  cat("\n")
  print(table(gg.dat %>% select_(v_)))
  cat("\n")
}

best.simul.id <- unique(gg.dat$simul.id)

gg.dat <- res.df %>% filter(is.element(simul.id, best.simul.id)) %>%
  select(simul.id, pfg, year, sens, spec, one_of(params.to.test)) %>% 
  mutate(sens_spec = sens + spec, pfg_short = sub("_.*$", "", pfg)) %>%
  gather(stat.name, stat.val, sens, spec, sens_spec) 

## which pfgs are upgraded in the best simulations
# pdf(file.path(output.dir, paste0("eval_fhd_pfg_imporve_top_20simulbyarea_sens_spec.pdf")), width = 7 , height = 21)

#+ eval_fhd_pfg_imporve_top_20simulbyarea_sens_spec, dev = c('png'), fig.width = 7, fig.height = 21
  gg <- ggplot(gg.dat%>% ungroup, aes(stat.name, stat.val, fill = factor(stat.name))) + geom_boxplot(varwidth = TRUE) + facet_grid(pfg_short~.) + 
    scale_fill_discrete(name = "stat") + xlab("")
  print(gg)
# dev.off()

##' According to this exploration I think that the best conbination of varaibles to use are:
##' - "envsuit.option" = 2 (or 1)
##' - "seeding.params" = 2
##' - "dispers.mode" = 3
##' - "global.abund" = 3
##' - "global.resource.thresh" = 2
##' - "max.by.cohort" = 1 (or 3)

params %>% filter(area == 1, seeding.params == 2, dispers.mode == 3, global.abund == 3, global.resource.thresh == 2, max.by.cohort == 1)

# envsuit.option seeding.params dispers.mode global.abund global.resource.thresh max.by.cohort area simul.id
# 1              1              2            3            3                      2             1    1      163
# 2              2              2            3            3                      2             1    1      164

#######################################################################

##' Integrate the pfg's cohexistance in performances calculations


gg.dat <- eval.df %>% filter(is.element(year, c("year830", "year840", "year850"))) %>% ## keep only year of interest
  dplyr::select(simul.id, pfg, year, nb.occ, nb.abs, nb.pred.occ, nb.pred.abs, sens, spec, one_of(params.to.test)) %>%  ## remove useless columns
  mutate(sens_spec = sens + spec, pfg_short = sub("_.*$", "", pfg)) %>% ## create the sens + spec stat
  gather(stat.name, stat.value, sens, spec, sens_spec) %>% ## reshape statistiqueq
  distinct() ## remove duplicated rows
gg.dat$stat.val[is.na(gg.dat$stat.val) & gg.dat$nb.occ > 0 & gg.dat$nb.abs > 0] <- 0 ## set to 0 the stats that are not calculated because of species disparition

gg.dat.all <- gg.dat %>% group_by(simul.id, stat.name, year) %>% 
  summarize(pfg = "all", nb.pfg = sum(nb.pred.occ > 0, na.rm = TRUE),
            stat.value = mean(stat.value, na.rm = TRUE)) %>% left_join(params) %>% ungroup

gg.dat.all <- gg.dat.all %>% filter(stat.name == "sens_spec") %>% group_by(area, stat.name) %>%
  filter( stat.value >= quantile(stat.value, 0.90), nb.pfg >= quantile(nb.pfg, 0.90)) %>%
  gather_("param.name", "param.value", params.to.test) %>%
  distinct()

summary(gg.dat.all)
sort(table(gg.dat.all$simul.id))

# pdf(file.path(output.dir, paste0("eval_fhd_barplot_top_10pc_nbpfg_sens_spec_sensspec.pdf")), width = 7 , height = 14)

#+ eval_fhd_barplot_top_10pc_nbpfg_sens_spec_sensspec, dev = c('png'), fig.width = 7, fig.height = 14
gg <- ggplot(gg.dat.all, aes(x = param.value, fill = factor(param.value))) + geom_bar() + facet_grid(param.name~stat.name) + 
  scale_fill_discrete(name = "parameter value") + xlab("")
print(gg)
# dev.off()

# pdf(file.path(output.dir, paste0("eval_fhd_boxplot_top_5pc_sens.pdf")), width = 7 , height = 14)

#+ eval_fhd_boxplot_top_5pc_sens, dev = c('png'), fig.width = 7, fig.height = 14
gg <- ggplot(gg.dat.all, aes(param.value, stat.value, fill = factor(param.value))) + geom_boxplot(varwidth = TRUE) + facet_grid(param.name~stat.name) + 
  scale_fill_discrete(name = "parameter value") + xlab("")
print(gg)
# dev.off()

##' identificate the parameters list that lead to the pest simuls

##' According to this exploration I think that the best conbination of varaibles to use are:
##' - "envsuit.option" = 1
##' - "seeding.params" = 2/3
##' - "dispers.mode" = 1
##' - "global.abund" = 1/3
##' - "global.resource.thresh" = 2
##' - "max.by.cohort" = 2/3
##' 

params %>% filter(area == 1, envsuit.option == 1, is.element(seeding.params, c(2,3)), dispers.mode == 1, 
                  is.element(global.abund, c(1, 3)), global.resource.thresh == 2, is.element(max.by.cohort, c(2,3)))

# envsuit.option seeding.params dispers.mode global.abund global.resource.thresh max.by.cohort area simul.id
# 1              1              2            1            1                      2             2    1      483
# 2              1              3            1            1                      2             2    1      485
# 3              1              2            1            3                      2             2    1      531
# 4              1              3            1            3                      2             2    1      533
# 5              1              2            1            1                      2             3    1      867
# 6              1              3            1            1                      2             3    1      869
# 7              1              2            1            3                      2             3    1      915
# 8              1              3            1            3                      2             3    1      917


##' As an extra test we can buid a linear model to check the influence of parameters on our stetistics

gg.dat <- eval.df %>% filter(is.element(year, c("year830", "year840", "year850"))) %>% ## keep only year of interest
  dplyr::select(simul.id, pfg, year, nb.occ, nb.abs, nb.pred.occ, nb.pred.abs, sens, spec, one_of(params.to.test)) %>%  ## remove useless columns
  mutate(tss = sens + spec - 1, pfg_short = sub("_.*$", "", pfg)) 

gg.dat$sens[is.na(gg.dat$sens) & gg.dat$nb.occ > 0 & gg.dat$nb.abs > 0] <- 0 ## set to 0 the stats that are not calculated because of species disparition
gg.dat$spec[is.na(gg.dat$spec) & gg.dat$nb.occ > 0 & gg.dat$nb.abs > 0] <- 0 ## set to 0 the stats that are not calculated because of species disparition
gg.dat$sens_spec[is.na(gg.dat$tss) & gg.dat$nb.occ > 0 & gg.dat$nb.abs > 0] <- 0 ## set to 0 the stats that are not calculated because of species disparition





gg.dat.all <- gg.dat %>% group_by(simul.id, year) %>% 
  summarize(pfg = "all", nb.pfg = sum(nb.pred.occ > 0, na.rm = TRUE),
            sens = mean(sens, na.rm = TRUE),
            spec = mean(spec, na.rm = TRUE),
            tss = mean(tss, na.rm = TRUE)) %>% left_join(params) %>% ungroup %>% data.frame
for(ptt_ in params.to.test){
  gg.dat.all[, ptt_] <- factor(gg.dat.all[, ptt_])
}

head(gg.dat.all)

library(lme4)

lm.spec <- lmer(spec ~ envsuit.option + seeding.params + dispers.mode + global.abund + global.resource.thresh + max.by.cohort + (1|year) + (1|area), data = gg.dat.all)
summary(lm.spec)

lm.sens <- lmer(sens ~ envsuit.option + seeding.params + dispers.mode + global.abund + global.resource.thresh + max.by.cohort + (1|year) + (1|area), data = gg.dat.all)
summary(lm.sens)

lm.tss <- lmer(tss ~ envsuit.option + seeding.params + dispers.mode + global.abund + global.resource.thresh + max.by.cohort + (1|year) + (1|area), data = gg.dat.all)
summary(lm.tss)


lm.nb.pfg <- lmer(nb.pfg ~ envsuit.option + seeding.params + dispers.mode + global.abund + global.resource.thresh + max.by.cohort  + (1|year) + (1|area), data = gg.dat.all)
summary(lm.nb.pfg)

nb.rand <- 5
fixed.var <- c("envsuit.option", "seeding.params", "dispers.mode", "global.abund", "global.resource.thresh", "max.by.cohort")
mod.names <- c("lm.sens", "lm.spec", "lm.tss", "lm.nb.pfg")
names(mod.names) <- c("sensitivity", "specificity", "TSS", "pfg co-exitance")

gg.dat <- NULL

for(mod.id_ in 1:length(mod.names)){
  cat("\n> mod:", mod.id_)
  vi.tab <- matrix(NA, nb.rand, length(fixed.var), dimnames = list(paste0("rand_", 1:nb.rand), fixed.var))
  for(fixed.var_ in fixed.var){
    pred.ref_ <- predict(get(mod.names[mod.id_]), type = 'response', re.form=NA)
    for(nb.rand_ in 1:nb.rand){
      dat_ <- gg.dat.all %>% select_(.dots = fixed.var) %>% 
        mutate_(.dots=setNames(list(lazyeval::interp(~ sample(varname), varname = as.name(fixed.var_))), fixed.var_))
      pred_ <- predict(get(mod.names[mod.id_]), type = 'response', dat_, re.form=NA)
      vi_ <- 1 - abs(cor(pred.ref_, pred_, method = "spearman"))
      vi.tab[nb.rand_, fixed.var_] <- vi_
    }
  }
  
  vi.mean <- colMeans(vi.tab)
  gg.dat <- bind_rows(gg.dat, data.frame(vi = vi.mean, var = factor(names(vi.mean)), group = "models.params", stat = names(mod.names)[mod.id_], stringsAsFactors = FALSE))
}

#+ varimport_radar, dev = c('png'), fig.width = 14, fig.height = 11
gg <- ggplot(gg.dat, aes(x = factor(var), y = vi, group = group)) +
  geom_point() + coord_cartesian(ylim = c(0,1)) +
  facet_wrap(~ stat) +
  # geom_line() + 
  coord_polar() + labs(title = paste0("variable importance"), x = "", y ="") 

print(gg)