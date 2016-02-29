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

rm(list = ls())

output.dir <- "~/fhdmpt_simul_comparaisons_merged"
dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)

##' execute the scp command to retrieave results from all cigri clusters

## copy results from luke
system("cp -r ~/fhdmpt_simul_comparaisons/* ~/fhdmpt_simul_comparaisons_merged/")

## copy results from other clusters
clust.names <- c("froggy", "gofree", "ceciccluster", "fontaine")
for(cn_ in clust.names){
  cat("\n> getting", cn_, "results...")
  system(paste0("scp -r ", cn_, ":fhdmpt_simul_comparaisons/* ~/fhdmpt_simul_comparaisons_merged/"))
}

## get input parameters file
params <- read.table("~/params_fhdmpt.txt", header = FALSE, sep = " ", stringsAsFactors = FALSE)
head(params)
