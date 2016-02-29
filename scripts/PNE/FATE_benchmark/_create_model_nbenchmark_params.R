##' ----------------------------------------------------------------------------
##' @title Create the parameter files to check the influence of fatehd models
##'   on the simulation outputs
##' @description Here we will create ~ 3 x 5000 parameter files (globalparams, 
##'   simulparams and namespaceparams) that will be used to run the simualation
##'   parameters.
##' @author damien g.
##' @date 23/02/2016
##' @licence GPL-2
##' ----------------------------------------------------------------------------

library(dplyr)

simul.dir <- "/nfs_scratch2/dgeorges/FATE_newHS_devel/SIMUL_6STRATA_TEST/"


## global simul parameters
envsuit.option <- 1:2
seeding.params <- 1:4
dispers.mode <- 1:3
## test global namespace parameters
global.abund <- 1:4
global.resource.thresh <- 1:4
max.by.cohort <- 1:4
## studied area
area <- 1:3

params.df <- expand.grid(envsuit.option = envsuit.option,
                         seeding.params = seeding.params,
                         dispers.mode = dispers.mode,
                         global.abund = global.abund,
                         global.resource.thresh = global.resource.thresh,
                         max.by.cohort = max.by.cohort,
                         area = area)
params.df$simul.id <- 1:nrow(params.df)

## save a copy of the parameters on the hard drive
write.table(params.df, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t", file = "~/params_formal_fhdmpt.txt")

p_ <- params.df[1, , drop = FALSE]
create.params.file <- function(p_){    
  ## deal with global simul param
  if(p_$envsuit.option == 1){
    ENVSUIT_OPTION <- 1
  } else if(p_$envsuit.option == 2){
    ENVSUIT_OPTION <- 2
  } else stop("unknown envsuit.option param")
  
  if(p_$seeding.params == 1){
    SEEDING_TIMESTEP <- 1
    SEEDING_DURATION <- 300
  } else if(p_$seeding.params == 2){
    SEEDING_TIMESTEP <- 10
    SEEDING_DURATION <- 300
  } else if(p_$seeding.params == 3){
    SEEDING_TIMESTEP <- 10
    SEEDING_DURATION <- 500
  } else if(p_$seeding.params == 4){
    SEEDING_TIMESTEP <- 1
    SEEDING_DURATION <- 100
  } else stop("unknown seeding.params param")
  
  if(p_$dispers.mode == 1){
    MODE_DISPERS <- 1
  } else if(p_$dispers.mode == 2){
    MODE_DISPERS <- 2
  } else if(p_$dispers.mode == 3){
    MODE_DISPERS <- 3
  } else stop("unknown dispers.mode param")
  
  if(p_$global.abund == 1){
    GLOBAL_LOW_ABUND <- 3000000
    GLOBAL_MEDIUM_ABUND <- 7000000
    GLOBAL_HIGH_ABUND <- 10000000
  } else if(p_$global.abund == 2){
    GLOBAL_LOW_ABUND <- 6000000
    GLOBAL_MEDIUM_ABUND <- 14000000
    GLOBAL_HIGH_ABUND <- 20000000
  } else if(p_$global.abund == 3){
    GLOBAL_LOW_ABUND <- 1000000
    GLOBAL_MEDIUM_ABUND <- 5000000
    GLOBAL_HIGH_ABUND <- 8000000
  } else if(p_$global.abund == 4){
    GLOBAL_LOW_ABUND <- 5000000
    GLOBAL_MEDIUM_ABUND <- 9000000
    GLOBAL_HIGH_ABUND <- 12000000
  } else stop("unknown global.abund param")
  
  if(p_$global.resource.thresh == 1){
    GLOBAL_LOW_RESOURCES_TRESH <- 9000000
    GLOBAL_MEDIUM_RESOURCES_TRESH <- 6000000
  } else if(p_$global.resource.thresh == 2){
    GLOBAL_LOW_RESOURCES_TRESH <- 19000000
    GLOBAL_MEDIUM_RESOURCES_TRESH <- 13000000
  } else if(p_$global.resource.thresh == 3){
    GLOBAL_LOW_RESOURCES_TRESH <- 7000000
    GLOBAL_MEDIUM_RESOURCES_TRESH <- 4000000
  } else if(p_$global.resource.thresh == 4){
    GLOBAL_LOW_RESOURCES_TRESH <- 11000000
    GLOBAL_MEDIUM_RESOURCES_TRESH <- 8000000
  } else stop("unknown global.resource.thresh param")
  
  if(p_$max.by.cohort == 1){
    GLOBAL_MAX_BY_COHORT <- 7000000
  } else if(p_$max.by.cohort == 2){
    GLOBAL_MAX_BY_COHORT <- 5000000
  } else if(p_$max.by.cohort == 3){
    GLOBAL_MAX_BY_COHORT <- 9000000
  } else if(p_$max.by.cohort == 4){
    GLOBAL_MAX_BY_COHORT <- 14000000
  } else stop("unknown max.by.cohort param")
  
  if(p_$area == 1){
    zone <- 1
  } else if(p_$area == 2){
    zone <- 2
  } else if(p_$area == 3){
    zone <- 3
  } else stop("unknown dispers.mode param")
  
  ## create the params files
  params.dir <- "/nfs_scratch2/dgeorges/FATE_newHS_devel/SIMUL_6STRATA_TEST/DATA"
  simul.rel.dir <- "/home/dgeorges/fhdmpt/SIMUL_6STRATA_TEST"

  ### create global parameters files
  glob.params.dir <- file.path(params.dir, "GLOBAL_PARAMS")
  dir.create(glob.params.dir, showWarnings = FALSE, recursive = TRUE)
  
  cat("## FATEHD global parameters file (automatically generated)\n",
      "## date : ", date(), "\n",
      "NB_CPUS 1", "\n",
      "SUCC_MOD 2", "\n",
      "NB_FG 24", "\n",
      "NB_STRATUM 7", "\n",
      "SEEDING_TIMESTEP ", SEEDING_TIMESTEP, "\n",
      "SEEDING_DURATION ", SEEDING_DURATION, "\n",
      "SIMULATION_TIME 850", "\n",
      "DO_DISTURBANCES 1", "\n",
      "NB_DISTURBANCES 4", "\n",
      "NB_SUBDISTURBANCES 4", "\n",
      "FREQ_DISTURBANCES 1 1 1 1", "\n",
      "DO_SOIL_COMPETITION 0", "\n",
      "DO_FIRE_DISTURBANCES 0", "\n",
      "DO_DROUGHT_DISTURBANCES 0", "\n",
      "DO_HAB_STABILITY 0", "\n",
      "DO_ALIENS_DISTURBANCE 0", "\n",
      "FREQ_ALIENS 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0", "\n",
      "ENVSUIT_OPTION ", ENVSUIT_OPTION, sep = "", 
      file = file.path(glob.params.dir, paste0("Global_parameters_", p_$simul.id, ".txt")))
  
  ### create global parameters files
  namespace.params.dir <- file.path(params.dir, "NAMESPACE_CONSTANTS")
  dir.create(namespace.params.dir, showWarnings = FALSE, recursive = TRUE)
  
  cat("## FATEHD namespace constants parameters file (automatically generated)\n",
      "## date : ", date(), "\n",
      "GLOBAL_LOW_ABUND ", sprintf("%d", GLOBAL_LOW_ABUND), "\n",
      "GLOBAL_MEDIUM_ABUND ", sprintf("%d", GLOBAL_MEDIUM_ABUND), "\n",
      "GLOBAL_HIGH_ABUND ", sprintf("%d", GLOBAL_HIGH_ABUND), "\n",
      "GLOBAL_LOW_RESOURCES_TRESH ", sprintf("%d", GLOBAL_LOW_RESOURCES_TRESH), "\n",
      "GLOBAL_MEDIUM_RESOURCES_TRESH ", sprintf("%d", GLOBAL_MEDIUM_RESOURCES_TRESH), "\n",
      "GLOBAL_FULL_SOIL_COVERAGE 9000000", "\n",
      "GLOBAL_MAX_BY_COHORT ", sprintf("%d", GLOBAL_MAX_BY_COHORT), sep = "", 
      file = file.path(namespace.params.dir, paste0("Namespace_constants_", p_$simul.id, ".txt")))
  
  ### create param simul file
  paramsimul.params.dir <- file.path(params.dir, "../PARAM_SIMUL")
  dir.create(paramsimul.params.dir, showWarnings = FALSE, recursive = TRUE)
  
  cat("## FATEHD simulation parrameter file (automatically generated)\n",
      "## date : ", date(), "\n",
      "--GLOBAL_PARAMS--", "\n",
      file.path(simul.rel.dir, "DATA", "GLOBAL_PARAMS", paste0("Global_parameters_", p_$simul.id, ".txt")), "\n",
      "--NAMESPACE_CONSTANTS--", "\n",
      file.path(simul.rel.dir, "DATA", "NAMESPACE_CONSTANTS", paste0("Namespace_constants_", p_$simul.id, ".txt")), "\n",
      "--SAVE_DIR--", "\n", 
      file.path("RESULTS", p_$simul.id), "\n",
      "--ARRAYS_SAVING_YEARS--", "\n",
      file.path(simul.rel.dir, "DATA", "SAVE", "arraysSavingYears.txt"), "\n",
      "--MASK--", "\n",
      file.path(simul.rel.dir, "DATA", "MASK", paste0("zone", zone), "maskEcrins.asc"), "\n",
      "--PFG_LIFE_HISTORY_PARAMS--", "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "SUCC", paste0("dispmod", MODE_DISPERS), "SUCC_C1_ruderal.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "SUCC", paste0("dispmod", MODE_DISPERS), "SUCC_C2_VaccUli.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "SUCC", paste0("dispmod", MODE_DISPERS), "SUCC_C3_cushion.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "SUCC", paste0("dispmod", MODE_DISPERS), "SUCC_C4_shrubs.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "SUCC", paste0("dispmod", MODE_DISPERS), "SUCC_C5_CallVulg.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "SUCC", paste0("dispmod", MODE_DISPERS), "SUCC_C6_VaccMyrt.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "SUCC", paste0("dispmod", MODE_DISPERS), "SUCC_H1_alpines.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "SUCC", paste0("dispmod", MODE_DISPERS), "SUCC_H10_megaph.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "SUCC", paste0("dispmod", MODE_DISPERS), "SUCC_H2_dryGrass.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "SUCC", paste0("dispmod", MODE_DISPERS), "SUCC_H3_coolGrass.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "SUCC", paste0("dispmod", MODE_DISPERS), "SUCC_H4_undergrowth.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "SUCC", paste0("dispmod", MODE_DISPERS), "SUCC_H5_poorGrass.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "SUCC", paste0("dispmod", MODE_DISPERS), "SUCC_H6_megaph.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "SUCC", paste0("dispmod", MODE_DISPERS), "SUCC_H7_rocks.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "SUCC", paste0("dispmod", MODE_DISPERS), "SUCC_H8_notgraz.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "SUCC", paste0("dispmod", MODE_DISPERS), "SUCC_H9_NardStri.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "SUCC", paste0("dispmod", MODE_DISPERS), "SUCC_P1_ubacPion.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "SUCC", paste0("dispmod", MODE_DISPERS), "SUCC_P2_river.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "SUCC", paste0("dispmod", MODE_DISPERS), "SUCC_P3_ubac2.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "SUCC", paste0("dispmod", MODE_DISPERS), "SUCC_P4_larix.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "SUCC", paste0("dispmod", MODE_DISPERS), "SUCC_P5_alpext.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "SUCC", paste0("dispmod", MODE_DISPERS), "SUCC_P6_alpint.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "SUCC", paste0("dispmod", MODE_DISPERS), "SUCC_P7_acer.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "SUCC", paste0("dispmod", MODE_DISPERS), "SUCC_P8_BetAlb.txt"), "\n",
      "--PFG_DISPERSAL_PARAMS--", "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DISP", "DISP_C1_ruderal.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DISP", "DISP_C2_VaccUli.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DISP", "DISP_C3_cushion.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DISP", "DISP_C4_shrubs.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DISP", "DISP_C5_CallVulg.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DISP", "DISP_C6_VaccMyrt.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DISP", "DISP_H1_alpines.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DISP", "DISP_H10_megaph.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DISP", "DISP_H2_dryGrass.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DISP", "DISP_H3_coolGrass.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DISP", "DISP_H4_undergrowth.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DISP", "DISP_H5_poorGrass.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DISP", "DISP_H6_megaph.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DISP", "DISP_H7_rocks.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DISP", "DISP_H8_notgraz.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DISP", "DISP_H9_NardStri.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DISP", "DISP_P1_ubacPion.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DISP", "DISP_P2_river.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DISP", "DISP_P3_ubac2.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DISP", "DISP_P4_larix.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DISP", "DISP_P5_alpext.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DISP", "DISP_P6_alpint.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DISP", "DISP_P7_acer.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DISP", "DISP_P8_BetAlb.txt"), "\n",
      "--PFG_ENVSUIT--", "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "ENVSUIT", "_SCALED_SENS_95", paste0("zone", zone), "HS_f0_C1.asc"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "ENVSUIT", "_SCALED_SENS_95", paste0("zone", zone), "HS_f0_C2.asc"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "ENVSUIT", "_SCALED_SENS_95", paste0("zone", zone), "HS_f0_C3.asc"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "ENVSUIT", "_SCALED_SENS_95", paste0("zone", zone), "HS_f0_C4.asc"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "ENVSUIT", "_SCALED_SENS_95", paste0("zone", zone), "HS_f0_C5.asc"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "ENVSUIT", "_SCALED_SENS_95", paste0("zone", zone), "HS_f0_C6.asc"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "ENVSUIT", "_SCALED_SENS_95", paste0("zone", zone), "HS_f0_H1.asc"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "ENVSUIT", "_SCALED_SENS_95", paste0("zone", zone), "HS_f0_H10.asc"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "ENVSUIT", "_SCALED_SENS_95", paste0("zone", zone), "HS_f0_H2.asc"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "ENVSUIT", "_SCALED_SENS_95", paste0("zone", zone), "HS_f0_H3.asc"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "ENVSUIT", "_SCALED_SENS_95", paste0("zone", zone), "HS_f0_H4.asc"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "ENVSUIT", "_SCALED_SENS_95", paste0("zone", zone), "HS_f0_H5.asc"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "ENVSUIT", "_SCALED_SENS_95", paste0("zone", zone), "HS_f0_H6.asc"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "ENVSUIT", "_SCALED_SENS_95", paste0("zone", zone), "HS_f0_H7.asc"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "ENVSUIT", "_SCALED_SENS_95", paste0("zone", zone), "HS_f0_H8.asc"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "ENVSUIT", "_SCALED_SENS_95", paste0("zone", zone), "HS_f0_H9.asc"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "ENVSUIT", "_SCALED_SENS_95", paste0("zone", zone), "HS_f0_P1.asc"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "ENVSUIT", "_SCALED_SENS_95", paste0("zone", zone), "HS_f0_P2.asc"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "ENVSUIT", "_SCALED_SENS_95", paste0("zone", zone), "HS_f0_P3.asc"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "ENVSUIT", "_SCALED_SENS_95", paste0("zone", zone), "HS_f0_P4.asc"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "ENVSUIT", "_SCALED_SENS_95", paste0("zone", zone), "HS_f0_P5.asc"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "ENVSUIT", "_SCALED_SENS_95", paste0("zone", zone), "HS_f0_P6.asc"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "ENVSUIT", "_SCALED_SENS_95", paste0("zone", zone), "HS_f0_P7.asc"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "ENVSUIT", "_SCALED_SENS_95", paste0("zone", zone), "HS_f0_P8.asc"), "\n",
      "--PFG_DISTURBANCES_PARAMS--", "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DIST", "DIST_C1_ruderal.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DIST", "DIST_C2_VaccUli.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DIST", "DIST_C3_cushion.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DIST", "DIST_C4_shrubs.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DIST", "DIST_C5_CallVulg.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DIST", "DIST_C6_VaccMyrt.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DIST", "DIST_H1_alpines.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DIST", "DIST_H10_megaph.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DIST", "DIST_H2_dryGrass.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DIST", "DIST_H3_coolGrass.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DIST", "DIST_H4_undergrowth.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DIST", "DIST_H5_poorGrass.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DIST", "DIST_H6_megaph.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DIST", "DIST_H7_rocks.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DIST", "DIST_H8_notgraz.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DIST", "DIST_H9_NardStri.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DIST", "DIST_P1_ubacPion.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DIST", "DIST_P2_river.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DIST", "DIST_P3_ubac2.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DIST", "DIST_P4_larix.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DIST", "DIST_P5_alpext.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DIST", "DIST_P6_alpint.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DIST", "DIST_P7_acer.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "PFGS", "DIST", "DIST_P8_BetAlb.txt"), "\n",
      "--DIST_MASK--", "\n",
      file.path(simul.rel.dir, "DATA", "MASK", paste0("zone", zone), "noDisturb.asc"), "\n",
      file.path(simul.rel.dir, "DATA", "MASK", paste0("zone", zone), "noDisturb.asc"), "\n",
      file.path(simul.rel.dir, "DATA", "MASK", paste0("zone", zone), "noDisturb.asc"), "\n",
      file.path(simul.rel.dir, "DATA", "MASK", paste0("zone", zone), "noDisturb.asc"), "\n",
      "--DIST_CHANGE_TIME--", "\n",
      file.path(simul.rel.dir, "DATA", "SCENARIO", "disturbances_changing_times.txt"), "\n",
      "--DIST_CHANGE_MASK--", "\n",
      file.path(simul.rel.dir, "DATA", "SCENARIO", paste0("zone", zone), "disturbances_changing_masks_t600.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "SCENARIO", paste0("zone", zone), "disturbances_changing_masks_t601.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "SCENARIO", paste0("zone", zone), "disturbances_changing_masks_t800.txt"), "\n",
      file.path(simul.rel.dir, "DATA", "SCENARIO", paste0("zone", zone), "disturbances_changing_masks_t801.txt"), "\n",
      "--END_OF_FILE--", sep = "", 
      file = file.path(paramsimul.params.dir, paste0("ParamsSimul_", p_$simul.id, ".txt")))
  
  return(TRUE)
}



cpf <- params.df %>% rowwise %>% do(status = data.frame(create.params.file(.)))

##' @note trick to change a param in a list of files (e.g. the dispersal mode)
##' find dispmod2 -name 'SUCC_*.txt' -type f -exec sed -i 's/MODE_DISPERS 1/MODE_DISPERS 2/g' {} \;

##' create the array saving time param file
saving.times <- c(100, 150, 300, 350, 600, 650, 800, 810, 820, 830, 840, 850)
cat(saving.times, sep = "\n",
    file = file.path(simul.dir, "DATA", "SAVE", "arraysSavingYears.txt"))

##' find ref parameters ids
params.df.ref <- params.df %>% filter(envsuit.option == 1,
                                      seeding.params == 1,
                                      dispers.mode == 1,
                                      global.abund == 1,
                                      global.resource.thresh == 1,
                                      max.by.cohort == 1) 
params.df.ref

##' create the grid parameter files
grid.params <- params.df %>% transmute(job.id = simul.id, 
                                       path.to.param.test = paste0("SIMUL_6STRATA_TEST/PARAM_SIMUL/ParamsSimul_", simul.id, ".txt"),
                                       path.to.param.ref = paste0("SIMUL_6STRATA_REF/PARAM_SIMUL/ParamsSimul_ref", area, ".txt"))

write.table(grid.params, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = " ", file = "~/params_fhdmpt.txt")
