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

simul.basedir <- "~/Work/FATEHD/benchmarking/fatehd_params_test/grid/model_full_params_test/fhdmfpt" 
params.file <- file.path(simul.basedir, "Data/params.fhdmfpt.txt")
dir.create(dirname(params.file), showWarnings = FALSE, recursive = TRUE)

##' # Create the parameter file

create.param.line <- function(){
  ## global_xx_abund
  var.temp <- round(sort(runif(3, 0, 20000000)))
  global.low.abund <- var.temp[1]
  global.medium.abund <- var.temp[2]
  global.high.abund <- var.temp[3]
  
  ## global_xx_resource_treshold
  var.temp <- round(sort(runif(2, 0, 20000000)))
  global.low.resources.tresh <- var.temp[1]
  global.medium.resources.tresh <- var.temp[2]
  
  ## global_max_by_cohort
  var.temp <- round(sort(runif(1, 0, 20000000)))
  global.max.by.cohort <- var.temp[1]
  
  ## global_full_soil_coverage
  var.temp <- round(sort(runif(1, 0, 20000000)))
  global.full.soil.coverage <- var.temp[1]
  
  ## envsuit_option
  var.temp <- sample(1:2, 1)
  envsuit.option <- var.temp[1]
  
  ## seeding_timestep
  var.temp <- sample(1:50, 1)
  seeding.timestep <- var.temp[1]
  
  ## seeding_duration
  var.temp <- sample(0:850, 1)
  seeding.duration <- var.temp[1]
  
  ## soil_default_value
  var.temp <- round(sort(runif(1, 0, 50)), digits = 2)
  soil.default.value <- var.temp[1]
  
  ## soil_percent
  var.temp <- round(sort(runif(1, 0, 1)), digits = 2)
  soil.percent <- var.temp[1]
  
  ## soil_categories_threshold
  var.temp <- round(sort(runif(2, 0, 100)), digits = 2)
  soil.categories.threshold.low <- var.temp[1]
  soil.categories.threshold.medium <- var.temp[2]
  
  ## disturbance_changing_times
  var.temp <- sort(sample(0:850, 2))
  disturbance.changing.times.1 <- var.temp[1]
  disturbance.changing.times.2 <- var.temp[2]
  
  ## hs_type
  var.temp <- sample(1:4, 1)
  hs.type <- var.temp[1]
  
  ## area
  var.temp <- sample(1:3, 1)
  area <- var.temp[1]
  
  return(data.frame(
    global.low.abund = global.low.abund,
    global.medium.abund = global.medium.abund,
    global.high.abund = global.high.abund,  
    global.low.resources.tresh = global.low.resources.tresh,
    global.medium.resources.tresh = global.medium.resources.tresh,
    global.max.by.cohort = global.max.by.cohort,
    global.full.soil.coverage = global.full.soil.coverage,
    envsuit.option = envsuit.option,
    seeding.timestep = seeding.timestep,
    seeding.duration = seeding.duration,
    soil.default.value = soil.default.value,
    soil.percent = soil.percent,
    soil.categories.threshold.low = soil.categories.threshold.low,
    soil.categories.threshold.medium = soil.categories.threshold.medium,
    disturbance.changing.times.1 = disturbance.changing.times.1,
    disturbance.changing.times.2 = disturbance.changing.times.2,
    hs.type = hs.type,
    area = area))
}

params <- bind_rows(lapply(1:10000, function(i) create.param.line() ))
params$simul.id <- 1:nrow(params)

write.table(params, params.file, sep ="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
