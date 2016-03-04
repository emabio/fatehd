# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
# file : create_disturbances_parameters.R
# authors : Ceres B. Isabelle B. Damien G. and Maya G. 
# date : 17-12-2014
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #

# Description =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
#   this script is an example to show how we can produce
#   FATE-HD disturbances parameters files from a subset of 
#   parameters and some systematic rules. It comes whith
#   the 2nd FATE-HD example "Succession, Dispersal and Disturbances".
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #

## clear the workspace and set the right working directory 
rm(list=ls())

(load("~/Work/FATEHD/benchmarking/workdir/pfg_benchmark_obj.RData"))

params.list <- commandArgs(trailingOnly=TRUE)
params.val <- c("REF", "TEST")[as.numeric(params.list[1:lenght(params.names)])]
names(params.val) <- params.names

simul_folder = "~/Work/FATEHD/benchmarking/workdir/FATE_6STRATA_TEST_PFG/"
dir.create(simul_folder)
simul_name = basename(simul_folder)
setwd(simul_folder)

#################################################################################################
# PARAMETERS OF THE SUCCESSION (FATE)
#################################################################################################

#################################################################################################
### INPUT FILE
( dat = read.table("editParamsC24_newH8.csv", sep = ";", h = T) )

### SUCCESSION FILE
( sucParams = read.table(file.path(simul_folder,simul_name,"/DATA/PFGS/SUCC/successionParamsTable.csv"),sep='\t', h=T, row.names=1, skip=1) )

#################################################################################################
### NB OF PFG, NAMES and OUTPUT DIR CREATION
n <- nrow(dat) # number of groups
sp <- dat$group # groups' names


## change the parameters value according to case we study
## REF: no change
## TEST: palabity cat <- palabity cat +/- 1

for(p_ in grep("DIST", names(params.val), value = TRUE)){
  if(params.val[p_] == "TEST"){
    elt_ <- which(is.element(dat$group, names(pfg.group)[pfg.group == sub("^.*[[:punct:]]", "", p_)]))
    for(elt__ in elt_){
      if(dat$palatability[elt_] <= 2) dat$palatability[elt_] <- dat$palatability[elt_] + 1 else dat$palatability[elt_] <- dat$palatability[elt_] - 1
    }
  }
}

dist_path = paste(simul_name,"/DATA/PFGS/DIST",sep="")
dir.create(dist_path, showWarnings = F, recursive = T)  # create the directory if does not exist

## initialize a list where all parameters will be stored
param_store_df <- data.frame(matrix(NA,0,n))
colnames(param_store_df) <- sp
 
###############################################################################################
### NUMBER OF DISTURBANCES
## 1: mowing
## 2: grazing intensity 1
## 3: grazing intensity 2
## 4: grazing intensity 3
dist_names <- c("mowing", "graz1", "graz2", "graz3")
( NDIST = length(dist_names) )

###############################################################################################
### PROPORTION OF KILLED PROPAGULES
PROP_KILLED <- param_store_df
PROP_KILLED[1:NDIST,] <- 0 # here 0 for all perturbations and all PFG
rownames(PROP_KILLED) <- paste("PROP_KILLED", "_", dist_names, sep="")

###############################################################################################
### RESPONSES' STAGES
nStages <- 4 # number of responses' stages for each PFG
BREAK_AGE <- param_store_df
BREAK_AGE[1: (NDIST * (nStages -1)), ] <- 0
rownames(BREAK_AGE) <- paste("BREAK_AGE", "_", rep(dist_names, each = nStages-1), "_", paste(1:(nStages-1), "to", 2:(nStages), sep=""), sep="" )

## A12 = maturity-2 for herbaceous, 1 year otherwise
## A23 = min(age12, maturity)
## A34 = min(age12, longevity-2)

### in that case :
## herbaceous: 0 ---> maturity-2  ---> maturity   ---> longevity-2  ---> longevity
## others:     0 ---> 1 year      --->  age12     ---> age12        ---> longevity
##             0 ---> 1 year      ---> maturity   ---> longevity-2  ---> longevity

brk_ages_tmp <- NULL
brk_ages_tmp <- rbind( brk_ages_tmp, 
                       ifelse(dat$type=='H', dat$maturity-2, 1),
                       apply(rbind(as.numeric(sucParams['CHANG_STR_AGES_to_str2',]), as.numeric(sucParams['MATURITY',])),2 , min),
                       apply(rbind(as.numeric(sucParams['CHANG_STR_AGES_to_str2',]), as.numeric(sucParams['LONGEVITY',])-2),2 , min) )
for( i in 1:NDIST) {
  BREAK_AGE[(1 + ((i-1)*(nStages-1))):((nStages-1) + ((i-1)*(nStages-1))) ,] <- brk_ages_tmp
}


###############################################################################################
### RESPROUT AGES
RESPR_AGE <- param_store_df
RESPR_AGE[1: (NDIST * nStages ), ] <- 0
rownames(RESPR_AGE) <- paste("RESPR_AGE", "_", rep(dist_names, each = nStages), "_", paste(1:nStages, sep=""), sep="" )

## stage 1
RESPR_AGE[seq(1,nrow(RESPR_AGE), by=nStages),] <- 0
## stage 2 
RESPR_AGE[seq(2,nrow(RESPR_AGE), by=nStages),] <- 0
## stage 3 : matures resprout at maturity-2; juveniles not affected
RESPR_AGE[seq(3,nrow(RESPR_AGE), by=nStages),] <- rep( apply( rbind(sucParams['MATURITY',]-2, sucParams['CHANG_STR_AGES_to_str2',]), 2, min), each = NDIST)
## stage 4 : resprout at longevity-2
RESPR_AGE[seq(4,nrow(RESPR_AGE), by=nStages),] <- unlist(rep( sucParams['LONGEVITY',] - 2, each = NDIST))

###############################################################################################
### FATES : PROPORTION OF KILLED or RESPROUTING INDIVIDUALS
## 0 : 0%
## 1 : 10%
## 2 : 50%
## 3 : 90%
## 4 : 100%
## 5 : 40%
## 6 : 80%

FATES <- param_store_df
FATES[1: (NDIST * nStages * 2 ), ] <- 0 # lines are pairs (prop of individuals killed and resprouting) for each stage of each disturbance
rownames(FATES) <- paste("FATES", "_", rep(dist_names, each = nStages*2),  "_", rep(1:nStages, each=2), "_" ,c("kill", "respr"), sep="" )

## THIS PARAMETER HAS TO BE FILLED BY HAND !

## MOWING #################################
### AGE 1 : trees 80% killed 0% resprout
FATES["FATES_mowing_1_kill", dat$type=="P"] <- 6
### AGE 2 : trees and bushes 100% killed 0% resprout
FATES["FATES_mowing_2_kill", dat$type=="P" | dat$type=="C"] <- 4
### AGE 3 : trees 100% killed 0% resprout / bushes 50% killed 50% resprout / herbs specific rules
FATES["FATES_mowing_3_kill", dat$type=="P"] <- 4
FATES["FATES_mowing_3_kill", dat$type=="C"] <- 2
FATES["FATES_mowing_3_respr", dat$type=="C"] <- 2
FATES["FATES_mowing_3_kill", "H2_dryGrass"] <- 5
FATES["FATES_mowing_3_respr", "H2_dryGrass"] <- 2
FATES["FATES_mowing_3_kill", "H9_NardStri"] <- 0
FATES["FATES_mowing_3_respr", "H9_NardStri"] <- 3
### AGE 4 : everybody is killed
FATES["FATES_mowing_4_kill", ] <- 4

## GRAZING INTENSITY 1 ####################
### AGE 1 : pfgs are killed from 0% to 10% according to their palability index / trees are all killed
FATES["FATES_graz1_1_kill", ] <- c(0, 0, 1, 1)[dat$palatability + 1]  
FATES["FATES_graz1_1_kill", dat$type=="P"] <- 4
### AGE 2 : pfgs are killed from 0% to 10% according to their palability index / trees are unaffected
FATES["FATES_graz1_2_kill", ] <- c(0, 0, 1, 1)[dat$palatability + 1]  
FATES["FATES_graz1_2_kill", dat$type=="P"] <- 0
### AGE 3 : pfgs resprout from 0% to 50% according to their palability index / trees are unaffected
FATES["FATES_graz1_3_respr", ] <- c(0, 0, 1, 2)[dat$palatability + 1]  
FATES["FATES_graz1_3_respr", dat$type=="P"] <- 0
### AGE 4 : pfgs resprout from 0% to 50% according to their palability index / trees are unaffected
FATES["FATES_graz1_4_respr", ] <- c(0, 0, 1, 1)[dat$palatability + 1]  
FATES["FATES_graz1_4_respr", dat$type=="P"] <- 0

## GRAZING INTENSITY 2 ####################
### AGE 1 : pfgs are killed from 0% to 50% according to their palability index / trees are all killed
FATES["FATES_graz2_1_kill", ] <- c(0, 0, 1, 2)[dat$palatability + 1]  
FATES["FATES_graz2_1_kill", dat$type=="P"] <- 4
### AGE 2 : pfgs are killed from 0% to 50% according to their palability index / P4 are 50% killed / other trees are unaffected
FATES["FATES_graz2_2_kill", ] <- c(0, 0, 1, 2)[dat$palatability + 1]  
FATES["FATES_graz2_2_kill", dat$type=="P"] <- 0
FATES["FATES_graz2_2_kill", "P4_larix"] <- 2
### AGE 3 : pfgs resprout from 0% to 100% according to their palability index / trees are unaffected
FATES["FATES_graz2_3_respr", ] <- c(0, 1, 3, 4)[dat$palatability + 1]  
FATES["FATES_graz2_3_respr", dat$type=="P"] <- 0
### AGE 4 : pfgs resprout from 0% to 50% and are killed from 0% to 10% according to their palability index / trees are unaffected
FATES["FATES_graz2_4_respr", ] <- c(0, 0, 2, 2)[dat$palatability + 1]  
FATES["FATES_graz2_4_respr", dat$type=="P"] <- 0
FATES["FATES_graz2_4_kill", ] <- c(0, 0, 1, 1)[dat$palatability + 1]  
FATES["FATES_graz2_4_kill", dat$type=="P"] <- 0

## GRAZING INTENSITY 3 ####################
### AGE 1 : pfgs are killed from 0% to 90% according to their palability index / trees are all killed
FATES["FATES_graz3_1_kill", ] <- c(0, 0, 2, 3)[dat$palatability + 1]  
FATES["FATES_graz3_1_kill", dat$type=="P"] <- 4
### AGE 2 : pfgs are killed from 0% to 90% according to their palability index / trees are killed at 40%
FATES["FATES_graz3_2_kill", ] <- c(0, 0, 2, 3)[dat$palatability + 1]  
FATES["FATES_graz3_2_kill", dat$type=="P"] <- 5
### AGE 3 : pfgs resprout from 0% to 100% and are killed from 0% to 10% according to their palability index / trees are unaffected
FATES["FATES_graz3_3_respr", ] <- c(0, 2, 4, 3)[dat$palatability + 1]  
FATES["FATES_graz3_3_respr", dat$type=="P"] <- 0
FATES["FATES_graz3_3_kill", ] <- c(0, 0, 0, 1)[dat$palatability + 1]  
FATES["FATES_graz3_3_kill", dat$type=="P"] <- 0
### AGE 4 : pfgs resprout from 0% to 50% and are killed from 0% to 50% according to their palability index / trees are unaffected
FATES["FATES_graz3_4_respr", ] <- c(0, 1, 2, 2)[dat$palatability + 1]  
FATES["FATES_graz3_4_respr", dat$type=="P"] <- 0
FATES["FATES_graz3_4_kill", ] <- c(0, 0, 1, 2)[dat$palatability + 1]  
FATES["FATES_graz3_4_kill", dat$type=="P"] <- 0

###############################################################################################
### END OF SEED DORMANCY : % of seeds activated by the perturbation
## 0 : 0%
## 1 : 10%
## 2 : 50%
## 3 : 90%
## 4 : 100%
## 5 : 40%
## 6 : 80%
ACTIVATED_SEED <- param_store_df
ACTIVATED_SEED[1:NDIST, ] <- 0
rownames(ACTIVATED_SEED) <- paste("ACTIVATED_SEED", "_", dist_names, sep="")

#################################################################################################
# writing parameter files
#################################################################################################

## list of all supported parameters

available_param_list <- c("PROP_KILLED",
                          "BREAK_AGE",
                          "RESPR_AGE",
                          "FATES",
                          "ACTIVATED_SEED")

## create empty params files for all pfgs
dist_file_names <- file.path(dist_path, paste("DIST_",sp,".txt", sep=""))
names(dist_file_names) <- sp
## remove old files if exists
file.remove(dist_file_names)
## create the new ones
file.create(dist_file_names, showWarnings = F)

## create a dataframe where all parameters will be stored
all_params_table <- param_store_df

for(parm in available_param_list){
  if(exists(parm)){
    ## add params to global param file
    all_params_table <- rbind(all_params_table, get(parm))
    ## add params in each pfg param file
    for(spp in sp){
      line_to_add <- paste(parm, " ", paste( get(parm)[,spp], collapse = " " ), "\n", sep="")
      cat(line_to_add , 
          file = dist_file_names[spp],
          append = TRUE)
    } ## end loop on species
  }
} ## end loop on parameters

## write also the summary parameters table
outFile= file.path(dist_path,"disturbancesParamsTable.csv")
write.table(all_params_table,file=outFile, quote=FALSE, row.names=TRUE, col.names= NA, sep="\t")

