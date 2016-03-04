# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
# file : create_succession_parameters.R
# authors : Ceres B. Isabelle B. Damien G. and Maya G. 
# date : 06-01-2015
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #

# Description =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
#   Ptoduce FATE-HD succession parameters files from a subset of 
#   parameters and some systematic rules. 
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #


## clear the workspace and set the right working directory 
rm(list=ls())
dat.dir <- "~/fhdpfgpt/Data"
(load(file.path(dat.dir, "pfg_benchmark_obj.RData")))

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
( dat = read.table(file.path(dat.dir, "editParamsC24_newH8.csv"), sep = ";", h = T) )

#################################################################################################
### GLOBAL PARAMETERS DEFINITION

## HEIGHT STRATA DEFINITION (for light competition and PFG growth)
# 5 strata (+ germinants = 0)
# ! this parameter has been set up in fate-hd global parameters file
nb_height_strat = 6

## STRATA HEIGHT LIMITS
# Hs01 = 0  ## must always be 0
# Hs12 = 150
# Hs23 = 400
# Hs34 = 1000
# Hs45 = 2000
#height_strat_limits <- c(0,150,400,1000,2000)
height_strat_limits <- c(0,70,150,400,1000,2000)

## SOIL STRATA DEFINITION (for soil resources competition) ==> OPTIONAL
nb_soil_strat = 3 ## n + 2 soil_strat_ treshold because ve consider 
                  ## also the strat before first soil threshold and 
                  ## the one after last soil threshold
# ! this parameter has been set up in fate-hd global parameters file

# strata limits
soil_strat_limits <- c(17, 28) ## starta 'less than 22' and 'more than 27' will be added automatically within FATEHD

#################################################################################################
### NB OF PFG, NAMES and OUTPUT DIR CREATION
n <- nrow(dat) # number of groups
sp <- dat$group # groups' names

succ_path = paste(simul_name,"/DATA/PFGS/SUCC",sep="")
dir.create(succ_path, showWarnings = F, recursive = T)  # create the directory if does not exist

## initialize a list where all parameters will be stored
param_store_df <- data.frame(matrix(NA,0,n))
colnames(param_store_df) <- sp

## add max achieved stratum info in dat according to height_strat_limits
dat$strata <- sapply(dat$height, function(h){ max( c(0,which( height_strat_limits < h ), na.rm=T) ) })

#################################################################################################
### PFG NAMES
NAME <- param_store_df
NAME[1,] <-  dat$group
rownames(NAME) <- "NAME"

#################################################################################################
### MATURITY AGE
MATURITY <- param_store_df
MATURITY[1,] <-  dat$maturity
rownames(MATURITY) <- "MATURITY"

#################################################################################################
### LONGEVITY
## death precedes seed productivity in the model thus longevity param = longevity+1
LONGEVITY <- param_store_df
LONGEVITY[1,] <-  dat$longevity + 1
rownames(LONGEVITY) <- "LONGEVITY"

#################################################################################################
### MAX ABUNDANCE
## can be understood as the maximum shade a PFT can make. 
## we decided to put more for species spreading in several strata
## 1: 3ex ; 2: 7ex ; 3: 10ex 
MAX_ABUNDANCE <- param_store_df 
MAX_ABUNDANCE[1,] <- 3
rownames(MAX_ABUNDANCE) <- "MAX_ABUNDANCE"

MAX_ABUNDANCE[1, which(dat$type=="H")] <- 1 # also because of the phenology
MAX_ABUNDANCE[1, which(dat$type=="C" & dat$strata==1)] <- 1 # herbaceous Chamaephyts same as herbaceous
MAX_ABUNDANCE[1, which(dat$type=="C" & dat$strata>1)] <- 2 # Chamaephyts type shrub: intermediate

## change the parameters value according to case we study
## REF: no change
## TEST: H -> MAX_ABUNDANCE + 1
##       C -> MAX_ABUNDANCE + 1
##       P -> MAX_ABUNDANCE - 1

for(p_ in grep("TAS", names(params.val), value = TRUE)){
  if(params.val[p_] == "TEST"){
    elt_ <- is.element(dat$group, names(pfg.group)[pfg.group == sub("^.*[[:punct:]]", "", p_)])
    ## H -> MAX_ABUNDANCE + 1; C -> MAX_ABUNDANCE + 1
    if(grepl("[[:punct:]][H,C]", p_)) {
      MAX_ABUNDANCE[1, elt_] <- MAX_ABUNDANCE[1, elt_] + 1
    }
    ## P -> MAX_ABUNDANCE - 1
    if(grepl("[[:punct:]][P]", p_)) {
      MAX_ABUNDANCE[1, elt_] <- MAX_ABUNDANCE[1, elt_] - 1
    }
  }
}



#################################################################################################
### RATIO IMMATURE SIZE / MATURE SIZE
## type H: 100% (4)
## type C et Strata==1: 100% (4)
## type C et Strata>1: 50% (2)
## type P et Hmax<10m: 50% (2)
## type P et Hmax>10m: 10% (1)
## 80% (3) ?

IMM_SIZE <- param_store_df 
IMM_SIZE[1,] <- 4
rownames(IMM_SIZE) <- "IMM_SIZE"

IMM_SIZE[1, dat$type == "H"] <- 4
IMM_SIZE[1, dat$type == "H" & dat$strata > 1 ] <- 3
IMM_SIZE[1, dat$type == "C" & dat$strata == 1 ] <- 4
IMM_SIZE[1, dat$type == "C" & dat$strata > 1 ] <- 2
IMM_SIZE[1, dat$type == "P" & dat$height < 1000 ] <- 2
IMM_SIZE[1, dat$type == "P" & dat$height >= 1000 ] <- 1

## change the parameters value according to case we study
## REF: no change
## TEST: H -> IMM_SIZE - 1
##       C -> IMM_SIZE - 1
##       P -> IMM_SIZE + 1

for(p_ in grep("TAS", names(params.val), value = TRUE)){
  if(params.val[p_] == "TEST"){
    elt_ <- is.element(dat$group, names(pfg.group)[pfg.group == sub("^.*[[:punct:]]", "", p_)])
    ## H -> IMM_SIZE - 1; C -> IMM_SIZE - 1
    if(grepl("[[:punct:]][H,C]", p_)) {
      IMM_SIZE[1, elt_] <- IMM_SIZE[1, elt_] - 1
    }
    ## P -> IMM_SIZE + 1
    if(grepl("[[:punct:]][P]", p_)) {
      IMM_SIZE[1, elt_] <- IMM_SIZE[1, elt_] + 1
    }
  }
}


#################################################################################################
### GROWTH : AGES TO MOVE FROM ONE STRATUM TO ANOTHER
CHANG_STR_AGES <- param_store_df 
## initialize transition ages between strata 0 and 1 (end of "germinants")
CHANG_STR_AGES[1,] <- 0 # 0 for every PFT
## other transitions (10000 = sometimes never reach in simulations)
CHANG_STR_AGES[2:nb_height_strat,] <- 10000
rownames(CHANG_STR_AGES) <- paste("CHANG_STR_AGES", "_to_str", 1:nb_height_strat, sep="")

## equation (logistic growth curve with 2 points to parameterize it)
for (i in 1:n) { 
	if(IMM_SIZE[1,i]==1) k = 0.1
	if(IMM_SIZE[1,i]==2) k = 0.5
	if(IMM_SIZE[1,i]==3) k = 0.9
	
# 	if(!(IMM_SIZE[1,i]==4)){ # no curve for herbaceous 
	if(dat$strata[i]>1){ # no curve for herbaceous 
	  # at age=maturity/2, height= IMM_SIZE*height	
	  # at age=longevity, height= height
	  k = -log(1-k)/(dat$maturity[i]/2) 
	  A=1:dat$longevity[i]
	  
	  # negative binomiale curve
	  H=dat$height[i]*(1-exp(-k*A))

	  # calculation of transition ages depending on strata heights
	  for(str in 2:nb_height_strat ){
	    CHANG_STR_AGES[str,i] <- ifelse(is.na(A[which(H>=height_strat_limits[str])][1]), CHANG_STR_AGES[str,i], A[which(H>=height_strat_limits[str])][1])
	  }
	} else if(dat$height[i]>height_strat_limits[2]) CHANG_STR_AGES[2,i] <- 0 # for herbaceous, direct into strata max
}

#################################################################################################
### DISPERSAL : is PFG widely dispersed ?
## 0 = no
## 1 = yes
WIDE_DISPERS <- param_store_df 
WIDE_DISPERS[1,] <- 0 # dispersal module is separated
rownames(WIDE_DISPERS) <- "WIDE_DISPERS"

### DISPERSAL MODULE
## 0 = no dispersal
## 1 = dispersal everywhere
## 2 = negative exponential kernel
## 3 = negative exponential kernel + probability decreasing with distance
MODE_DISPERS <- param_store_df 
MODE_DISPERS[1,] <- 1
rownames(MODE_DISPERS) <- "MODE_DISPERS"

#################################################################################################
### GERMINATION RATE DEPENDING ON LIGHT CONDITIONS
## 3 light conditions => 1 = Low ; 2 = Medium ; 3 = High
## thus 3 parameters
## these rates should express a deviation from the germination rate in optimal conditions (=100%)
ACTIVE_GERM <- param_store_df 
ACTIVE_GERM[1:3,] <- NA
rownames(ACTIVE_GERM) <- paste("ACTIVE_GERM", "_with_", c("low", "medium", "high"), sep="")

## N.B. 2= 50% ; 3= 90% ; 4= 100%
## other possibilities:  
## 0 : 0%
## 1 : 10%
## 2 : 50%
## 3 : 90%
## 4 : 100%
## 5 : 40%
## 6 : 80%

## woody species have little variation in germination rate depending on light conditions
ACTIVE_GERM[1, dat$type=="P" | dat$type=="C"] = 3 ## low light conditions
ACTIVE_GERM[2, dat$type=="P" | dat$type=="C"] = 3 ## medium light conditions
ACTIVE_GERM[3, dat$type=="P" | dat$type=="C"] = 3 ## high light conditions

## herbaceous germinate less in the shadow
ACTIVE_GERM[1, dat$type=="H"] = 2 ## low light conditions
ACTIVE_GERM[2, dat$type=="H"] = 6 ## medium light conditions
ACTIVE_GERM[3, dat$type=="H"] = 3 ## high light conditions

## change the parameters value according to case we study
## REF: no change
## TEST: H -> ACTIVE_GERM  <- 3,3,3
##       C -> ACTIVE_GERM <- 2,6,3
##       P -> ACTIVE_GERM <- 2,6,3

for(p_ in grep("LIGHT", names(params.val), value = TRUE)){
  if(params.val[p_] == "TEST"){
    elt_ <- is.element(dat$group, names(pfg.group)[pfg.group == sub("^.*[[:punct:]]", "", p_)])
    ## H -> ACTIVE_GERM  <- 3,3,3
    if(grepl("[[:punct:]][H]", p_)) {
      ACTIVE_GERM[1, elt_] <- ACTIVE_GERM[2, elt_] <- ACTIVE_GERM[3, elt_] <- 3
    }
    ## P -> MAX_ABUNDANCE - 1
    if(grepl("[[:punct:]][C,P]", p_)) {
      ACTIVE_GERM[1, elt_] <- 2
      ACTIVE_GERM[2, elt_] <- 6
      ACTIVE_GERM[3, elt_] <- 3
    }
  }
}

#################################################################################################
### SHADE TOLERANCE
## 0 = non tolerant
## 1 = tolerant
## 3 light conditions (L = Low ; M = Medium ; H = High)
## x 3 life stages (Ge = Greminant ; Im = Immature ; Ma = Mature)
## thus 9 parameters

SHADE_TOL <- param_store_df 
SHADE_TOL[1:9,] <- NA
rownames(SHADE_TOL) <- c("GeL", "GeM", "GeH", "ImL", "ImM", "ImH", "MaL", "MaM", "MaH")

## change the parameters value according to case we study
## REF: no change
## TEST: light cat <- light cat + 1

for(p_ in grep("LIGHT", names(params.val), value = TRUE)){
  if(params.val[p_] == "TEST"){
    elt_ <- is.element(dat$group, names(pfg.group)[pfg.group == sub("^.*[[:punct:]]", "", p_)])
    dat$light[elt_] <- dat$light[elt_] + 1
  }
}


## parameterisation according to Landolt classes and parcimony
for (i in 1:n){
	if(dat$light[i] ==4) {
	  SHADE_TOL["GeL",i] <- SHADE_TOL["ImL",i] <- SHADE_TOL["MaL",i] <- 1
	  SHADE_TOL["GeM",i] <- SHADE_TOL["ImM",i] <- SHADE_TOL["MaM",i] <- 1
	  SHADE_TOL["GeH",i] <- SHADE_TOL["ImH",i] <- SHADE_TOL["MaH",i] <- 0
		}
	if(dat$light[i] ==5) {
	  SHADE_TOL["GeL",i] <- SHADE_TOL["ImL",i] <- SHADE_TOL["MaL",i] <- 1
	  SHADE_TOL["GeM",i] <- SHADE_TOL["ImM",i] <- SHADE_TOL["MaM",i] <- 1
	  SHADE_TOL["GeH",i] <- SHADE_TOL["ImH",i] <- SHADE_TOL["MaH",i] <- 0
		}
	if(dat$light[i] ==6) {
	  SHADE_TOL["GeL",i] <- SHADE_TOL["ImL",i] <- SHADE_TOL["MaL",i] <- 1
	  SHADE_TOL["GeM",i] <- SHADE_TOL["ImM",i] <- SHADE_TOL["MaM",i] <- 1
	  SHADE_TOL["GeH",i] <- SHADE_TOL["ImH",i] <- SHADE_TOL["MaH",i] <- 1
		}
	if(dat$light[i] ==7) {
	  SHADE_TOL["GeL",i] <- SHADE_TOL["ImL",i] <- SHADE_TOL["MaL",i] <- 0
	  SHADE_TOL["GeM",i] <- SHADE_TOL["ImM",i] <- SHADE_TOL["MaM",i] <- 1
	  SHADE_TOL["GeH",i] <- SHADE_TOL["ImH",i] <- SHADE_TOL["MaH",i] <- 1
		}
	if(dat$light[i] ==8) {
	  SHADE_TOL["GeL",i] <- SHADE_TOL["ImL",i] <- SHADE_TOL["MaL",i] <- 0
	  SHADE_TOL["GeM",i] <- SHADE_TOL["ImM",i] <- SHADE_TOL["MaM",i] <- 1
	  SHADE_TOL["GeH",i] <- SHADE_TOL["ImH",i] <- SHADE_TOL["MaH",i] <- 1
		}
	if(dat$light[i] ==9) {
	  SHADE_TOL["GeL",i] <- SHADE_TOL["ImL",i] <- SHADE_TOL["MaL",i] <- 0
	  SHADE_TOL["GeM",i] <- SHADE_TOL["ImM",i] <- SHADE_TOL["MaM",i] <- 0
	  SHADE_TOL["GeH",i] <- SHADE_TOL["ImH",i] <- SHADE_TOL["MaH",i] <- 1
		}
}

## we assum that all germinants are tolerant to LOW and MEDIUM light
SHADE_TOL["GeL",] <- SHADE_TOL["GeM",] <- 1

## we have to adjust so that big trees continue growing when they are in the upper strata
SHADE_TOL["MaH",which(dat$type=="P")] <- 1 # last stratum in which all trees must tolerate the light
# the stratum before the last stratum must also tolerate the light (for the highest trees) because there will never be enough shadow in this stratum
SHADE_TOL["ImH",which(dat$type=="P" & CHANG_STR_AGES[3,] < dat$maturity )] <- 1 

## update SHADE_TOL rownames
rownames(SHADE_TOL) <- paste("SHADE_TOL", "_for_", rownames(SHADE_TOL), sep="")


#################################################################################################
### SEED POOLS (active and dormant) LIFE SPAN
## available seeds will exponentially decrease according to seed pool life span parameter
SEED_POOL_LIFE <- param_store_df
SEED_POOL_LIFE[1:2,] <- 0 ## 2 lines because of 2 seed pool (1 = active ;  2 = dormant)
rownames(SEED_POOL_LIFE) <- paste("SEED_POOL_LIFE", "_", c("active", "dormant"), sep="")

### DORMANCY
## 0 = no
## 1 = yes
SEED_DORMANCY <- param_store_df
SEED_DORMANCY[1,] <- 0
rownames(SEED_DORMANCY) <- "SEED_DORMANCY"

#################################################################################################
### PFG CONTRIBUTION TO SOIL (retroaction)
## when the PFG is covering the whole area, "amount" of nutrient given back to soil
## ~ nutrient richness
SOIL_CONTRIB<- param_store_df
SOIL_CONTRIB[1,] <- dat$Soil_contrib
rownames(SOIL_CONTRIB) <- "SOIL_CONTRIB"

#################################################################################################
### SOIL TOLERANCE
## 0 = non tolerant
## 1 = tolerant
## 3 categories of soil (0, 1, 2)
## x 3 life stages (Ge = Germinant ; Im = Immature ; Ma = Mature)
## thus 9 parameters

SOIL_TOL <- param_store_df 
SOIL_TOL[1:9,] <- 0
rownames(SOIL_TOL) <- c(paste("Ge",seq(3),sep=""),paste("Im",seq(3),sep=""),paste("Ma",seq(3),sep=""))

#Marta:
## 1a: filled according to tolerance class determined beforehand
dat$tol
min_tol_class <- sapply(dat$tol, function(x){ifelse(x<=2,1,2)})
max_tol_class <- sapply(dat$tol, function(x){ifelse(x>=2,3,2)})

for(sp_id in 1:ncol(SOIL_TOL)){ ## loop over species
  for(ls in 1:3){ ## loop over life stages
    SOIL_TOL[ ( (ls-1) * nb_soil_strat ) + (min_tol_class[sp_id]:max_tol_class[sp_id]),sp_id] <- 1
  }
}


# ## 1b : filled according to input parameters

# #Marta: in this case needs to be based on thresholds defined based on landolt, not on LNC (as specified in soil_strat_limits), because they are not on same scale!
# soil_strat_limits_tol<- c(2,3.5) ?

# min_tol_class <- sapply(dat$Soil_tol_min, function(x){head(which(soil_strat_limits_tol>=x),1)})
# max_tol_class <- sapply(dat$Soil_tol_max, function(x){tail(which(soil_strat_limits_tol <x),1)})

# for(sp_id in 1:ncol(SOIL_TOL)){ ## loop over species
  # for(ls in 1:3){ ## loop over life stages
    # SOIL_TOL[ ( (ls-1) * nb_soil_strat ) + (min_tol_class[sp_id]:max_tol_class[sp_id]),sp_id] <- 1
  # }
# }


## 2 : filled by hand, for users who want to define more finely soil constraints
### All woody PFGs are tolerant to rich soils when they are mature
#SOIL_TOL["Ma3",which(dat$type=="P")] <- 1
#SOIL_TOL["Ma3",which(dat$type=="C")] <- 1

### All woody PFGs are ALSO tolerant to rich soils when they are Immature
#SOIL_TOL["Im3",which(dat$type=="P")] <- 1
#SOIL_TOL["Im3",which(dat$type=="C")] <- 1

### All woody PFGs are ALWAYS tolerant to rich soils 
#SOIL_TOL["Ge3",which(dat$type=="P")] <- 1
#SOIL_TOL["Ge3",which(dat$type=="C")] <- 1


### Dryas
#i=2
#SOIL_TOL["Ma1",i] <- 1
#SOIL_TOL["Ma4",i] <- 1
#SOIL_TOL["Ma5",i] <- 1

## update SOIL_TOL rownames
rownames(SOIL_TOL) <- paste("SOIL_TOL", "_for_", rownames(SOIL_TOL), sep="")


#################################################################################################
### IS THE PFG AN ALIEN THAT WILL BE INTRODUCED LATER
IS_ALIEN <- param_store_df
#IS_ALIEN[1,] <- dat$is_alien
IS_ALIEN[1,] <- 0
rownames(IS_ALIEN) <- "IS_ALIEN"

#################################################################################################
# writing parameter files
#################################################################################################

## list of all supported parameters
available_param_list <- c("NAME",
                          "MATURITY", 
                          "LONGEVITY",
                          "MAX_ABUNDANCE", 
                          "IMM_SIZE", 
                          "CHANG_STR_AGES",
                          "IS_ALIEN",
                          "WIDE_DISPERS",
                          "MODE_DISPERS",
                          "ACTIVE_GERM",
                          "SHADE_TOL",
                          "SEED_POOL_LIFE",
                          "SEED_DORMANCY",
                          "SOIL_CONTRIB",
                          "SOIL_TOL")

## create empty params files for all pfgs
succ_file_names <- file.path(succ_path, paste("SUCC_",sp,".txt", sep=""))
#succ_file_names <- file.path(succ_path, paste("SUCC_ImCTol_",sp,".txt", sep=""))
#succ_file_names <- file.path(succ_path, paste("SUCC_AllCTol_",sp,".txt", sep=""))
names(succ_file_names) <- sp
## remove old files if exist
file.remove(succ_file_names)
## create the new ones
file.create(succ_file_names, showWarnings = F)

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
          file = succ_file_names[spp],
          append = TRUE)
    } ## end loop on species
  }
} ## end loop on parameters

## write also the summary parameters table
outFile= file.path(succ_path,"successionParamsTable.csv")
#outFile= file.path(succ_path,"successionParamsTable_ImCTol.csv")
#outFile= file.path(succ_path,"successionParamsTable_AllCTol.csv")
write.table(all_params_table,file=outFile, quote=FALSE, row.names=TRUE, col.names= NA, sep="\t") 



