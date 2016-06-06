# library(sensitivity)
# 
# pfg.list <- c(paste0("C", 1:6), paste0("C", 1:10), paste0("P", 1:8))
# moris.factors <- c("hs.type", "area", "soil.tol.thresh.min", "soil.tol.thresh.max", paste0("soil.contrib.", pfg.list), paste0("soil.tol.", pfg.list))
# moris.factors.levels <- c(3, 3, rep(10, 2 + 2 * length(pfg.list)))
# 
# morris.plan <-  morris(model = NULL, 
#                        factors = moris.factors, 
#                        r = 1,
#                        design = list(type = "oat", levels = moris.factors.levels, grid.jump = 5))
# 
# summary(morris.plan$X)
# dim(morris.plan$X)

library(pse)
library(dplyr)

dat <- read.csv2("~/Work/FATEHD/benchmarking/fatehd_params_test/grid/soil_params_test/fhdsoilpt/Data/editParamsC24_newH8.csv", 
                  stringsAsFactors = FALSE)
pfg.list <- dat$group

lhs.factors <- c("hs.type", "area", "glob.soil.tol.soil.default.value", "glob.soil.tol.soil.percent", "glob.soil.tol.thresh.min", "glob.soil.tol.thresh.max", paste0("pfg.soil.contrib.", pfg.list), paste0("pfg.soil.tol.", pfg.list))

lhs.factors.df <- data.frame(factors = lhs.factors)
# lhs.factors.df <- data.frame(factors = lhs.factors,
#                              quant.funct = c(rep("qunif", 2), rep("qnorm", 2 + 2 * length(pfg.list))),
#                              quant.params = NA)

f.quant.params <- function(x){
  sd.fact <- 4
  if(x$factors == "hs.type") return(data.frame( quant.funct = "qunif", quant.params = list(min = 0, max = 2)))
  if(x$factors == "area") return(data.frame( quant.funct = "qunif", quant.params = list(min = 0, max = 3)))
  if(x$factors == "glob.soil.tol.soil.default.value") return(data.frame( quant.funct = "qunif", quant.params = list(min = 0, max = 30)))
  if(x$factors == "glob.soil.tol.soil.percent") return(data.frame( quant.funct = "qunif", quant.params = list(min = 0, max = 1)))
  if(x$factors == "glob.soil.tol.thresh.min") return(data.frame( quant.funct = "qnorm", quant.params = list(mean = 20, sd = 20/sd.fact)))
  if(x$factors == "glob.soil.tol.thresh.max") return(data.frame( quant.funct = "qnorm", quant.params = list(mean = 25, sd = 25/sd.fact)))
  if(grepl("pfg.soil.contrib.", x$factors)){ mean_ <- as.numeric(dat$Soil_contrib[dat$group == sub("pfg.soil.contrib.", "", x$factors)]); return(data.frame( quant.funct = "qnorm", quant.params = list(mean = mean_, sd = mean_/sd.fact)))}
  if(grepl("pfg.soil.tol.", x$factors)) return(list( quant.funct = "qunif", quant.params = list(min = 0, max = 3)))
}

lhs.factors.df2 <- lhs.factors.df %>% rowwise %>% do(data.frame(f.quant.params(x = .)))

## test with a optimized lhs
lhs.q <- sapply(1:nrow(lhs.factors.df2), function(i) lhs.factors.df2[i, "quant.funct"] %>% as.character) 
lhs.q.arg <- lapply(1:nrow(lhs.factors.df2), function(i){ if(lhs.factors.df2[i, "quant.funct"] %>% as.character == "qunif") return(list(min = lhs.factors.df2[i, "quant.params.min"] %>% as.numeric, max = lhs.factors.df2[i, "quant.params.max"] %>% as.numeric)) else return(list(mean = lhs.factors.df2[i, "quant.params.mean"] %>% as.numeric, sd = lhs.factors.df2[i, "quant.params.sd"] %>% as.numeric))}) 

# new.cluster <- parallel::makeForkCluster(nnodes = 8)
# uncoupledLHS <- LHS(model = NULL,
#                     factors = lhs.factors,
#                     N = 10000,
#                     repetitions = 1,
#                     q = lhs.q,
#                     q.arg  = lhs.q.arg,
#                     cl = new.cluster)
# write.csv(get.data(uncoupledLHS), file="mydata.par.csv")
# stopCluster(new.cluster)

uncoupledLHS <- LHS(model = NULL,
                    factors = lhs.factors,
                    N = 10000,
                    repetitions = 1,
                    q = lhs.q,
                    q.arg  = lhs.q.arg,
                    cl = new.cluster,
                    method = 'random')

uncoupledLHS.dat <- get.data(uncoupledLHS)
head(uncoupledLHS.dat)

write.csv(uncoupledLHS.dat, file = "~/Work/FATEHD/benchmarking/fatehd_params_test/grid/soil_params_test/fhdsoilpt/Data/uncoupledLHS.csv")


