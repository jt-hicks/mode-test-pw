#Run loop function
setwd('Q:/mode-test-pw')
# drat:::add("mrc-ide") # install.packages("drat") if this errors, then try again
# install.packages("didehpc")
library("odin.dust")
library("odin")
library("patchwork")
library('mcstate')
library(didehpc)
library(pkgdepends)

##Data generation##
source("MiP-given/model_parameters.R")
source("MiP-given/equilibrium-init-create-stripped.R")
source("data_gen.R")

data_raw <- data_gen(EIR_volatility = 0.4)

root <- "contexts"
sources <- c("loop_function_main.R","MiP-given/model_parameters.R",
             "MiP-given/equilibrium-init-create-stripped.R")
config <- didehpc::didehpc_config(cores = 4, parallel = TRUE)


ctx <- context::context_save("contexts", sources = sources,
                             packages = c("statmod"),
                             package_sources = conan::conan_sources(c("mrc-ide/odin.dust",'mrc-ide/mcstate')))

obj <- didehpc::queue_didehpc(ctx,config = config)

obj$cluster_load(TRUE)
obj$config

t <- obj$enqueue(pf_loop(data=data_raw,n_loop = 100))
obj$unsubmit(t$id)
start.time <- Sys.time()
test <- pf_loop(data_file=data_raw,n_loop = 3,volatility=0.4)
test2 <- pf_loop(data_file=data_raw,n_loop = 3,volatility=0.8)
print(Sys.time()-start.time)

t$status()
t$result()
t$wait(100)

t$times()
t$log()

#Group 1
pars <- expand.grid(volatility = c(0.2, 0.4), freq = c(1, 2))
group1 <- obj$enqueue_bulk(pars, pf_loop, data=data_raw, n_loop=100)
obj$unsubmit(group1$ids[4])
group1$wait(100)
group1$results()
group1$times()
group1$log()
results_group1 <- group1$results()
names(results_group1) <- c('vol=0.2,freq=1','vol=0.4,freq=1','vol=0.2,freq=2','vol=0.4,freq=2')
saveRDS(results_group1,'C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PMCMC/Development/results_group1_lik_020822.rds')

#Group 2
pars2 <- expand.grid(volatility = c(0.04, 0.4), freq = c(1, 10))
group2 <- obj$enqueue_bulk(pars2, pf_loop, data=data_raw, n_loop=100)
obj$unsubmit(group2$ids[3])
group2$wait(100)
group2$results()
group2$times()
group2$log()
results_group2 <- group2$results()
names(results_group2) <- c('vol=0.04,freq=1','vol=0.4,freq=1','vol=0.04,freq=10','vol=0.4,freq=10')
saveRDS(results_group2,'C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PMCMC/Development/results_group2_lik_020822.rds')

#Group 3
pars <- expand.grid(volatility = c(0.2, 0.4), freq = c(2))
group3 <- obj$enqueue_bulk(pars, pf_loop, data=data_raw, n_loop=100)
obj$unsubmit(group1$ids[4])
group1$wait(100)
group3$results()
group3$times()
group1$log()
results_group3 <- group3$results()
names(results_group3) <- c('vol=0.2,freq=2','vol=0.4,freq=2','vol=0.2,freq=2','vol=0.4,freq=2')
saveRDS(results_group3,'C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PMCMC/Development/results_group3_lik_020822.rds')

#Group 4
pars2 <- expand.grid(volatility = c(0.04, 0.4), freq = c(10))
group4 <- obj$enqueue_bulk(pars2, pf_loop, data=data_raw, n_loop=100)
obj$unsubmit(group4$ids[1])
group2$wait(100)
group4$results()
group4$times()
group2$log()
results_group2 <- group2$results()
names(results_group2) <- c('vol=0.04,freq=1','vol=0.4,freq=1','vol=0.04,freq=10','vol=0.4,freq=10')
saveRDS(results_group2,'C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PMCMC/Development/results_group2_lik_020822.rds')
