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
library(dplyr)

##Data generation##
source("MiP-given/model_parameters.R")
source("MiP-given/equilibrium-init-create-stripped.R")
source("data_gen.R")
source("loop_function_main.R")

data_raw <- data_gen(EIR_volatility = 0.3)
saveRDS(data_raw,'data_sim2.RDS')
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
obj$login()
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

#Experiment Cae Adda:
# Generate 5 simulated dataset at volatility = 0.3
# For each, fun a 200 particle pf with varying volatilities (0.1, 0.2, 0.3, 0.4 and 0.5) - do not vary stochastic schedule
vols <- rep(0.3,5)
data_vol0.3 <- lapply(vols, FUN=data_gen)
saveRDS(data_vol0.3,'C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PMCMC/Development/cae_adda_data.rds')

ca_grid <- expand.grid(volatility = seq(0.1,0.5,0.1), data = data_vol0.3)
ca_grid_test <- expand.grid(volatility = 0.1, data = as.data.frame(data_vol0.3[1]))
ca_test <- obj$enqueue_bulk(ca_grid_test, pf_loop, n_particles=200, n_loop=100)
ca <- obj$enqueue_bulk(ca_grid, pf_loop, n_particles=200, n_loop=100)
ca_test$times()
ca_test$results()
ca_1 <- obj$lapply(seq(0.1,0.5,0.1), FUN=pf_loop, n_particles=200, n_loop=100,data = as.data.frame(data_vol0.3[1]),name='cae_adda1')
ca_1$times()
# ca_2$times()
# ca_2$status()
ca_3$times()
## ca_4$times()
# ca_5$times()

ca_3$status()
ca_1$tasks$`61700ab0e536c632291e81112ec72949`$result()
ca_1$tasks$d682bbb8253949272ba3ffb96a35322f$result()
ca_1$tasks$`77cb8e923c13e44f84ed6dae0997868a`$result()
ca_1$tasks$`3a02dac89051fef8a3adc2930a6479f2`$log()
ca_2$tasks$`07031d1edea13c8d2292350cddcc63c3`$result()
ca_2$tasks$`8090b6dfd2038fdf39c3abf3e58ec085`$log()
ca_5$tasks$`0f5bc2396b4954951604cb1bb356b69d`$result()

ca_2 <- obj$lapply(seq(0.1,0.5,0.1), FUN=pf_loop, n_particles=200, n_loop=100,data = as.data.frame(data_vol0.3[2]),name='cae_adda2')
ca_3 <- obj$lapply(seq(0.1,0.5,0.1), FUN=pf_loop, n_particles=200, n_loop=100,data = as.data.frame(data_vol0.3[3]),name='cae_adda3')
ca_4 <- obj$lapply(seq(0.1,0.5,0.1), FUN=pf_loop, n_particles=200, n_loop=100,data = as.data.frame(data_vol0.3[4]),name='cae_adda4')
ca_5 <- obj$lapply(seq(0.1,0.5,0.1), FUN=pf_loop, n_particles=200, n_loop=100,data = as.data.frame(data_vol0.3[5]),name='cae_adda5')

ca_1_results <- ca_1$results()
names(ca_1_results) <- c('0.1','0.2','0.3','0.4','0.5')
saveRDS(ca_1_results,'C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PMCMC/Development/results_ca1_100822.rds')


ca_2_results <- ca_2$results()
names(ca_2_results) <- c('0.1','0.2','0.3','0.4','0.5')
saveRDS(ca_2_results,'C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PMCMC/Development/results_ca2_100822.rds')

ca_3_results <- ca_3$results()
names(ca_3_results) <- c('0.1','0.2','0.3','0.4','0.5')
saveRDS(ca_3_results,'C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PMCMC/Development/results_ca3_100822.rds')


ca_4_results <- ca_4$results()
names(ca_4_results) <- c('0.1','0.2','0.3','0.4','0.5')
saveRDS(ca_4_results,'C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PMCMC/Development/results_ca4_100822.rds')

ca_5_results <- ca_5$results()
names(ca_5_results) <- c('0.1','0.2','0.3','0.4','0.5')
saveRDS(ca_5_results,'C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PMCMC/Development/results_ca5_100822.rds')

ca1_times <- ca_1$times()
ca1_times$run <- 1

ca2_times <- ca_2$times()
ca2_times$run <- 2

ca3_times <- ca_3$times()
ca3_times$run <- 3

ca4_times <- ca_4$times()
ca4_times$run <- 4

ca5_times <- ca_5$times()
ca5_times$run <- 5

ca_times <- rbind(ca1_times,ca2_times,ca3_times,ca4_times,ca5_times)
saveRDS(ca_times,'C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PMCMC/Development/times_ca_100822.rds')

##Error=NA test
start.time <- Sys.time()
test <- pf_loop(volatility = 0.9,data_file = data_raw, n_particles = 200, n_loop=10)
difftime(Sys.time(),start.time,units='secs')
test

#Tolerance experiments
vols <- c(0.3,0.3,0.9,0.9)
data_tol <- lapply(vols, FUN=data_gen)
saveRDS(data_tol,'C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PMCMC/Development/tol_data.rds')
data_tol <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PMCMC/Development/tol_data.rds')

tol_grid <- expand.grid(volatility = 0.9, max_steps = c(1e4,1e6), atol = c(10e-6,10e-4), rtol = c(10e-6,10e-4), data_file = data_tol)
tol_grid <- expand.grid(volatility = 0.9, max_steps = c(1e4,1e6), atol = c(10e-6,10e-4), rtol = c(10e-6,10e-4), data_file = c(1:4))
tol <- obj$enqueue_bulk(tol_grid, pf_loop, n_particles=200, n_loop=20)
tol <- obj$task_bundle_get('emigrational_bobwhite')
tol$times()
tol$status()
tol$tasks$`0c0b6b3a5e4c12bfa6b82803e36dd76a`$log()

tol$tasks$`66c7738901d8fba83b2a845c4908309d`$result()
tol$tasks$b56d10ccba23989240215c5b27e1356b$result()
tol_grid[1,]
tol$tasks$de1b43a7909d761bbda80612c2f3d21c$result()
tol_test$tasks$`74d906af9bbee57179b8b59bfc890f37`$log()
tol_test$tasks$`85b2919f23656d6d3e5f152f9e05b0e0`$result()

tol_results <- tol$results()

df_tol_results <- data.frame(volatility=numeric(),
                             max_steps=numeric(),
                             atol=numeric(),
                             rtol=numeric(),
                             data_file=factor(),
                             value=character(),
                             variable=character()
                             )
for(i in 1:length(tol_results)){
  temp <- tol_results[[i]] %>% 
    rename(value=`200`)
  temp <- cbind(temp,tol_grid[i,])
  df_tol_results<- rbind(df_tol_results,temp)
}

df_tol_results$data_vol <- ifelse(df_tol_results$data_file <=2,0.3,0.9)
saveRDS(df_tol_results,'C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PMCMC/Development/tol_results.rds')
saveRDS(tol_grid,'C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PMCMC/Development/tol_params.rds')
