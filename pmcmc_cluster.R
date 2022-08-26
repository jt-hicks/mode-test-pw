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

data_raw <- readRDS('data_raw_120822.RDS')
root <- "contexts"
sources <- c("run_pmcmc.R",
             "MiP-given/model_parameters.R","MiP-given/equilibrium-init-create-stripped.R")
config <- didehpc::didehpc_config(cores = 4, parallel = TRUE)


ctx <- context::context_save("contexts", sources = sources,
                             package_sources = conan::conan_sources(c("mrc-ide/odin.dust",'mrc-ide/mcstate')))

obj <- didehpc::queue_didehpc(ctx,config = config)

obj$cluster_load(TRUE)
obj$config
obj$login()
pmcmc_sim2 <- obj$enqueue(run_pmcmc(data=data_raw,EIR_vol=0.3,n_particles = 100))
pmcmc_sim <- obj$task_get("8681bf6c1d3fde4486ffd630f0229273")
pmcmc_sim2 <- obj$task_get("77442db45ba2aa8ca1c44dc8028d5cfe")
pmcmc_sim$id #"8681bf6c1d3fde4486ffd630f0229273"
pmcmc_sim$status()
pmcmc_sim$times()
pmcmc_sim$log()
pmcmc_sim2$id #"77442db45ba2aa8ca1c44dc8028d5cfe"
pmcmc_sim2$status()
pmcmc_sim2$times()
pmcmc_sim2$log()
obj$unsubmit("8681bf6c1d3fde4486ffd630f0229273")

data_raw_cmis <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/WesternKenya/data_raw_cmis.RDS')
data_raw_cmis_all <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/WesternKenya/data_raw_cmis_all.RDS')
pmcmc_cmis2 <- obj$enqueue(run_pmcmc(data=data_raw_cmis,EIR_vol=0.3,init_EIR_user = 10,n_particles = 100))
pmcmc_desktop_3 <- run_pmcmc(data=data_raw_cmis,EIR_vol=0.3,init_EIR_user = 20,n_particles = 100)
pmcmc_cmis_first <- obj$task_get('7e4d21482a7cf2bb1803414ca0606710')
pmcmc_cmis1 <- obj$task_get("817207f4d579cf9ab67e72bda6beffc6")
pmcmc_cmis1$id #"817207f4d579cf9ab67e72bda6beffc6"
pmcmc_cmis1$status()
pmcmc_cmis1$times()
pmcmc_cmis1$log()
pmcmc_cmis2 <- obj$task_get("381937e3430b312e5f7c27fc2394b8f3")
pmcmc_cmis2$id #"381937e3430b312e5f7c27fc2394b8f3"
pmcmc_cmis2$status()
pmcmc_cmis2$times()
pmcmc_cmis2$log()

pars <- data.frame(max_steps = c(2e4,4e4,1e5), atol = c(5e-6,10e-5,10e-4), rtol = c(5e-6,10e-5,10e-4))
cmis_bulk <- obj$enqueue_bulk(pars, run_pmcmc, data=data_raw_cmis,EIR_vol=0.3,init_EIR_user = 10,n_particles = 100)
cmis_bulk$status()

na <- 21
age = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)
fd <- vector(length=20L)
for(i in 1:20){
  fd[i] <- 1-(1-state$fD0)/(1+((age[i]+age[i+1])/2/state$aD)^state$gammaD)
}
p_det <- state$d1 + (1-state$d1)/(1 + fd[]*(state$init_ID[,]/state$ID0)^state$kD)
prev0to59 <- state$init_T[1:state$age59,] + state$init_D[1:state$age59,]  + state$init_A[1:state$age59,]*p_det[1:state$age59,]
prev <- sum(prev0to59[,])/sum(state$den[1:state$age59])


source("run_pmcmc.R")
source("MiP-given/model_parameters.R")
source("MiP-given/equilibrium-init-create-stripped.R")
pmcmc_desktop <- run_pmcmc(data=data_raw_cmis,EIR_vol=0.3,n_particles = 100)
saveRDS(pmcmc_desktop,'./pmcmc-run/pmcmc_desktop_1.RDS')

proposal_dist <- cov(pmcmc_desktop_1$pars)
pmcmc_desktop_2 <- run_pmcmc(data=data_raw_cmis,EIR_vol=0.3,proposal_dist=proposal_dist,n_particles = 100)
saveRDS(pmcmc_desktop_2,'./pmcmc-run/pmcmc_desktop_2.RDS')

pmcmc_desktop_3 <- run_pmcmc(data=data_raw_cmis,EIR_vol=0.3,init_EIR_user = 10,n_particles = 100, max_steps = 1e7,
                             atol = 1e-3,
                             rtol = 1e-3)
saveRDS(pmcmc_desktop_3,'./pmcmc-run/pmcmc_desktop_3.RDS')

pmcmc_desktop_4 <- run_pmcmc(data=data_raw_cmis,EIR_vol=0.3,init_EIR_user = 5,n_particles = 100,
                             max_steps = 1e8,
                             atol = 1e-2,
                             rtol = 1e-2)
data_raw_cmis_all <- readRDS('./pmcmc-run/data_raw_cmis_all.RDS')
pmcmc_desktop_all <- run_pmcmc(data=data_raw_cmis_all,EIR_vol=0.3,init_EIR = 5,n_particles = 100,
                             max_steps = 1e8,
                             atol = 1e-2,
                             rtol = 1e-2)
saveRDS(pmcmc_desktop_all,'./pmcmc-run/pmcmc_desktop_all_EIR.RDS')
pmcmc_desktop_exp <- run_pmcmc(data=data_raw_cmis_all,EIR_vol=0.3,init_EIR = 5,n_particles = 100,
                               max_steps = 1e8,
                               atol = 1e-2,
                               rtol = 1e-2)
saveRDS(pmcmc_desktop_exp,'./pmcmc-run/pmcmc_desktop_all_exp.RDS')

pmcmc_desktop_exp <- run_pmcmc(data=data_raw_cmis_all,EIR_vol=0.1,init_EIR = 5,n_particles = 200,
                               max_steps = 1e6,
                               atol = 1e-6,
                               rtol = 1e-6)
saveRDS(pmcmc_desktop_exp,'./pmcmc-run/pmcmc_desktop_all_exp.RDS')
saveRDS(pmcmc_desktop_exp,'./pmcmc-run/cmis_good1.RDS')


data_raw_ng <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/data_raw_NG_anc.RDS')
ng_test <- run_pmcmc(data=data_raw_ng,EIR_vol=0.3,init_EIR = 40,n_particles = 200,
                     max_steps = 1e6,
                     atol = 1e-4,
                     rtol = 1e-4)
saveRDS(ng_test,'ng_test.RDS')

data_raw_bf <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/data_raw_BF_anc.RDS')
bf_test <- run_pmcmc(data=data_raw_bf,EIR_vol=0.1,init_EIR = 10,n_particles = 100,
                     max_steps = 1e6,
                     atol = 1e-4,
                     rtol = 1e-4)
saveRDS(bf_test,'bf_test.RDS')

data_sim <- readRDS('data_sim_180822.RDS')

data_raw <- data_gen(EIR_volatility = 0.6, init_EIR = 20)
saveRDS(data_raw,'data_sim2.RDS')

sim_test <- run_pmcmc(data=data_raw,EIR_vol=0.1,init_EIR = 20,n_particles = 100,
                     max_steps = 1e6,
                     atol = 1e-4,
                     rtol = 1e-4)
saveRDS(sim_test,'sim_test.RDS')
