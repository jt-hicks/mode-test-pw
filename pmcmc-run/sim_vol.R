##Generate various simulated data sets at various initial EIR and EIR volatility
##values, then fit with PMCMC on cluster
library("odin.dust")
library("odin")
library("patchwork")
library('mcstate')
library(didehpc)
library(pkgdepends)
library(dplyr)
library("coda")
library(binom)
library(ggplot2)

##Function to generate simulated data and run pmcmc
source('data_gen.R')
init_EIR_list <- c(1,10,100)
EIR_vol_list <- c(0.3,0.6,1.0,1.5,2.0)
pars <- expand.grid(EIR_volatility = EIR_vol_list,init_EIR = init_EIR_list)
proposal_dist <- matrix(c(0.0336,-0.000589,-0.000589,0.049420))

root <- "contexts"
sources <- c("run_pmcmc.R",
             "data_gen.R",
             "gen_run.R",
             "MiP-given/model_parameters.R",
             "MiP-given/equilibrium-init-create-stripped.R")
config <- didehpc::didehpc_config(cores = 4, parallel = TRUE)


ctx <- context::context_save("contexts", sources = sources,
                             packages = c('odin','dde'),
                             package_sources = conan::conan_sources(c("mrc-ide/odin.dust",'mrc-ide/mcstate')))

obj <- didehpc::queue_didehpc(ctx,config = config)

obj$cluster_load(TRUE)
obj$config
obj$login()
proposal_dist
sim_vol_1 <- obj$enqueue_bulk(pars, gen_run, proposal=matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2))
bundles <- obj$task_bundle_list()
obj$task_bundle_info()

sim_vol_1$status()
sim_vol_1$tasks$`080b7cb622066786f9f093de5f994612`$log()
sim_vol_1.1 <- sim_vol_1$tasks$`4f5d3410769b1a817b32f709295775a0`$result()
sim_vol_1.1_times <- sim_vol_1$tasks$`4f5d3410769b1a817b32f709295775a0`$times()
sim_vol_1.1_mcmc <- coda::as.mcmc(cbind(sim_vol_1.1$probabilities, sim_vol_1.1$pars))
1 - coda::rejectionRate(sim_vol_1.1_mcmc)
coda::effectiveSize(sim_vol_1.1_mcmc)
plot(sim_vol_1.1_mcmc)
raftery.diag(sim_vol_1.1_mcmc)
accept.rate <- function(submission){
  if(submission$status() == 'ERROR'){
    if(grepl('particles reported errors',submission$log()[4])){
      return('PMCMC ERROR')
    }
    else if(grepl('Error in dde',submission$log()[4])){
      return('SIMULATION ERROR')
    }
    else{
      return('OTHER ERROR')
    }
  }
  else if(submission$status() == 'CANCELLED'){
    return('RUN CANCELLED')
  }
  else if(submission$status() == 'RUNNING'){
    return('STILL RUNNING')
  }
  result <- submission$result()
  time <- submission$times()$running
  mcmc <- coda::as.mcmc(cbind(result$probabilities, result$pars))
  ar <- 1 - coda::rejectionRate(mcmc)
  return(c(acceptance = ar[[1]],run_time = time))
}
sim_vol_1$tasks$`4f5d3410769b1a817b32f709295775a0`$status()
grepl('particles reported errors',sim_vol_1$tasks[[7]]$log()[4])
grepl('Error in dde',sim_vol_1$tasks[[7]]$log()[4])
grepl('particles reported errors',sim_vol_1$tasks[[5]]$log()[4])
grepl('Error in dde',sim_vol_1$tasks[[5]]$log()[4])

log <- sim_vol_1$tasks[[7]]$log()
tasks <- sim_vol_1$tasks
result <- tasks[[1]]$result()
lapply(tasks,accept.rate)
lapply(sim_vol_2$tasks,accept.rate)
lapply(sim_vol_3$tasks,accept.rate)
lapply(sim_vol_4$tasks,accept.rate)
lapply(sim_vol_5$tasks,accept.rate)
lapply(sim_vol_6$tasks,accept.rate)

##Get bundles
bundles <- obj$task_bundle_list()
bundle_info <- obj$task_bundle_info()
bundle_info_sim <- bundle_info[bundle_info$`function` == 'gen_run' & bundle_info$length==15,]
bundle_info_sim <- bundle_info_sim[8:13,]

##Take list of bulk submissions and summarize
df.all <- data.frame(log_prior = numeric(),
                     log_likelihood = numeric(),
                     log_posterior = numeric(),
                     EIR_SD = numeric(),
                     log_init_EIR = numeric(),
                     step = integer(),
                     true_EIR_SD = numeric(),
                     true_log_init_EIR = numeric(),
                     group = integer()
                     )
ar.all <- data.frame(acceptance = numeric(),
                     run_time = numeric(),
                     true_EIR_SD = numeric(),
                     true_log_init_EIR = numeric(),
                     group = integer()
)
init_EIR_list <- c(1,10,100)
EIR_vol_list <- c(0.3,0.6,1.0,1.5,2.0)
pars <- expand.grid(EIR_volatility = EIR_vol_list,init_EIR = init_EIR_list)
task <- obj$task_get('6e555cdf471ef2a29eb4542f31b63aee')
task$log()
for(i in 1:nrow(bundle_info_sim)){
  group <- obj$task_bundle_get(bundle_info_sim[i,1])
  ar <- as.data.frame(t(as.data.frame(lapply(group$tasks,accept.rate))))
  ar$true_EIR_SD = pars$EIR_volatility
  ar$true_log_init_EIR = pars$init_EIR
  ar$group = i
  ar.all <- rbind(ar.all,ar)
  for(j in 1:15){
    if(group$tasks[[j]]$status() != 'COMPLETE'){
      next
    }
    result <- group$tasks[[j]]$result()
    mcmc <- as.data.frame(coda::as.mcmc(cbind(result$probabilities,result$pars)))
    mcmc$true_EIR_SD = pars[j,1]
    mcmc$true_log_init_EIR = pars[j,2]
    mcmc$group = i
    mcmc <- cbind(step = 1:nrow(mcmc), mcmc) 
    df.all <- rbind(df.all,mcmc)
  }
}

saveRDS(ar.all,'bulk_ar_080922.rds')
saveRDS(df.all,'bulk_mcmc_080922.rds')

library(reshape2)
library(ggplot2)
library(ggpubr)
library(zoo)
library(patchwork)
library(RColorBrewer)
RColorBrewer::display.brewer.all()
theme_set(theme_minimal()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.border = element_rect(colour = "black",fill=NA)))

windows(80,50)
ggplot(df.all,aes(x=EIR_SD,color=as.factor(true_log_init_EIR),linetype=as.factor(group)))+
  geom_density()+
  facet_wrap(~as.factor(true_EIR_SD))+
  geom_vline(aes(xintercept = true_EIR_SD))+
  scale_color_brewer(palette='Set1')+
  theme(legend.position = 'bottom')+
  labs(x='Volatility',y='Density',color='EIR',linetype='Run')
windows(80,30)
ggplot(df.all,aes(x=log_init_EIR,color=as.factor(true_EIR_SD),linetype=as.factor(group)))+
  geom_density()+
  facet_wrap(~as.factor(true_log_init_EIR))+
  geom_vline(aes(xintercept = log(true_log_init_EIR)))+
  scale_color_brewer(palette='Set1')+
  theme(legend.position = 'bottom')+
  labs(x='EIR',y='Density',color='Volatility',linetype='Run')
windows(80,50)
ggplot(df.all,aes(x=step,y=log_posterior,color=as.factor(group)))+
  geom_line()+
  facet_grid(as.factor(true_EIR_SD)~as.factor(true_log_init_EIR))+
  scale_color_brewer(palette='Paired')+
  theme(legend.position = 'bottom')+
  labs(x='Step',y='Log Posterior',color='Run')
ggplot(df.all,aes(x=step,y=EIR_SD,color=as.factor(group)))+
  geom_line()+
  facet_grid(as.factor(true_EIR_SD)~as.factor(true_log_init_EIR))+
  scale_color_brewer(palette='Paired')+
  theme(legend.position = 'bottom')+
  labs(x='Step',y='Volatility',color='Run')
windows(30,30)
ggplot(df.all,aes(x=step,y=log_init_EIR,color=as.factor(group)))+
  geom_line()+
  facet_grid(as.factor(true_EIR_SD)~as.factor(true_log_init_EIR))+
  scale_color_brewer(palette='Paired')+
  theme(legend.position = 'bottom')+
  labs(x='Step',y='Log EIR',color='Run')

ar_eir <- ggplot(ar.all,aes(x=as.factor(true_log_init_EIR),y=as.numeric(acceptance)))+
  geom_point()+
  # facet_grid(as.factor(true_EIR_SD)~as.factor(true_log_init_EIR))+
  scale_y_continuous(limits=c(0,1))+
  labs(x='EIR',y='Acceptance Rate',color='Run')

ar_vol <- ggplot(ar.all,aes(x=as.factor(true_EIR_SD),y=as.numeric(acceptance)))+
  geom_point()+
  # facet_grid(as.factor(true_EIR_SD)~as.factor(true_log_init_EIR))+
  scale_y_continuous(limits=c(0,1))+
  labs(x='Volatility',y='Acceptance Rate',color='Run')
windows(80,50)
ar.plot <- ar_vol + ar_eir

time_eir <- ggplot(ar.all,aes(x=as.factor(true_log_init_EIR),y=as.numeric(run_time)/3600))+
  geom_point()+
  # facet_grid(as.factor(true_EIR_SD)~as.factor(true_log_init_EIR))+
  # scale_y_continuous(limits=c(0,1))+
  labs(x='EIR',y='Run Time (Hours)',color='Run')
time_eir
time_ar <- ggplot(ar.all,aes(x=as.factor(true_EIR_SD),y=as.numeric(run_time)/3600))+
  geom_point()+
  # facet_grid(as.factor(true_EIR_SD)~as.factor(true_log_init_EIR))+
  # scale_y_continuous(limits=c(0,1))+
  labs(x='Volatility',y='Run Time (Hours)',color='Run')
time_ar
time.plot <- time_ar + time_eir
windows(60,50)
ar.plot/time.plot

ar.all.errors <- ar.all%>%
  mutate(error=ifelse(acceptance %in% c('PMCMC ERROR','STILL RUNNING'),'PMCMC Error',
                      ifelse(acceptance=='SIMULATION ERROR','Simulation Error','Success')))

table(ar.all.errors$error)
table(ar.all.errors$true_EIR_SD,ar.all.errors$true_log_init_EIR)
table(ar.all.errors[ar.all.errors$error=='Success',]$true_EIR_SD,ar.all.errors[ar.all.errors$error=='Success',]$true_log_init_EIR)
table(ar.all.errors[ar.all.errors$error=='Simulation Error',]$true_EIR_SD,ar.all.errors[ar.all.errors$error=='Simulation Error',]$true_log_init_EIR)
table(ar.all.errors[ar.all.errors$error!='Simulation Error',]$true_EIR_SD,ar.all.errors[ar.all.errors$error!='Simulation Error',]$true_log_init_EIR)
table(ar.all.errors[ar.all.errors$error!='Simulation Error'&ar.all.errors$error=='PMCMC Error',]$true_EIR_SD,ar.all.errors[ar.all.errors$error!='Simulation Error'&ar.all.errors$error=='PMCMC Error',]$true_log_init_EIR)

#obj$unsubmit(sim_vol_1$ids)

sim_vol_2 <- obj$enqueue_bulk(pars, gen_run, proposal=matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2))
sim_vol_2$status()


sim_vol_3 <- obj$enqueue_bulk(pars, gen_run, proposal=matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2))
sim_vol_3$status()

sim_vol_4 <- obj$enqueue_bulk(pars, gen_run, proposal=matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2))
sim_vol_4$status()

sim_vol_5 <- obj$enqueue_bulk(pars, gen_run, proposal=matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2))
sim_vol_5$status()

sim_vol_6 <- obj$enqueue_bulk(pars, gen_run, proposal=matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2))
sim_vol_6$status()

##Get bundles
bundles <- obj$task_bundle_list()
bundle_info <- obj$task_bundle_info()


##Bulk simulated data creation
source('data_gen.R')
source("run_pmcmc.R")
source("MiP-given/model_parameters.R")
source("MiP-given/equilibrium-init-create-stripped.R")
init_EIR_list <- c(1,10,100)
EIR_vol_list <- c(0.3,0.6,0.9,1.5,2.0)
#EIR_vol_list <- c(0.3,0.6,0.9,1.2)
pars <- expand.grid(EIR_volatility = EIR_vol_list,init_EIR = init_EIR_list)
proposal_dist <- matrix(c(0.0336,-0.000589,-0.000589,0.049420))
i <- rep(1:nrow(pars),10)
data_sim_1 <- with(pars, lapply(i, function(j){data_gen(EIR_volatility[j],init_EIR[j])}))
sim_bulk_1 <- obj$enqueue_bulk(1:length(data_sim_1), function(i,data_sim){
  run_pmcmc(data_sim[[i]],proposal_matrix=matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
            n_particles=200,
            max_EIR=1000,
            max_steps = 1e7,
            atol = 1e-5,
            rtol = 1e-6,
            n_steps = 10000)
},data_sim=data_sim_1)
sim_bulk_1$status()
sim_bulk_1$name #bundle name: 'paralyzing_bettong'
saveRDS(data_sim_1,'data_sim_max1000.RDS')
test_result <- sim_bulk_1$tasks$e860007297d850fd28981c8f134d9516$result()
sim_test_mcmc <- coda::as.mcmc(cbind(test_result$probabilities, test_result$pars))
1 - coda::rejectionRate(sim_test_mcmc)
coda::effectiveSize(sim_test_mcmc)
windows(60,40)
plot(sim_test_mcmc)
raftery.diag(sim_test_mcmc)


# obj$unsubmit(sim_bulk_1$ids)
obj$cluster_load(TRUE)
obj$login()

# saveRDS(data_sim_1,'data_sim_test.RDS')
# sim_bulk_test <- obj$enqueue_bulk(1:2, function(i,data_sim){
#   run_pmcmc(data_sim[[i]],proposal_matrix=matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
#                                                                          n_particles=200,
#                                                                          max_EIR=1000,
#                                                                          max_steps = 1e7,
#                                                                          atol = 1e-5,
#                                                                          rtol = 1e-6,
#                                                                          n_steps = 1000)
# },data_sim=data_sim_1)
# sim_bulk_test$status()
# sim_bulk_test$tasks$`12953245ea2f640ef650f9e7ba85cc24`$log()
# obj$unsubmit(sim_bulk_test$ids)
##newspaperish_queenbee
sim_bulk_1$tasks$e71f17d6b23c0c0cdf5098800dbc21f5$log()
data_sim_2 <- with(pars, lapply(i, function(j){data_gen(EIR_volatility[j],init_EIR[j],max_EIR = 10000)}))
saveRDS(data_sim_2,'data_sim_max10000.RDS')

check <- data_sim_2[[length(data_sim_2)]]
data_sim_3 <- with(pars, lapply(i, function(j){data_gen(EIR_volatility[j],init_EIR[j],max_EIR = 100000)}))
saveRDS(data_sim_3,'data_sim_max100000.RDS')

data_sim_4 <- with(pars, lapply(i, function(j){data_gen(EIR_volatility[j],init_EIR[j],max_EIR = 100000)}))

check.prev <- function(df){
  prev_neg <- sum(ifelse(df$prev_true<0,1,0))
  eir_neg <- sum(ifelse(df$EIR_true<0,1,0))
  inc_neg <- sum(ifelse(df$inc_true<0,1,0))
  return(c(prev=prev_neg,eir=eir_neg,inc=inc_neg))
}
negs <- lapply(1:length(data_sim_2),function(i) check.prev(data_sim_2[[i]]))
df.negs <- as.data.frame(t(as.data.frame(negs)))
plot.truth <- function(df){
  ggplot(df,aes(x=t,y=prev_true))+
    geom_line(color='blue',size=1)
}
plot.truth(data_sim[[6]])

plot.eir <- function(df){
  ggplot(df,aes(x=t,y=EIR_true))+
    geom_line(color='blue',size=1)
}
windows(40,50)
eir <- plot.eir(data_sim_3[[10]])
prev <- plot.truth(data_sim_3[[10]])
eir/prev

eir1 <- plot.eir(data_sim_1[[20]])
prev1 <- plot.truth(data_sim_1[[20]])
eir1/prev1

eir2 <- plot.eir(data_sim_2[[30]])
prev2 <- plot.truth(data_sim_2[[30]])
eir2/prev2
test <- data_sim_2[[30]]

data_gen(EIR_volatility = 0.8, init_EIR=1000, max_EIR = 100000)
