library("odin.dust")
library("odin")
library("patchwork")
library('mcstate')
library(RColorBrewer)
source("MiP-given/model_parameters.R")
source("MiP-given/equilibrium-init-create-stripped.R")
library("coda")
library(binom)


########################################
####### Function to plot particle trajectories
addCIs<-function(df,Ys,Ns){
  df$mean<-NA
  df$upper<-NA
  df$lower<-NA
  CIs<-ifelse(is.na(Ns) | is.na(Ys),NA,binom.confint(Ys,Ns,method="exact"))
  df$mean[Ns>0]<-CIs$mean[Ns>0]
  df$upper[Ns>0]<-CIs$upper[Ns>0]
  df$lower[Ns>0]<-CIs$lower[Ns>0]
  return(df)
}
plot_particle_filter <- function(history, true_history, times, obs_end = NULL) {
  cis <- addCIs(true_history,true_history$positive,true_history$tested)
  if (obs_end>max(times)) {
    times <- c(times,seq(max(times)+30,obs_end,30))
    cis[(nrow(cis)+1):length(times),] <- NA
  }
  par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
  matplot(times, t(history[1, , -1]), type = "l",
          xlab = "Time", ylab = "Prevalence",
          col = "#A6CEE3", lty = 1)
  matpoints(times, cis$mean, pch = 19,
            col = "#1F78B4")
  arrows(cis$t, cis$lower, cis$t, cis$upper, length=0.05, angle=90, code=3, col = "#1F78B4")
}


pmcmc_desktop_1 <- as.list(readRDS('./pmcmc-run/pmcmc_desktop_1.RDS'))
mcmc_1 <- coda::as.mcmc(cbind(pmcmc_run$probabilities, pmcmc_run$pars),
                        )
1 - coda::rejectionRate(mcmc_1)
coda::effectiveSize(mcmc_1)
cov(pmcmc_desktop_1$pars)
summary(mcmc_1)
windows(60,50)
plot(mcmc_1)

pmcmc_desktop_2 <- readRDS('./pmcmc-run/pmcmc_desktop_2.RDS')
mcmc_2 <- coda::as.mcmc(cbind(pmcmc_desktop_2$probabilities, pmcmc_desktop_2$pars),
)
1 - coda::rejectionRate(mcmc_2)
coda::effectiveSize(mcmc_2)
cov(pmcmc_desktop_2$pars)
summary(mcmc_2)
windows(60,50)
plot(mcmc_2)
min(pmcmc_desktop_2$pars$init_EIR)
history <- pmcmc_desktop_2$trajectories$state
plot_particle_filter(history,true_history=data_raw_cmis,times=data_raw_cmis$t)

pmcmc_desktop_3 <- readRDS('./pmcmc-run/pmcmc_desktop_3.RDS')
mcmc_3 <- coda::as.mcmc(cbind(pmcmc_desktop_3$probabilities, pmcmc_desktop_3$pars),
)
1 - coda::rejectionRate(mcmc_3)
coda::effectiveSize(mcmc_3)
cov(pmcmc_desktop_3$pars)
summary(mcmc_3)
windows(60,50)
plot(mcmc_3)
min(pmcmc_desktop_3$pars$init_EIR)
history <- pmcmc_desktop_3$trajectories$state

plot_particle_filter(history,true_history=data_raw_cmis,times=data_raw_cmis$t)

pmcmc_desktop_all <- readRDS('./pmcmc-run/pmcmc_desktop_all.RDS')
mcmc_all <- coda::as.mcmc(cbind(pmcmc_desktop_all$probabilities, pmcmc_desktop_all$pars),
)
1 - coda::rejectionRate(mcmc_all)
coda::effectiveSize(mcmc_all)
cov(pmcmc_desktop_all$pars)
summary(mcmc_all)
windows(60,50)
plot(mcmc_all)
min(pmcmc_desktop_all$pars$init_EIR)
history <- pmcmc_desktop_all$trajectories$state
plot_particle_filter(history,true_history=data_raw_cmis_all,times=data_raw_cmis_all$t)
matplot(data_raw_cmis$t, t(history['EIR', , -1]), type = "l",
        xlab = "Time", ylab = "EIR",
        col = "#A6CEE3", lty = 1, ylim = range(history['EIR', , -1]))
matplot(data_raw_cmis$t, t(history['inc', , -1]), type = "l",
        xlab = "Time", ylab = "Clinical Incidence",
        col = "#A6CEE3", lty = 1, ylim = range(history['inc', , -1]))


pmcmc_desktop_exp <- readRDS('./pmcmc-run/pmcmc_desktop_exp.RDS')
mcmc_exp <- coda::as.mcmc(cbind(pmcmc_desktop_exp$probabilities, pmcmc_desktop_exp$pars),
)
1 - coda::rejectionRate(mcmc_exp)
coda::effectiveSize(mcmc_exp)
cov(pmcmc_desktop_exp$pars)
summary(mcmc_exp)
windows(60,50)
plot(mcmc_exp)
min(pmcmc_desktop_exp$pars)
history <- pmcmc_desktop_exp$trajectories$state
plot_particle_filter(history,true_history=data_raw_cmis_all,times=data_raw_cmis_all$t)
matplot(data_raw_cmis$t, t(history['EIR', , -1]), type = "l",
        xlab = "Time", ylab = "EIR",
        col = "#A6CEE3", lty = 1, ylim = range(history['EIR', , -1]))
matplot(data_raw_cmis$t, t(history['inc', , -1]), type = "l",
        xlab = "Time", ylab = "Clinical Incidence",
        col = "#A6CEE3", lty = 1, ylim = range(history['inc', , -1]))

ng_test <- readRDS('./pmcmc-run/ng_test.RDS')
ng_mcmc <- coda::as.mcmc(cbind(ng_test$probabilities, ng_test$pars),
)
1 - coda::rejectionRate(ng_mcmc)
coda::effectiveSize(ng_mcmc)
cov(ng_test$pars)
summary(ng_mcmc)
windows(60,50)
plot(ng_mcmc)
min(ng_test$pars$init_EIR)
history <- ng_test$trajectories$state
plot_particle_filter(history,true_history=data_raw_ng,times=data_raw_ng$t)

bf_test <- readRDS('./pmcmc-run/bf_test.RDS')
bf_mcmc <- coda::as.mcmc(cbind(bf_test$probabilities, bf_test$pars),
)
1 - coda::rejectionRate(bf_mcmc)
coda::effectiveSize(bf_mcmc)
cov(bf_test$pars)
summary(bf_mcmc)
windows(60,50)
plot(bf_mcmc)
min(bf_test$pars$init_EIR)
history <- bf_test$trajectories$state
plot_particle_filter(history,true_history=data_raw_bf,times=data_raw_bf$t)

matplot(times, t(history[1, , -1]), type = "l",
        xlab = "Time", ylab = "Prevalence",
        col = "#A6CEE3", lty = 1, ylim = range(cis[,4:6]))
##Figures##
library(reshape2)
library(ggplot2)
library(ggpubr)
library(zoo)
library(patchwork)

theme_set(theme_minimal()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.border = element_rect(colour = "black",fill=NA)))
setwd('../pmcmc-run')
cmis_good1 <- readRDS('cmis_good1.RDS')
data_raw_cmis_all <- readRDS('./data_raw_cmis_all.RDS')
cmis <- cmis_good1$trajectories$state
cmis_prev_traj <- cbind(as.data.frame(t(cmis[1, , -1])),data_raw_cmis_all$t)%>%
  rename(t=`data_raw_cmis_all$t`)%>%
  melt(id='t')
cmis_data_cis <- addCIs(data_raw_cmis_all,data_raw_cmis_all$positive,data_raw_cmis_all$tested)

cmis_data_cis$date <- (cmis_data_cis$t-30)+as.Date("2015-04-08")
cmis_prev_traj$date <- (cmis_prev_traj$t-30)+as.Date("2015-04-08")

cmis_prev <- ggplot(cmis_prev_traj)+
  geom_line(aes(x=date,y=value,group=variable),col = "#A6CEE3",alpha=0.2)+
  geom_point(data=cmis_data_cis,aes(x=date,y=mean),pch = 19,
             col = "#666666")+
  geom_errorbar(data=cmis_data_cis,aes(x=date,ymin=lower,ymax=upper),width = 0,
                col = "#666666")+
  scale_y_continuous(limits=c(0,max(cmis_data_cis$upper)))+
  scale_x_date(date_labels = "%b %y")+
  labs(title='Western Kenya',x='Date',y='RDT Prevalence')
cmis_prev

cmis_eir_traj <- cbind(as.data.frame(t(cmis[2, , -1])),data_raw_cmis_all$t)%>%
  rename(t=`data_raw_cmis_all$t`)%>%
  melt(id='t')
cmis_eir_traj$date <- (cmis_eir_traj$t-30)+as.Date("2015-04-08")
cmis_eir_med <- cmis_eir_traj %>%
  group_by(date) %>%
  summarise(median(value))%>%
  rename(median=`median(value)`)

cmis_eir <- ggplot(cmis_eir_traj)+
  geom_line(aes(x=date,y=value,group=variable),col = "#A6CEE3",alpha=0.2)+
  geom_line(data=cmis_eir_med,aes(x=date,y=median),col = "#1F78B4")+
  scale_y_continuous(limits=c(0,max(cmis_eir_med$median)))+
  scale_x_date(date_labels = "%b %y")+
  labs(x='Date',y='EIR')
cmis_eir

cmis_inc_traj <- cbind(as.data.frame(t(cmis[3, , -1])),data_raw_cmis_all$t)%>%
  rename(t=`data_raw_cmis_all$t`)%>%
  melt(id='t')
cmis_inc_traj$date <- (cmis_inc_traj$t-30)+as.Date("2015-04-08")
cmis_inc_med <- cmis_inc_traj %>%
  group_by(date) %>%
  summarise(median(value))%>%
  rename(median=`median(value)`)
cmis_inc <- ggplot(cmis_inc_traj)+
  geom_line(aes(x=date,y=value,group=variable),col = "#A6CEE3",alpha=0.2)+
  geom_line(data=cmis_inc_med,aes(x=date,y=median),col = "#1F78B4")+
  scale_y_continuous(limits=c(0,max(cmis_inc_med$median)))+
  scale_x_date(date_labels = "%b %y")+
  labs(x='Date',y='Clinical Incidence')
cmis_inc

cmis <- ggarrange(cmis_prev,cmis_eir,cmis_inc,
          ncol = 1, nrow = 3)
ng_good <- readRDS('../ng_good_1.RDS')
ng <- ng_good$trajectories$state
data_raw_ng <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/data_raw_NG_anc.RDS')

ng_prev_traj <- cbind(as.data.frame(t(ng[1, , -1])),data_raw_ng$t)%>%
  rename(t=`data_raw_ng$t`)%>%
  melt(id='t')
ng_data_cis <- addCIs(data_raw_ng,data_raw_ng$positive,data_raw_ng$tested)
ng_data_cis$date <- (ng_data_cis$t-30)+as.Date("2020-11-01")
ng_prev_traj$date <- (ng_prev_traj$t-30)+as.Date("2020-11-01")

ng_prev <- ggplot(ng_prev_traj)+
  geom_line(aes(x=date,y=value,group=variable),col = "#A6CEE3",alpha=0.2)+
  geom_point(data=ng_data_cis,aes(x=date,y=mean),pch = 19,
             col = "#666666")+
  geom_errorbar(data=ng_data_cis,aes(x=date,ymin=lower,ymax=upper),width = 0,
                col = "#666666")+
  scale_y_continuous(limits=c(0,max(c(ng_data_cis$upper,ng_prev_traj$value))))+
  scale_x_date(date_labels = "%b %y")+
  labs(title='Nigeria',x='Date',y='RDT Prevalence')
ng_prev

ng_eir_traj <- cbind(as.data.frame(t(ng[2, , -1])),data_raw_ng$t)%>%
  rename(t=`data_raw_ng$t`)%>%
  melt(id='t')
ng_eir_traj$date <- (ng_eir_traj$t-30)+as.Date("2020-11-01")
ng_eir_med <- ng_eir_traj %>%
  group_by(date) %>%
  summarise(median(value))%>%
  rename(median=`median(value)`)
ng_eir <- ggplot(ng_eir_traj)+
  geom_line(aes(x=date,y=value,group=variable),col = "#A6CEE3",alpha=0.2)+
  geom_line(data=ng_eir_med,aes(x=date,y=median),col = "#1F78B4")+
  scale_y_continuous(limits=c(0,max(c(ng_eir_traj$value,ng_eir_med$median))))+
  scale_x_date(date_labels = "%b %y")+
  labs(x='Date',y='EIR')
ng_eir

ng_inc_traj <- cbind(as.data.frame(t(ng[3, , -1])),data_raw_ng$t)%>%
  rename(t=`data_raw_ng$t`)%>%
  melt(id='t')
ng_inc_traj$date <- (ng_inc_traj$t-30)+as.Date("2020-11-01")
ng_inc_med <- ng_inc_traj %>%
  group_by(date) %>%
  summarise(median(value))%>%
  rename(median=`median(value)`)
ng_inc <- ggplot(ng_inc_traj)+
  geom_line(aes(x=date,y=value,group=variable),col = "#A6CEE3",alpha=0.2)+
  geom_line(data=ng_inc_med,aes(x=date,y=median),col = "#1F78B4")+
  scale_y_continuous(limits=c(0,max(c(ng_inc_traj$value,ng_inc_med$median))))+
  scale_x_date(date_labels = "%b %y")+
  labs(x='Date',y='Clinical Incidence')
ng_inc

ng_plot <- ggarrange(ng_prev,ng_eir,ng_inc,
                  ncol = 1, nrow = 3)
ng_plot

bf_good <- readRDS('../bf_good1.RDS')
bf <- bf_good$trajectories$state
data_raw_bf <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/data_raw_BF_anc.RDS')

bf_prev_traj <- cbind(as.data.frame(t(bf[1, , -1])),data_raw_bf$t)%>%
  rename(t=`data_raw_bf$t`)%>%
  melt(id='t')
bf_data_cis <- addCIs(data_raw_bf,data_raw_bf$positive,data_raw_bf$tested)
bf_data_cis$date <- (bf_data_cis$t-30)+as.Date("2020-09-01")
bf_prev_traj$date <- (bf_prev_traj$t-30)+as.Date("2020-09-01")

bf_prev <- ggplot(bf_prev_traj)+
  annotate("rect", xmin = min(bf_prev_traj$date), xmax = as.Date('2020-10-31'), ymin = 0, ymax = max(c(bf_data_cis$upper,bf_prev_traj$value)),alpha = .1,fill = "#0072B2")+
  annotate("rect", xmin = as.Date('2021-06-01'), xmax = as.Date('2021-10-31'), ymin = 0, ymax = max(c(bf_data_cis$upper,bf_prev_traj$value)),alpha = .1,fill = "#0072B2")+
  geom_line(aes(x=date,y=value,group=variable),col = "#A6CEE3",alpha=0.2)+
  geom_point(data=bf_data_cis,aes(x=date,y=mean),pch = 19,
             col = "#666666")+
  geom_errorbar(data=bf_data_cis,aes(x=date,ymin=lower,ymax=upper),width = 0,
                col = "#666666")+
  scale_y_continuous(limits=c(0,max(c(bf_data_cis$upper,bf_prev_traj$value))))+
  scale_x_date(date_labels = "%b %y")+
  labs(title='Burkina Faso',x='Date',y='RDT Prevalence')
bf_prev

bf_eir_traj <- cbind(as.data.frame(t(bf[2, , -1])),data_raw_bf$t)%>%
  rename(t=`data_raw_bf$t`)%>%
  melt(id='t')
bf_eir_traj$date <- (bf_eir_traj$t-30)+as.Date("2020-09-01")
bf_eir_med <- bf_eir_traj %>%
  group_by(date) %>%
  summarise(median(value))%>%
  rename(median=`median(value)`)
bf_eir <- ggplot(bf_eir_traj)+
  annotate("rect", xmin = min(bf_eir_traj$date), xmax = as.Date('2020-10-31'), ymin = 0, ymax = max(c(bf_eir_traj$value,bf_eir_med$median)),alpha = .1,fill = "#0072B2")+
  annotate("rect", xmin = as.Date('2021-06-01'), xmax = as.Date('2021-10-31'), ymin = 0, ymax = max(c(bf_eir_traj$value,bf_eir_med$median)),alpha = .1,fill = "#0072B2")+
  geom_line(aes(x=date,y=value,group=variable),col = "#A6CEE3",alpha=0.2)+
  geom_line(data=bf_eir_med,aes(x=date,y=median),col = "#1F78B4")+
  scale_y_continuous(limits=c(0,max(c(bf_eir_traj$value,bf_eir_med$median))))+
  scale_x_date(date_labels = "%b %y")+
  labs(x='Date',y='EIR')
bf_eir

bf_inc_traj <- cbind(as.data.frame(t(bf[3, , -1])),data_raw_bf$t)%>%
  rename(t=`data_raw_bf$t`)%>%
  melt(id='t')
bf_inc_traj$date <- (bf_inc_traj$t-30)+as.Date("2020-09-01")
bf_inc_med <- bf_inc_traj %>%
  group_by(date) %>%
  summarise(median(value))%>%
  rename(median=`median(value)`)
bf_inc <- ggplot(bf_inc_traj)+
  annotate("rect", xmin = min(bf_inc_traj$date), xmax = as.Date('2020-10-31'), ymin = 0, ymax = max(c(bf_inc_traj$value,bf_inc_med$median)),alpha = .1,fill = "#0072B2")+
  annotate("rect", xmin = as.Date('2021-06-01'), xmax = as.Date('2021-10-31'), ymin = 0, ymax = max(c(bf_inc_traj$value,bf_inc_med$median)),alpha = .1,fill = "#0072B2")+
  geom_line(aes(x=date,y=value,group=variable),col = "#A6CEE3",alpha=0.2)+
  geom_line(data=bf_inc_med,aes(x=date,y=median),col = "#1F78B4")+
  scale_y_continuous(limits=c(0,max(c(bf_inc_traj$value,bf_inc_med$median))))+
  labs(x='Date',y='Clinical Incidence')+
  scale_x_date(date_labels = "%b %y")
bf_inc

bf_plot <- ggarrange(bf_prev,bf_eir,bf_inc,
                     ncol = 1, nrow = 3)
bf_plot
windows(90,50)
(cmis_prev/cmis_eir/cmis_inc) + (ng_prev/ng_eir/ng_inc) + (bf_prev/bf_eir/bf_inc) + plot_layout(ncol=3,nrow=3)
cmis_prev+ng_prev+bf_prev + cmis_eir+ng_eir+bf_eir + cmis_inc+ng_inc+bf_inc + plot_layout(ncol=3,nrow=3)
ggarrange(cmis,ng_plot,bf_plot,
          ncol = 3, nrow = 1)
windows(90,50)
ggarrange(cmis,ng_plot,bf_plot,
          ncol = 3, nrow = 1)
sim_mcmc <- coda::as.mcmc(cbind(sim_test$probabilities, sim_test$pars),
)
1 - coda::rejectionRate(sim_mcmc)
coda::effectiveSize(sim_mcmc)
cov(sim_test$pars)
summary(bf_mcmc)
windows(60,50)
plot(sim_mcmc)
history <- sim_test$trajectories$state
plot_particle_filter(history,true_history=data_raw,times=data_raw$t)

sim <- sim_test$trajectories$state
simprev_traj <- cbind(as.data.frame(t(sim[1, , -1])),data_raw$t)%>%
  rename(t=`data_raw$t`)%>%
  melt(id='t')
simdata_cis <- addCIs(data_raw,data_raw$positive,data_raw$tested)
simdata_cis$date <- (simdata_cis$t-30)+as.Date("2020-09-01")
simprev_traj$date <- (simprev_traj$t-30)+as.Date("2020-09-01")

simprev <- ggplot(simprev_traj)+
  geom_line(aes(x=t,y=value,group=variable),col = "#A6CEE3",alpha=0.2)+
  geom_point(data=simdata_cis,aes(x=t,y=mean),pch = 19,
             col = "black")+
  geom_errorbar(data=simdata_cis,aes(x=t,ymin=lower,ymax=upper),width = 0,
                col = "black")+
  labs(x='Day',y='RDT Prevalence')
simprev

simeir_traj <- cbind(as.data.frame(t(sim[2, , -1])),data_raw$t)%>%
  rename(t=`data_raw$t`)%>%
  melt(id='t')
simeir_traj$date <- (simeir_traj$t-30)+as.Date("2020-09-01")
simeir_med <- simeir_traj %>%
  group_by(t) %>%
  summarise(median(value))%>%
  rename(median=`median(value)`)
simeir <- ggplot(simeir_traj)+
  geom_line(aes(x=t,y=value,group=variable),col = "#FB9A99",alpha=0.2)+
  geom_line(data=simeir_med,aes(x=t,y=median),col = "#E31A1C")+
  scale_y_continuous(limits=c(0,max(simeir_med$median)))+
  labs(x='Day',y='EIR')
simeir

siminc_traj <- cbind(as.data.frame(t(sim[3, , -1])),data_raw$t)%>%
  rename(t=`data_raw$t`)%>%
  melt(id='t')
siminc_traj$date <- (siminc_traj$t-30)+as.Date("2020-09-01")
siminc_med <- siminc_traj %>%
  group_by(t) %>%
  summarise(median(value))%>%
  rename(median=`median(value)`)
siminc <- ggplot(siminc_traj)+
  geom_line(aes(x=t,y=value,group=variable),col = "#FB9A99",alpha=0.2)+
  geom_line(data=siminc_med,aes(x=t,y=median),col = "#E31A1C")+
  scale_y_continuous(limits=c(0,max(siminc_med$median)))+
  labs(x='Day',y='Clinical Incidence')
siminc
windows(90,30)
simplot <- ggarrange(simprev,simeir,siminc,
                     ncol = 3, nrow = 1)
simplot


##Cluster runs - prevalence by age
pmcmc_run1 <- readRDS('cmis_all_run1.RDS')
pmcmc_run3 <- readRDS('cmis_all_run3.RDS')
history_run1 <- pmcmc_run1$trajectories$state
history_run3 <- pmcmc_run3$trajectories$state

#[agegroup,mcmc sample,time]
df_ages_run1 <- data.frame(agegroup=numeric(),
                      month=numeric(),
                      value=numeric())
for(age in 4:24){
  for(time in 2:62){
    temp <- data.frame(agegroup=age-3,
                       month=time,
                       value=history_run1[age,101:1000,time])
    df_ages_run1 <- rbind(df_ages_run1,temp)
  }
}
df_ages_run3 <- data.frame(agegroup=numeric(),
                           month=numeric(),
                           value=numeric())
for(age in 4:24){
  for(time in 2:62){
    temp <- data.frame(agegroup=age-3,
                       month=time,
                       value=history_run3[age,101:1000,time])
    df_ages_run3 <- rbind(df_ages,temp)
  }
}
df_ages_both <- rbind(df_ages_run1,df_ages_run3)
rm(df_ages_run1,df_ages_run3)

windows(160,90)
ggplot(df_ages_both,aes(x=as.factor(agegroup),y=value))+
  geom_boxplot()+
  facet_wrap(~month)

ggplot(df_ages[df_ages$month<=13,],aes(x=as.factor(agegroup),y=value))+
  geom_boxplot()+
  facet_wrap(~month)

##By Age group fitting
mcmc_1 <- coda::as.mcmc(cbind(pmcmc_run$probabilities, pmcmc_run$pars),
)
1 - coda::rejectionRate(mcmc_1)
coda::effectiveSize(mcmc_1)
cov(pmcmc_run$pars)
summary(mcmc_1)
windows(60,50)
plot(mcmc_1)
history <- pmcmc_run$trajectories$state

df_ages_run <- data.frame(agegroup=numeric(),
                           month=numeric(),
                           value=numeric())
for(age in 4:28){
  for(time in 2:62){
    temp <- data.frame(agegroup=age-3,
                       month=time,
                       value=history[age,21:100,time])
    df_ages_run <- rbind(df_ages_run,temp)
  }
}

windows(160,90)
theme_set(theme_minimal()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.border = element_rect(colour = "black",fill=NA)))


ggplot(df_ages_run,aes(x=as.factor(agegroup),y=value))+
  geom_boxplot()+
  facet_wrap(~month)
ggplot(df_ages_run[df_ages_run$month==2,],aes(x=as.factor(agegroup),y=value))+
  geom_boxplot()
