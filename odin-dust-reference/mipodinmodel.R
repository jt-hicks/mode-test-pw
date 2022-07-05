
## MODEL VARIABLES
##------------------------------------------------------------------------------

na <- user() # number of age categories
nh <- user() # number of biting heterogeneity categories
ft <- user() # proportion of cases treated

##------------------------------------------------------------------------------
##################
## HUMAN STATES ##
##################
##------------------------------------------------------------------------------


# Human states as specified in full transmission model
# http://journals.plos.org/plosmedicine/article?id=10.1371%2Fjournal.pmed.1000324
# http://www.nature.com/articles/ncomms4136

# fitted parameters for human compartments:
eta <- user() # death rate for exponential population distribution
age_rate[] <- user() # rate at which humans move through age categories
dim(age_rate) <- na
het_wt[] <- user() # weights of the different heterogeneous biting categories
dim(het_wt) <- nh
rA <- user() # rate of movement from A -> U
rT <- user() # rate of treatment working: T -> P
rD <- user() #  rate from D -> A
rU <- user() # rate of clearance of subpatent infection U -> S
rP <- user() # rate at which prophylaxis wears off P -> S

# S - SUSCEPTIBLE
init_S[,] <- user()
dim(init_S) <- c(na,nh)
initial(S[,]) <- init_S[i,j]
dim(S) <- c(na,nh)

deriv(S[1, 1:nh]) <- -(lag_rates/dE)*FOI[i,j,lag_rates]*S[i,j] + rP*P[i,j] + rU*U[i,j] +
  eta*H*het_wt[j] - (eta+age_rate[i])*S[i,j]
deriv(S[2:na, 1:nh]) <- -(lag_rates/dE)*FOI[i,j,lag_rates]*S[i,j] + rP*P[i,j] + rU*U[i,j] -
  (eta+age_rate[i])*S[i,j] + age_rate[i-1]*S[i-1,j]

# T- SUCCESSFULLY TREATED
init_T[,] <- user()
dim(init_T) <- c(na,nh)
initial(T[,]) <- init_T[i,j]
dim(T) <- c(na,nh)

deriv(T[1, 1:nh]) <- ft*clin_inc[i,j] - rT*T[i,j] -
  (eta+age_rate[i])*T[i,j]
deriv(T[2:na, 1:nh]) <- ft*clin_inc[i,j] - rT*T[i,j] -
  (eta+age_rate[i])*T[i,j] + age_rate[i-1]*T[i-1,j]

# D - CLEAR DISEASE
init_D[,] <- user()
dim(init_D) <- c(na,nh)
initial(D[,]) <- init_D[i,j]
dim(D) <- c(na,nh)

deriv(D[1, 1:nh]) <- (1-ft)*clin_inc[i,j] - rD*D[i,j] -
  (eta+age_rate[i])*D[i,j]
deriv(D[2:na, 1:nh]) <- (1-ft)*clin_inc[i,j] - rD*D[i,j] -
  (eta+age_rate[i])*D[i,j] + age_rate[i-1]*D[i-1,j]

# A - ASYMPTOMATIC DISEASE
init_A[,] <- user()
dim(init_A) <- c(na,nh)
initial(A[,]) <- init_A[i,j]
dim(A) <- c(na,nh)

deriv(A[1, 1:nh]) <- (1-phi[i,j])*(lag_rates/dE)*FOI[i,j,lag_rates]*Y[i,j] - (lag_rates/dE)*FOI[i,j,lag_rates]*A[i,j] +
  rD*D[i,j] - rA*A[i,j] - (eta+age_rate[i])*A[i,j]
deriv(A[2:na, 1:nh]) <- (1-phi[i,j])*(lag_rates/dE)*FOI[i,j,lag_rates]*Y[i,j] - (lag_rates/dE)*FOI[i,j,lag_rates]*A[i,j] +
  rD*D[i,j] - rA*A[i,j] - (eta+age_rate[i])*A[i,j] + age_rate[i-1]*A[i-1,j]

# U - SUBPATENT DISEASE
init_U[,] <- user()
dim(init_U) <- c(na,nh)
initial(U[,]) <- init_U[i,j]
dim(U) <- c(na,nh)

deriv(U[1, 1:nh]) <- rA*A[i,j] - (lag_rates/dE)*FOI[i,j,lag_rates]*U[i,j] - rU*U[i,j] -
  (eta+age_rate[i])*U[i,j]
deriv(U[2:na, 1:nh]) <- rA*A[i,j] - (lag_rates/dE)*FOI[i,j,lag_rates]*U[i,j] - rU*U[i,j] -
  (eta+age_rate[i])*U[i,j] + age_rate[i-1]*U[i-1,j]

# P - PROPHYLAXIS
init_P[,] <- user()
dim(init_P) <- c(na,nh)
initial(P[,]) <- init_P[i,j]
dim(P) <- c(na,nh)

deriv(P[1, 1:nh]) <- rT*T[i,j] - rP*P[i,j] - (eta+age_rate[i])*P[i,j]
deriv(P[2:na, 1:nh]) <- rT*T[i,j] - rP*P[i,j] - (eta+age_rate[i])*P[i,j] +
  age_rate[i-1]*P[i-1,j]

# The number of individuals able to acquire clinical malaria
dim(Y) <- c(na,nh)
Y[1:na, 1:nh] <- S[i,j]+A[i,j]+U[i,j]

# The number of new cases at this timestep
dim(clin_inc) <- c(na,nh)
clin_inc[1:na, 1:nh] <- phi[i,j]*(lag_rates/dE)*FOI[i,j,lag_rates]*Y[i,j]

# Sum compartments over all age, heterogeneity and intervention categories
Sh <- sum(S[,])
Th <- sum(T[,])
Dh <- sum(D[,])
Ah <- sum(A[,])
Uh <- sum(U[,])
Ph <- sum(P[,])
H <- Sh + Th + Dh + Ah + Uh + Ph


##------------------------------------------------------------------------------
#####################
## PREGNANCY STATES ##
######################
##------------------------------------------------------------------------------

#Gravidity states: calculates risk of previous exposure to MiP
#Will be used to calculate parasite detectability depending on a woman's pregnancy history
nrates<-user()
sample_rates[]<-user()
sample_transition_rates[]<-user()
wane_rates <- user()
dim(sample_rates)<-nrates
dim(sample_transition_rates)<-nrates

dim(pregs)<-nrates
init_pregs[] <- user()
dim(init_pregs) <- nrates
initial(pregs[])<-init_pregs[i]

deriv(pregs[1])<-prev_pcr*sample_rates[i] -
  (sample_transition_rates[i]+wane_rates)*pregs[i]
deriv(pregs[2:nrates])<-prev_pcr*sample_rates[i] +
  sample_transition_rates[i-1]*pregs[i-1] - (sample_transition_rates[i]+wane_rates)*pregs[i]
previous_pregs<-pregs[nrates]
#output(previous_pregs)<-previous_pregs

# index of the requested age vector min
agestart <- user()
# index of the requested age vector just prior to the maximum value for the prevalence calculation
ageend <- user()
# Define pregnancy-specific A and U recovery rates
rA_preg <- user()
rU_preg <- user()
# Specify number of age categories that are child-bearing
preg_range <- (ageend - agestart) + 1

# The number of pregnant individuals able to acquire clinical malaria
dim(Y_preg) <- c(preg_range,nh)
Y_preg[1:preg_range, 1:nh] <- S_preg[i,j]+A_preg[i,j]+U_preg[i,j]

# The number of new pregnant cases at this timestep
dim(clin_inc_preg) <- c(preg_range,nh)
clin_inc_preg[1:preg_range, 1:nh] <- phi[as.integer((i+agestart)-1),j]*(lag_rates/dE)*FOI[as.integer((i+agestart)-1),j,lag_rates]*Y_preg[i,j]

init_S_preg[,] <- user()
dim(S_preg) <- c(preg_range,nh)
initial(S_preg[1:preg_range,1:nh]) <- init_S_preg[i,j]
dim(init_S_preg) <- c(preg_range,nh)

deriv(S_preg[1:preg_range,1:nh]) <- - (lag_rates/dE)*FOI[as.integer((i+agestart)-1),j,lag_rates]*S_preg[i,j] +
  rP*P_preg[i,j] + rU_preg*U_preg[i,j]

init_T_preg[,] <- user()
dim(T_preg) <- c(preg_range,nh)
initial(T_preg[1:preg_range,1:nh]) <- init_T_preg[i,j]
dim(init_T_preg) <- c(preg_range,nh)

deriv(T_preg[1:preg_range,1:nh]) <- ft*clin_inc_preg[i,j] - rT*T_preg[i,j]

init_D_preg[,] <- user()
dim(D_preg) <- c(preg_range,nh)
initial(D_preg[1:preg_range,1:nh]) <- init_D_preg[i,j]
dim(init_D_preg) <- c(preg_range,nh)

deriv(D_preg[1:preg_range,1:nh]) <- (1-ft)*clin_inc_preg[i,j] - rD*D_preg[i,j]

init_A_preg[,] <- user()
dim(A_preg) <- c(preg_range,nh)
initial(A_preg[1:preg_range,1:nh]) <- init_A_preg[i,j]
dim(init_A_preg) <- c(preg_range,nh)

deriv(A_preg[1:preg_range,1:nh]) <- (1-phi[as.integer((i+agestart)-1),j])*(lag_rates/dE)*FOI[as.integer((i+agestart)-1),j,lag_rates]*Y_preg[i,j] -
  (lag_rates/dE)*FOI[as.integer((i+agestart)-1),j,lag_rates]*A_preg[i,j] + rD*D_preg[i,j] - rA_preg*A_preg[i,j]

init_U_preg[,] <- user()
dim(U_preg) <- c(preg_range,nh)
initial(U_preg[1:preg_range,1:nh]) <- init_U_preg[i,j]
dim(init_U_preg) <- c(preg_range,nh)

deriv(U_preg[1:preg_range,1:nh]) <- rA_preg*A_preg[i,j] -
  (lag_rates/dE)*FOI[as.integer((i+agestart)-1),j,lag_rates]*U_preg[i,j] - rU_preg*U_preg[i,j]

init_P_preg[,] <- user()
dim(P_preg) <- c(preg_range,nh)
initial(P_preg[1:preg_range,1:nh]) <- init_P_preg[i,j]
dim(init_P_preg) <- c(preg_range,nh)

deriv(P_preg[1:preg_range,1:nh]) <- rT*T_preg[i,j] - rP*P_preg[i,j]

##------------------------------------------------------------------------------
#####################
## IMMUNITY STATES ##
#####################
##------------------------------------------------------------------------------

# See supplementary materials S1 from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6

# ICM - Maternally acquired immunity acquired by babies from mothers (assumed to be proportional to the immunity of a 15 to 30 year old woman)
# ICA - Immunity acquired due to exposure to previous infection, increases with age
# IC - Clinical immunity. Upon infection, immunity against clinical case. IC = ICA + ICM
# IB - Infection blocking immunity, chances of preventing infection upon receiving infectious bite
# ID - Detection immunity, when immunity suppresses parasite densities this makes it less likely that diagnostics will detect parasite infection

# fitted immunity parameters:
dCM <- user() # decay of maternal immunity
uCA <- user() # scale parameter (see Supplementary mats. 3.1.2)
dCA <- user() # decay for clinical immunity
dB <- user() # decay for infection blocking immunity
uB <- user() # scale param for IB immunity
dID <- user() # decay for detection immunity
uD <- user() # scale param for ID immunity
x_I[] <- user() # intermediate variable for calculating immunity functions
dim(x_I) <- na
age20l <- user(integer=TRUE) # lower index of age 20 age compartment
age20u <- user(integer=TRUE) # upper index of age 20 age compartment
age_20_factor <- user() # factor calculated in equilibrium solution
PM <- user() # immunity constant

# ICM - maternally acquired immunity
init_ICM[,] <- user()
dim(init_ICM) <- c(na,nh)
initial(ICM[,]) <- init_ICM[i,j]
dim(ICM) <- c(na,nh)
dim(init_ICM_pre) <- c(nh)
init_ICM_pre[1:nh] <- PM*(ICA[age20l,i] + age_20_factor*(ICA[age20u,i]-ICA[age20l,i]))
output(ICM_pre) <- init_ICM_pre[1]

deriv(ICM[1, 1:nh]) <- -1/dCM*ICM[i,j] + (init_ICM_pre[j]-ICM[i,j])/x_I[i]
deriv(ICM[2:na, 1:nh]) <- -1/dCM*ICM[i,j] - (ICM[i,j]-ICM[i-1,j])/x_I[i]

# ICA - exposure driven immunity
init_ICA[,] <- user()
dim(init_ICA) <- c(na,nh)
initial(ICA[,]) <- init_ICA[i,j]
dim(ICA) <- c(na,nh)

deriv(ICA[1, 1:nh]) <- (lag_rates/dE)*FOI[i,j,lag_rates]/((lag_rates/dE)*FOI[i,j,lag_rates] * uCA + 1) - 1/dCA*ICA[i,j] -ICA[i,j]/x_I[i]
deriv(ICA[2:na, 1:nh]) <- (lag_rates/dE)*FOI[i,j,lag_rates]/((lag_rates/dE)*FOI[i,j,lag_rates] * uCA + 1) - 1/dCA*ICA[i,j] - (ICA[i,j]-ICA[i-1,j])/x_I[i]

# clinical immunity is a combination of maternal and exposure-driven immunity
dim(IC) <- c(na,nh)
IC[,] <- ICM[i,j] + ICA[i,j]

# phi - probability of clinical disease, dependent on clinical immunity
phi0 <- user()
phi1 <- user() # these parameters characterise the hill function
IC0 <- user() # for probability of clinical disease
kC <- user() # See supplementary materials 1.1.3
dim(phi) <- c(na,nh)
phi[1:na,1:nh] <- phi0*((1-phi1)/(1+(IC[i,j]/IC0)^kC) + phi1)

# IB - infection blocking immunity
init_IB[,] <- user()
dim(init_IB) <- c(na,nh)
initial(IB[,]) <- init_IB[i,j]
dim(IB) <- c(na,nh)

deriv(IB[1, 1:nh]) <- EIR[i,j]/(EIR[i,j]* uB + 1) - IB[i,j]/dB - IB[i,j]/x_I[i]
deriv(IB[2:na, 1:nh]) <- EIR[i,j]/(EIR[i,j]* uB + 1) - IB[i,j]/dB - (IB[i,j]-IB[i-1,j])/x_I[i]

# b - probability of disease from infectious bite, depends on infection blocking immunity
b0 <- user() # these parameters characterise the hill function for b
b1 <- user() # prob of infection from bite with zero immunity
kB <- user() #
IB0 <- user()
dim(b) <- c(na,nh)
b[1:na, 1:nh] <- b0 * ((1-b1)/(1+(IB[i,j]/IB0)^kB)+b1)

# detection immunity
init_ID[,] <- user()
dim(init_ID) <- c(na,nh)
initial(ID[,]) <- init_ID[i,j]
dim(ID) <- c(na,nh)

deriv(ID[1, 1:nh]) <- (lag_rates/dE)*FOI[i,j,lag_rates]/((lag_rates/dE)*FOI[i,j,lag_rates]*uD + 1) - ID[i,j]/dID - ID[i,j]/x_I[i]
deriv(ID[2:na, 1:nh]) <- (lag_rates/dE)*FOI[i,j,lag_rates]/((lag_rates/dE)*FOI[i,j,lag_rates]*uD + 1) - ID[i,j]/dID - (ID[i,j]-ID[i-1,j])/x_I[i]

# p_det - probability of detection by microscopy, immunity decreases chances of
# infection because it pushes parasite density down
aD <- user()
fD0 <- user()
gammaD <- user()
d1 <- user()
ID0 <- user()
kD <- user()
dim(age) <- na
age[] <- user() # vector of age categories supplied by user

dim(fd) <- na
fd[1:na] <- 1-(1-fD0)/(1+(age[i]/aD)^gammaD)
dim(p_det) <- c(na,nh)
p_det[,] <- d1 + (1-d1)/(1 + fd[i]*(ID[i,j]/ID0)^kD)

# Force of infection, depends on level of infection blocking immunity
dim(FOI_lag) <- c(na,nh)
FOI_lag[1:na, 1:nh] <- EIR[i,j] * (if(IB[i,j]==0) b0 else b[i,j])

# Current FOI depends on humans that have been through the latent period
dE <- user() # latent period of human infection.
lag_rates <- user()

FOI_eq[,] <- user()
dim(FOI_eq) <- c(na,nh)
init_FOI[,,] <- FOI_eq[i,j]*dE/lag_rates
dim(init_FOI) <- c(na,nh,lag_rates)
initial(FOI[,,]) <- init_FOI[i,j,k]
dim(FOI) <- c(na,nh,lag_rates)

deriv(FOI[,,1]) <- FOI_lag[i,j] - (lag_rates/dE)*FOI[i,j,1]
deriv(FOI[,,2:lag_rates]) <- (lag_rates/dE)*FOI[i,j,k-1] - (lag_rates/dE)*FOI[i,j,k]

# EIR -rate at which each age/het/int group is bitten
# rate for age group * rate for biting category * FOI for age group * prop of
# infectious mosquitoes
dim(foi_age) <- na
foi_age[] <- user()
dim(rel_foi) <- nh
rel_foi[] <- user()
av0 <- user()
#####SECTION BELOW NEEDED IF FEEDING IN BETAA VALUES######
DY <- user()
dim(EIR) <- c(na,nh)
EIR[,] <- rel_foi[j] * foi_age[i] * Iv*av0/omega
#output(EIR_annual) <- Iv*DY*av0/omega

#####Feed in EIR values; COMMENT OUT IF FEEDING IN BETAA VALS#########
# DY <- user()
# dim(EIR) <- c(na,nh)
# EIR_times[]<-user()
# EIR_seq[]<-user()
# EIR_valsd[]<-EIR_seq[i]/DY
# EIR_td<-interpolate(EIR_times, EIR_valsd, "constant")
#
# dim(EIR_times)<-user()
# dim(EIR_seq)<-user()
# dim(EIR_valsd)<-length(EIR_seq)
#
# EIR[,] <- EIR_td * rel_foi[j] * foi_age[i]
# output(EIR) <- EIR
# output(EIR_td) <- EIR_td
# output(Ivout) <- Iv
#
# output(omega) <- omega
#

# ##------------------------------------------------------------------------------
# #####################
# ## MOSQUITO STATES ##
# #####################
# ##------------------------------------------------------------------------------

# See supplementary materials S1 from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6
# fitted entomological parameters:
mv0 <- user() # initial mosquito density
mu0 <- user() # baseline mosquito death rate
tau1 <- user() # duration of host-seeking behaviour
tau2 <- user() # duration of resting behaviour
p10 <- user() # prob of surviving 1 feeding cycle
p2 <- user() #prob of surviving one resting cycle
fv <- 1/( tau1 + tau2 ) # mosquito feeding rate (zbar from intervention param.)
mu <- -fv*log(p10*p2) # mosquito death rate
omega <- user() #normalising constant for biting rates
# cA is the infectiousness to mosquitoes of humans in the asmyptomatic compartment broken down
# by age/het/int category, infectiousness depends on p_det which depends on detection immunity
cU <- user() # infectiousness U -> mosq
cD <- user() # infectiousness D -> mosq
cT <- user() # T -> mosq
gamma1 <- user() # fitted value of gamma1 characterises cA function
dim(cA) <- c(na,nh)
cA[,] <- cU + (cD-cU)*p_det[i,j]^gamma1
# Current hum->mos FOI depends on the number of individuals now producing gametocytes (12 day lag)
delayGam <- user() # Lag from parasites to infectious gametocytes
delayMos <- user() # Extrinsic incubation period.

# Number of mosquitoes that become infected at each time point
surv <- exp(-mu*delayMos)

# Force of infection from humans to mosquitoes
lag_ratesMos <- user()

FOIv_eq <- user()
initial(FOIv[]) <- FOIv_eq*delayGam/lag_ratesMos
dim(FOIv) <- lag_ratesMos

FOIvijk[1:na, 1:nh] <- av0 * (cT*T[i,j] + cD*D[i,j] + cA[i,j]*A[i,j] + cU*U[i,j]) * rel_foi[j] *foi_age[i]/omega
dim(FOIvijk) <- c(na,nh)
lag_FOIv <- sum(FOIvijk)

deriv(FOIv[1]) <- lag_FOIv - (lag_ratesMos/delayGam)*FOIv[1]
deriv(FOIv[2:lag_ratesMos]) <- (lag_ratesMos/delayGam)*FOIv[i-1] -
  (lag_ratesMos/delayGam)*FOIv[i]

ince <- FOIv[lag_ratesMos] * lag_ratesMos/delayGam * Sv

initial(ince_delay[]) <- FOIv_eq*init_Sv*mv0*delayMos/lag_ratesMos
dim(ince_delay) <- lag_ratesMos

deriv(ince_delay[1]) <- ince - (lag_ratesMos/delayMos)*ince_delay[1]
deriv(ince_delay[2:lag_ratesMos]) <- (lag_ratesMos/delayMos)*ince_delay[i-1] -
  (lag_ratesMos/delayMos)*ince_delay[i]

incv <- ince_delay[lag_ratesMos]*lag_ratesMos/delayMos *surv

# feed in betaa values from a random walk
####Comment in if reading in betaa values#####
initial(betaa_td) <- 0.65
deriv(betaa_td) <- 0
#update(betaa_td) <- betaa_td + 0.15

#####Comment next line out if reading in betaa values#####
# betaa_td <- mv * mu

# Sv - Susceptible mosquitoes
# Ev - latently infected (exposed) mosquitoes. Number of compartments used to simulate delay in becoming infectious
# Iv - Infectious mosquitoes

# initial state values:
init_Sv <- user()
init_Ev <- user()
init_Iv <- user()
initial(Sv) <- init_Sv * mv0
initial(Ev) <- init_Ev * mv0
initial(Iv) <- init_Iv * mv0

deriv(Sv) <- -ince - mu*Sv + betaa_td
deriv(Ev) <- ince - incv - mu*Ev
deriv(Iv) <- incv - mu*Iv

# Total mosquito population
mv = Sv+Ev+Iv

##------------------------------------------------------------------------------
###################
## MODEL OUTPUTS ##
###################
##------------------------------------------------------------------------------

# Outputs for each compartment across the sum across all ages, biting heterogeneities and intervention categories
output(Sout) <- sum(S[,])
output(Tout) <- sum(T[,])
output(Dout) <- sum(D[,])
output(Aout) <- sum(A[,])
output(Uout) <- sum(U[,])
output(Pout) <- sum(P[,])

# Outputs for each pregnancy compartment across the sum across all ages, biting heterogeneities and intervention categories
output(S_preg_out) <- sum(S_preg[,])
output(T_preg_out) <- sum(T_preg[,])
output(D_preg_out) <- sum(D_preg[,])
output(A_preg_out) <- sum(A_preg[,])
output(U_preg_out) <- sum(U_preg[,])
output(P_preg_out) <- sum(P_preg[,])
output(Y_preg_out) <- sum(Y_preg[,])

# Outputs for clinical incidence and prevalence on a given day
# population densities for each age category
den[] <- user()
dim(den) <- na
# index of the age vector above 59 months
age59 <- user(integer=TRUE)
# index of the age vector above 5 years
age05 <- user(integer=TRUE)

dim(prev0to59) <- c(age59,nh)
prev0to59[1:age59,] <- T[i,j] + D[i,j]  + A[i,j]*p_det[i,j]
output(prev) <- sum(prev0to59[,])/sum(den[1:age59])

# slide positivity in 0 -5 year age bracket
dim(clin_inc0to5) <- c(age05,nh)
clin_inc0to5[1:age05,] <- clin_inc[i,j]
output(inc05) <- sum(clin_inc0to5)/sum(den[1:age05])
output(inc) <- sum(clin_inc[,])

## Pregnancy Prevalence Outputs##
alphaA <- user()
alphaU <- user()
dim(prev_cba) <- c(preg_range,nh)
dim(prev_rdt) <- c(preg_range,nh)

prev_cba[1:preg_range,] <- D_preg[i,j] + A_preg[i,j]*p_det[as.integer((i+agestart)-1),j]^alphaA +
  U_preg[i,j]*p_det[as.integer((i+agestart)-1),j]^alphaU
prev_pcr <- sum(prev_cba[,])/sum(den[agestart:ageend])
prev_rdt[1:preg_range,] <- D_preg[i,j]  + A_preg[i,j]*p_det[as.integer((i+agestart)-1),j]

#output(prev_pcr) <- prev_pcr
output(prev_rdt_primi) <- sum(prev_rdt[,])/sum(den[agestart:ageend])
