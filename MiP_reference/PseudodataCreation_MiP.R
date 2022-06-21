# Create a simulated data set to test binomial likelihood comparison function.
# load the model package (or have it built and reloaded as described above)
library(ICDMM)
library(odin)
library(dde)
library(ggplot2)

#Call functions housed in other files for simplicity
#Function that calculates the initial equilibrium state
source("MiP_model/equilibrium-init-create.R")
#Function that defines model parameters
source("MiP_model/model_parameters.R")
#Function that actually runs the odin code and generates model output
source("MiP_model/run_model_function.R")
#Odin object that specifies the model in the odin language
source("MiP_model/MiP_odin_model_nodelay.R")

# create a vector of age categories
init_age <- c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,7.5,10,15,20,30,40,50,60,70,80)

# provide a value of the annual EIR for this model run
init_EIR <- 10
init_betaa <- 0.65

# provide the length of time (in days) that you want to run the model for
time_period <- 365*5

# provide a value for the proportion of cases that are treated (referred to as ft in the paper)
prop_treated <- 0.4

# pregnancy category initial values
rA_preg <- 0.00512821 #same as main compartments
#rA_preg <- 0.1
rU_preg <- 0.00906627 #same as main compartments


# Create random walk beta sequence to pass to model
init_beta <- -log(0.9)
genRandWalk <- function(x,vol,randWalk) {
  if (x == 0)    return (randWalk)
  else return(genRandWalk(x-1,vol,c(randWalk,randWalk[length(randWalk)]*exp(rnorm(1)*vol))))
}
beta_times<-seq(0,time_period,by=30)
beta_volatility<-0.5
beta_vals=genRandWalk(length(beta_times)-1,beta_volatility,init_beta)
keep <- beta_vals
save(keep,file='MiP_model/beta_vals_mip.rData')

# Run model

model_run <- run_model(age=init_age, init_EIR=init_EIR, init_ft=prop_treated,
                       time=time_period,
                       betaa_vals=beta_vals,
                       betaa_times=beta_times,
                       rA_preg=rA_preg, rU_preg=rU_preg,
                       comm_age_min = 0,
                       comm_age_max = 59/12,
                       anc_age_min = 15,
                       anc_age_max = 20,
                       lag_rates = 10,
                       lag_ratesMos = 10
)
output <- as.data.frame(model_run[c('t','inc05','prev')])
plot(model_run$t, model_run$inc05)
plot(model_run$t, model_run$prev)

#Generate pseudo data from run instance
#Select proportions at whole days
days <- output[,c('t','prev')]
days$tested <- round(rnorm(n = length(days$t), mean = 10000, sd = 100))
days$positive <- rbinom(n = length(days$t), size = days$tested, p = days$prev)
test_data <- days[seq(1,nrow(days),30),c('t','tested','positive')]
write.csv(test_data, 'MiP_model/casedata_monthly.csv', row.names = FALSE)

ggplot(test_data, aes(positive)) +
  geom_histogram()
matplot(test_data$t, test_data$positive/test_data$tested, type = "l", lty = 1, #col = cols_prev,
        xlab = "time", ylab = "number positive")
plot(beta_vals)
