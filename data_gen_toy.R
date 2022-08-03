# Create a simulated data set to test binomial likelihood comparison function.
# Read in odin model
data_gen_toy <- function(volatility,init_beta=-log(0.9),model_file='MiP-given/RM_RW_odin_model.R'){
# Create random walk beta sequence to pass to model
  source(model_file)
  genRandWalk <- function(x,vol,randWalk) {
    if (x == 0)    return (randWalk)
    else return(genRandWalk(x-1,vol,c(randWalk,randWalk[length(randWalk)]*exp(rnorm(1)*vol))))
  }
  beta_times<-seq(0,1800,by=30)
  beta_volatility<-volatility
  beta_vals <- genRandWalk(length(beta_times)-1,beta_volatility,init_beta)
  keep <- beta_vals

  # Run model
  init_Ih <- 0.8
  init_Sv <- 100
  init_Iv <- 1
  nrates <- 15
  
  model_run <- odin_model(init_Ih=init_Ih,
                          init_Sv=init_Sv,
                          init_Iv=init_Iv,
                          beta_times=beta_times,
                          beta_vals=beta_vals,
                          nrates=nrates)
  t <- seq(30, 1830, by = 30)
  output <- as.data.frame(model_run$run(t))
  
  #Generate pseudo data from run instance
  #Select proportions at whole days
  days <- output[,c('t','Ih')]
  days$tested <- round(rnorm(n = length(days$t), mean = 10000, sd = 100))
  days$positive <- rbinom(n = length(days$t), size = days$tested, p = days$Ih)
  test_data <- days[,c('t','tested','positive')]

  
  matplot(days$t, days$positive/days$tested, type = "l", lty = 1, #col = cols_prev,
          xlab = "time", ylab = "number positive")
  return(test_data)
}