eir_gen <- function(init_EIR=100,EIR_volatility=0.5,max_EIR=1000){
################## generate the data ######################
# generate random walk of EIR (recursive fn)
genRandWalk <- function(x,vol,randWalk) {
  if (x == 0)    return (randWalk)
  else return(genRandWalk(x-1,vol,c(randWalk,min(exp(log(randWalk[length(randWalk)])+rnorm(1)*vol),max_EIR))))
}
EIR_times<-seq(0,1800,by=30)

### just a random walk on logscale
EIR_vals=genRandWalk(length(EIR_times)-1,EIR_volatility,init_EIR)

return(EIR_vals)
}
