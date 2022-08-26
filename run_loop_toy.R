#Run loop function
setwd('Q:/mode-test-pw')
# drat:::add("mrc-ide") # install.packages("drat") if this errors, then try again
# install.packages("didehpc")
library(didehpc)
library(pkgdepends)

source('data_gen_toy.R')
data_raw <- data_gen_toy(volatility = 0.3)
root <- "contexts"
sources <- "loop_function_toy.R"
config <- didehpc::didehpc_config(cores = 4, parallel = TRUE)


ctx <- context::context_save("contexts", sources = sources,
                             #packages = c("parallelly"),
                             package_sources = conan::conan_sources(c("mrc-ide/odin.dust",'mrc-ide/mcstate')))

obj <- didehpc::queue_didehpc(ctx,config = config)

obj$cluster_load(TRUE)
obj$config

t <- obj$enqueue(pf_loop(data_raw,n_loop = 500,volatility=0.3))
start.time <- Sys.time()
pf_loop(data_raw,n_loop = 10,volatility=0.3,freq=1)
print(Sys.time()-start.time)

t$status()
t$result()
t$wait(500)

t$times()
t$log()

df <- as.data.frame(t$result())
colnames(df) <- c(10,50,100,200,500)

saveRDS(df,'C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PMCMC/Development/toymodel_lik_030822-2.rds')
lapply(df,var)
