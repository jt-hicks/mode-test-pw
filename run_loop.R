#Run loop function
setwd('Q:/mode-test-pw')
drat:::add("mrc-ide") # install.packages("drat") if this errors, then try again
install.packages("didehpc")
library(didehpc)
library(pkgdepends)

root <- "contexts"
sources <- "loop_function.R"
config <- didehpc::didehpc_config(cores = 4, parallel = TRUE)


ctx <- context::context_save("contexts", sources = sources,
                             #packages = c("parallelly"),
                             package_sources = conan::conan_sources(c("mrc-ide/odin.dust",'mrc-ide/mcstate')))

obj <- didehpc::queue_didehpc(ctx,config = config)

obj$cluster_load(TRUE)
obj$config

t <- obj$enqueue(pf_loop(data='casedata_monthly.csv',n_loop = 500))
start.time <- Sys.time()
pf_loop(data='casedata_monthly.csv',n_loop = 10)
print(Sys.time()-start.time)

t$status()
t$result()
t$wait(500)

t$times()
t$log()

df <- as.data.frame(t$result())
colnames(df) <- c(10,50,100,200,500)
