#Run loop function
setwd('Q:/mode-test-main')

library(didehpc)
library(pkgdepends)

root <- "contexts"
sources <- "loop_function.R"

ctx <- context::context_save("contexts", sources = sources,
                             package_sources = conan::conan_sources(c("mrc-ide/odin.dust",'mrc-ide/mcstate')))
obj <- didehpc::queue_didehpc(ctx)

obj$cluster_load(TRUE)

t <- obj$enqueue(pf_loop(data='casedata_monthly.csv',n_loop = 10))
t$status()
t$result()
t$wait(100)

t$times()
t$log()

df <- as.data.frame(t$result())
colnames(df) <- c(10,50,100,200,500)
