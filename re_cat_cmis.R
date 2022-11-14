pmcmc_history <- history_byage3
grav_code <- 2
anc_index_min <- 5
anc_index_max <- 24
burnin <- 50
data <- ANC_cMIS_for_pmcmc
re_cat_cmis <- function(pmcmc_history,anc_index_min=5,
                       anc_index_max=24,
                       data=ANC_cMIS_for_pmcmc,
                       burnin=50){
  
  ##Calculate new totals from observed data
  obs_tested <- as.matrix(data[,grep('tested_cmis',names(data))])
  tested_sums <- rowSums(obs_tested[,anc_index_min:anc_index_max],na.rm = TRUE)
  
  ##Calculate new positives
  obs_positives <- as.matrix(data[,grep('positive_cmis',names(data))])
  positive_sums <- rowSums(obs_positives[,anc_index_min:anc_index_max],na.rm = TRUE)
  
  obs_prev_grav <- positive_sums/tested_sums
  
  obs_cis <- addCIs(df=data.frame(positive_sums/tested_sums),Ys=positive_sums,Ns=tested_sums)
                     
  hist_trunc <- pmcmc_history[(anc_index_min+3):(anc_index_max+3),(burnin+1):1000,-1]
  i<-1
  hist_newprev <- sapply(1:(1000-burnin), function(i) {
    exp <- data.frame(matrix(mapply(`*`,t(hist_trunc[,i,1:61]),obs_tested[1:61,anc_index_min:anc_index_max]),nrow = 61))
    sum_exp <- rowSums(exp,na.rm=TRUE)
    newprev <- sum_exp/tested_sums[1:61]
    newprev[62:83] <- rowMeans(data.frame(t(hist_trunc[,i,62:83])))
    return(newprev)
  })
  out <- as.data.frame(hist_newprev) %>%
    mutate(t=1:83)%>%
    melt(id='t')%>%
    mutate(time=t)
  return(list(model=out,observed=obs_cis))
}
hist_trunc[1,1,1]
hist_grav[1,1,1]
plot_particle_filter(hist_newprev,true_history = data.frame(tested=total_sums,positive=positive_sums),times =1:83)                     

prev <- positive_sums/total_sums
par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
matplot(1:83, hist_newprev, type = "l",
        xlab = "Time", ylab = "Prevalence",
        col = "#A6CEE3",
        lty = 1, ylim = c(0,1))
matpoints(1:83, obs_prev_grav, pch = 19,
          col = "#1F78B4")
