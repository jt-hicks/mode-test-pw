pmcmc_history <- history_byage3
grav_code <- 2
anc_index_min <- 8
anc_index_max <- 27
burnin <- 50
data <- ANC_cMIS_for_pmcmc
re_cat_anc <- function(pmcmc_history,anc_index_min=8,
                       anc_index_max=27,
                       data=ANC_cMIS_for_pmcmc,
                       burnin=50,
                       grav_code=c(1,2,3)){
  
  change_prev_by_OR<-function(prev,OR){
    get_odds<-prev/(1-prev)*OR
    return(get_odds/(1+get_odds))
  }
  
  ##Get correct OR by gravidity
  or_list <- c(1.1422,0.8599,0.9137)
  or <- or_list[grav_code]
  
  ##Calculate new totals from observed data
  obs_tested <- as.matrix(data[,grep(paste0('tested_g',grav_code),names(data))])
  tested_sums <- rowSums(obs_tested[,1:20],na.rm = TRUE)
  
  ##Calculate new positives
  obs_positives <- as.matrix(data[,grep(paste0('positive_g',grav_code),names(data))])
  positive_sums <- rowSums(obs_positives[,1:20],na.rm = TRUE)
  
  obs_prev_grav <- positive_sums/tested_sums
  
  obs_cis <- addCIs(df=data.frame(positive_sums/tested_sums),Ys=positive_sums,Ns=tested_sums)%>%
    mutate(grav=grav_code)
                     
  hist_trunc <- pmcmc_history[anc_index_min:anc_index_max,(burnin+1):1000,-1]
  
  hist_grav <- change_prev_by_OR(hist_trunc,or)
  i<-1
  hist_newprev <- sapply(1:(1000-burnin), function(i) {
    exp <- data.frame(matrix(mapply(`*`,t(hist_grav[,i,]),obs_tested),nrow = 83))
    sum_exp <- rowSums(exp,na.rm=TRUE)
    sum_exp/tested_sums
  })
  out <- as.data.frame(hist_newprev) %>%
    mutate(t=1:83)%>%
    melt(id='t')%>%
    mutate(time=t,
           grav=grav_code)
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
