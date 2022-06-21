#Modifications from original model:
# Removed intervention components (formerly num_int, indexed as 'k')
# Seasonality removed
# Allow an interpolated series of values for mosquito birth rate
#     (betaa)
# Removed larval compartments
# Added pregnancy components, including:
#     Pregnancy compartments
#     Calculation of risk of previous exposure to MiP
#     PCR and RDT prevalence for pregnancy compartments

#------------------------------------------------
#' Equilibrium initialisation list creation
#'
#' \code{equilibrium_init_create} creates an equilibrium initialisation state to be
#' used within later model runs
#' @param age_vector Vector of age brackets.
#' @param het_brackets Integer number of biting heteogenity compartments.
#' @param country String for country of interest. If NULL the seasonal parameters
#' will attempt to be loaded using just the admin unit, however if there is ambiguity
#' in the admin unit an error will be thrown. If both NULL then no seasonality is
#' assumed. Default = NULL.
#' @param admin_unit String for admin unit with country for loading seasonal
#' parameters. If country is NULL, the admin unit will attempt to be located,however
#' if there is ambiguity in the admin unit an error will be thrown. If both country
#' and admin_unit are NULL then no seasonality is assumed. Default = NULL.
#' @param ft Numeric for the frequency of people seeking treatment.
#' @param EIR Numeric for desired annual EIR.
#' @param model_param_list List of epidemiological parameters created by
#'
#' @importFrom stringi stri_trans_general
#' @importFrom statmod gauss.quad.prob
#'
#'
#' @export

equilibrium_init_create <- function(age_vector, het_brackets,
                                    country = NULL, admin_unit = NULL, ft,
                                    EIR, model_param_list, rA_preg
                                    )
{
  # mpl is shorter :)
  mpl <- model_param_list

  ## Check Parameters
  if(!is.numeric(age_vector)) stop("age_vector provided is not numeric")
  if(!is.numeric(het_brackets)) stop("het_brackets provided is not numeric")
  if(!(is.null(country) | is.character(country))) stop("country specified is not character string")
  if(!(is.null(admin_unit) | is.character(admin_unit))) stop("admin_unit specified is not character string")
  if(!is.numeric(ft)) stop("ft provided is not numeric")
  if(!is.numeric(EIR)) stop("EIR provided is not numeric")

  ## Handle parameters
  # database for admin units is all in Latin-ASCII for CRAN reasons so must
  # encode parameters accordingly
  if(!is.null(country)) country <- stringi::stri_trans_general(country,"Latin-ASCII")
  if(!is.null(admin_unit)) admin_unit <- stringi::stri_trans_general(admin_unit, "Latin-ASCII")

  ## population demographics
  age <- age_vector * mpl$DY
  na <- as.integer(length(age))  # number of age groups
  nh <- as.integer(het_brackets)  # number of heterogeneity groups
  h <- statmod::gauss.quad.prob(nh, dist = "normal")
  age0 <- 2
  age1 <- 10
  num_int <- mpl$num_int

  age_rate <- age_width <- age_mid_point <- den <- c()
  for (i in 1:(na-1))
  {
    age_width[i] <- age[i+1] - age[i]
    age_rate[i] <- 1/(age[i + 1] - age[i])  # vector of rates at which people leave each age group (1/age group width)
    age_mid_point[i] <- 0.5 * (age[i] + age[i + 1])  # set age group vector to the midpoint of the group

  }
  age_rate[na] = 0


  den <- 1/(1 + age_rate[1]/mpl$eta)
  for (i in 1:(na-1))
  {
    den[i+1] <- age_rate[i] * den[i]/(age_rate[i+1] + mpl$eta)  # proportion in each age_vector group
  }

  age59 <- which(age_vector * 12 > 59)[1] - 1  # index of age vector before age is >59 months
  age05 <- which(age_vector > 5)[1] - 1  # index of age vector before age is 5 years

  # index of the requested age vector min
  agestart <- as.integer(which(age_vector > 15)[1] - 1)  # index of age vector before age is >59 months
  # index of the requested age vector just prior to the maximum value for the prevalence calculation
  ageend <- as.integer(which(age_vector > 49)[1] - 1)

  ## force of infection
  foi_age <- c()
  for (i in 1:na)
  {
    foi_age[i] <- 1 - (mpl$rho * exp(-age[i]/mpl$a0))  #force of infection for each age group
  }
  fden <- foi_age * den
  omega <- sum(fden)  #normalising constant

  ## heterogeneity
  het_x <- h$nodes
  het_wt <- h$weights
  den_het <- outer(den, het_wt)
  rel_foi <- exp(-mpl$sigma2/2 + sqrt(mpl$sigma2) * het_x)/sum(het_wt * exp(-mpl$sigma2/2 + sqrt(mpl$sigma2) * het_x))

  ## EIR
  EIRY_eq <- EIR  # initial annual EIR

  EIRd_eq <- EIRY_eq/mpl$DY
  EIR_eq <- outer(foi_age, rel_foi) * EIRd_eq

  ## Immunity and FOI
  x_I <- den[1]/mpl$eta
  for (i in 2:na)
  {
    x_I[i] <- den[i]/(den[i - 1] * age_rate[i - 1])  #temporary variables
  }
  fd <- 1 - (1 - mpl$fD0)/(1 + (age/mpl$aD)^mpl$gammaD)

  # maternal immunity begins at a level proportional to the clinical
  # immunity of a 20 year old, this code finds that level
  age20i <- rep(0, na)
  for (i in 2:na)
  {
    age20i[i] <- ifelse(age[i] >= (20 * mpl$DY) & age[i - 1] < (20 * mpl$DY),
                        i, age20i[i - 1])
  }
  age20u <- as.integer(age20i[na])
  age20l <- as.integer(age20u - 1)
  age_20_factor <- (20 * mpl$DY - age[age20l] - 0.5 * age_width[age20l]) *
    2/(age_width[age20l] + age_width[age20u])

  # finding initial values for all immunity states
  IB_eq <- matrix(0, na, nh)
  FOI_eq <- matrix(0, na, nh)
  ID_eq <- matrix(0, na, nh)
  ICA_eq <- matrix(0, na, nh)
  ICM_init_eq <- vector(length = nh, mode = "numeric")
  ICM_eq <- matrix(0, na, nh)
  cA_eq <- matrix(0, na, nh)
  FOIvij_eq <- matrix(0, na, nh)
  p_det_eq <- matrix(0, na, nh)
  for (j in 1:nh)
  {
    for (i in 1:na)
    {
      IB_eq[i, j] <- (ifelse(i == 1, 0, IB_eq[i - 1, j]) +
                        EIR_eq[i,j]/(EIR_eq[i, j] * mpl$uB + 1) * x_I[i])/(1 + x_I[i]/mpl$dB)
      FOI_eq[i, j] <- EIR_eq[i, j] * ifelse(IB_eq[i, j] == 0, mpl$b0,
                                            mpl$b0 * ((1 - mpl$b1)/(1 + (IB_eq[i, j]/mpl$IB0)^mpl$kB) + mpl$b1))
      ID_eq[i, j] <- (ifelse(i == 1, 0, ID_eq[i - 1, j]) +
                        FOI_eq[i, j]/(FOI_eq[i, j] * mpl$uD + 1) * x_I[i])/(1 + x_I[i]/mpl$dID)
      ICA_eq[i, j] <- (ifelse(i == 1, 0, ICA_eq[i - 1, j]) +
                         FOI_eq[i,j]/(FOI_eq[i, j] * mpl$uCA + 1) * x_I[i])/(1 + x_I[i]/mpl$dCA)
      p_det_eq[i, j] <- mpl$d1 + (1 - mpl$d1)/(1 + fd[i] * (ID_eq[i, j]/mpl$ID0)^mpl$kD)
      cA_eq[i, j] <- mpl$cU + (mpl$cD - mpl$cU) * p_det_eq[i, j]^mpl$gamma1
    }
  }
  # needs to be calculated after because it references ICA
  for (j in 1:nh)
  {
    for (i in 1:na)
    {
      ICM_init_eq[j] <- mpl$PM * (ICA_eq[age20l, j] + age_20_factor *
                                    (ICA_eq[age20u, j] - ICA_eq[age20l, j]))
      ICM_eq[i, j] <- ifelse(i == 1,
                             ICM_init_eq[j], ICM_eq[i - 1,j])/(1 + x_I[i]/mpl$dCM)
    }
  }

  IC_eq <- ICM_eq + ICA_eq
  phi_eq <- mpl$phi0 * ((1 - mpl$phi1)/(1 + (IC_eq/mpl$IC0)^mpl$kC) + mpl$phi1)


  # human states
  gamma <- mpl$eta + c(age_rate[1:(na - 1)], 0)
  delta <- c(mpl$eta, age_rate[1:(na - 1)])

  betaT <- matrix(rep(mpl$rT + gamma, rep(nh, na)), ncol = nh, byrow = TRUE)
  betaD <- matrix(rep(mpl$rD + gamma, rep(nh, na)), ncol = nh, byrow = TRUE)
  betaP <- matrix(rep(mpl$rP + gamma, rep(nh, na)), ncol = nh, byrow = TRUE)

  aT <- FOI_eq * phi_eq * ft/betaT
  aD <- FOI_eq * phi_eq * (1 - ft)/betaD
  aP <- mpl$rT * aT/betaP

  Z_eq <- array(dim = c(na, nh, 4))
  Z_eq[1, , 1] <- den_het[1, ]/(1 + aT[1, ] + aD[1, ] + aP[1, ])
  Z_eq[1, , 2] <- aT[1, ] * Z_eq[1, , 1]
  Z_eq[1, , 3] <- aD[1, ] * Z_eq[1, , 1]
  Z_eq[1, , 4] <- aP[1, ] * Z_eq[1, , 1]

  for (j in 1:nh)
  {
    for (i in 2:na)
      {
      Z_eq[i, j, 1] <- (den_het[i, j] - delta[i] * (Z_eq[i - 1, j, 2]/betaT[i, j] +
                                                      Z_eq[i - 1, j, 3]/betaD[i, j] +
                                                      (mpl$rT *  Z_eq[i - 1, j, 2]/betaT[i, j]
                                                       + Z_eq[i - 1, j, 4])/betaP[i, j]))/(1 + aT[i, j] + aD[i, j] + aP[i, j])
      Z_eq[i, j, 2] <- aT[i, j] * Z_eq[i, j, 1] + delta[i] * Z_eq[i -
                                                                    1, j, 2]/betaT[i, j]
      Z_eq[i, j, 3] <- aD[i, j] * Z_eq[i, j, 1] + delta[i] * Z_eq[i -
                                                                    1, j, 3]/betaD[i, j]
      Z_eq[i, j, 4] <- aP[i, j] * Z_eq[i, j, 1] + delta[i] * (mpl$rT *
                                                                Z_eq[i - 1, j, 2]/betaT[i, j] + Z_eq[i - 1, j, 4])/betaP[i,j]

    }
  }

  Y_eq <- matrix(Z_eq[, , 1], nrow = na, ncol=nh)
  T_eq <- matrix(Z_eq[, , 2], nrow = na, ncol=nh)
  D_eq <- matrix(Z_eq[, , 3], nrow = na, ncol=nh)
  P_eq <- matrix(Z_eq[, , 4], nrow = na, ncol=nh)

  betaS <- apply(FOI_eq, MARGIN = 2, FUN = function(x, y)
  {
    x + y
  }, y = gamma)
  betaA <- apply(FOI_eq * phi_eq
                 + mpl$rA, MARGIN = 2, FUN = function(x, y)
  {
    x + y
  }, y = gamma)
  betaU <- apply(FOI_eq + mpl$rU, MARGIN = 2, FUN = function(x, y)
  {
    x + y
  }, y = gamma)

  A_eq <- matrix(ncol = nh, nrow = na)
  U_eq <- matrix(ncol = nh, nrow = na)
  S_eq <- matrix(ncol = nh, nrow = na)

  for (i in 1:na)
  {
    for (j in 1:nh)
    {
      A_eq[i, j] <- (delta[i] * ifelse(i == 1, 0, A_eq[i - 1, j]) +
                       FOI_eq[i, j] * (1 - phi_eq[i, j]) * Y_eq[i, j] +
                       mpl$rD * D_eq[i,j])/(betaA[i, j] + FOI_eq[i, j] * (1 - phi_eq[i, j]))
      U_eq[i, j] <- (mpl$rA * A_eq[i, j] + delta[i] * ifelse(i == 1,
                                                             0, U_eq[i - 1, j]))/betaU[i, j]
      S_eq[i, j] <- Y_eq[i, j] - A_eq[i, j] - U_eq[i, j]
      FOIvij_eq[i, j] <- foi_age[i] * mpl$av0 * (mpl$cT * T_eq[i, j] + mpl$cD *
                                                   D_eq[i, j] + cA_eq[i, j] * A_eq[i, j] + mpl$cU * U_eq[i, j]) * rel_foi[j]/omega
    }
  }

  # pregnancy categories
  # human pregnancy states
  # Specify number of age categories that are child-bearing
  preg_range <- (ageend - agestart) + 1

  betaT_preg <- matrix(rep(rep(mpl$rT,na), rep(nh, na)), ncol = nh, byrow = TRUE)
  betaD_preg <- matrix(rep(rep(mpl$rD,na), rep(nh, na)), ncol = nh, byrow = TRUE)
  betaP_preg <- matrix(rep(rep(mpl$rP,na), rep(nh, na)), ncol = nh, byrow = TRUE)

  aT_preg <- FOI_eq * phi_eq * ft/betaT_preg
  aD_preg <- FOI_eq * phi_eq * (1 - ft)/betaD_preg
  aP_preg <- mpl$rT * aT_preg/betaP_preg

  Z_eq_preg <- array(dim = c(preg_range, nh, 4))

  for (j in 1:nh)
  {
    for (i in 1:preg_range)
    {
      Z_eq_preg[i, j, 1] <- den_het[as.integer((i+agestart)-1),j]/
                            (1 + aT_preg[as.integer((i+agestart)-1),j] + aD_preg[as.integer((i+agestart)-1),j] + aP_preg[as.integer((i+agestart)-1),j])
      Z_eq_preg[i, j, 2] <- aT_preg[as.integer((i+agestart)-1), j] * Z_eq_preg[i, j, 1]
      Z_eq_preg[i, j, 3] <- aD_preg[as.integer((i+agestart)-1), j] * Z_eq_preg[i, j, 1]
      Z_eq_preg[i, j, 4] <- mpl$rT * Z_eq_preg[i, j, 2]/betaP_preg[as.integer((i+agestart)-1), j]

    }
  }

  Y_eq_preg <- matrix(Z_eq_preg[, , 1], nrow = preg_range, ncol=nh)
  T_eq_preg <- matrix(Z_eq_preg[, , 2], nrow = preg_range, ncol=nh)
  D_eq_preg <- matrix(Z_eq_preg[, , 3], nrow = preg_range, ncol=nh)
  P_eq_preg <- matrix(Z_eq_preg[, , 4], nrow = preg_range, ncol=nh)

  betaS_preg <- FOI_eq
  betaA_preg <- FOI_eq * phi_eq + rA_preg
  betaU_preg <- FOI_eq + rU_preg

  A_eq_preg <- matrix(ncol = nh, nrow = preg_range)
  U_eq_preg <- matrix(ncol = nh, nrow = preg_range)
  S_eq_preg <- matrix(ncol = nh, nrow = preg_range)
  prev_cba_eq <- matrix(ncol = nh, nrow = preg_range)

  for (i in 1:preg_range)
  {
    for (j in 1:nh)
    {
      A_eq_preg[i, j] <- (FOI_eq[as.integer((i+agestart)-1), j] * (1 - phi_eq[as.integer((i+agestart)-1), j]) * Y_eq_preg[i, j] +
                       mpl$rD * D_eq_preg[i,j])/(betaA_preg[as.integer((i+agestart)-1), j] + FOI_eq[as.integer((i+agestart)-1), j] * (1 - phi_eq[as.integer((i+agestart)-1), j]))
      U_eq_preg[i, j] <- (rA_preg * A_eq_preg[i, j])/betaU_preg[as.integer((i+agestart)-1), j]
      S_eq_preg[i, j] <- Y_eq_preg[i, j] - A_eq_preg[i, j] - U_eq_preg[i, j]
      prev_cba_eq[i, j] <- D_eq_preg[i,j] + A_eq_preg[i,j]*p_det_eq[as.integer((i+agestart)-1),j]^mpl$alphaA +
        U_eq_preg[i,j]*p_det_eq[as.integer((i+agestart)-1),j]^mpl$alphaU
    }
  }
  prev_pcr_eq <- sum(prev_cba_eq[,])/sum(den[agestart:ageend])

  #gravidity
  nrates <- mpl$nrates
  sample_rates <- mpl$sample_rates
  sample_transition_rates <- mpl$sample_transition_rates
  wane_rates <- mpl$wane_rates
  pregs_eq <- vector(length = nrates)
  for (i in 1:nrates){
    pregs_eq[i] <- (prev_pcr_eq*sample_rates[i] +
                      ifelse(i==1, 0, sample_transition_rates[i-1]*pregs_eq[i-1]))/
      (sample_transition_rates[i] + wane_rates)
  }
  previous_pregs_eq <- pregs_eq[nrates]

  # mosquito states
  FOIv_eq <- sum(FOIvij_eq)
  Iv_eq <- FOIv_eq * mpl$Surv0/(FOIv_eq + mpl$mu0)
  Sv_eq <- mpl$mu0 * Iv_eq/(FOIv_eq * mpl$Surv0)
  Ev_eq <- 1 - Sv_eq - Iv_eq
  # mosquito density needed to give this EIR
  mv0 <- omega * EIRd_eq/(Iv_eq * mpl$av0)
  beta_eq <- mpl$mu0 * mv0

  # larval states
  # K0 <- 2 * mv0 * mpl$dLL * mpl$mu0 * (1 + mpl$dPL * mpl$muPL) * mpl$gammaL * (mpl$lambda + 1)/(mpl$lambda/(mpl$muLL *
  #                                                                                                             mpl$dEL) - 1/(mpl$muLL * mpl$dLL) - 1)
  # PL_eq <- 2 * mpl$dPL * mpl$mu0 * mv0
  # LL_eq <- mpl$dLL * (mpl$muPL + 1/mpl$dPL) * PL_eq
  # EL_eq <- (LL_eq/mpl$dLL + mpl$muLL* LL_eq * (1 + mpl$gammaL * LL_eq/K0))/(1/mpl$dEL - mpl$muLL * mpl$gammaL * LL_eq/K0)


  # better het bounds for equilbirum initialisation in individual model
  zetas <- rlnorm(n = 1e5,meanlog = -mpl$sigma2/2, sdlog = sqrt(mpl$sigma2))
  while(sum(zetas>100)>0){
    zetas[zetas>100] <- rlnorm(n = sum(zetas>100),meanlog = -mpl$sigma2/2, sdlog = sqrt(mpl$sigma2))
  }

  wt_cuts <- round(cumsum(het_wt)*1e5)
  zeros <- which(wt_cuts==0)
  wt_cuts[zeros] <- 1:length(zeros)
  larges <- which(wt_cuts==1e5)
  wt_cuts[larges] <- (1e5 - (length(larges)-1)):1e5
  wt_cuts <- c(0,wt_cuts)
  het_bounds <- sort(zetas)[wt_cuts]
  het_bounds[length(het_bounds)] <- (mpl$max_age/365)+1

  #prevalence for fitting
  #Prevalence in general population for user-specified age range
  age_min <- as.integer(which(age_vector > mpl$comm_age_min)[1] - 1)
  age_max <- as.integer(which(age_vector > mpl$comm_age_max)[1] - 1)
  range_comm_age <- (age_max - age_min) + 1

  prev_comm_rdt <- matrix(nrow = range_comm_age, ncol = nh)
  prev_comm_pcr <- matrix(nrow = range_comm_age, ncol = nh)

  for (i in 1:range_comm_age)
  {
    for (j in 1:nh)
    {
      prev_comm_rdt[i,j] <- D_eq[as.integer((i+age_min)-1),j] +
        A_eq[as.integer((i+age_min)-1),j]*p_det_eq[as.integer((i+age_min)-1),j]
      prev_comm_pcr[i,j] <- D_eq[as.integer((i+age_min)-1),j] +
        A_eq[as.integer((i+age_min)-1),j]*p_det_eq[as.integer((i+age_min)-1),j]^mpl$alphaA +
        U_eq[as.integer((i+age_min)-1),j]*p_det_eq[as.integer((i+age_min)-1),j]^mpl$alphaU
    }
  }
  prev_comm_rdt_eq <- sum(prev_comm_rdt[,])/sum(den[age_min:age_max])
  prev_comm_pcr_eq <- sum(prev_comm_pcr[,])/sum(den[age_min:age_max])
  #weighted average of p_det over the age and het categories of interest
  p_det_comm <- sum(p_det_eq[age_min:age_max,]*den_het[age_min:age_max,])/sum(den_het[age_min:age_max,])
  S_inf <- sum(S_eq[age_min:age_max,])
  A_inf <- sum(A_eq[age_min:age_max,])
  U_inf <- sum(U_eq[age_min:age_max,])
  D_inf <- sum(D_eq[age_min:age_max,])
  P_inf <- sum(P_eq[age_min:age_max,])
  T_inf <- sum(T_eq[age_min:age_max,])

  #Prevalence in pregnanty women aged 15-20 years
  anc_min <- as.integer(which(age_vector > mpl$anc_age_min)[1] - 1)  # index of age vector before age is 2 years
  anc_max <- as.integer(which(age_vector > mpl$anc_age_max)[1] - 1)  # index of age vector before age is 10 years
  range_anc_age <- (anc_max - anc_min) + 1

  prev_anc_rdt <- matrix(ncol = nh, nrow = range_anc_age)
  prev_anc_pcr <- matrix(ncol = nh, nrow = range_anc_age)
  for (i in 1:range_anc_age)
  {
    for (j in 1:nh)
    {
      prev_anc_rdt[i,j] <- D_eq_preg[i,j] + A_eq_preg[i,j]*p_det_eq[as.integer((i+anc_min)-1),j]
      prev_anc_pcr[i,j] <- D_eq_preg[i,j] + A_eq_preg[i,j]*p_det_eq[as.integer((i+anc_min)-1),j]^mpl$alphaA +
        U_eq_preg[i,j]*p_det_eq[as.integer((i+anc_min)-1),j]^mpl$alphaU
    }
  }
  prev_anc_rdt_eq <- sum(prev_anc_rdt[,])/sum(den[anc_min:anc_max])
  prev_anc_pcr_eq <- sum(prev_anc_pcr[,])/sum(den[anc_min:anc_max])
  #weighted average of p_det over the age and het categories of interest
  p_det_anc <- sum(p_det_eq[anc_min:anc_max,]*den_het[anc_min:anc_max,])/sum(den_het[anc_min:anc_max,])
  S_eq_anc <- sum(S_eq_preg[1:range_anc_age,])
  A_eq_anc <- sum(A_eq_preg[1:range_anc_age,])
  U_eq_anc <- sum(U_eq_preg[1:range_anc_age,])
  D_eq_anc <- sum(D_eq_preg[1:range_anc_age,])
  P_eq_anc <- sum(P_eq_preg[1:range_anc_age,])
  T_eq_anc <- sum(T_eq_preg[1:range_anc_age,])
  #print(beta_eq)
  # cat('in eq function\nprev_anc_pcr: ',prev_anc_pcr,'\nden[agestart]: ',den[agestart],
  #     '\nsum prev_anc_pcr: ',sum(prev_anc_pcr[,]),'\nagestart: ',agestart,
  #     '\nrA_preg: ',rA_preg,'\nEIR: ',EIRY_eq,'\n')
  ## collate init
  res <- list(init_S = S_eq, init_T = T_eq, init_D = D_eq, init_A = A_eq, init_U = U_eq,
              init_P = P_eq, init_Y = Y_eq, init_IB = IB_eq, init_ID = ID_eq, init_ICA = ICA_eq,
              init_ICM = ICM_eq, ICM_init_eq = ICM_init_eq, init_Iv = Iv_eq, init_Sv = Sv_eq,
              init_Ev = Ev_eq, 
              age_width = age_width, age_rate = age_rate, het_wt = het_wt, het_x = het_x,
              den_het = den_het,
              omega = omega, foi_age = foi_age, rel_foi = rel_foi,
              mv0 = mv0, na = na, nh = nh, ni = num_int, x_I = x_I,
              FOI_eq = FOI_eq, EIR_eq = EIR_eq, cA_eq = cA_eq,
              den = den, age59 = age59, age05 = age05,
              agestart = agestart,
              ageend = ageend,
              age = age_vector*mpl$DY, ft = ft, FOIv_eq = FOIv_eq,
              betaS = betaS, betaA = betaA, betaU = betaU, FOIvij_eq=FOIvij_eq,
              age_mid_point = age_mid_point, het_bounds = het_bounds, pi = pi,
              age20l = age20l, age20u = age20u, age_20_factor = age_20_factor,
              init_pregs = pregs_eq,init_S_preg = S_eq_preg,
              init_T_preg = T_eq_preg, init_D_preg = D_eq_preg,
              init_A_preg = A_eq_preg, init_U_preg = U_eq_preg, init_P_preg = P_eq_preg,
              prev_pcr_eq = prev_pcr_eq, previous_pregs_eq = previous_pregs_eq,
              EIRd_eq = EIRd_eq,
              p_det_eq = p_det_eq,
              prev_cba_eq = prev_cba_eq,
              prev_cba_eq_sum = sum(prev_cba_eq[,]),
              den_cba = sum(den[agestart:ageend]),
              prev_comm_rdt_eq = prev_comm_rdt_eq,
              prev_comm_pcr_eq = prev_comm_pcr_eq,
              prev_anc_rdt_eq = prev_anc_rdt_eq,
              prev_anc_pcr_eq = prev_anc_pcr_eq,
              p_det_comm =p_det_comm,
              p_det_anc = p_det_anc,
              S_inf = S_inf, A_inf = A_inf, U_inf = U_inf,
              D_inf = D_inf, P_inf = P_inf, T_inf = T_inf,
              S_eq_anc = S_eq_anc, A_eq_anc = A_eq_anc, U_eq_anc = U_eq_anc,
              D_eq_anc = D_eq_anc, P_eq_anc = P_eq_anc, T_eq_anc = T_eq_anc
              )

  res <- append(res,mpl)

  return(res)
}




