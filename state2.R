library(tidyverse)
library(glue)

set.seed(20200601)
i <- Sys.getenv("SLURM_ARRAY_TASK_ID")
print(paste0("Entering loop number:", i))
i=as.numeric(i)
##########hospital  https://www.cdc.gov/coronavirus/2019-ncov/covid-data/covidview/index.html
US=read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv")
h_r=98.4*326687501/100000/sum(US$X6.26.20)
###################
data_repo <- "https://raw.githubusercontent.com/umich-cphds/cov-ind-19-data/master/"
data_state=read.csv(paste0(data_repo, Sys.Date()-4, "/covid19india_data.csv"),sep='\t')

top20case=c("Andhra Pradesh","Assam","Bihar","Delhi","Madhya Pradesh",
            "Gujarat","Haryana","Jharkhand","Jammu and Kashmir","Karnataka","Kerala",
            "Maharashtra","Odisha","Punjab","Rajasthan","Tamil Nadu","Telangana","Uttar Pradesh",
            "Uttarakhand","West Bengal")
N_state=read.csv("indian states population.csv",sep=',',header = T)
N_state1=N_state %>% filter(State.or.union.territory %in% top20case) %>%
  select(State.or.union.territory,Population)
N_state1$Population=as.numeric(gsub(",", "", N_state1$Population))

data_state1 =
  data_state %>%
  filter(Name %in% top20case)
data_state1$Date=as.Date(data_state1$Date)
indicator=unique(data_state1$State)
data_state2 =
  data_state1 %>%
  filter(State %in% indicator[i])
nameindicator=data_state1 %>%filter(State %in% indicator[i])%>%select(Name)
N_state=as.numeric(N_state1 %>%filter(State.or.union.territory %in% nameindicator[1,1])%>%select(Population))
cases_state=diff(c(0,data_state2$Cases))
if(sum(cases_state<0)>0){cases_state=ifelse(cases_state<0,0,cases_state)}
recovered_state=diff(c(0,data_state2$Recovered))
if(sum(recovered_state<0)>0){recovered_state=ifelse(recovered_state<0,0,recovered_state)}
deaths_state=diff(c(0,data_state2$Deaths))
if(sum(deaths_state<0)>0){deaths_state=ifelse(deaths_state<0,0,deaths_state)}

start_point=which(cases_state>0)[1]
H0 <- data_state2[start_point,3]*0.5
R0 <- data_state2[start_point,4]
D0 <- data_state2[start_point,5]
I0 <- data_state2[start_point,3]-H0-R0-D0
if(I0<0){I0=0}
propA0 <- 1
A0 <- I0 * propA0
E0 <- (sum(cases_state[start_point:(start_point+5)])) * (1 + propA0)
S0 <- N_state - E0 - I0 - A0 - H0 - R0 - D0
yinit <- c(S = S0, E = E0, I = I0, R = R0, H = H0, A = A0, D = D0)
time=seq.Date(from=min(data_state2$Date),to=as.Date("2020/06/18"),by='days')
state=cbind.data.frame(cases_state,recovered_state,deaths_state)[c(start_point:length(time)),]
myfix_pars <- c(a = 1, De = 5.1, delta=h_r/7 , gamma = (1-0.01)/22.3, eta = 0.85/22.3, omega=0.01/15.4, nu=0.15/15.4  , N = N_state)
############
fun_pred <- function(init_obs, times, pars, myfix_pars) {
  ## myode function
  myode <- function(stage_pars, fix_pars, old_values) {
    ## stage pars
    b <- stage_pars[1]
    r <- stage_pars[2]
    n <- stage_pars[3]
    ## fixed pars
    a <- fix_pars[1]
    De <- fix_pars[2]
    delta <- fix_pars[3]
    gamma <- fix_pars[4]
    eta <- fix_pars[5]
    omega <- fix_pars[6]
    nu <- fix_pars[7]
    N <- fix_pars[8]
    ## old values
    S <- old_values[1]
    E <- old_values[2]
    I <- old_values[3]
    R <- old_values[4]
    H <- old_values[5]
    A <- old_values[6]
    D <- old_values[7]
    N1 <- S+E+A+R
    ## new values
    S_new <- S - b * S * (I + a * A) / (N-H-D) + n - n * S / N1
    E_new <- E + b * S * (I + a * A) / (N-H-D) - E / De - n * E / N1
    I_new <- I + E * r / De  - (gamma+omega)*(1-delta*7)*I-(delta)*I
    R_new <- R + (A + I*(1-delta*7)) *gamma  + H *eta - n * R / N1
    H_new <- H + I *delta - H *(nu+eta)
    A_new <- A + E * (1 - r) / De - A *(omega+gamma) - n * A / N1
    D_new <- D+omega*(A+I*(1-delta*7))+nu*H
    estN_case <- E * r / De
    estN_recovred <- (I*(1-delta*7)) * gamma +eta*H
    estN_death <- (I*(1-delta*7)) * omega +nu*H
    ##
    return(c(S_new, E_new, I_new, R_new, H_new, A_new,D_new, estN_case, estN_recovred, estN_death))
  }
  ## matrix for indiave the results
  ymat <- matrix(0, length(times), length(init_obs) + 4)
  ymat[, 1] <- times
  colnames(ymat) <- c("time", names(init_obs), "estN_case", "estN_recovred", "estN_death")
  ## stage 1
  mystage_pars <- c(b = pars[1], r = pars[2], n = 0)
  for(i in 1:length(times)) {
    if(i == 1) {
      myold_values <- init_obs
    } else {
      myold_values <- ymat[i - 1, 2:8]
    }
    ymat[i, 2:11] <- myode(stage_pars = mystage_pars, fix_pars = myfix_pars, old_values = myold_values)
  }
  return(ymat)
}
fun_mle <- function(init_pars, init_state_num, observed_num, opt_num) {
  ## negtive likelihood
  ## parinit <- (b12, b3, b4, b5, r12, r3, r4, r5)
  ftime <- 1:nrow(observed_num)
  negLL_func <- function(pars){
    ypred <- fun_pred(init_obs = init_state_num, times = ftime, pars = pars, myfix_pars=myfix_pars )
    ypred_I <- ypred[, 9]
    ypred_R <- ypred[, 10]
    ypred_D <- ypred[, 11]
    #cat(ypred)
    #cat(ypred)
    try(p <- dpois(observed_num[,1], ypred_I, log=T)+dpois(observed_num[,2], ypred_R, log=T)+ dpois(observed_num[,3], ypred_D, log=T))
    if(any(is.nan(p))){
      logL <- -Inf
    }else{
      logL <- sum(p)
    }
    return(-logL)
  }
  result_mat <- matrix(NA, opt_num, length(init_pars) * 2 + 2)
  colnames(result_mat) <- c("likelihood", "convergence","b1","r1","b1_sd","r1_sd")
  result_mat[, 1] <- -Inf
  for(i in 1:opt_num) {
    passed <- FALSE
    while (!passed) {
      #try(mle_opt <- optim(init_pars, negLL_func,lower=c(0,0,0,0,0,0),upper=c(Inf,Inf,Inf,1,1,1),method="L-BFGS-B", hessian = T))
      try(mle_opt <- optim(init_pars, negLL_func, method="SANN", hessian = T))
      #try(mle_opt <- optim(init_pars, negLL_func, method = "Nelder-Mead", hessian = F))
      #if(all(mle_opt$par[6:10]<=1)){
      passed <- exists("mle_opt")
      #}
    }
    result_mat[i, 1] <- -mle_opt$value
    result_mat[i, 2] <- mle_opt$convergence
    result_mat[i, 3:4] <- mle_opt$par
    if(any(result_mat[i, 4]>1)){
      for (j in 4){
        result_mat[i, j] <- ifelse(result_mat[i, j]>1,1,result_mat[i, j])
      }
    }
    result_mat[i, 5:6] <- sqrt(diag(solve(mle_opt$hessian)))
    #result_mat[i, 11:18] <- 0.01
    cat(i, "MLE run finished!", fill = T)
    rm(mle_opt)
  }
  result <- result_mat[which(result_mat[, 1] == max(result_mat[, 1])), ]
  return(result)
}
fun_seir <- function(init_obs, times, pars,myfix_pars) {
  ## mysample function
  mysample <- function(stage_pars, fix_pars, old_values) {

    ## stage pars
    b <- stage_pars[1]
    r <- stage_pars[2]
    n <- stage_pars[3]
    ## fixed pars
    a <- fix_pars[1]
    De <- fix_pars[2]
    delta <- fix_pars[3]
    gamma <- fix_pars[4]
    eta <- fix_pars[5]
    omega <- fix_pars[6]
    nu <- fix_pars[7]
    N <- fix_pars[8]
    ## old values
    S <- old_values[1]
    E <- old_values[2]
    I <- old_values[3]
    R <- old_values[4]
    H <- old_values[5]
    A <- old_values[6]
    D <- old_values[7]
    N1 <- S+E+A+R
    ##
    pS_vec <- c(b * (I + a * A) / (N-H-D), n / N1, 1 - b * (I + a * A) / (N-H-D) - n / N1)
    sample_S <- rmultinom(1, size = S, prob = pS_vec)
    ##
    pE_vec <- c(r / De, (1-r) / De, n / N1, 1 - 1 / De - n / N1)
    sample_E <- rmultinom(1, size = E, prob = pE_vec)
    ##
    pI_vec <- c(gamma*(1-delta*7), omega*(1-delta*7), delta, 1 - (gamma*(1-delta*7)+omega*(1-delta*7)+delta))
    if((gamma*(1-delta*7)+omega*(1-delta*7)+delta) > 1) {
      pI_vec <- c(gamma*(1-delta*7), omega*(1-delta*7),1-(gamma*(1-delta*7)+omega*(1-delta*7)), 0)
    }
    sample_I <- rmultinom(1, size = I, prob = pI_vec)
    ##
    pA_vec <- c(gamma, omega, n / N1, 1 - (gamma+omega) - n / N1)
    sample_A <- rmultinom(1, size = A, prob = pA_vec)
    ##
    pH_vec <- c(eta, nu, 1 - (eta+nu))
    sample_H <- rmultinom(1, size = H, prob = pH_vec)
    ##
    pR_vec <- c(n / N1, 1 - n / N1)
    sample_R <- rmultinom(1, size = R, prob = pR_vec)
    sample_D <- D
    ## new values
    S_new <- sample_S[3] + n
    E_new <- sample_E[4] + sample_S[1]
    I_new <- sample_I[4] + sample_E[1]
    R_new <- sample_R[2] + sample_I[1] + sample_A[1] + sample_H[1]
    H_new <- sample_H[3] + sample_I[3]
    A_new <- sample_A[4] + sample_E[2]
    D_new <- sample_D + sample_I[2] + sample_A[2] + sample_H[2]

    estN_case <- sample_E[1]
    estN_recovred <- sample_I[1] + sample_H[1]
    estN_death <- sample_I[2] + sample_H[2]
    ##
    return(c(S_new, E_new, I_new, R_new, H_new, A_new, D_new, estN_case,estN_recovred, estN_death))
  }
  ## matrix for save the results
  ymat <- matrix(0, length(times), length(init_obs) + 4)
  ymat[, 1] <- times
  colnames(ymat) <- c("time", names(init_obs), "estN_case","estN_recovred", "estN_death")
  ## stage 1
  mystage_pars <- c(b = pars[1], r = pars[2], n = 0)
  for(i in 1:length(times)) {
    if(i == 1) {
      myold_values <- init_obs
    } else {
      myold_values <- ymat[i - 1, 2:8]
    }
    ymat[i, 2:11] <- mysample(stage_pars = mystage_pars, fix_pars = myfix_pars, old_values = myold_values)
  }
  return(ymat)
}
fun_R0 <- function(estpar,myfix_pars) {
  a <- myfix_pars[1]
  De <- myfix_pars[2]
  delta <- myfix_pars[3]
  gamma <- myfix_pars[4]
  eta <- myfix_pars[5]
  omega <- myfix_pars[6]
  nu <- myfix_pars[7]
  N <- myfix_pars[8]
  Dh=7
  b <- estpar[1]
  r <- estpar[2]
  R01 <- (1-r)*(22.3^2*gamma+15.4^2*omega)*a*b+r*b*(h_r*(+eta*22.3^2+nu*15.4^2)+(1-h_r)*(gamma*22.3^2+omega*15.4^2))
  R0 <- c(mean(R01))
  return(R0)
}
fun_mcmc <- function(init_pars, step_pars, init_state_num, observed_num, niter = 1000000, BurnIn = 200000, trace_num = 100, runMLE = T){
  ## use mle to set the initial pars for mcmc
  if(runMLE) {
    mle_optim <- fun_mle(init_pars = init_pars, init_state_num = init_state_num, observed_num = observed_num, opt_num = 3)
    mle_pars <- mle_optim
    init_pars <- mle_pars[3:4]
    cat("The MLE estimates: ", round(init_pars, digits=4), fill = T)
  } else {
    mle_pars <- NULL
    init_pars <- init_pars
  }

  ## function for pars
  ## pars = c(b12, b3, b3, b5, r12, r3, r4, r5)
  Prior_func <- function(pars){
    if(any(pars < 0) || any(pars[2] > 1)){
      return (0)
    }else{
      return (1) ## Set a non-informative prior
    }
  }
  ## function for log likelihood
  LL_func <- function(pars){
    ypred <- fun_pred(init_obs = init_state_num, times = ftime, pars = pars, myfix_pars=myfix_pars)
    ypred_I <- ypred[, 9]
    ypred_R <- ypred[, 10]
    ypred_D <- ypred[, 11]
    #cat(ypred)
    #cat(ypred)
    try(p <- dpois(observed_num[,1], ypred_I, log=T)+dpois(observed_num[,2], ypred_R, log=T)+ dpois(observed_num[,3], ypred_D, log=T))
    if(any(is.nan(p))){
      logL <- -Inf
    }else{
      logL <- sum(p)
    }
    return(logL)
  }
  ## R0 estimation
  ftime <- 1:nrow(observed_num)
  predat <- fun_seir(init_obs = init_state_num, times = ftime, pars = init_pars,myfix_pars=myfix_pars)
  R0_est <- fun_R0(estpar = init_pars,myfix_pars=myfix_pars)
  ## buil the matrix to store the results
  np <- length(init_pars)
  pmat <- matrix(0, ((niter + BurnIn) / trace_num) + 1, np+1) ## parameters + R0 for five periods
  colnames(pmat) <- c("b1","r1","R01")
  pmat[1, 1:np] <- init_pars
  pmat[1, (np+1):(np+1)] <- R0_est
  rm(R0_est, predat)
  ##
  pars_now <- init_pars

  cat("MCMC:", fill = T)
  for(i in 2:(niter + BurnIn)){
    #pars_now <- pmat[i-1, ]
    #pars_new <- rep(0, np)
    pars_new <- rep(0, np)
    for(j in 1:np){
      pars_new[j] <- rnorm(1, mean = pars_now[j], sd = step_pars[j])
    }
    A <- 0
    #cat(pars_new)
    if(Prior_func(pars_new) > 0){
      ll_pars_new <- LL_func(pars = pars_new)
      if(ll_pars_new != -Inf) {
        ll_pars_now <- LL_func(pars = pars_now)
        A <-  10^(ll_pars_new - ll_pars_now) # Prior_func(pars_new) / Prior_func(pars_now) * 10^(ll_pars_new - ll_pars_now)
      }else{
        cat(i, " logL =", ll_pars_new, fill = T)
      }
    }

    if(runif(1) < A){
      pars_now <- pars_new
    }
    if(i %% trace_num == 0) {
      predat <- fun_seir(init_obs = init_state_num, times = ftime, pars = pars_now, myfix_pars=myfix_pars)
      R0_est <- fun_R0(estpar = pars_now,myfix_pars=myfix_pars)
      pmat[(i / trace_num) + 1, 1:np] <- pars_now
      pmat[(i / trace_num) + 1, (np+1):(np+1)] <- R0_est
      rm(R0_est, predat)
    }
    if(i%%10000 == 0) cat("Iter", i, " A =", round(A, digits=5), " : ", round(pars_now, digits=4), fill = T)

  }
  est_list <- list(mle_estimates = mle_pars, mcmc_estimates = pmat[-c(1:(BurnIn / trace_num + 1)), ])
  ## print output
  mcmc_estimates = pmat[-c(1:(BurnIn / trace_num + 1)), ]
  cat("########################   Brief Summary    #####################################", fill = T)
  cat("##              b1  r1  R01  ", fill = T)
  cat("##  Mean :", round(apply(mcmc_estimates, 2, mean), 3), fill = T)
  cat("##  2.5% :", round(apply(mcmc_estimates, 2, function(x) quantile(x, 0.025)), 3), fill = T)
  cat("## 97.5% :", round(apply(mcmc_estimates, 2, function(x) quantile(x, 0.975)), 3), fill = T)
  cat("###############################   Done    ########################################", fill = T)
  return(est_list)
}
pars_estimate_state <- fun_mcmc(init_pars = c(1.6,0.2), step_pars =  c(1.6,0.2)/ 100, init_state_num = yinit, observed_num = state, niter = 1000000, BurnIn = 200000, trace_num = 100, runMLE = T)
mcmc_pars_estimate_state <- round(pars_estimate_state$mcmc_estimates, 3)
write.table(mcmc_pars_estimate_state, paste0("mcmc_pars_estimate_2_",indicator[i],".txt"), quote = F, row.names = F, sep = "\t")

