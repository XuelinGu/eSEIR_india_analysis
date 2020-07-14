library(tidyverse)
library(glue)
set.seed(20200601)
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
a=data.frame()
for (i in 1:20){
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
  yinit=as.data.frame(t(yinit))
  yinit$group=indicator[i]
  a=rbind.data.frame(a,yinit)
}
#########################
#####################2 time periods
###########################
##############combine all state pars results
c=data.frame()
for (i in 1:20){
  a1=read.csv(paste0('india/mcmc_pars_estimate_2_',indicator[i],'.txt'),header= T,sep = '')
  b=as.data.frame(t(round(apply(a1,2, mean),3)))
  b1=as.data.frame(t(round(apply(a1, 2, function(x) quantile(x, 0.975)), 3)))
  b2=as.data.frame(t(round(apply(a1, 2, function(x) quantile(x, 0.025)), 3)))
  b=cbind.data.frame(b,b1,b2)
  b$group=indicator[i]
  c=rbind.data.frame(c,b)
}
d=merge(a,c,by='group')
write.table(d, "india/mcmc_pars_estimate_of_all_india_states_2.txt", quote = F, row.names = F, sep = "\t")
#################estimate and predict using parameters estimated
myfix_pars <- c(a = 1, De = 5.1, delta=h_r/7 , gamma = (1-0.01)/22.3, eta = 0.85/22.3, omega=0.01/15.4, nu=0.15/15.4  , N = N_state)
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
mydate <- c(paste("Mar", 15:31), paste("Apr", 1:30), paste("May", 1:31),paste("Jun", 1:30),paste("Jul", 1:31),paste("Aug", 1:30))
ptime <- 1:length(mydate)
prediction_state=data.frame()
for(i in 1:20){
  prediction_state_tmp=data.frame()
  data_state2 =
    data_state1 %>%
    filter(State %in% indicator[i])
  nameindicator=data_state1 %>%filter(State %in% indicator[i])%>%select(Name)
  cases_state=diff(c(0,data_state2$Cases))
  if(sum(cases_state<0)>0){cases_state=ifelse(cases_state<0,0,cases_state)}
  recovered_state=diff(c(0,data_state2$Recovered))
  if(sum(recovered_state<0)>0){recovered_state=ifelse(recovered_state<0,0,recovered_state)}
  deaths_state=diff(c(0,data_state2$Deaths))
  if(sum(deaths_state<0)>0){deaths_state=ifelse(deaths_state<0,0,deaths_state)}
  start_point=which(cases_state>0)[1]
  H0 <- d[i,6]
  R0 <- d[i,5]
  D0 <- d[i,8]
  I0 <- d[i,4]
  A0 <- d[i,7]
  E0 <- d[i,3]
  S0 <- d[i,2]
  yinit <- c(S = S0, E = E0, I = I0, R = R0, H = H0, A = A0, D = D0)
  time=seq.Date(from=min(data_state2$Date),to=as.Date("2020/06/26"),by='days')
  prediction=as.data.frame(fun_pred(init_obs=yinit,times=ptime,pars=as.numeric(d[i,9:10]),myfix_pars=myfix_pars))
  prediction$group=indicator[i]
  prediction$daily_confirmed=c(cases_state[start_point:length(time)],rep(NA, 169-104+start_point-1))
  prediction$daily_death=c(deaths_state[start_point:length(time)],rep(NA, 169-104+start_point-1))
  prediction$daily_recovered=c(recovered_state[start_point:length(time)],rep(NA, 169-104+start_point-1))
  prediction$confirmed=c(data_state2$Cases[start_point:length(time)],rep(NA, 169-104+start_point-1))
  prediction$death=c(data_state2$Deaths[start_point:length(time)],rep(NA, 169-104+start_point-1))
  prediction$recovered=c(data_state2$Recovered[start_point:length(time)],rep(NA, 169-104+start_point-1))
  prediction$Date=seq.Date(from=data_state2$Date[start_point],to=data_state2$Date[start_point]+168,by='days')
  prediction_state_tmp=prediction
  
  prediction=as.data.frame(fun_pred(init_obs=yinit,times=ptime,pars=as.numeric(d[i,12:13]),myfix_pars=myfix_pars))
  prediction_state_tmp=cbind.data.frame(prediction_state_tmp,prediction)
  
  prediction=as.data.frame(fun_pred(init_obs=yinit,times=ptime,pars=as.numeric(d[i,15:16]),myfix_pars=myfix_pars))
  prediction_state_tmp=cbind.data.frame(prediction_state_tmp,prediction)
  
  prediction_state=rbind.data.frame(prediction_state,prediction_state_tmp)
}

colnames(prediction_state)[21:30] <- c("S_lower","E_lower","I_lower","R_lower","H_lower","A_lower","D_lower","estN_case_lower","estN_recovred_lower","estN_death_lower" )
colnames(prediction_state)[32:41] <- c("S_upper","E_upper","I_upper","R_upper","H_upper","A_upper","D_upper","estN_case_upper","estN_recovred_upper","estN_death_upper" )

write.table(prediction_state, "india/estimate_and_prediction_of_all_india_states_2.txt", quote = F, row.names = F, sep = "\t")
#####################
#################4 time periods
########################
#################combine all state pars results
c=data.frame()
for (i in 1:20){
  a1=read.csv(paste0('india/mcmc_pars_estimate_',indicator[i],'.txt'),header= T,sep = '')
  b=as.data.frame(t(round(apply(a1,2, mean),3)))
  b1=as.data.frame(t(round(apply(a1, 2, function(x) quantile(x, 0.975)), 3)))
  b2=as.data.frame(t(round(apply(a1, 2, function(x) quantile(x, 0.025)), 3)))
  b=cbind.data.frame(b,b1,b2)
  b$group=indicator[i]
  c=rbind.data.frame(c,b)
}
d=merge(a,c,by='group')
write.table(d, "india/mcmc_pars_estimate_of_all_india_states_4.txt", quote = F, row.names = F, sep = "\t")
################estimate and predict using parameters estimated
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
  mystage_pars <- c(b = pars[1], r = pars[5], n = N_state*0.00002)
  for(i in 1:17) {
    if(i == 1) {
      myold_values <- init_obs
    } else {
      myold_values <- ymat[i - 1, 2:8]
    }
    ymat[i, 2:11] <- myode(stage_pars = mystage_pars, fix_pars = myfix_pars, old_values = myold_values)
  }
  mystage_pars <- c(b = pars[2], r = pars[6], n = 0)
  for(i in 18:31) {
    myold_values <- ymat[i - 1, 2:8]
    ymat[i, 2:11] <- myode(stage_pars = mystage_pars, fix_pars = myfix_pars, old_values = myold_values)
  }
  mystage_pars <- c(b = pars[3], r = pars[7], n = 0)
  for(i in 32:50) {
    myold_values <- ymat[i - 1, 2:8]
    ymat[i, 2:11] <- myode(stage_pars = mystage_pars, fix_pars = myfix_pars, old_values = myold_values)
  }
  mystage_pars <- c(b = pars[4], r = pars[8], n = 0)
  for(i in 51:length(times)) {
    myold_values <- ymat[i - 1, 2:8]
    ymat[i, 2:11] <- myode(stage_pars = mystage_pars, fix_pars = myfix_pars, old_values = myold_values)
  }
  return(ymat)
}
prediction_state=data.frame()
for(i in 1:20){
  prediction_state_tmp=data.frame()
  data_state2 =
    data_state1 %>%
    filter(State %in% indicator[i])
  nameindicator=data_state1 %>%filter(State %in% indicator[i])%>%select(Name)
  cases_state=diff(c(0,data_state2$Cases))
  if(sum(cases_state<0)>0){cases_state=ifelse(cases_state<0,0,cases_state)}
  recovered_state=diff(c(0,data_state2$Recovered))
  if(sum(recovered_state<0)>0){recovered_state=ifelse(recovered_state<0,0,recovered_state)}
  deaths_state=diff(c(0,data_state2$Deaths))
  if(sum(deaths_state<0)>0){deaths_state=ifelse(deaths_state<0,0,deaths_state)}
  start_point=which(cases_state>0)[1]
  H0 <- d[i,6]
  R0 <- d[i,5]
  D0 <- d[i,8]
  I0 <- d[i,4]
  A0 <- d[i,7]
  E0 <- d[i,3]
  S0 <- d[i,2]
  yinit <- c(S = S0, E = E0, I = I0, R = R0, H = H0, A = A0, D = D0)
  time=seq.Date(from=min(data_state2$Date),to=as.Date("2020/06/26"),by='days')
  prediction=as.data.frame(fun_pred(init_obs=yinit,times=ptime,pars=as.numeric(d[i,9:16]),myfix_pars=myfix_pars))
  prediction$group=indicator[i]
  prediction$daily_confirmed=c(cases_state[start_point:length(time)],rep(NA, 169-104+start_point-1))
  prediction$daily_death=c(deaths_state[start_point:length(time)],rep(NA, 169-104+start_point-1))
  prediction$daily_recovered=c(recovered_state[start_point:length(time)],rep(NA, 169-104+start_point-1))
  prediction$confirmed=c(data_state2$Cases[start_point:length(time)],rep(NA, 169-104+start_point-1))
  prediction$death=c(data_state2$Deaths[start_point:length(time)],rep(NA, 169-104+start_point-1))
  prediction$recovered=c(data_state2$Recovered[start_point:length(time)],rep(NA, 169-104+start_point-1))
  prediction$Date=seq.Date(from=data_state2$Date[start_point],to=data_state2$Date[start_point]+168,by='days')
  prediction_state_tmp=prediction
  
  prediction=as.data.frame(fun_pred(init_obs=yinit,times=ptime,pars=as.numeric(d[i,21:28]),myfix_pars=myfix_pars))
  prediction_state_tmp=cbind.data.frame(prediction_state_tmp,prediction)
  
  prediction=as.data.frame(fun_pred(init_obs=yinit,times=ptime,pars=as.numeric(d[i,33:40]),myfix_pars=myfix_pars))
  prediction_state_tmp=cbind.data.frame(prediction_state_tmp,prediction)
  
  prediction_state=rbind.data.frame(prediction_state,prediction_state_tmp)
}

colnames(prediction_state)[21:30] <- c("S_lower","E_lower","I_lower","R_lower","H_lower","A_lower","D_lower","estN_case_lower","estN_recovred_lower","estN_death_lower" )
colnames(prediction_state)[32:41] <- c("S_upper","E_upper","I_upper","R_upper","H_upper","A_upper","D_upper","estN_case_upper","estN_recovred_upper","estN_death_upper" )
write.table(prediction_state, "india/estimate_and_prediction_of_all_india_states_4.txt", quote = F, row.names = F, sep = "\t")