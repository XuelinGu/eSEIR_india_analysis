library(tidyverse)
library(glue)
set.seed(20200601)
##########hospital  https://www.cdc.gov/coronavirus/2019-ncov/covid-data/covidview/index.html
US=read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv")
h_r=98.4*326687501/100000/sum(US$X6.26.20)
############eSEIR model analysis for indian data
data_repo <- "https://raw.githubusercontent.com/umich-cphds/cov-ind-19-data/master/"
data_india=read.csv(paste0(data_repo, Sys.Date()-4 , "/jhu_data_mod.csv"),sep='\t')%>%
  filter(Country == "India" )
data_india$Date=as.Date(data_india$Date)
colnames(data_india)
cases_india=diff(c(0,data_india$Cases))
recovered_india=diff(c(0,data_india$Recovered))
deaths_india=diff(c(0,data_india$Deaths))
N_india=1340000000
##################
start_point=which(data_india$Date==as.Date('2020/03/15'))
H0 <- data_india[start_point-7,3]*0.5
R0 <- data_india[start_point,4]
D0 <- data_india[start_point,5]
I0 <- data_india[start_point,3]-H0-R0-D0
propA0 <- 1   
A0 <- I0 * propA0
E0 <- (sum(cases_india[start_point:(start_point+5)])) * (1 + propA0)
S0 <- N_india - E0 - I0 - A0 - H0 - R0 - D0
yinit <- c(S = S0, E = E0, I = I0, R = R0, H = H0, A = A0, D = D0)
pars_start <- c(1.6,0.4,0.4,0.4,0.2,0.2,0.2,0.2)  
time=seq.Date(from=min(data_india$Date),to=as.Date("2020/06/18"),by='days')
india=cbind.data.frame(cases_india,recovered_india,deaths_india)[c(start_point:length(time)),]
############
####estimate
#####there is no uncertainty/95%CI
#a1=read.csv('mcmc_pars_estimate_india_2.txt',header= T,sep = '')
#b=as.data.frame(t(round(apply(a1,2, mean),3)))
#b1=as.data.frame(t(round(apply(a1, 2, function(x) quantile(x, 0.975)), 3)))
#b2=as.data.frame(t(round(apply(a1, 2, function(x) quantile(x, 0.025)), 3)))
#b=cbind.data.frame(b,b1,b2)
#b$group='india'
#d=b
myfix_pars <- c(a = 1, De = 5.1, delta=h_r/7 , gamma = (1-0.0138)/22.3, eta = 0.85/22.3, omega=0.0138/15.4, nu=0.15/15.4  , N = N_india)
############
mydate <- c(paste("Mar", 15:31), paste("Apr", 1:30), paste("May", 1:31),paste("Jun", 1:30),paste("Jul", 1:31),paste("Aug", 1:30))
ptime <- 1:length(mydate)
prediction_state=data.frame()

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
time=seq.Date(from=min(data_india$Date),to=as.Date("2020/06/26"),by='days')
prediction=as.data.frame(fun_pred(init_obs=yinit,times=ptime,pars=c(0.18,1),myfix_pars=myfix_pars))
prediction$group='India'
prediction$daily_confirmed=c(cases_india[start_point:length(time)],rep(NA, 169-104))
prediction$daily_death=c(deaths_india[start_point:length(time)],rep(NA,169-104))
prediction$daily_recovered=c(recovered_india[start_point:length(time)],rep(NA, 169-104))
prediction$confirmed=c(data_india$Cases[start_point:length(time)],rep(NA, 169-104))
prediction$death=c(data_india$Deaths[start_point:length(time)],rep(NA, 169-104))
prediction$recovered=c(data_india$Recovered[start_point:length(time)],rep(NA, 169-104))
prediction$Date=seq.Date(from=data_india$Date[start_point],to=data_india$Date[start_point]+168,by='days')
write.table(prediction, "india/estimate_and_prediction_of_india_2.txt", quote = F, row.names = F, sep = "\t")

############
####estimate_4_time period
#####
a1=read.csv('mcmc_pars_estimate_india_4.txt',header= T,sep = '')
b=as.data.frame(t(round(apply(a1,2, mean),3)))
b1=as.data.frame(t(round(apply(a1, 2, function(x) quantile(x, 0.975)), 3)))
b2=as.data.frame(t(round(apply(a1, 2, function(x) quantile(x, 0.025)), 3)))
b=cbind.data.frame(b,b1,b2)
b$group='india'
d=b
myfix_pars <- c(a = 1, De = 5.1, delta=h_r/7 , gamma = (1-0.0138)/22.3, eta = 0.85/22.3, omega=0.0138/15.4, nu=0.15/15.4  , N = N_india)
############
mydate <- c(paste("Mar", 15:31), paste("Apr", 1:30), paste("May", 1:31),paste("Jun", 1:30),paste("Jul", 1:31),paste("Aug", 1:30))
ptime <- 1:length(mydate)
prediction_state=data.frame()

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
  mystage_pars <- c(b = pars[1], r = pars[5], n = N_india*0.00002)
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

time=seq.Date(from=min(data_india$Date),to=as.Date("2020/06/26"),by='days')
prediction=as.data.frame(fun_pred(init_obs=yinit,times=ptime,pars=as.numeric(d[1:8]),myfix_pars=myfix_pars))
prediction$group='India'
prediction$daily_confirmed=c(cases_india[start_point:length(time)],rep(NA, 169-104))
prediction$daily_death=c(deaths_india[start_point:length(time)],rep(NA,169-104))
prediction$daily_recovered=c(recovered_india[start_point:length(time)],rep(NA, 169-104))
prediction$confirmed=c(data_india$Cases[start_point:length(time)],rep(NA, 169-104))
prediction$death=c(data_india$Deaths[start_point:length(time)],rep(NA, 169-104))
prediction$recovered=c(data_india$Recovered[start_point:length(time)],rep(NA, 169-104))
prediction$Date=seq.Date(from=data_india$Date[start_point],to=data_india$Date[start_point]+168,by='days')
write.table(prediction, "india/estimate_and_prediction_of_india_4.txt", quote = F, row.names = F, sep = "\t")

