library("stochvol")
dat <- read.table("AABA.csv",sep=",",header=TRUE)
dat$Date <- strptime(dat$Date,"%Y-%m-%d")
dat$year <- as.numeric(format(dat$Date, format="%Y"))
dat$month <- as.numeric(format(dat$Date, format="%m"))
dat$day <- as.numeric(format(dat$Date, format="%d"))
time <- dat$year + dat$month/12 + dat$day/365
N <- nrow(dat)
AABA <- dat$Close[1:N]
log_returns <- diff(log(AABA),lag=1)
demeaned_returns <- logret(AABA,demean=TRUE)

require(pomp)
AABA_statenames <- c("H","Y_state")
AABA_rp_names <- c("rho","mu_h","phi","sigma_eta")
AABA_ivp_names <- c("H_0")
AABA_paramnames <- c(AABA_rp_names,AABA_ivp_names)
AABA_covarnames <- "covaryt"

expit<-function(real){1/(1+exp(-real))}
logit<-function(p.arg){log(p.arg/(1-p.arg))}

rproc1 <- "
double beta,omega;
omega = rnorm(0,sigma_eta * sqrt( 1- phi*phi ) * sqrt(1-rho*rho));
beta = Y_state * sigma_eta * sqrt( 1- phi*phi );
H = mu_h*(1 - phi) + phi*H + beta * rho * exp(-H/2) + omega;
"
rproc2.sim <- "
Y_state = rnorm( 0,exp(H/2) );
"

rproc2.filt <- "
Y_state = covaryt;
"
AABA_rproc.sim  <- paste(rproc1,rproc2.sim)
AABA_rproc.filt <- paste(rproc1,rproc2.filt)

AABA_initializer <- "
  H = H_0;
Y_state = rnorm( 0,exp(H/2) );
"

AABA_rmeasure <- "
   y=Y_state;
"

AABA_dmeasure <- "
lik=dnorm(y,0,exp(H/2),give_log);
"

AABA_toEstimationScale <- "
Tsigma_eta = log(sigma_eta);
Tphi = logit(phi);
Trho = logit(-rho);
"


AABA_fromEstimationScale <- "
Tsigma_eta = exp(sigma_eta);
Tphi = expit(phi);
Trho = -expit(rho);
"

AABA.filt <- pomp(data=data.frame(y=demeaned_returns,
                                  time=1:length(demeaned_returns)),
                  statenames=AABA_statenames,
                  paramnames=AABA_paramnames,
                  covarnames=AABA_covarnames,
                  times="time",
                  t0=0,
                  covar=data.frame(covaryt=c(0,demeaned_returns),
                                   time=0:length(demeaned_returns)),
                  tcovar="time",
                  rmeasure=Csnippet(AABA_rmeasure),
                  dmeasure=Csnippet(AABA_dmeasure),
                  rprocess=discrete.time.sim(step.fun=Csnippet(AABA_rproc.filt),delta.t=1),
                  initializer=Csnippet(AABA_initializer),
                  toEstimationScale=Csnippet(AABA_toEstimationScale), 
                  fromEstimationScale=Csnippet(AABA_fromEstimationScale)
)

params_test <- c(
  mu_h = -1,
  rho = -0.65,
  phi = expit(4),       
  sigma_eta = exp(-0.07),
  H_0 = 0)

run_level <- 3 
AABA_Np <-          c(100,1e3,2e3)
AABA_Nmif <-        c(10, 100,200)
AABA_Nreps_eval <-  c(4,  10,  20)
AABA_Nreps_local <- c(10, 20, 20)
AABA_Nreps_global <-c(10, 20, 100)

require(doParallel)
registerDoParallel()

require(doMC)
registerDoMC()

AABA_rw.sd_rp <- 0.02
AABA_rw.sd_ivp <- 0.1
AABA_cooling.fraction.50 <- 0.5

stew("local_test_fixed.rda",{
  t.if1 <- system.time({
    if1 <- foreach(i=1:AABA_Nreps_local[run_level],
                   .packages='pomp', .combine=c,
                   .options.multicore=list(set.seed=TRUE)) %dopar% try(
                     mif2(AABA.filt,
                          start=params_test,
                          Np=AABA_Np[run_level],
                          Nmif=AABA_Nmif[run_level],
                          cooling.type="geometric",
                          cooling.fraction.50=AABA_cooling.fraction.50,
                          transform=TRUE,
                          rw.sd = rw.sd(
                            rho       = -AABA_rw.sd_rp,
                            mu_h      = AABA_rw.sd_rp,
                            phi       = AABA_rw.sd_rp,
                            sigma_eta = AABA_rw.sd_rp,
                            
                            H_0       = ivp(AABA_rw.sd_ivp)
                          )
                     )
                   )
    
    L.if1 <- foreach(i=1:AABA_Nreps_local[run_level],.packages='pomp',
                     .combine=rbind,.options.multicore=list(set.seed=TRUE)) %dopar% 
                     {
                       logmeanexp(
                         replicate(AABA_Nreps_eval[run_level],
                                   logLik(pfilter(AABA.filt,params=coef(if1[[i]]),Np=AABA_Np[run_level]))
                         ),
                         se=TRUE)
                     }
  })
},seed=318817883,kind="L'Ecuyer")


run_level<-1

AABA_box <- rbind(
  rho=c(-0.2,-0.01),
  mu_h=c(-8,-0.1),
  phi = c(0.6,0.99),
  sigma_eta = c(0.001,15),
  H_0 = c(-1,1)
)

stew(file="global_test_fixed.rda",{
  t.box <- system.time({
    if.box <- foreach(i=1:AABA_Nreps_global[run_level],.packages='pomp',.combine=c,
                      .options.multicore=list(set.seed=TRUE)) %dopar%  
      mif2(
        if1[[1]],
        start=apply(AABA_box,1,function(x)runif(1,x))
      )
    
    L.box <- foreach(i=1:AABA_Nreps_global[run_level],.packages='pomp',.combine=rbind,
                     .options.multicore=list(set.seed=TRUE)) %dopar% {
                       set.seed(87933+i)
                       logmeanexp(
                         replicate(AABA_Nreps_eval[run_level],
                                   logLik(pfilter(AABA.filt,params=coef(if.box[[i]]),Np=AABA_Np[run_level]))
                         ), 
                         se=TRUE)
                     }
  })
},seed=290860873,kind="L'Ecuyer")


r.global_test_fixed <- data.frame(logLik=L.box[,1],logLik_se=L.box[,2],t(sapply(if.box,coef)))
if(run_level>1) write.table(r.global_test_fixed,file="globalfixed_params.csv",append=TRUE,col.names=FALSE,row.names=FALSE)
summary(r.global_test$logLik,digits=5)
