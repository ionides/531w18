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
AABA_statenames <- c("H","G","Y_state")
AABA_rp_names <- c("sigma_nu","mu_h","phi","sigma_eta")
AABA_ivp_names <- c("G_0","H_0")
AABA_paramnames <- c(AABA_rp_names,AABA_ivp_names)
AABA_covarnames <- "covaryt"

expit<-function(real){1/(1+exp(-real))}
logit<-function(p.arg){log(p.arg/(1-p.arg))}

rproc1 <- "
double beta,omega,nu;
omega = rnorm(0,sigma_eta * sqrt( 1- phi*phi ) * sqrt(1-tanh(G)*tanh(G)));
nu = rnorm(0, sigma_nu);
G += nu;
beta = Y_state * sigma_eta * sqrt( 1- phi*phi );
H = mu_h*(1 - phi) + phi*H + beta * tanh( G ) * exp(-H/2) + omega;
"
rproc2.sim <- "
Y_state = rnorm( 0,exp(H/2) );
"

rproc2.filt <- "
Y_state = covaryt;
"
AABA_rproc.sim <- paste(rproc1,rproc2.sim)
AABA_rproc.filt <- paste(rproc1,rproc2.filt)

AABA_initializer <- "
G = G_0;
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
Tsigma_nu = log(sigma_nu);
Tphi = logit(phi);
"

AABA_fromEstimationScale <- "
Tsigma_eta = exp(sigma_eta);
Tsigma_nu = exp(sigma_nu);
Tphi = expit(phi);
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
  sigma_nu = exp(-4.5),  
  mu_h = -0.25,       
  phi = expit(4),     
  sigma_eta = exp(-0.07),
  G_0 = 0,
  H_0=0
)

run_level <- 4 
AABA_Np <-          c(100,1e3,2e3,1e4)
AABA_Nmif <-        c(10, 100,200,300)
AABA_Nreps_eval <-  c(4,  10,  20,20)
AABA_Nreps_local <- c(10, 20, 20,20)
AABA_Nreps_global <-c(10, 20, 100,100)

AABA_rw.sd_rp <- 0.02
AABA_rw.sd_ivp <- 0.1
AABA_cooling.fraction.50 <- 0.5

require(doParallel)
cores <-20
registerDoParallel(cores)
mcopts <- list(set.seed=TRUE)

stew("local_test_unfixed1.rda",{
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
                            sigma_nu  = AABA_rw.sd_rp,
                            mu_h      = AABA_rw.sd_rp,
                            phi       = AABA_rw.sd_rp,
                            sigma_eta = AABA_rw.sd_rp,
                            G_0       = ivp(AABA_rw.sd_ivp),
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
},seed=318817884,kind="L'Ecuyer")

It=30
nprof=40
profile.box <- profileDesign(  
  sigma_nu=seq(0.00001,0.5,length.out=It),
  lower=c(mu_h=-7,phi=0.6,sigma_eta=0.01,H_0=-4,G_0=-0.4),
  upper=c(mu_h=-6,phi=0.99,sigma_eta=1,H_0=1,G_0=1.5),
  nprof=nprof
)

stew(file=sprintf("profile_sigmanu-30.rda",It),{
  
  t_global.4 <- system.time({
    prof.llh<- foreach(i=1:(It*nprof),.packages='pomp', .combine=rbind, .options.multicore=mcopts) %dopar%{
      # Find MLE
      mif2(
        if1[[1]],
        start=unlist(profile.box[i,]),
        Np=3000,Nmif=100,
        rw.sd=rw.sd(
          mu_h      = AABA_rw.sd_rp,
          phi       = AABA_rw.sd_rp,
          sigma_eta = AABA_rw.sd_rp,
          G_0       = ivp(AABA_rw.sd_ivp),
          H_0       = ivp(AABA_rw.sd_ivp)
        )
      )->mifs_global.4
      
      evals = replicate(10, logLik(pfilter(mifs_global.4,Np=5000)))
      ll=logmeanexp(evals, se=TRUE)        
      
      data.frame(as.list(coef(mifs_global.4)),
                 loglik = ll[1],
                 loglik.se = ll[2])
    }
  })
},seed=931129,kind="L'Ecuyer")

require(ggplot2)
require(plyr)
require(reshape2)
require(magrittr)
prof.llh %>% 
  ddply(~sigma_nu,subset,rank(-loglik)<=10) %>%
  subset(select=AABA_paramnames) -> pars


## mif2 again
stew(file=sprintf("profile_sigmanu1.rda",It),{
  
  t_global.5 <- system.time({
    prof.llh<- foreach(i=1:(nrow(pars)),.packages='pomp', .combine=rbind, .options.multicore=mcopts) %dopar%{
      # Find MLE
      mif2(
        if1[[1]],
        start=unlist(pars[i,]),
        Np=5000,Nmif=100,
        rw.sd=rw.sd(
          mu_h      = AABA_rw.sd_rp,
          phi       = AABA_rw.sd_rp,
          sigma_eta = AABA_rw.sd_rp,
          G_0       = ivp(AABA_rw.sd_ivp),
          H_0       = ivp(AABA_rw.sd_ivp)
        )
      )->mifs_global.5
      # evaluate llh 
      pf= replicate(10,pfilter(mifs_global.5,Np=5000))
      evals=sapply(pf,logLik)
      ll=logmeanexp(evals, se=TRUE)  
      nfail=sapply(pf,getElement,"nfail")
      
      data.frame(as.list(coef(mifs_global.5)),
                 loglik = ll[1],
                 loglik.se = ll[2],
                 nfail.max=max(nfail))
    }
  })
},seed=931129,kind="L'Ecuyer")

prof.llh %<>%
  subset(nfail.max==0) %>%
  mutate(sigma_nu=signif(sigma_nu,digits=6)) %>%
  ddply(~sigma_nu,subset,rank(-loglik)<=1)

a=max(prof.llh$loglik)
b=a-1.92
CI=which(prof.llh$loglik>=b)
c=prof.llh$sigma_nu[min(CI)]
d=prof.llh$sigma_nu[max(CI)]


prof.llh %>%
  ggplot(aes(x=sigma_nu,y=loglik))+
  geom_point()+
  geom_smooth(method="loess")+
  geom_hline(aes(yintercept=a),linetype="dashed")+
  geom_hline(aes(yintercept=b),linetype="dashed")+
  geom_vline(aes(xintercept=c),linetype="dashed")+
  geom_vline(aes(xintercept=d),linetype="dashed")

c(lower=c,upper=d)