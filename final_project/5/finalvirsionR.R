## ------------------------------------------------------------------------
require(ggplot2)
#sample <- read.csv("hep-a-mi.csv",header = TRUE)
#head(sample)
hepa <- read.csv("hep-a-mi2.csv",header = TRUE)
head(hepa)
hep <- hepa$population
num <- hepa$number
plot(hep~num,type='l')


## ------------------------------------------------------------------------
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  encoding="UTF-8"
  )
set.seed(594709947L)
require(ggplot2)
theme_set(theme_bw())
require(plyr)
require(reshape2)
require(magrittr)
require(foreach)
require(pomp)
stopifnot(packageVersion("pomp")>="0.69-1")
bsflu_data <- read.csv("hep-a-mi2.csv")

## ------------------------------------------------------------------------
bsflu_statenames <- c("S","E","I","C","R")
bsflu_paramnames <- c("Beta","mu_I","rho","mu_R1","mu_R2","mu_R3","mu_R4","mu_R5")
bsflu_obsnames <- colnames(bsflu_data)[2]
bsflu_dmeasure <- "
  lik = dpois(population,rho*R+1e-6,give_log);
"

bsflu_rmeasure <- "
  population= rpois(rho*R+1e-6);

"

bsflu_rprocess <- "
  double t1 = rbinom(S,1-exp(-Beta*I*dt));
  double t2 = rbinom(E,1-exp(-dt*mu_I));
  double t3 = rbinom(I,1-exp(-dt*mu_R1));
  double t4 = rbinom(C,1-exp(-dt*mu_R2));
  double t5 = rbinom(R,1-exp(-dt*mu_R3));
  double t6 = rbinom(I,1-exp(-dt*mu_R4));
  double t7 = rbinom(S,1-exp(-dt*mu_R5));
  S -= t1-t5+t7;
  E += t1 - t2;
  I += t2 - t3-t6;
  C += t3 - t4;
  R += t4 - t5+t6+t7;
"

bsflu_fromEstimationScale <- "
 TBeta = exp(Beta);
 Tmu_I = exp(mu_I);
 Tmu_R1 = exp(mu_R1);
 Tmu_R2 = exp(mu_R2);
 Tmu_R3 = exp(mu_R3);
 Tmu_R4 = exp(mu_R4);
 Tmu_R5 = exp(mu_R5);
 Trho = expit(rho);
"

bsflu_toEstimationScale <- "
 TBeta = log(Beta);
 Tmu_I = log(mu_I);
 Tmu_R1 = log(mu_R1);
 Tmu_R2 = log(mu_R2);
 Tmu_R3 = log(mu_R3);
 Tmu_R4 = log(mu_R4);
 Tmu_R5 = log(mu_R5);
 Trho = logit(rho);
"

bsflu_initializer <- "
 S=15000;
 E=0;
 I=0;
 C=0;
 R=0;
"

require(pomp)
stopifnot(packageVersion("pomp")>="0.75-1")
bsflu2 <- pomp(
  data=bsflu_data,
  times="number",
  t0=0,
  rprocess=euler.sim(
    step.fun=Csnippet(bsflu_rprocess),
    delta.t=1
  ),
  rmeasure=Csnippet(bsflu_rmeasure),
  dmeasure=Csnippet(bsflu_dmeasure),
  fromEstimationScale=Csnippet(bsflu_fromEstimationScale),
  toEstimationScale=Csnippet(bsflu_toEstimationScale),
  obsnames = bsflu_obsnames,
  statenames=bsflu_statenames,
  paramnames=bsflu_paramnames,
  initializer=Csnippet(bsflu_initializer)
)

## ------------------------------------------------------------------------
plot(bsflu2)

## ------------------------------------------------------------------------
run_level <- 1
switch(run_level,
       {bsflu_Np=200; bsflu_Nmif=10; bsflu_Neval=10; bsflu_Nglobal=10; bsflu_Nlocal=10},
       {bsflu_Np=2000; bsflu_Nmif=100; bsflu_Neval=10; bsflu_Nglobal=10; bsflu_Nlocal=10},
       {bsflu_Np=20000; bsflu_Nmif=300; bsflu_Neval=10; bsflu_Nglobal=10; bsflu_Nlocal=10},
       
)



require(doParallel)
cores <- 20  # The number of cores on this machine 
registerDoParallel(cores)
mcopts <- list(set.seed=TRUE)

set.seed(396658101,kind="L'Ecuyer")
bsflu_box <- rbind(
  Beta=c(0.001,0.01),
  mu_I=c(0.5,2),
  mu_R1=c(0.5,2),
  mu_R2=c(0.5,2),
  mu_R3=c(0.5,2),
  mu_R4=c(0.5,2),
  mu_R5=c(0.5,2),
  rho = c(0.5,1)
)

bsflu_rw.sd <- 0.02
bsflu_cooling.fraction.50 <- 0.5
hep_rw.sd <-0.2
stew(file=sprintf("hepeefdinalso_eval-%d.rda",run_level),{
  
  t_global <- system.time({
    mifs_global <- foreach(i=1:bsflu_Nglobal,.packages='pomp', .combine=c, .options.multicore=mcopts) %dopar% 
      mif2(     
      bsflu2,
      start=c(apply(bsflu_box,1,function(x)runif(1,x[1],x[2]))),
      Np=bsflu_Np,
      Nmif=bsflu_Nmif,
      cooling.type="geometric",
      cooling.fraction.50=bsflu_cooling.fraction.50,
      transform=TRUE,
      rw.sd=rw.sd(
        Beta=bsflu_rw.sd,
        mu_I=bsflu_rw.sd,
        rho=bsflu_rw.sd,
        mu_R1=hep_rw.sd,
        mu_R2=hep_rw.sd,
        mu_R3=hep_rw.sd,
        mu_R4=hep_rw.sd,
        mu_R5=hep_rw.sd
      )
    )
  })
},seed=1270401374,kind="L'Ecuyer")

## ------------------------------------------------------------------------




## ------------------------------------------------------------------------
stew(file=sprintf("hepeefinall-lik_global_eval-%d.rda",run_level),{
  t_global_eval <- system.time({
    liks_global <- foreach(i=1:bsflu_Nglobal,.packages='pomp',.combine=rbind, .options.multicore=mcopts) %dopar% {
      evals <- replicate(bsflu_Neval, logLik(pfilter(bsflu2,params=coef(mifs_global[[i]]),Np=bsflu_Np)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=442141592,kind="L'Ecuyer")

results_global <- data.frame(logLik=liks_global[,1],logLik_se=liks_global[,2],t(sapply(mifs_global,coef)))
summary(results_global$logLik,digits=5)


## ------------------------------------------------------------------------
t_global

## ------------------------------------------------------------------------
pairs(~logLik+Beta+rho+mu_I+mu_R1+mu_R2+mu_R3+mu_R4+mu_R5,data=subset(results_global,logLik>max(logLik)-3500))

## ------------------------------------------------------------------------
idx <- which.max(results_global$logLik)
results_global[idx,]

## ------------------------------------------------------------------------
sims <- simulate(bsflu2,params=coef(mifs_global[[idx]]),nsim=1,as.data.frame=TRUE,include.data=TRUE)


