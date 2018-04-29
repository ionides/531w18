require(ggplot2)
require(plyr)
require(reshape2)
require(magrittr)
require(foreach)
require(doParallel)
require(pomp)
require(doMC)
require(tseries)

##switch
cores <- 20
run_level <- 3


ex_data <- read.table("exchange.txt")
ex=log(ex_data$Rate)
exx=diff(ex)
mu0=mean(exx)
delta0=sqrt(var(exx))

ex_statenames <- c("N","e")
ex_paramnames <- c("mu","delta")
ex_dmeasure <- Csnippet("lik = dnorm(Rate,0,N,give_log);")
ex_rmeasure <- Csnippet("Rate = rnorm(0,N);")
ex_rprocess <- Csnippet(" e = rnorm(0,1);N = N*exp((mu-delta*delta/2)+delta*e);")
ex_initializer <-"N=rnorm(6.5,0.3);e=rnorm(1,1);"
stopifnot(packageVersion("pomp")>="0.75-1")
ex_fromEstimationScale <- "Tmu = exp(mu);Tdelta = expit(delta);"
ex_toEstimationScale <- "Tmu = log(mu);Tdelta = logit(delta);"

ex <<- pomp(
  data=ex_data,
  times="Day",
  t0=0,
  rprocess=discrete.time.sim(step.fun=ex_rprocess,delta.t=1),
  rmeasure=ex_rmeasure,
  dmeasure=ex_dmeasure,
  obsnames ="Rate",
  statenames=ex_statenames,
  paramnames=ex_paramnames,
  initializer=Csnippet(ex_initializer))

simulate(ex,params=c(mu=mu0,delta=delta0),nsim=20,states=TRUE) -> x
matplot(time(ex),t(x["N",1:20,]),type='l',lty=1, xlab="time",ylab="Rate",bty='l',col='blue')
lines(time(ex),obs(ex,"Rate"),lwd=2,col='black')


switch(run_level,
       {ex_Np=100; ex_Nmif=10; ex_Neval=5; ex_Nglobal=5; ex_Nlocal=5},
       {ex_Np=1000; ex_Nmif=100; ex_Neval=10; ex_Nglobal=10; ex_Nlocal=10},
       {ex_Np=6000; ex_Nmif=300; ex_Neval=10; ex_Nglobal=100; ex_Nlocal=20})


pf <- pfilter(ex,Np=ex_Np,params=c(mu=-0.0001,delta=0.002))
logLik(pf)

pf <- replicate(10,pfilter(ex,Np=ex_Np,params=c(mu=-0.0001,delta=0.002)))
ll <- sapply(pf,logLik)
logmeanexp(ll,se=TRUE)

cat(c("logLik","logLik_se","mu", "delta"),file="ex_params.txt",append=F,sep="\t")
cat("\n",file="ex_params.txt",append=TRUE,sep="\t")
for (mu in seq(from = -0.05, to = 0.05, length.out = 10)){
  for (delta in seq(from =0.0005, to = 0.005, length.out = 10)){
      pf <- replicate(10,pfilter(ex,Np=ex_Np,params=c(mu=mu,delta=delta)))
      ll <- sapply(pf,logLik)
      cat(c(logmeanexp(ll,se=TRUE)[1],logmeanexp(ll,se=TRUE)[2],mu,delta),file="ex_params.txt",append=TRUE,sep="\t")
      cat("\n",file="ex_params.txt",append=TRUE,sep="\t")
  }
}


sliceDesign(
  c(mu=mu0,delta=delta0),
  mu=rep(seq(from=-0.5,to=0.5,length=40),each=3),
  delta=rep(seq(from=0,to=0.5,length=40),each=3)) -> p


registerDoParallel(cores)
set.seed(998468235L,kind="L'Ecuyer")
mcopts <- list(preschedule=FALSE,set.seed=TRUE)
foreach (theta=iter(p,"row"),.combine=rbind,.inorder=FALSE,.options.multicore=mcopts,.export=ls(globalenv())) %dopar% 
         {
           require(pomp)
           pfilter(ex,params=unlist(theta),Np=ex_Np) -> pf
           theta$loglik <- logLik(pf)
           theta
         } -> p

foreach (v=c("mu","delta")) %do% 
{
  x <- subset(p,slice==v)
  plot(x[[v]],x$loglik,xlab=v,ylab="loglik")
}


expand.grid(mu=seq(from=-0.5,to=0.5,length=50),delta=seq(from=0,to=0.5,length=50)) -> p

foreach (theta=iter(p,"row"),.combine=rbind,
         .inorder=FALSE,.options.multicore=mcopts) %dopar% 
         {
           pfilter(ex,params=unlist(theta),Np=ex_Np) -> pf
           theta$loglik <- logLik(pf)
           theta
         } -> p
pp <- mutate(p,loglik=ifelse(loglik>max(loglik)-10000,loglik,NA))
ggplot(data=pp,mapping=aes(x=mu,y=delta,z=loglik,fill=loglik))+
  geom_tile(color=NA)+
  geom_contour(color='black',binwidth=3)+
  scale_fill_gradient()+
  labs(x=expression(mu),y=expression(delta))


ex2 <<- pomp(
  data=ex_data,
  times="Day",
  t0=0,
  rprocess=euler.sim(step.fun=ex_rprocess,delta.t=1),
  rmeasure=ex_rmeasure,
  dmeasure=ex_dmeasure,
  fromEstimationScale=Csnippet(ex_fromEstimationScale),
  toEstimationScale=Csnippet(ex_toEstimationScale),
  obsnames ="Rate",
  statenames=ex_statenames,
  paramnames=ex_paramnames,
  initializer=Csnippet(ex_initializer))
plot(ex2)

ex_params <- data.matrix(read.table("ex_params.txt",row.names=NULL,header=TRUE))
ex_mle <- ex_params[which.max(ex_params[,"logLik"]),][ex_paramnames]
ex_fixed_params <- c(mu=mu0,delta=delta0)


registerDoParallel(cores)
mcopts <- list(set.seed=TRUE)
set.seed(396658101,kind="L'Ecuyer")
stew(file=sprintf("pf-%d.rda",run_level),{
  t_pf <- system.time(
    pf <- foreach(i=1:20,.packages='pomp', .options.multicore=mcopts,.export=ls(globalenv())) %dopar% try(
                    pfilter(ex2,params=ex_mle,Np=ex_Np)) )
},seed=1320290398,kind="L'Ecuyer")

(L_pf <- logmeanexp(sapply(pf,logLik),se=TRUE))

ex_rw.sd <- 0.002
ex_cooling.fraction.50 <- 0.5

stew(file=sprintf("local_search-%d.rda",run_level),{
  t_local <- system.time({
    mifs_local <- foreach(i=1:ex_Nlocal,.packages='pomp', .combine=c, .options.multicore=mcopts,.export=ls(globalenv())) %dopar%  {
      mif2(ex2,start=ex_mle,Np=ex_Np,
        Nmif=ex_Nmif,cooling.type="geometric",cooling.fraction.50=ex_cooling.fraction.50,transform=TRUE,
        rw.sd=rw.sd(
          mu=ex_rw.sd,
          delta=ex_rw.sd))}})},seed=900242057,kind="L'Ecuyer")

stew(file=sprintf("lik_local-%d.rda",run_level),{
  t_local_eval <- system.time({
    liks_local <- foreach(i=1:ex_Nlocal,.packages='pomp',.combine=rbind,.export=ls(globalenv())) %dopar% {
      evals <- replicate(ex_Neval, logLik(pfilter(ex2,params=coef(mifs_local[[i]]),Np=ex_Np)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=900242057,kind="L'Ecuyer")

results_local <- data.frame(logLik=liks_local[,1],logLik_se=liks_local[,2],t(sapply(mifs_local,coef)))
summary(results_local$logLik,digits=5)
pairs(~logLik+mu+delta,data=subset(results_local,logLik>max(logLik)-500))

