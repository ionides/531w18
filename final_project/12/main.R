setwd("~/Desktop/531w18/final_project_me")
source('functions.R')


# Read data
dat = read.csv('data//mass_shootings.csv')

dat$Date = as.Date(dat$Date,format="%m/%d/%Y")

dat = dat[,c('Fatalities','Injured','Date')]

dat$Year = year(dat$Date)

dat$Month = month(dat$Date)

# Sum the fatalities per month
dat = sql("
select
  Year
  ,Month
  ,sum(Fatalities) as Fatalities
from dat
group by 1,2
having Year <> 2018
order by 1,2
")

#Fill in the missing month-years
dat0 = expand.grid(c(min(dat$Year):max(dat$Year)),
            c(1:12),0)
names(dat0) = c('Year', 'Month', 'Fatalities')

# Data at the year-month level
dat_month = sql("
select
  d0.Year as Year
  ,d0.Month as Month
  ,d0.Fatalities + d.Fatalities as Fatalities
from dat0 d0 left join dat d on d.Year=d0.Year and d.Month=d0.Month
order by 1,2
")

dat_month[is.na(dat_month$Fatalities),'Fatalities'] = 0

plot(dat_month$Fatalities,type = 'l',bty='n')
plot(diff(dat_month$Fatalities),type='l',bty='n')
acf(dat_month$Fatalities)

# Data at the year level
dat_year = sql("
select
  Year
  ,sum(Fatalities) as Fatalities
from dat_month
group by 1
")

plot(dat_year$Fatalities,type = 'l',bty='n')
plot(diff(dat_year$Fatalities),type='l',bty='n')
acf(dat_year$Fatalities)

gg_time_plot = ggplot(dat_year, aes(Year, Fatalities)) + 
  geom_line(size=.8) +
  ggtitle('US Mass Shooting Fatalities 1982 - 2017')
gg_time_plot = formatGG(gg_time_plot)

gg_acf = autoplot(acf(dat_year$Fatalities, plot = FALSE))+
  ggtitle('AutoCorrelation Plot')
gg_acf = formatGG(gg_acf)

dat_demeaned = diff(dat_year$Fatalities) - mean(diff(dat_year$Fatalities))
gg_time_plot_demeaned = ggplot(data.frame(Year = dat_year$Year[2:nrow(dat_year)], Fatalities = dat_demeaned), aes(Year, Fatalities)) + 
  geom_line(size=.8) +
  ggtitle('US Mass Shooting Fatalities\n Differenced/Demeaned')
gg_time_plot_demeaned = formatGG(gg_time_plot_demeaned)

multiplot(gg_time_plot, gg_time_plot_demeaned, gg_acf)

###############################################################################
# RICKER MODEL
###############################################################################

# Create pomp object for year mass shootings

skel = "DN = r*N*exp(-N);" # add in determnistic skeleton

# Add stochastic Ricker model 
stochStep = "
e = rnorm(0,sigma);
N = r*N*exp(-N+e);
"

# Add in the (Poisson) measurement model  
ms_rmeasure = "Fatalities = rpois(phi*N);"
ms_dmeasure = "lik = dpois(Fatalities,phi*N, FALSE);"
ms_statenames = c("N","e")
ms_paramnames = c("phi","r","sigma")
ms_initializer = "
  N = 1;
  e = 6;
"

# Build pomp object
ms_pomp = pomp(
  data = dat_year,
  times = "Year",
  t0 = 1981,
  skeleton = map(Csnippet(skel)),
  rprocess = discrete.time.sim(
    step.fun=Csnippet(stochStep),
    delta.t=1
    ),
  rmeasure = Csnippet(ms_rmeasure),
  dmeasure = Csnippet(ms_dmeasure),
  statenames = ms_statenames,
  paramnames = ms_paramnames,
  initializer = Csnippet(ms_initializer),
  cdir=getwd()
  )

#pomp(ms_pomp,rmeasure=rmeas,dmeasure=dmeas,statenames=c("N"),paramnames=c("phi")) -> ms_pomp


########################################################
# Construct likelihood slices for phi and r parameters
########################################################

sliceDesign(
  c(N.0=1,e.0=6,r=7,sigma=0.6,phi=10),
  phi=rep(seq(from=1,to=20,length=50),each=3),
  r=rep(seq(from=0.1,to=15,length=50),each=3)) -> p

set.seed(998468235L,kind="L'Ecuyer")
mcopts = list(preschedule=FALSE,set.seed=TRUE)

foreach (theta=iter(p,"row"),.combine=rbind,
         .inorder=FALSE,.options.multicore=mcopts) %dopar% 
 {
   pfilter(ms_pomp,params=unlist(theta),Np=5000) -> pf
   theta$loglik <- logLik(pf)
   theta
 } -> p


foreach (v=c("phi","r")) %do% 
{
  x <- subset(p,slice==v)
  plot(x[[v]],x$loglik,xlab=v,ylab="loglik")
}

####################################################################
# Visualize Likelihood Surface for Ricker Model
####################################################################
expand.grid(phi=seq(from=5,to=50,length=50),
            r=seq(from=0,to=15,length=50),
            sigma=0.8,
            N.0 = 1,
            e.0 = 6) -> p

foreach (theta=iter(p,"row"),.combine=rbind,
         .inorder=FALSE,.options.multicore=mcopts) %dopar% 
 {
   pfilter(ms_pomp,params=unlist(theta),Np=5000) -> pf
   theta$loglik <- logLik(pf)
   theta
 } -> p

pp = mutate(p,loglik=ifelse(loglik>max(loglik)-100,loglik,NA))

ggplot(data=pp,mapping=aes(x=r,y=phi,z=loglik,fill=loglik))+
  geom_tile(color=NA)+
  geom_contour(color='black',binwidth=3)+
  scale_fill_gradient()+
  labs(x=expression(r),y=expression(phi))


#####################################################################################
# Run iterated filtering algorithm to maximize likelihood surface
#####################################################################################

run_level = 2
switch(run_level,
       {ms_Np=100; ms_Nmif=10; ms_Neval=10; ms_Nglobal=10; ms_Nlocal=10}, 
       {ms_Np=2000; ms_Nmif=200; ms_Neval=10; ms_Nglobal=10; ms_Nlocal=10}, 
       {ms_Np=60000; ms_Nmif=300; ms_Neval=10; ms_Nglobal=100; ms_Nlocal=20}
)

# From looking at likelihood surface
ms_mle = c(phi = 40, r = 2.5, sigma = 0.8)

ms_rw.sd = 0.02
ms_cooling.fraction.50 = 0.5

######################################################################################
# Conduct LOCAL SEARCH of the likelihood surface using IF2 algorithm to maximize
# likelihood over parameter space
######################################################################################
stew(file=sprintf("local_search_ricker-%d.rda",run_level),{
  
  t_local <- system.time({
    mifs_local <- foreach(i=1:ms_Nlocal,.packages='pomp', .combine=c, .options.multicore=mcopts) %dopar%  {
      mif2(
        ms_pomp,
        start=ms_mle,
        Np=ms_Np,
        Nmif=ms_Nmif,
        cooling.type="geometric",
        cooling.fraction.50=ms_cooling.fraction.50,
        transform=TRUE,
        rw.sd=rw.sd(phi=ms_rw.sd, r=ms_rw.sd, sigma=ms_rw.sd)
      )
    }
  })
  
},seed=900242057,kind="L'Ecuyer")


################################################################################################
# Some parameter perturbations remain in last filtering iteration in mif2 call above, so inference
# is not entirely reliable. Thus, we evaluate the likelihood using replicated particle filters
# at each point estimate
################################################################################################

stew(file=sprintf("lik_local_ricker-%d.rda",run_level),{
    t_local_eval <- system.time({
    liks_local <- foreach(i=1:ms_Nlocal,.packages='pomp',.combine=rbind) %dopar% {
      evals <- replicate(ms_Neval, logLik(pfilter(ms_pomp,params=coef(mifs_local[[i]]),Np=ms_Np)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=900242057,kind="L'Ecuyer")

results_local <- data.frame(logLik=liks_local[,1],logLik_se=liks_local[,2],t(sapply(mifs_local,coef)))
summary(results_local$logLik,digits=5)

# Show geometry of the likelihood surface in a neighborhood of the point estimates
pairs(~logLik+phi+r+sigma,data=subset(results_local,logLik>max(logLik)-50))


######################################################################################
# Conduct GLOBAL SEARCH of the likelihood surface using IF2 algorithm to maximize
# likelihood over parameter space
######################################################################################

# One can specify a large box in parameter space that contains all parameter vectors
# which seem remotely sensible. If an estimation method gives stable conclusions with
# starting values drawn randomly from this box, this gives some confidence that an adequate
# global search has been carried out. 
ms_box = rbind(
  phi=c(30,40),
  r=c(1,4),
  sigma = c(.9,3)
)

# Carry out the likelihood maximizations from diverse starting points
# We reset only the starting parameters mifs_global[[1]] since the rest of the 
# call to mif2 can be read in from mifs_global[[1]]    
stew(file=sprintf("box_eval_ricker-%d.rda",run_level),{
  
  t_global <- system.time({
    mifs_global <- foreach(i=1:ms_Nglobal,.packages='pomp', .combine=c, .options.multicore=mcopts) %dopar%  
      mif2(
        mifs_local[[1]],
        start=c(apply(ms_box,1,function(x)runif(1,x[1],x[2]))),
        transform=TRUE
    )
  })
},seed=1270401374,kind="L'Ecuyer")

# Again, evaluate the likelihood and standard error
# by using replicated particle filters at each point estimate
stew(file=sprintf("lik_global_eval_ricker-%d.rda",run_level),{
  t_global_eval <- system.time({
    liks_global <- foreach(i=1:ms_Nglobal,.packages='pomp',.combine=rbind, .options.multicore=mcopts) %dopar% {
      evals <- replicate(ms_Neval, logLik(pfilter(ms_pomp,params=coef(mifs_global[[i]]),Np=ms_Np)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=442141592,kind="L'Ecuyer")

results_global <- data.frame(logLik=liks_global[,1],logLik_se=liks_global[,2],t(sapply(mifs_global,coef)))
summary(results_global$logLik,digits=5)


plot(mifs_global)


conv.rec(mifs_global)

########################################################################################
# Simulate from POMP model using MLE estimates from above
########################################################################################

# Simulate from whole POMP
# Parameters below appear to maximize likelihood surface
coef(ms_pomp) = c(N.0=1,e.0=6,r=1.8,sigma=0.93,phi=35.3)

sims = simulate(ms_pomp,nsim=3,as.data.frame=TRUE,include.data=TRUE)

dat_sim = sims

dat_sim$Year = dat_year$Year

levels(dat_sim$sim) = c('Data', 'Simulation 1', 'Simulation 2', 'Simulation 3')
gg_sim = ggplot(data=dat_sim,mapping=aes(x=time,y=Fatalities))+
  geom_line()+
  ggtitle('Ricker Simulations')+
  ylab('Mass Shooting Fatalities')+
  xlab('Year')+
  facet_wrap(~sim)
formatGG(gg_sim)

########################################################################################
# Build and simulate from financial stochastic volatility model (notes 14)
########################################################################################

ms_statenames = c("H","G","Y_state")
ms_rp_names = c("sigma_nu","mu_h","phi","sigma_eta")
ms_ivp_names = c("G_0","H_0")
ms_paramnames = c(ms_rp_names,ms_ivp_names)
ms_covarnames = "covaryt"

rproc1 = "
  double beta,omega,nu;
  omega = rnorm(0,sigma_eta * sqrt( 1- phi*phi ) * sqrt(1-tanh(G)*tanh(G)));
  nu = rnorm(0, sigma_nu);
  G += nu;
  beta = Y_state * sigma_eta * sqrt( 1- phi*phi );
  H = mu_h*(1 - phi) + phi*H + beta * tanh( G ) * exp(-H/2) + omega;
"
rproc2.sim = "
  Y_state = rnorm( 0,exp(H/2) );
 "

rproc2.filt = "
  Y_state = covaryt;
 "
ms_rproc.sim = paste(rproc1,rproc2.sim)
ms_rproc.filt = paste(rproc1,rproc2.filt)

ms_initializer = "
  G = G_0;
  H = H_0;
  Y_state = rnorm( 0,exp(H/2) );
"

ms_rmeasure = "
   y=Y_state;
"

ms_dmeasure = "
   lik=dnorm(y,0,exp(H/2),give_log);
"

# Transform parameters to be defined on whole real line

ms_toEstimationScale = "
  Tsigma_eta = log(sigma_eta);
  Tsigma_nu = log(sigma_nu);
  Tphi = logit(phi);
"

ms_fromEstimationScale = "
  Tsigma_eta = exp(sigma_eta);
  Tsigma_nu = exp(sigma_nu);
  Tphi = expit(phi);
"

# Build pomp model for filtering and parameter estimation
ms.filt = pomp(data=data.frame(y=dat_demeaned,
                     time=1:length(dat_demeaned)),
              statenames=ms_statenames,
              paramnames=ms_paramnames,
              covarnames=ms_covarnames,
              times="time",
              t0=0,
              covar=data.frame(covaryt=c(0,dat_demeaned),
                     time=0:length(dat_demeaned)),
              tcovar="time",
              rmeasure=Csnippet(ms_rmeasure),
              dmeasure=Csnippet(ms_dmeasure),
              rprocess=discrete.time.sim(step.fun=Csnippet(ms_rproc.filt),delta.t=1),
              initializer=Csnippet(ms_initializer),
              toEstimationScale=Csnippet(ms_toEstimationScale), 
              fromEstimationScale=Csnippet(ms_fromEstimationScale)
)


########################################################################################
# Simlate from model
########################################################################################
expit<-function(real){1/(1+exp(-real))}
logit<-function(p.arg){log(p.arg/(1-p.arg))}
params_test <- c(
     sigma_nu = exp(-4.5),  
     mu_h = 5,  	 
     phi = expit(3),	 
     sigma_eta = exp(-0.07),
     G_0 = 6,
     H_0=0
  )


sim1.sim <- pomp(ms.filt, 
               statenames=ms_statenames,
               paramnames=ms_paramnames,
               covarnames=ms_covarnames,
               rprocess=discrete.time.sim(step.fun=Csnippet(ms_rproc.sim),delta.t=1)
)


sim1.sim = simulate(sim1.sim,seed=1,params=params_test)

# Build filtering object by copying new simulated data into 
# the covariate slot, and put back apropriate version of rprocess
sim1.filt = pomp(sim1.sim, 
  covar=data.frame(
    covaryt=c(obs(sim1.sim),NA),
    time=c(timezero(sim1.sim),time(sim1.sim))),
  tcovar="time",
  statenames=ms_statenames,
  paramnames=ms_paramnames,
  covarnames=ms_covarnames,
  rprocess=discrete.time.sim(step.fun=Csnippet(ms_rproc.filt),delta.t=1)
)

run_level = 3
ms_Np =          c(100,1e3,5000)
ms_Nmif =        c(10, 100,200)
ms_Nreps_eval =  c(4,  10,  20)
ms_Nreps_local = c(10, 20, 20)
ms_Nreps_global =c(10, 20, 100)

###########################################################################
# Run iterated filtering on the mass shooting data with stochastic volatility model
###########################################################################

ms_rw.sd_rp = 0.02
ms_rw.sd_ivp = 0.1
ms_cooling.fraction.50 = 0.5

stew(file=sprintf("mif2_stochastic_vol-%d.rda",run_level),{
   t.if1 <- system.time({
   if1 <- foreach(i=1:ms_Nreps_local[run_level],
                  .packages='pomp', .combine=c,
                  .options.multicore=list(set.seed=TRUE)) %dopar% try(
                    mif2(ms.filt,
                         start=params_test,
                         Np=ms_Np[run_level],
                         Nmif=ms_Nmif[run_level],
                         cooling.type="geometric",
                         cooling.fraction.50=ms_cooling.fraction.50,
                         transform=TRUE,
                         rw.sd = rw.sd(
                            sigma_nu  = ms_rw.sd_rp,
                            mu_h      = ms_rw.sd_rp,
                            phi       = ms_rw.sd_rp,
                            sigma_eta = ms_rw.sd_rp,
                            G_0       = ivp(ms_rw.sd_ivp),
                            H_0       = ivp(ms_rw.sd_ivp)
                         )
                    )
                  )
    
    L.if1 <- foreach(i=1:ms_Nreps_local[run_level],.packages='pomp',
                      .combine=rbind,.options.multicore=list(set.seed=TRUE)) %dopar% 
                      {
                        logmeanexp(
                          replicate(ms_Nreps_eval[run_level],
                                    logLik(pfilter(ms.filt,params=coef(if1[[i]]),Np=ms_Np[run_level]))
                          ),
                          se=TRUE)
                      }
  })
},seed=318817883,kind="L'Ecuyer")

r.if1 <- data.frame(logLik=L.if1[,1],logLik_se=L.if1[,2],t(sapply(if1,coef)))
summary(r.if1$logLik,digits=5)

pairs(~logLik+sigma_nu+mu_h+phi+sigma_eta,data=subset(r.if1,logLik>max(logLik)-20))

######################################################################################
# Conduct GLOBAL SEARCH of the likelihood surface using IF2 algorithm to maximize
# likelihood over parameter space
######################################################################################

ms_box = rbind(
 sigma_nu=c(.04,0.09),
 mu_h    =c(6,9),
 phi = c(0.1,0.4),
 sigma_eta = c(.05,.5),
 G_0 = c(8,10),
 H_0 = c(1,3)
)

stew(file=sprintf("box_eval_stochastic_vol-%d.rda",run_level),{
  t.box <- system.time({
    
    if.box <- foreach(i=1:ms_Nreps_global[run_level],.packages='pomp',.combine=c,
                  .options.multicore=list(set.seed=TRUE)) %dopar%
      
    mif2(
      if1[[1]],
      start=apply(ms_box,1,function(x)runif(1,x[1],x[2]))
      )
                  
    
    L.box <- foreach(i=1:ms_Nreps_global[run_level],.packages='pomp',.combine=rbind,
                      .options.multicore=list(set.seed=TRUE)) %dopar% {
                        set.seed(87932+i)
                        logmeanexp(
                          replicate(
                            ms_Nreps_eval[run_level],
                            logLik(pfilter(ms.filt,params=coef(if.box[[i]]),Np=ms_Np[run_level]))
                          ), 
                          se=TRUE)
                      }
  })
},seed=290860873,kind="L'Ecuyer")

r.box = data.frame(logLik=L.box[,1],logLik_se=L.box[,2],t(sapply(if.box,coef)))
summary(r.box$logLik,digits=5)

pairs(~logLik+log(sigma_nu)+mu_h+phi+sigma_eta+H_0,data=subset(r.box,logLik>max(logLik)-10))

# Look at diagnostic plots
plot(if.box)

################################################################################################
# Simulate from stochastic volatility model
################################################################################################
sim1.sim <- pomp(ms.filt, 
               statenames=ms_statenames,
               paramnames=ms_paramnames,
               covarnames=ms_covarnames,
               rprocess=discrete.time.sim(step.fun=Csnippet(ms_rproc.sim),delta.t=1)
)

params_test <- c(
     sigma_nu = .06,  
     mu_h = 6.27,  	 
     phi = .257,	 
     sigma_eta = .035,
     G_0 = 9.76,
     H_0=1.1
  )

sim1.sim = simulate(sim1.sim,seed=1,params=params_test, as.data.frame = T, include.data = T, nsim = 3)

dat_sim = sim1.sim

dat_sim$Year = dat_year$Year[2:nrow(dat_year)]

levels(dat_sim$sim) = c('Data', 'Simulation 1', 'Simulation 2', 'Simulation 3')

gg_sim = ggplot(data=dat_sim,mapping=aes(x=Year,y=y))+
  geom_line()+
  ggtitle('Stochastic Volatility Simulations')+
  ylab('Mass Shooting Fatalities')+
  xlab('Year')+
  facet_wrap(~sim)
formatGG(gg_sim)


conv.rec(if.box)
################################################################################################
# Simulate from stochastic volatility model
################################################################################################
fit.garch.benchmark <- garch(dat_demeaned,grad = "numerical", trace = FALSE)
L.garch.benchmark <- logLik(fit.garch.benchmark)
AIC(L.garch.benchmark)

# AIC for previous stochastic volatility model:
2*6 - 2*median(r.box$logLik)







