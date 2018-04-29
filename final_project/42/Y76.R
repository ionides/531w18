getwd()
    set.seed(1234567890)
    library(ggplot2)
    theme_set(theme_bw())
    library(plyr)
    library(reshape2)
    library(magrittr)
    library(foreach)
    library(doMC)
    library(pomp)
    library(doParallel)
    
    
    SARS_data <- read.table("SARS_Beijing.csv", sep=',', header = T)
    colnames(SARS_data) <- c('B','C','K','day')
    pdf(file ='1_new.pdf')
    plot(x= SARS_data$day, y=SARS_data$B , xlab='day', ylab='Number of Infection')
    plot(x= SARS_data$day, y=SARS_data$C , xlab='day', ylab='Number of Recovered')
    plot(x= SARS_data$day, y=SARS_data$K , xlab='day', ylab='Number of Death')
    dev.off()
####################
    SARS_b_data <- SARS_data[1:60,]
    
    sir_step <- Csnippet("
                         double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
                         double dN_IR = rbinom(I,1-exp(-gamma1*dt));
                         double dN_ID = rbinom(I,1-exp(-gamma2*dt));
                         S -= dN_SI;
                         I += dN_SI - dN_IR - dN_ID;
                         R += dN_IR;
                         D += dN_ID;
                         H += dN_IR;
                         ")
    
    sir_init <- Csnippet("
                         S = nearbyint(N)-1;
                         I = 93;
                         R = 43;
                         D = 25;
                         H = 43;
                         ")
    
    dmeas <- Csnippet("lik = dbinom(B,H,rho,give_log);")
    rmeas <- Csnippet("B = rbinom(H,rho);")
    
    pomp(SARS_b_data,times="day",t0=0,
         rprocess=euler.sim(sir_step,delta.t=1/5),
         initializer=sir_init,rmeasure=rmeas,dmeasure=dmeas,
         zeronames="H",statenames=c("H","S","I","R","D"),
         paramnames=c("Beta","gamma1","gamma2","rho","N")) -> sir
    
    
    sims <- simulate(sir,params=c(Beta=1.2,gamma1=0.8, gamma2 = 0.2, rho=0.5,N=200000),nsim=100,
                     as.data.frame=TRUE,include.data=TRUE)
    pdf(file='2_new.pdf')
    ggplot(sims,mapping=aes(x=time,y=B,group=sim,color=sim=="data"))+
      geom_line()+guides(color=FALSE)
    dev.off()
    
####################
    SARS_b_statenames <- c("S","I","R1","R2","R3")
    SARS_b_paramnames <- c("Beta","mu_I","rho","mu_R1","mu_R2","mu_R3")
    (SARS_b_obsnames <- colnames(SARS_b_data)[1:3])
    
    SARS_b_dmeasure <- "
    lik = dpois(B,rho*R1+1e-6,give_log);
    "
    
    SARS_b_rmeasure <- "
    B = rpois(rho*R1+1e-6);
    C = rpois(rho*R2);
    K = rpois(rho*R3);
    "
    
    SARS_b_rprocess <- "
    double t1 = rbinom(S,1-exp(-Beta*I*dt));
    double t2 = rbinom(I,1-exp(-dt*mu_I));
    double t3 = rbinom(R1,1-exp(-dt*mu_R1));
    double t4 = rbinom(R2,1-exp(-dt*mu_R2));
    double t5 = rbinom(R3,1-exp(-dt*mu_R3));
    S -= t1;
    I += t1 - t2;
    R1 += t2 - t3;
    R2 += t3 - t4;
    R3 += t4 - t5;
    "
    
    SARS_b_fromEstimationScale <- "
    TBeta = exp(Beta);
    Tmu_I = exp(mu_I);
    Trho = expit(rho);
    "
    
    SARS_b_toEstimationScale <- "
    TBeta = log(Beta);
    Tmu_I = log(mu_I);
    Trho = logit(rho);
    "
    
    SARS_b_initializer <- "
    S=20000;
    I=1;
    R1=0;
    R2=0;
    R3=0;
    "
    
    stopifnot(packageVersion("pomp")>="0.75-1")
    SARS_b2 <- pomp(
      data=SARS_b_data,
      times="day",
      t0=0,
      rprocess=euler.sim(
        step.fun=Csnippet(SARS_b_rprocess),
        delta.t=1/12
      ),
      rmeasure=Csnippet(SARS_b_rmeasure),
      dmeasure=Csnippet(SARS_b_dmeasure),
      fromEstimationScale=Csnippet(SARS_b_fromEstimationScale),
      toEstimationScale=Csnippet(SARS_b_toEstimationScale),
      obsnames = SARS_b_obsnames,
      statenames=SARS_b_statenames,
      paramnames=SARS_b_paramnames,
      initializer=Csnippet(SARS_b_initializer)
    )
    pdf(file='3_new.pdf')
    plot(SARS_b2)
    dev.off()
    
    sims <- simulate(SARS_b2,params=c(Beta=2, mu_I=0.04, rho =0.3, mu_R1=0.05, mu_R2 =0.1, mu_R3= 0.1, n= 1000000),nsim=10,
                     as.data.frame=TRUE,include.data=TRUE)
    
    
    pdf(file='2_new.pdf') 
    ggplot(sims,mapping=aes(x=time,y=B,group=sim,color=sim=="data"))+
      geom_line()+guides(color=FALSE)
   dev.off()
    
###################

    run_level <- 2
    switch(run_level,
           {SARS_b_Np=100; SARS_b_Nmif=10; SARS_b_Neval=10; SARS_b_Nglobal=10; SARS_b_Nlocal=10}, 
           {SARS_b_Np=20000; SARS_b_Nmif=100; SARS_b_Neval=10; SARS_b_Nglobal=10; SARS_b_Nlocal=10}, 
           {SARS_b_Np=60000; SARS_b_Nmif=300; SARS_b_Neval=10; SARS_b_Nglobal=100; SARS_b_Nlocal=20}
    )
    
    require(doParallel)
    cores <- 80
    registerDoParallel(cores)
    mcopts <- list(set.seed=TRUE)
    
    SARS_b_rw.sd <- 0.02
    SARS_b_cooling.fraction.50 <- 0.5
    
    SARS_b_params <- data.matrix(read.table("mif_SARS_b_params.csv",row.names=NULL,header=TRUE,sep=''))
    SARS_b_mle <- SARS_b_params[which.max(SARS_b_params[,"logLik"]),][SARS_b_paramnames]
    SARS_b_fixed_params <- c(mu_R1=1/(sum(SARS_b_data$B)/3000),mu_R2=1/(sum(SARS_b_data$C)/3000), mu_R3=1/(sum(SARS_b_data$K)/3000))
    
    
    mcopts <- list(preschedule=FALSE,set.seed=TRUE)
    
    library(foreach)
    library(doMC)
    library(pomp)
    library(doParallel)
    registerDoMC(cores=80)  
    stew(file=sprintf("local_search-%d.rda",run_level),{
      
      t_local <- system.time({
        mifs_local <- foreach(i=1:SARS_b_Nlocal,.packages='pomp', .combine=c, .options.multicore=mcopts) %dopar%  {
          mif2(
            SARS_b2,
            start=SARS_b_mle,
            Np=SARS_b_Np,
            Nmif=SARS_b_Nmif,
            cooling.type="geometric",
            cooling.fraction.50=SARS_b_cooling.fraction.50,
            transform=TRUE,
            rw.sd=rw.sd(
              Beta=SARS_b_rw.sd,
              mu_I=SARS_b_rw.sd,
              rho=SARS_b_rw.sd
            )
          )
          
        }
      })
      
    },seed=1234567890,kind="L'Ecuyer")
    
    stew(file=sprintf("lik_local-%d.rda",run_level),{
      t_local_eval <- system.time({
        liks_local <- foreach(i=1:SARS_b_Nlocal,.packages='pomp',.combine=rbind) %dopar% {
          evals <- replicate(SARS_b_Neval, logLik(pfilter(SARS_b2,params=coef(mifs_local[[i]]),Np=SARS_b_Np)))
          logmeanexp(evals, se=TRUE)
        }
      })
    },seed=1234567890,kind="L'Ecuyer")
    
    results_local <- data.frame(logLik=liks_local[,1],logLik_se=liks_local[,2],t(sapply(mifs_local,coef)))
    
    summary(results_local$logLik,digits=10)
    
    local_loglik_max=max(results_local$logLik)
    results_local_max= subset(results_local,results_local$logLik==local_loglik_max)
    results_local_max
    
    results_local %>% 
      subset(logLik==max(logLik)) %>% unlist() -> mle
  
    
    pdf(file='4_new.pdf')
    pairs(~logLik+Beta+mu_I+rho,data=subset(results_local,logLik>max(logLik)-500))
    dev.off()

pdf(file='3.1_new.pdf')
SARS_b2 %>%
simulate(params=mle,nsim=50,as.data.frame=TRUE,include.data=TRUE) %>%
ggplot(mapping=aes(x=time,y=B,group=sim,alpha=(sim=="data")))+
scale_alpha_manual(name="",values=c(`TRUE`=1,`FALSE`=0.2),
labels=c(`FALSE`="simulation",`TRUE`="data"))+
geom_line()+
theme_bw()
dev.off()
####################
    SARS_b_box <- rbind(
      Beta=c(0.001,0.01),
      mu_I=c(0.5,10),
      rho = c(0.5,1)
    )
    
    stew(file=sprintf("box_eval-%d.rda",run_level),{
      
      t_global <- system.time({
        mifs_global <- foreach(i=1:SARS_b_Nglobal,.packages='pomp', .combine=c, .options.multicore=mcopts) %dopar%  mif2(
          mifs_local[[1]],
          start=c(apply(SARS_b_box,1,function(x)runif(1,x[1],x[2])),SARS_b_fixed_params)
        )
      })
    },seed=1234567,kind="L'Ecuyer")
    
    stew(file=sprintf("lik_global_eval-%d.rda",run_level),{
      t_global_eval <- system.time({
        liks_global <- foreach(i=1:SARS_b_Nglobal,.packages='pomp',.combine=rbind, .options.multicore=mcopts) %dopar% {
          evals <- replicate(SARS_b_Neval, logLik(pfilter(SARS_b2,params=coef(mifs_global[[i]]),Np=SARS_b_Np)))
          logmeanexp(evals, se=TRUE)
        }
      })
    },seed=98765432,kind="L'Ecuyer")
    
    results_global <- data.frame(logLik=liks_global[,1],logLik_se=liks_global[,2],t(sapply(mifs_global,coef)))
    
    global_loglik_max=max(results_global$logLik)
    results_global_max= subset(results_global,results_global$logLik==global_loglik_max)
    results_global_max
    
    results_global %>% 
      subset(logLik==max(logLik)) %>% unlist() -> mle
    
    pdf(file='5.1_new.pdf')
    SARS_b2 %>%
      simulate(params=mle,nsim=50,as.data.frame=TRUE,include.data=TRUE) %>%
      ggplot(mapping=aes(x=time,y=B,group=sim,alpha=(sim=="data")))+
      scale_alpha_manual(name="",values=c(`TRUE`=1,`FALSE`=0.2),
                         labels=c(`FALSE`="simulation",`TRUE`="data"))+
      geom_line()+
      theme_bw()
    dev.off()
    
    
    summary(results_global$logLik,digits=8)
    
    if (run_level>2) 
      write.table(rbind(results_local,results_global),
                  file="mif_SARS_b_params.csv",append=TRUE,col.names=FALSE,row.names=FALSE)
    
    
    pdf(file= '5_new.pdf')
    pairs(~logLik+Beta+mu_I+rho,data=subset(results_global,logLik>max(logLik)-2500))
    dev.off()
    

    
    pdf(file='6_new.pdf')
    plot(mifs_global)
    dev.off()
###################################    

    
########################  
    
    
    
    
    
