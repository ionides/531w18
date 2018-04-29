#' ---
#' title: "Analyzing the Volatility of Bitcoin Market Price"
#' author: "Jun Luo"
#' date: "4/23/2018"
#' output: html_document
#' ---
#' 
## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
#' 
#' ### 1. Introduction
#' Cryptocurrency is without any doubt one of the most discussed techniques in the recent years. As the very first cryptocurrency, Bitcoin has been in the center of discussion and the market value of Bitcoin went from merely \$0.3 in the early 2011 to almost \$20000 at the end of 2017, which is illustrated in the left figure below.     
#' 
#' Lots of faith and money have been put into the market. Therefore, it is always relevant to study the market price and look for some insights, which might help in making an investment decision. One type of study on market prices focuses on the correlation between volatility and return. For example, one empirical observation tells us that negative shocks to a stockmarket index are associated with a subsequent increase in volatility, which may also be referred to as financial leverage [1].
#' 
#' In this project, I retrieved the market history for Bitcoin from 2011-01-01 to 2018-04-22 [2]. The plots below give us a basic idea how the market price for Bitcoin is changing from 2011 to 2018. The figure to the left is the plot for original price, the figure in the middle plots the market price in the log space, and the figure to the right illustrates the fluctuations of the demeaned return.
#' 
## ----echo=FALSE----------------------------------------------------------
data = read.csv("https://raw.githubusercontent.com/joebeav/531midterm/master/market_price_1.csv")
date <- data$Date
price <- data$Price
return <- diff(log(price))
demeaned <- return - mean(return)

#' 
#' ### 2. Garch Model
#' #### 2.1. Model Definition
#' We start with fitting a Garch model to our data.    
#' Assuming $Y_n$ is the return and $V_n$ is the volatility at time n, a GARCH(p,q) is defined as the follows:
#' $$Y_n = \epsilon_n\sqrt{V_n},$$
#' where
#' $$V_n = \alpha_0 + \sum_{j=1}^{p}{\alpha_jY_{n-j}^2} + \sum_{k=1}^{q}{\beta_k}V_{n-k}$$
#' and $\epsilon$ denotes white noise.
#' 
#' #### 2.2. Check Model Assumptions
#' Given the model definition, one could expect the volatility at time n to be correlated with previous returns and volatilities. Meanwhile, one could also expect the return at time n to be less correlated with previous returns as white noise is a major component for the return at time n. 
#' 
#' The auto-correlation for return and returns-squared are shown below. It is demonstrated that the auto-correlation for the original returns is weak and the auto-correlation for returns-squared is stronger, though the difference is not very large. These plots do comply to the assumption entailed by the model definition.
#' 
## ----echo=FALSE----------------------------------------------------------
par(mfrow = c(1,2))
acf(return)
acf(return^2)

#' 
#' #### 2.3. Select Model Based on AIC
## ----echo=TRUE, message=FALSE, warning=FALSE-----------------------------
require(tseries)
require(knitr)
Table_For_GARCH_AIC <- function(data,P,Q){
  table <- matrix(NA,(P),(Q))
  for(p in 1:P) {
    for(q in 1:Q) {
      temp.fit = garch(x = data, order = c(p,q), grad = "numerical", trace = FALSE)
      table[p,q] <- 2*length(temp.fit$coef) - 2*as.numeric(logLik(temp.fit))
    }
  }
  dimnames(table) <- list(paste("<b> p",1:P, "</b>", sep=""),paste("q",1:Q,sep=""))
  table
}
aic_table <- Table_For_GARCH_AIC(demeaned,6,6)
kable(aic_table,digits=2)


#' 
#' ### 3. POMP Model
#' Based on the work of Carles Bret?? [3], the model is defined as follows. 
#' $$Y_n=\mathrm{exp}(H_n/2)\epsilon_n,$$
#' $$H_n=\mu_h(1-\phi)+\phi H_{n-1}+\beta_{n-1} \mathrm{exp(-H_{n-1}/2)+\omega_n},$$
#' $$G_n=G_{n-1}+v_n$$
#' where $\beta_n=Y_{n}\sigma_{\eta}\sqrt{1-\phi^2}$, $\epsilon_{1:N}$, $v_{1:N}$ and $\omega_{1:N}$ are Gaussian white noise, and $R_n = \frac{\mathrm{exp}{(2G_n)}-1}{\mathrm{exp}{(2G_n)}+1}$, and {$G_n$} is Gaussian random walk.
#' 
#' In this model, $H_n$ is the log volatility, which is $H_n=log(\sigma^2_n)=2log(\sigma_n)$.
#' 
#' Below is how the model is built based on the instructions from class slides [4]. For computational convenience, we transformed the parameters to the whole real time scale using logit function, and then transform them back using expit function.
#' 
## ----echo=TRUE, warning=FALSE, message=FALSE-----------------------------
require(pomp)
btc_statenames <- c("H","G","Y_state")
btc_rp_names <- c("sigma_nu","mu_h","phi","sigma_eta")
btc_ivp_names <- c("G_0","H_0")
btc_paramnames <- c(btc_rp_names,btc_ivp_names)
btc_covarnames <- "covaryt"

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
btc_rproc.sim <- paste(rproc1,rproc2.sim)
btc_rproc.filt <- paste(rproc1,rproc2.filt)

btc_initializer <- "
  G = G_0;
  H = H_0;
  Y_state = rnorm( 0,exp(H/2) );
"
btc_rmeasure <- "
   y=Y_state;
"

btc_dmeasure <- "
   lik=dnorm(y,0,exp(H/2),give_log);
"

btc_toEstimationScale <- "
  Tsigma_eta = log(sigma_eta);
  Tsigma_nu = log(sigma_nu);
  Tphi = logit(phi);
"

btc_fromEstimationScale <- "
  Tsigma_eta = exp(sigma_eta);
  Tsigma_nu = exp(sigma_nu);
  Tphi = expit(phi);
"

btc.filt <- pomp(data=data.frame(y=demeaned,
                     time=1:length(demeaned)),
              statenames=btc_statenames,
              paramnames=btc_paramnames,
              covarnames=btc_covarnames,
              times="time",
              t0=0,
              covar=data.frame(covaryt=c(0,demeaned),
                     time=0:length(demeaned)),
              tcovar="time",
              rmeasure=Csnippet(btc_rmeasure),
              dmeasure=Csnippet(btc_dmeasure),
              rprocess=discrete.time.sim(step.fun=Csnippet(btc_rproc.filt),delta.t=1),
              initializer=Csnippet(btc_initializer),
              toEstimationScale=Csnippet(btc_toEstimationScale), 
              fromEstimationScale=Csnippet(btc_fromEstimationScale)
)

#' 
#' To fit the model, we use IF2 algorithm of Ionides et al. (2015) [5], implemented by mif2. We defined two different run_levels. run_level = 1 is for debugging, run_level = 2 and 3 are for finer simulations.
#' 
## ----echo=TRUE, message=FALSE--------------------------------------------
require(doParallel)
registerDoParallel()

run_level <- 3
btc_Np <-          c(100,1e3,1e4)
btc_Nmif <-        c(10,50,500)
btc_Nreps_eval <-  c(4,10,20)
btc_Nreps_local <- c(10,20,20)
btc_Nreps_global <-c(10,20,100)

expit<-function(real){1/(1+exp(-real))}
logit<-function(p.arg){log(p.arg/(1-p.arg))}

params_test <- c(
     sigma_nu = exp(-4.5),  
     mu_h = -0.25,       
     phi = expit(4),     
     sigma_eta = exp(-0.07),
     G_0 = 0,
     H_0=0
  )

btc_rw.sd_rp <- 0.02
btc_rw.sd_ivp <- 0.1
btc_cooling.fraction.50 <- 0.5

stew("mif1.rda",{
   t.if1 <- system.time({
   if1 <- foreach(i=1:btc_Nreps_local[run_level],
                  .packages='pomp', .combine=c,
                  .options.multicore=list(set.seed=TRUE)) %dopar% try(
                    mif2(btc.filt,
                         start=params_test,
                         Np=btc_Np[run_level],
                         Nmif=btc_Nmif[run_level],
                         cooling.type="geometric",
                         cooling.fraction.50=btc_cooling.fraction.50,
                         transform=TRUE,
                         rw.sd = rw.sd(
                            sigma_nu  = btc_rw.sd_rp,
                            mu_h      = btc_rw.sd_rp,
                            phi       = btc_rw.sd_rp,
                            sigma_eta = btc_rw.sd_rp,
                            G_0       = ivp(btc_rw.sd_ivp),
                            H_0       = ivp(btc_rw.sd_ivp)
                         )
                    )
                  )
    
    L.if1 <- foreach(i=1:btc_Nreps_local[run_level],.packages='pomp',
                      .combine=rbind,.options.multicore=list(set.seed=TRUE)) %dopar% 
                      {
                        logmeanexp(
                          replicate(btc_Nreps_eval[run_level],
                                    logLik(pfilter(btc.filt,params=coef(if1[[i]]),Np=btc_Np[run_level]))
                          ),
                          se=TRUE)
                      }
  })
},seed=318817883,kind="L'Ecuyer")