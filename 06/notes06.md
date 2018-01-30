---
title: "6. Extending the ARMA model: Seasonality and trend"
author: "Edward Ionides"
date: "1/30/2018"
output:
  html_document:
    theme: flatly
    toc: yes
    toc_depth: 2
    number_sections: true
    pandoc_args: [
      "--number-offset=6"
    ]
csl: ecology.csl
---


\newcommand\prob{\mathbb{P}}
\newcommand\E{\mathbb{E}}
\newcommand\var{\mathrm{Var}}
\newcommand\cov{\mathrm{Cov}}
\newcommand\loglik{\ell}
\newcommand\R{\mathbb{R}}
\newcommand\data[1]{#1^*}
\newcommand\params{\, ; \,}
\newcommand\transpose{\scriptsize{T}}
\newcommand\eqspace{\quad\quad\quad}
\newcommand\lik{\mathscr{L}}
\newcommand\loglik{\ell}
\newcommand\profileloglik[1]{\ell^\mathrm{profile}_#1}
\newcommand\ar{\phi}
\newcommand\ma{\psi}
\newcommand\AR{\Phi}
\newcommand\MA{\Psi}

Licensed under the Creative Commons attribution-noncommercial license, http://creativecommons.org/licenses/by-nc/3.0/.
Please share and remix noncommercially, mentioning its origin.  
![CC-BY_NC](cc-by-nc.png)




-------------------

------------------

<big><big><big>Objectives</big></big></big>

* Monthly time series often exhibit seasonal variation. January data are similar to observations at a different January, etc.

* Many time series exhibit a trend.

* We wish to extend the theoretical and practical elegance of the ARMA framework to cover these situations.

<br>

----------------------

---------------

## Seasonal autoregressive moving average (SARMA) models

* A general SARMA$(p,q)\times(P,Q)_{12}$ model for monthly data is
<br>
<br>
[S1] $\eqspace \ar(B)\AR(B^{12}) (Y_n-\mu) = \ma(B)\MA(B^{12}) \epsilon_n$,
<br>
<br>
where $\{\epsilon_n\}$ is a white noise process and
$$\begin{eqnarray}
\mu &=& \E[Y_n]
\\
\ar(x)&=&1-\ar_1 x-\dots -\ar_px^p,
\\ 
\ma(x)&=&1+\ma_1 x+\dots +\ma_qx^q, 
\\
\AR(x)&=&1-\AR_1 x-\dots -\AR_px^P,
\\ 
\MA(x)&=&1+\MA_1 x+\dots +\MA_qx^Q.
\end{eqnarray}$$

* We see that a SARMA model is a special case of an ARMA model, where the AR and MA polynomials are factored into a **monthly** polynomial in $B$ and an **annual** polynomial in $B^{12}$. The annual polynomial is also called the **seasonal** polynomial.

* Thus, everything we learned about ARMA models (including assessing causality, invertibility and reducibility) also applies to SARMA. 

* One could write a SARMA model for some **period** other than 12. For example, a  SARMA$(p,q)\times(P,Q)_{4}$ model could be appropriate for quarterly data. In principle, a SARMA$(p,q)\times(P,Q)_{52}$ model could be appropriate for weekly data, though in practice ARMA and SARMA may not work so well for higher frequency data. 

* Consider the following two models:
<br>
<br>
[S2] $\eqspace Y_n = 0.5 Y_{n-1} + 0.25 Y_{n-12} + \epsilon_n$,
<br>
<br>
[S3] $\eqspace Y_n = 0.5 Y_{n-1} + 0.25 Y_{n-12} - 0.125 Y_{n-13} + \epsilon_n$,
<br>
<br>

---------

--------

### Question: Which of [S2] and/or [S3] is a SARMA model? 

<br>

--------

-------

### Question: Why do we assume a multiplicative structure in [S1]?

* What theoretical and practical advantages (or disadvantages) arise from requiring that an ARMA model for seasonal behavior has polynomials that can be factored as a product of a monthly polynomial and an annual polynomial?

<br>

----------

---------

### Fitting a SARMA model

* Let's do this for the full, monthly, version of the Lake Huron data described in [Section 5.5](http://ionides.github.io/531w18/05/notes05.html#implementing-likelihood-based-inference-for-arma-models-in-r).

* First, we'll revisit reading in the data.















