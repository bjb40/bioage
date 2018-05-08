---
author: Bryce Bartlett
title: R package for Stochastic Window Sampling (using Bayesian Model Averaging) model to solve the Age-Period-Cohort identification problem.
---

# Introduction

## Overview

**What does it do?** Estimate hundreds of models using an MCMC method. Get a single estimate using Bayesian model averaging.

**Why?** To solve perfect correlation between age, period, and cohort, but eliminate sensitivity to arbitrary modeling decisions. In other words, to get consistent, sensible estimates despite exact correlation and unidentifiability. While this package is directed to a specific problem, the approach shows promise for any problem with high multicollinearity.

## Background

There is a long-standing problem in demography known as the age-period-cohort identification problem. There is a fairly straightforward issue: some social characteristic, say, happiness, may be correlated with age. People generally get happier as they get older. This is the age effect. Happiness can also be associated with current events: everyone is less happy during a recession. This is a period effect. Finally, different birth cohorts (sometimes called generations) have different average levels of happiness. Baby boomers are relatively unhappy compared to others. This is a cohort effect. (Note that mathematically these effects are separable if at least one is nonlinear).

Unfortunately, our measures for age, period, and cohort are perfectly correlated. If you know any two, you also know the third. If you know how old someone is (age) and when they were born (cohort), then you know what year (period) it is (age+cohort=period). Because these measures are perfectly correlated, they are not identified, and inestimable. 

There are a number of proposals to address this problem, but this package deploys a novel Bayesian Model Averaging approach I developed called the stochastic window sampling method. (The manuscript is still in preparation, but will be linked here once complete). This works by estimating age, period, and cohort effects using piecewise constant functions for hundreds of models. We have known for a long time that piecewise constant functions ("window" or "block" models in the APC literature) solve the identifiability problem, but these solutions are known to be sensitive to exactly where the "window" breaks are located. The stochastic window sampler uses MCMC methods to select locations for the window breaks, solving the sensitivity issue, *i.e.* once converged, the stochastic window sampler provides a single, consistent, estimate. It also provides a sensible classification for which ages, periods, and cohorts are most alike or most different.

# Installation

This package is still in development, but you can install and use it from github using the R library devtools. Here is the code-block:

```
install.packages('devtools')
library(devtools)
install_github('bjb40/apcwin')
```


# Example

This following codeblock uses simulated data to illustrate the arguments. This returns graphical displays and identifies best-fitting models and model-breaks. 

```
###
#load test data
#Note, you need to name your variables
#as follows a = age, p = period, and c= cohort
data(apcsim)
apcsim$c = apcsim$p-apcsim$a

#this draws 2500 samples on 4 cores for 10000 model samples
testsamp = apcsamp(dat=apcsim,
                   dv='y1',
                   method='ml',
                   samples=2500,
                   cores=4)

#this summarizes the results of the apc sampler
summarize(testsamp)


#this draws a posterior effect sample
#it takes a "sample" object (calculated by apcsamp)
testeff = draw_effs(testsamp,
                    tol=0.001)

#this plots the results of the dimensions
plot(testeff,alpha=0.05)
```




