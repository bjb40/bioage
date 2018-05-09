---
author: Bryce Bartlett
title: R package to estimate biological age training parameters using the Klemera-Doubal algorithm.
---

# Introduction

## Overview

**What does it do?** Estimate "biological age" from a range of biomarkers using the Klemera Doubal algorithm (2006).

**Why?** Estimating biological age gives better leverage on aging, senescence, and disease process than chronological age alone.

## Installation

This package is still in development, but you can install and use it from github using the R library devtools. Here is the code-block:

```
install.packages('devtools')
library(devtools)
install_github('bjb40/bioage')
```


## Example



```
#Train biological age parameters
train = kdm_calc(nhanes,agevar='age',
  biomarkers=c('sysbp','totchol','bun','cmv','mcv'))

#Use training data to calculate out-of-sample biological ages 
biocalc = kdm_calc(data,agevar='age',
  biomarkers=c('sysbp','totchol','bun','cmv','mcv'),
  fit=train$fit)

#combine biological ages and training data
data$bioage = extract_data(biocalc)[,'bioage']

```



Description of Algorithm:

Klemera P, Doubal S. 2006. A new approach to the concept and computation of biological age. *Mechanisms of Ageing and Development.* 127(3):240-48

