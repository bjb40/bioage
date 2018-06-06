#functions for building and manipulationg KDM algorithm

hd
#builds a fomula on the fly
form = function(y,x){
  return(as.formula(paste0(y,'~',x)))
}

#helper function to get effects from linear models on biomarkers
get_effs = function(mod){
  #input model object
  #output is list of
  #intercept, rmse, rsqaure

  m = summary(mod)
  res = data.frame(q=coef(mod)['(Intercept)'],k=coef(mod)[2],s=NA,r=NA)
  if('glm' %in% attr(mod,'class')){
    #sd of residuals is the RMSE
    res$s = sd(mod$residuals)
    #calculate s=pseudo r-squared
    res$r = 1-(mod$deviance/mod$null.deviance)
  }else{
    res$s = m$sigma
    res$r = m$r.squared
  }

  return(res)
}

####extend caclculations below to dataset...
weighted.var = function(vector,weights){
  #https://stats.stackexchange.com/questions/51442/weighted-variance-one-more-time
  #frequency weights... (vs. reliability weights?)
  v = vector[!is.na(vector)] #remove missing
  w = weights[!is.na(vector)]
  mu = weighted.mean(v,w)
  var = sum(w*(v-c(mu)))/(sum(w)-1)

  return(var)

}

#' Calculate biological ages in a dataset.
#'
#' @param data The dataset for calculating biological age.
#' @param agevar A character vector (length=1) indicating the name of the varialbe for age.
#' @param biomarkers A character vector indicating the names of the variables for the biomarkers to use in calculating biological age.
#' @param filter a list with biomarker names that identifies any restrictions in training data. See vignette or data description for example of use.
#' @param fit An S3 object for model fit. If the value is NULL, then the parameters to use for training biological age are calculated.
#' @param link "linear" is default and based on the original KDM algorithm; experimental use of log-linear link (use "log") is available for advanced users.
#' @param s_ba2 A particular fit parameter. Advanced users can modify this parameter to control the variance of biological age. If left NULL, defaults are used.
#' @param weightvar A character vector indicating survey weights. If supplied, a weighted regression is conducted. If not, weights are not used.
#' @param controls A character vector indicating control variables (if any) to be used for calculating biological age.
#' @return An object of class 'kdm'. This object is a list with two elements (data and fit),
#' and two methods (extract_data and extract_fit).
#' @examples
#' #(not run)
#' #Train biological age parameters
#' train = kdm_calc(nhanes,agevar='age',
#'   biomarkers=c('sysbp','totchol','bun','cmv','mcv'))
#'
#' #Use training data to calculate out-of-sample biological ages
#' biocalc = kdm_calc(data,agevar='age',
#'   biomarkers=c('sysbp','totchol','bun','cmv','mcv'),
#'   fit=train$fit)
#'
#' #combine biological ages calculated using training parameters
#' data$bioage = extract_data(biocalc)[,'bioage']
kdm_calc = function(data,
                    agevar,
                    biomarkers,
                    fit=NULL,
                    filter=NULL,
                    link='linear',
                    s_ba2=NULL,
                    weightvar=NULL,
                    controls=NULL){

  require(survey)
  require(dplyr)
  require(reshape2)

  #fit is a list an object from previously trained data used to fit
  #if null, it trains; otherwise calcs
  #link is 'log' or 'linear' ---- it specifies wehter to fit a log or linear fvalue
  #controls is vector of column names for controls

  train=data; rm(data)
  bm = biomarkers

  #identifying filter
  #need to add some tests to filter to capture errors for mis-specification
  if(is.null(filter)){
    filter = vector("list",length(bm))
    names(filter) = bm
  } else {
    leftout = bm[!bm %in% names(filter)]
    nl = vector('list',length(leftout))
    names(nl) = leftout
    filter = c(filter,nl); #rm(leftout,nl)
  }

  if(is.null(weightvar)){
    train$w = rep(1,nrow(train))
  } else{
    #cat('Using Weights ',weightvar,'.',sep='')
    train$w = unlist(train[,weightvar])
  }

  #identify filter conditions

  train2 = train
  for(l in seq_along(filter)){

    if(is.null(filter[[l]])){next} else{
      lim=with(train,eval(parse(text=filter[[l]])))
      train2[lim,names(filter)[l]] = NA
    }

  }

  design=svydesign(id=~1,weights=~w,data=train2)


  #pull the calculations out to make fit more general??

  iv = paste(c(agevar,controls),collapse='+')

  if(is.null(fit)){

    if(link=='linear'){
      lm_age = lapply(bm,function(marker){
        subset = eval(filter,1)[[marker]]
        svyglm(form(marker,iv),
               design=design,
               family=gaussian())
        })
    } else if(link=='log'){
      lm_age = lapply(bm,function(marker) glm(form(marker,iv), data=train,
                                              family = quasipoisson(link=log)))
    }

    agev = do.call(rbind,lapply(lm_age,get_effs))
    rm(lm_age)
    #agev$bm = row.names(agev)
    agev$bm = bm

    agev = agev %>%
      mutate(r1=abs((k/s)*sqrt(r)),
             r2=abs(k/s),
             n2=(k/s)^2)

    age.range=range(train[,agevar],na.rm=TRUE)

    rchar = sum(agev$r1)/sum(agev$r2)

    s_r = ((1-(rchar^2))/(rchar^2))*(((age.range[2]-age.range[1])^2)/(12*nrow(agev)))

  } else{
    agev = fit$agev
    s_r = fit$s_r
  }
  #end train conditional

  n1=train[,bm]

  for(m in colnames(n1)){
    row = which(agev$bm == m)
    if(link=='linear'){
      obs = (train[,m] - agev[row,'q'])*(agev[row,'k']/(agev[row,'s']^2))
    } else if(link=='log') {
      obs = (log(train[,m])-agev[row,'q'])*(agev[row,'k']/(agev[row,'s']^2))
    }

    n1[,m]=(obs)
  }

  #allow for proration of biological ages
  #based on missing data -- fill this out in text
  ba.nmiss = apply(n1,1,function(x) sum(is.na(x)))
  ba.obs = length(bm) - ba.nmiss
  ba.e_n = rowSums(n1,na.rm=TRUE)
  ba.e_d = sum(agev$n2,na.rm=TRUE)

  train = train %>%
    mutate(ba.eo = ba.e_n/ba.e_d,
           ba.e = (ba.eo/(ba.obs))*length(bm))

  #should calculated weighted means to use wieghted regressions; or delete altogether

  train$ba.ca = unlist(train$ba.e - train[,agevar])
  #s2 = var(train$ba.ca,na.rm=TRUE)
  t1 = (train$ba.ca - mean(train$ba.ca,na.rm=TRUE))^2
  s2=mean(t1,na.rm=TRUE)


  nobs = sum(!is.na(train$ba.ca))

  #s2 = c(attr(svymean(~ba.ca,design=design,na.rm=TRUE),'var'))*c(nobs)
  #s2 = weighted.var(train$ba.ca,train$w) #calculate using the sampling weight

  #levine provides s_ba2
  if(is.null(s_ba2)){
    s_ba2 = s2-s_r
  } else{
    s_ba2 = s_ba2
  }


  #use prevous sba2 --- can add an if above...
#  if(!is.null(fit)){
#    cat('\n\nold sba2=',s_ba2,'\n')
#    s_ba2 = fit[['s_ba2']]
#    cat('\n new sba2=',s_ba2,'\n')
#  }

  #cat('sba2',s_ba2,'\n')

  train$agevar = unlist(train[,agevar])

  train$bioage = unlist(
    (ba.e_n + (train$agevar/c(s_ba2)))/(ba.e_d+(1/c(s_ba2)))
  )

  #set to missing if more than 2 missing
  train$bioage = ifelse(ba.nmiss>2,NA,train$bioage)

  #  train = train %>%
  #    mutate(bioage = ((ba.e_n)+(agevar/s.ba2))/((ba.e_d)+ (1/s.ba2))) %>%
  #    select(-agevar)

  fit = list(agev=agev,s_r=s_r,link=link,s_ba2=s_ba2,s2=s2,nobs=nobs)

  #train$ba.e_n = ba.e_n
  #train$ba.e_d = ba.e_d

  #rename to bio__
  colnames(train)[which(colnames(train)=='bioage')] = paste0('bio',agevar)
  train$baaccel_diff = unlist(train[,paste0('bio',agevar)] - train[,agevar])
  train$baaccel_resid = residuals(
    lm(form(paste0('bio',agevar),agevar),
       data=train,
       na.action='na.exclude')
    )




  #train$s_r = s_r
  #train$s2 = s2
  #train$s_ba2 = s_ba2

  kdm = list(data = train, fit = fit)
  class(kdm) = append(class(kdm),'kdm')

  return(kdm)

}

#' Extracts and summarizes estimates for fitting biological age from a kdm object.
#'
#' @param kdmobj A kdm object estimated using the function kdm_calc.
#' @return A dataframe with trained parameters.
#' @examples
#' #(not run)
#' #Train biological age parameters
#' train = kdm_calc(nhanes,agevar='age',
#'   biomarkers=c('sysbp','totchol','bun','cmv','mcv'))
#'
#' myfit = extract_fit(train)
extract_fit = function(kdmobj){
  UseMethod('extract_fit',kdmobj)
}

#' @export
extract_fit.default = function(kdmobj){
  cat('\nError: Must use kdm object generated from kdm_calc.')
}

#' @export
extract_fit.kdm = function(kdmobj){

  fit = kdmobj$fit$agev
  fit$s_r = kdmobj$fit$s_r
  fit$s_ba2 = kdmobj$fit$s_ba2
  fit$link = kdmobj$fit$link
  fit$s2 = kdmobj$fit$s2
  fit$nobs = kdmobj$fit$nobs

  return(fit)
}


#' Extracts dataframe with age, controls, biomarkers, and biological age calculation from trianed data (using kdm_calc function).
#'
#' @param kdmobj A kdm object estimated using the function kdm_calc.
#' @return A dataframe with parameters and biological age.
#' @examples
#' #(not run)
#' #Train biological age parameters
#' train = kdm_calc(nhanes,agevar='age',
#'   biomarkers=c('sysbp','totchol','bun','cmv','mcv'))
#'
#' newdata = extract_data(train)
extract_data = function(kdmobj){
  UseMethod('extract_data',kdmobj)
}

extract_data.default = function(kdmobj){
  cat('\nError: Must use kdm object generated from kdm_calc.')
}

extract_data.hd = function(kdmobj){
  #kdmobj is  trained set
  #returns dataframe and prints to csv

  return(kdmobj$data)

}


extract_data.kdm = function(kdmobj){
  #kdmobj is  trained set
  #returns dataframe and prints to csv

  return(kdmobj$data)

}
