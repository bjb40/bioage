#functions for building and manipulationg KDM algorithm


#builds a fomula on the fly
form = function(y,x){
  return(as.formula(paste0(y,'~',x)))
}

#run linear regression on biomarkers; returns lists

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


kdm_calc = function(data,
                    agevar,
                    biomarkers,
                    fit=NULL,
                    link='linear',
                    s_ba2=NULL,
                    weightvar=NULL,
                    controls=NULL){

  require(survey)

  #fit is a list an object from previously trained data used to fit
  #if null, it trains; otherwise calcs
  #link is 'log' or 'linear' ---- it specifies wehter to fit a log or linear fvalue
  #controls is vector of column names for controls

  train=data; rm(data)
  bm = biomarkers

  if(is.null(weightvar)){
    train$w = rep(1,nrow(train))
  } else{
    #cat('Using Weights ',weightvar,'.',sep='')
    train$w = unlist(train[,weightvar])
  }

  design=svydesign(id=~1,weights=~w,data=train)


  #pull the calculations out to make fit more general??

  iv = paste(c(agevar,controls),collapse='+')

  if(is.null(fit)){

    if(link=='linear'){
      lm_age = lapply(bm,function(marker)
        svyglm(form(marker,iv),design=design,family=gaussian()))
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
  ba.nmiss = apply(n1,1,function(x) sum(is.na(x)))
  ba.obs = length(bm) - ba.nmiss
  ba.e_n = rowSums(n1,na.rm=TRUE)
  ba.e_d = sum(agev$n2,na.rm=TRUE)

  #prorate here or not...?

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
  if(!is.null(fit)){
    cat('\n\nold sba2=',s_ba2,'\n')
    s_ba2 = fit[['s_ba2']]
    cat('\n new sba2=',s_ba2,'\n')
  }

  #cat('sba2',s_ba2,'\n')

  train$agevar = unlist(train[,agevar])

  train$bioage = unlist(
    (ba.e_n + (train$agevar/c(s_ba2)))/(ba.e_d+(1/c(s_ba2)))
  )

  #  train = train %>%
  #    mutate(bioage = ((ba.e_n)+(agevar/s.ba2))/((ba.e_d)+ (1/s.ba2))) %>%
  #    select(-agevar)

  fit = list(agev=agev,s_r=s_r,link=link,s_ba2=s_ba2,s2=s2,nobs=nobs)

  train$ba.e_n = ba.e_n
  train$ba.e_d = ba.e_d

  #rename to bio__
  colnames(train)[which(colnames(train)=='bioage')] = paste0('bio',agevar)

  train$s_r = s_r
  train$s2 = s2
  train$s_ba2 = s_ba2

  return(
    list(data = train, fit = fit)
  )

}

#extraction functions for dataframe of fit object

extract_fit = function(kdmobj){
  #kdmobj is  trained set
  #returns dataframe

  fit = kdmobj$fit$agev
  fit$s_r = kdmobj$fit$s_r
  fit$s_ba2 = kdmobj$fit$s_ba2
  fit$link = kdmobj$fit$link
  fit$s2 = kdmobj$fit$s2
  fit$nobs = kdmobj$fit$nobs

  return(fit)
}

extract_data = function(kdmobj){
  #kdmobj is  trained set
  #returns dataframe and prints to csv

  return(kdmobj$data)

}
