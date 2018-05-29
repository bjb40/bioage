###
#these functions calculate homeostatic disregulation

#' Calculate homeostatic disregulation in a dataset.
#'
#' @param data The dataset for calculating homeostatic disregulation.
#' @param biomarkers A character vector indicating the names of the variables for the biomarkers to use in calculating homeostatic disregulation.
#' @param filter a list with biomarker names that identifies any restrictions in training data. These are generally a young, "healthy" reference group, i.e. a group within clinical guidelines.
#' @param fit An S3 object for model fit. If the value is NULL, then the parameters to use for training homeostatic disregulation are calculated.
#' @return An object of class 'hd'. This object is a list with two elements (data and fit),
#' and two methods (extract_data and extract_fit).
#' @examples
#' #(not run)
#' #Train homeostatic disregulation
#' train = hd(data=nhanes,
#'            biomarkers=c('sysbp','totchol','bun','cmv','mcv'),
#'            fit=list(sysbp='',
#'                     totchol='',
#'                     bun='',
#'                     cmv='',
#'                     mcv=''
#'            ))
#'
#' #Use training data to calculate out-of-sample biological ages
#' valudate = kdm_calc(data,fit=train$fit,
#'   biomarkers=c('sysbp','totchol','bun','cmv','mcv'))
#'
#' #combine biological ages calculated using training parameters
#' data$bioage = extract_data(biocalc)[,'bioage']
hd = function(data,biomarkers,fit=NULL,filter=NULL){

  #fit inlcudes vectors of means & the covariance matrix & filter
  #if null, should return observed range (min and max)

  if(is.null(fit)){
  #step 1: compose dataset based on filters
  #anyone who falls outside of any biomarkers is eliminated from analysis
  train = data

  for(l in seq_along(filter)){
    #print(filter[[l]])

    if(is.null(filter[[l]])){next} else{
      lim=with(train,eval(parse(text=filter[[l]])))
      train = train[lim & !is.na(lim),]
    }
  }

  cat(nrow(train), 'observations remaining after exclusions.\n')

  #step 2: standardize mean and covariance

  #export mean and covariance for each biomarker
    mcov = list(means=as.data.frame(apply(train[,biomarkers],2,mean,na.rm=TRUE)),
                cov = cov(train[,biomarkers],use='pairwise.complete.obs'))

    fit = list(mcov=mcov,nobs=nrow(train),filter=filter)

  }

  mcov = fit$mcov

  #step 3: transform projection dataset into z-scores based on standardized
  #        training data

  centered = data[,biomarkers]

  for(bm in biomarkers){
    #print(bm)
    #subtract mean
    centered[,bm] = centered[,bm]-mcov$means[row.names(mcov$means)==bm,]
    #divide by sd
    centered[,bm] = centered[,bm]/sqrt(diag(mcov$cov)[bm])
    #data[,biomarkers]

  }

  #step 4: calculate malhanobis distance, and return as "raw_dist"
  data$raw_dist = mahalanobis(centered,center=FALSE,cov=cov2cor(mcov$cov))

  #step 5: log distance & standardize
  data$hd = scale(log(data$raw_dist))


  hd = list(data=data,fit=fit)
  class(hd) = append(class(hd),'hd')

  return(hd)
}
