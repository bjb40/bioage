#functions and objects

#@@@@@@@@@@@@@@@@@@@@@
#Functions/objects
#@@@@@@@@@@@@@@@@@@@@@


rnd = function(db,rd){
  # rounds input to preserve leading zeros
  #
  # Args:
  #   db: an object with numeric types
  #   rd: length to round (including leading zeros, default=3)
  #
  # Returns:
  #   an object of of db translated to characters with leading zeros

  if(missing(rd)){rd=3}
  rdl=paste0('%.',rd,'f')
  return(sprintf(rdl,round(db,digits=rd)))
}

sig = function(pv){
  # returns stars based on pvalue
  #
  # Args:
  #   pv: a p-value
  #
  # Returns:
  #   a string with stars for values * <.05 **<.01 *** < .001
  s=' '
  if(length(pv)>0){
    if(pv<.001){s='***'} else if(pv<.01){s='**'} else if (pv<.05){s='*'} else if (pv<.1){s='+'}
  }
  return(s)

}


mnsd = function(x){
  #returns rounded mean and sd for vector
  return(
    c(rnd(mean(x,na.rm=TRUE)),
      paste0('(',rnd(sd(x,na.rm=TRUE)),')')
    )#end combine
  )#end return
}

