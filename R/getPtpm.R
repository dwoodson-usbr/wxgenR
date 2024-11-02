#' Precipitation transition probabilities
#'
#' Calculate seasonal precipitation transition probability matrix for each year.
#' Computes transition probabilities for wet/dry spells.
#' Calculates values for each year, then takes average of all years.
#'
#' @param dat.d Training data processed from prepData wrapper function.
#'
# @rawNamespace import(stats, except = filter)
#'
#'
#' @noRd
#'

"getPtpm" <- function(dat.d, traceThreshold){
  #calculate seasonal precipitation transition probability matrix for each year
  #
  yr.d = dat.d$year
  uyr = unique(yr.d)
  nyr = length(uyr)
  #precipitation data
  pcp.d = dat.d$prcp
  #transition probabilities for wet/dry spells
  #compute wet and dry transitions for each state (season) in a year, repeat for each year
  tpm.y2 = array(data = NA, dim = c(2,2,max(dat.d$states), length(uyr)))
  j=1
  for(j in 1:length(uyr)){
    k=1
    for(k in 1:max(dat.d$states)){
      if(sum(dat.d$states == k & yr.d == uyr[j]) == 0) next
      x = ts((pcp.d[dat.d$states == k & yr.d == uyr[j]] >= traceThreshold) + 0,1)
      tpm.tmp = transProbMatrix(x)
      tpm.y2[as.numeric(rownames(tpm.tmp))+1, as.numeric(colnames(tpm.tmp))+1,
             k,j] = tpm.tmp
    }#k
   }#j

  # collapse tpm.y2 over all years
  tpm.y = apply(tpm.y2, 1:3, mean, na.rm=T)
  tpm.y[is.na(tpm.y)] = 0
  #tpm.y2[is.na(tpm.y2)]=0

  for(j in 1:max(dat.d$states)){
    for(i in 1:nrow(tpm.y[,,j])){
      if(max(tpm.y[i,,j]) == 0) next
      tpm.y[i,,j] = tpm.y[i,,j]/sum(tpm.y[i,,j])
    }
  }

  #default
  olist=list("tpm.y2" = tpm.y2,"tpm.y" = tpm.y)
  return(olist)
} #end function

