#' N states probability
#'
#' Returns an integer vector corresponding to n states broken by equal
#' probability or equal distance.
#'
#' @noRd
#'
# @rawNamespace import(stats, except = filter)
#'
#'
nStatProb <-   function(x, n, limit.type = 'prob', limits = NULL, tie = 1, altobs = NULL ){
  # returns an integer vector corresponding to n states broken by equal
  # probability or equal distance
  #
  limit <-
    if(limit.type == 'prob')
      quantile(x,seq(0,1,1/n))
  else if(limit.type == 'equal')
    seq(min(x),max(x),by=diff(range(x))/n)
  else if(limit.type == 'manual')
    limits

  if(!is.null(altobs)) limit <- quantile(altobs,seq(0,1,1/n))

  b <- integer(length(x))

  for(i in 1:(n+1)){
    filter <-
      if(tie == 1)
        x >= limit[i] & x <= limit[i+1]
    else
      x > limit[i] & x <= limit[i+1]

    #only need to set the 1's because b is already 0's
    b[filter] <- as.integer(i-1)
  }

  # if(class(x) == 'ts')
  if(inherits(x, 'ts')){
    return(ts(b,start=start(x),end=end(x)))
  }else{
    return(b)
  }
} #end function

#' TPM
#'
#' Checks transition probability matrix.
#'
# @import msm
#'
#' @noRd
#'
#'
transProbMatrix <- function(x,ns=NULL,limits=NULL,tie=0){

  # require(msm)

  if(is.null(ns)){
    ns <- max(x)
    states <- x
    if(length(unique(states)) > 26) stop('Too many states, specify a smaller number.')
  }
  # else{
  #   states <- ntile.ts(x,n=ns,limit.type='manual',limits=limits,tie=tie)
  # }

  st <- statetable.msm(state,data=list(state=states))
  st/apply(st,1,sum)

} #end function

#Convert month strings to numeric
# Standardize input to match month names or abbreviations
match_month <- function(month) {
  month <- tolower(month)
  match <- match(tolower(substr(month, 1, 3)), tolower(month.abb))
  return(match)
}

#get days in month of any start and end month sequence
days_in_months <- function(sd, ed) {
  # Define days in each month for a leap year
  days_in_month_leap <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

  # Ensure sd and ed are numeric
  if(is.character(sd)) sd <- match_month(sd)
  if(is.character(ed)) ed <- match_month(ed)

  # Handle cases where ed < sd (e.g., water year October to September)
  if(ed < sd) ed <- ed + 12

  # Create a sequence of months from sd to ed
  months_seq <- (sd:ed) %% 12
  months_seq[months_seq == 0] <- 12

  # Extract the relevant days from the leap year vector
  days_in_month <- days_in_month_leap[months_seq]

  return(days_in_month)
}
