#' Select transition state
#'
#' Function selects and returns the transition state given a uniform random number between 0 and 1 and
#' the cumulative probability vector of the state sequence.
#'
#' @param uni Uniform random number between 0 and 1.
#' @param wt Cumulative probability vector of states.
#'
#'
#' @return Returns an object containing the transition state(s) based on the given cumulative probability vector and random numbers.
#'
#'
#' @examples
#'
#'  rand = runif(1)
#'
#'  print(rand)
#'
#'  selectState(uni = rand, wt = c(0.25, .55, 0.85, 1))
#'
#'
#' @export
#'

##
"selectState" <- function(uni, wt){
  #function returns state given a uniform random number and
  #the cumulative probability vector

  k = length(wt) #number of system states

  if (uni <= wt[1]){
    state = 1
    return(state)
  }

  if (uni == wt[k]){   #since max(wt[k])=1.0, so no need for (uni >= wt[k])
    state = k          #as we are using cumulative transition probability
    return(state)
  }

  jl = 1; ju = k

  while ((ju - jl) > 1 ){
    jn = ceiling((ju+jl)/2)
    if(wt[jn] > uni){
      ju = jn
    }
    else{
      jl = jn
    }
  } #end while

  state = ju
  return(state)

} #end function

