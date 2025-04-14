#' Spell length calculation
#'
#' Function to calculate the length (duration in years) of wet or dry periods.
#'
#' @param s A binary vector of 0 dry and 1 wet only.
#'
#' @return Returns a list object containing a vector of dry spell lengths and a vector of wet spell lengths.
#'
#' @examples
#'
#'  #use 0 for dry and 1 for wet years
#'  spells = c(0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0)
#'
#'  spellLengths(spells)
#'
#' @export
#'
#'
"spellLengths" <- function(s){
  #input s is a binary vector of 0 (dry) and 1 (wet) only
  n=length(s)
  l=list("0" = c(),"1" = c()) #"0": dry; "1":wet

  cnt=1

  for (i in 2:n){
    if (s[i-1] == s[i]){
      cnt = cnt+1
    }
    else{
      l[[as.character(s[i-1])]]=c(l[[as.character(s[i-1])]],cnt)
      cnt = 1
    }
  }#i
  l[[as.character(s[n])]]=c(l[[as.character(s[n])]],cnt)

  return(l)

} #end function

