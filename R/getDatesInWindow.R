#' Get dates in window
#'
#' Find grouping of dates around each Julian day of year (1-366) based on the window you set.
#' The start and end years for this function should include at least one leap year
#'  (i.e., the record should be at least 4-years in length), or else the function will
#'  return non-existing dates (February 29th during non-leap years).\cr
#'  \cr
#' Setting leapflag to true will set February 29th as NA for non-leap years.\cr
#' \cr
#' Setting leapflag to false will remove February 29th for non-leap years (recommended).\cr
#' \cr
#' The 'wwidth' variable is the semi-bandwidth that sets the window size to search
#' for adjacent days. Given a value of 'wwidth', the window size will be
#' 2*wwidth + 1. For example a 'wwidth' of 7 would give a window size of
#' 2*7+1 = 15.\cr
#' \cr
#' Other applications of this function might include a daily bias correction approach
#' where it is necessary to find N adjacent days for each day of year in order to train
#' the bias correction algorithm.
#'
#' @param syr Start year.
#' @param eyr End year.
#' @param smo Start month.
#' @param emo End month.
#' @param sdate Start date.
#' @param edate End date.
#' @param wwidth Window set for finding surrounding days (semi-bandwidth).
#' @param leapflag Set index for leap years (default = F).
#'
#' @return Returns a matrix with 366 rows (one for each Julian day of year, including leap days)
#'  and nCols; where nCols = (2 x wwidth + 1) x (eyr - syr + 1). Each row is specific to a certain
#'  Julian day (e.g., day 1) and contains the preceding and antecedent dates around that Julian day
#'  based on the window length you set. The dates will be fetched for each year in the range you set
#'  between the start and ending years (inclusive of the start and end years). Matrix values are either dates
#'  formatted as 'yyyymmdd' or NA values.
#'
#'
#' @examples
#'  getDatesInWindow(syr = 2000, eyr = 2005, smo = 10, emo = 09,
#'   sdate = 20001001, edate = 20050930, wwidth = 3, leapflag = FALSE)
#'
#' @export

"getDatesInWindow" <- function(syr, eyr, smo, emo, sdate, edate, wwidth, leapflag = FALSE){

#input(s)
#syr - starting year
#eyr - ending year
#wwidth - semi-bandwidth of window, e.g., wwidth=1 results in a 3-day window width, (day-1),day,(day+1)
#leapflag - boolean with default as TRUE, sets index for 02/29 (Julian day 60) for non-leap years to NA;
#           if FALSE, the window for Julian day 60 uses days from March.  For example, with wwidth=1,
#           and say, syr=1987 and eyr=1999, then for Julian day 60 (02/29), the days in the window will be,
#           19870228,19870301,19870302,19880228,19880229,19880301,19890228,19890301,10890302, so on and so forth.
#           With leapflag=TRUE, the dates will be, 19870228,NA,19870301,19880228,19880229,19880301,19890228,NA,19890301, etc.
#
#output(s)
#matrix with list of dates (columns)for a given Julian day (1 through 366; rows)

# if(leapflag == F | leapflag == FALSE){
#   leapflag == FALSE
# } else{leapflag == TRUE}

# Month days
month = c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# # Create a custom ordering function
# reorder_months <- function(months, start, end) {
#   # Create a sequence from start to end wrapping around
#   if (start <= end) {
#     return(months[start:end])
#   } else {
#     return(c(months[start:length(months)], months[1:end]))
#   }
# }
#
# # Reorder the months
# month = reorder_months(month, smo, emo)

# Print the result
# print(month)

# Function to generate the custom sequence of month numbers
generate_month_sequence <- function(start, end) {
  if (start <= end) {
    return(start:end)  # Simple case: no wrap around
  } else {
    return(c(start:12, 1:end))  # Wrap around case
  }
}

# Generate the sequence
month_sequence = generate_month_sequence(smo, emo)

# Print the result
# print(month_sequence)

yrlist = syr:eyr
nyrs = length(yrlist)

noleapdate = rep(NA,length(yrlist)) #list of years with no leap dates

i = 1
for (i in 1:length(yrlist)){
  yr = yrlist[i]
  if (yr%%4 != 0) date = (yr*10000+2*100+29)
    else if (yr%%100 == 0 && yr%%400 != 0) date = (yr*10000+2*100+29)
	  else date = NA
  noleapdate[i] = date
} #i

rindex = which(is.na(noleapdate) == TRUE) #get indices with NA's, to be removed
noleapdate = noleapdate[-rindex] #removes NA, needed for later in the code, intersect

datevec = rep(NA, wwidth)
juldayvec = rep(NA, wwidth)
mmddvec = rep(NA, wwidth)

yr = yrlist[1]
mm = month_sequence[1]
id = 1
for (yr in yrlist){
  jday = 0
  for (mm in month_sequence){
    # print(mm)
    ndays = month[mm]
    # print(ndays)
    for (id in 1:ndays){
      jday = jday+1
      #date=as.numeric(paste((yr*10000+mm*100+id),sprintf("%03d",jday),sep=""))
      date = (yr*10000+mm*100+id)
      if(date < sdate | date > edate) date = NA #set NA for dates not between sdate and edate (an issue for non-calendar year approaches)
      datevec = append(datevec, date)
      juldayvec = append(juldayvec, jday)
      mmdd = mm*100+id
      mmddvec = append(mmddvec, mmdd)
    } #id
  } #mm
} #yr

datevec = append(datevec, rep(NA, wwidth))
juldayvec = append(juldayvec, rep(NA, wwidth))
mmddvec = append(mmddvec, rep(NA, wwidth))

JMAT = matrix(NA, nrow = 366, ncol = nyrs*(2*wwidth+1)) #matrix with list of dates (columns)for a given Julian day (1 through 366)

#loop through all the 366 Julian days
julday=1
for (julday in 1:366){
  avec = which(juldayvec == julday)
  indexlist = vector()
  a = avec[1]
  for (a in avec){
    i1 = (a-wwidth)
    i2 = (a+wwidth)
    index = i1:i2

    if (leapflag){
      # print("leap = T")
      index[which((datevec[index]%in%noleapdate) == TRUE)] = NA
    }else{
      # print("leap = F")
      if (length(intersect(datevec[index],noleapdate)) > 0){
        rindex = which((datevec[index]%in%noleapdate) == TRUE) #remove index
        index = i1:(i2+1)
        index = index[-rindex]
      } #endif
    }
    indexlist = append(indexlist, index)
  } #a

  JMAT[julday,] = datevec[indexlist]

} #julday

#get month-day-julian day relationship
mdjday = vector()
jday = 0
for (mm in month_sequence){
  ndays = month[mm]
  for (id in 1:ndays){
    jday = jday+1
    mdjday[jday] = 100*mm+id
  } #id
} #mm


return(JMAT)

} #end function
