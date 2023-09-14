#' Read in data to train simulator
#'
#' Read training data and setup variables to facilitate simulation.
#'
#' @param trainingData A path to .csv file (or matrix, data frame in your environment)
#' with the following variables
#'  is required: year, month, day, prcp (daily precipitation),
#'   temp (daily temperature),
#'    and season (1, 2, ..., N, for N seasons - up to 20 seasons will work).
#'    Any units will work for precipitation and temperature, as long as they are
#'    consistent. Can be station data, basin averages, grid cells, etc.
#'
#' @param sdate Start date of training data (yyyymmdd).
#' If empty, the start date will be the beginning of your time series.
#' @param edate End date of training data (yyyymmdd).
#' If empty, the end date will be the end of your time series.
#'
#' @examples
#'
#' prepData(trainingData = "./MetData.csv",  sdate = 20000101, edate = 20201231)
#'
# @import lubridate
#'
#' @noRd


"prepData" <- function(trainingData, sdate, edate){

  # require("lubridate")
  # if(typeof(trainingData) == "character"){
  #   dat.d = read.table(trainingData, header=T, sep=",")
  # } else{
    dat.d = trainingData
  # }

  dat.d$date1 = dat.d$year*10000 + dat.d$month*100 + dat.d$day

  i1 = which(dat.d$date1 == sdate)
  i2 = which(dat.d$date1 == edate)

  dat.d = dat.d[i1:i2,]
  yr.d = dat.d$year
  mo.d = dat.d$month
  da.d = dat.d$day
  wk.d = week(ymd(paste(yr.d, mo.d, da.d, sep="-")))
  wk.d[wk.d == 53] = 52
  dat.d$week = wk.d
  dat.d$date = ymd(paste(yr.d, mo.d, da.d, sep="-"))
  dat.d$states = dat.d$season
  #
  #default
  return(dat.d)
} #end function


