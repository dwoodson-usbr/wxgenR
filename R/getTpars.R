#' Calculate temperature parameters
#'
#' Imports precipitation and temperature data,
#' then fits a linear model to predict daily
#' temperature based on the prior day’s temperature,
#' sine and cosine functions, monthly mean temperature,
#' and the occurrence of precipitation.
#'
#' @param dat.d Training data processed from prepData wrapper function.
#'
# @import lubridate
# @rawNamespace import(stats, except = filter)
# @importFrom dplyr group_by summarise left_join glimpse mutate relocate if_else filter

#'
#' @noRd
#'

"getTpars" <- function(dat.d, pcpOccFlag, traceThreshold, smo, emo){
  #calculate paramters for temperature simulation
  #
  # require("plyr")
  # require("dplyr")
  #temperature data, temp for time t
  tmp = dat.d$temp
  #temperature for time t-1
  ptmp = tmp
  ptmp[1] = NA
  ptmp[2:length(ptmp)] = tmp[1:(length(tmp)-1)]
  yr.d=dat.d$year
  uyr=unique(yr.d)
  nyr=length(uyr)
  mo.d=dat.d$month
  #define day of year
  tmp1 = tmp2 = c()

  if(smo == 1){
    k=1
    for(k in 1:nyr){
      origin.tmp = ymd(paste(uyr[k],smo,"01",sep="-"))
      start.tmp  = julian(ymd(subset(dat.d, year==uyr[k])$date)[1],origin=origin.tmp)
      end.tmp    = julian(ymd(subset(dat.d, year==uyr[k])$date)[nrow(subset(dat.d, year==uyr[k]))], origin=origin.tmp)
      tmp1 = c(tmp1, seq(from=start.tmp,to=end.tmp))
      tmp2 = c(tmp2, rep(end.tmp, nrow(subset(dat.d, year==uyr[k]))))
    }#k
  }else if(smo > 1){
    k=1
    for(k in 1:(nyr-1)){
      origin.tmp = ymd(paste(uyr[k],smo,"01",sep="-"))
      start.tmp  = julian(ymd(subset(dat.d, year==uyr[k] & month >= smo)$date)[1],origin=origin.tmp)
      end.tmp    = julian(ymd(subset(dat.d, year==uyr[k+1] & month <= emo)$date)[nrow(subset(dat.d, year==uyr[k+1] & month <= emo))], origin=origin.tmp)
      tmp1 = c(tmp1, seq(from=start.tmp,to=end.tmp))
      tmp2 = c(tmp2, rep(end.tmp, nrow(subset(dat.d, year==uyr[k] & month >= smo))+nrow(subset(dat.d, year==uyr[k+1] & month <= emo))))
    }#k
  }

  #change 0-364 to 1-365
  dat.d$jday = tmp1+1
  dat.d$tday = tmp2+1
  #define cos(t) and sin(t) for daily temp series
  ct <- cos((2*pi*dat.d$jday)/dat.d$tday)
  st <- sin((2*pi*dat.d$jday)/dat.d$tday)
  #monthly mean temperature by year
  montmp.obs = data.frame(year=yr.d, month=mo.d, temp=dat.d$temp) %>%
    group_by(year, month) %>%
    summarise(tavgm = mean(temp, na.rm = T))
  # montmp.obs = ddply(data.frame(year=yr.d, month=mo.d, temp=dat.d$temp),
  #                   .(year,month), summarise, tavgm=mean(temp))
  dat.d = left_join(dat.d, montmp.obs, by=c("year","month"))
  Rt <- dat.d$tavgm
  #define precip occurrence for daily temp series
  oc <- (dat.d$prcp >= traceThreshold) + 0
  #set NAs to 0 precipitation
  # oc[which(is.na(oc), T)] = 0
  dat.d$oc=oc
  dat.d$ct=ct   #cosine term
  dat.d$st=st   #sine term
  #define design matrix (covariates)
  #temp(t) is a function of:
  #[temp(t-1); cosine(t); sine(t); prec.occ(t); mon.mean.temp(t)]
  if(pcpOccFlag == TRUE){
    x.tmp <- cbind(ptmp, ct, st, oc, Rt)
  }else if(pcpOccFlag == FALSE){
    x.tmp <- cbind(ptmp, ct, st, Rt)
    }
  z.tmp <- lm(tmp ~ x.tmp)
  z.tmp.res <- z.tmp$residuals
  coeftmp <- z.tmp$coefficients
  tmp.sd <- numeric(12)
  for(i in 1:12) tmp.sd[i] <- sd(z.tmp.res[mo.d==i],na.rm=T)
  #
  #default
  olist=list("dat.d"=dat.d,"coeftmp"=coeftmp,"tmp.sd"=tmp.sd)
  return(olist)
} #end function

