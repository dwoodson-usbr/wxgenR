#' Precipitation simulator
#'
#' Simulate daily precipitation depth.
#'
#' @param dat.d Training data processed from prepData function.
#' @param syr Start year.
#' @param eyr End year.
#' @param wwidth Window set for finding surrounding days: +/- wwidth.
#' @param nsim Number of simulation years.
#' @param nrealz Number of realizations.
#' @param Xjday Julian day when prcp occurs, from simTPocc function
#' @param ekflag Simulate outside historical envelope?
#' @param parallelize Enable parallel computing for precip simulation, set T to enable
#'
# @import lubridate
# @import parallel
# @import doParallel
# @import foreach
# @import sm
# @rawNamespace import(stats, except = filter)
#'
#' @noRd
#'

"simPamt" <- function(dat.d,syr,eyr,wwidth,nsim,nrealz,Xjday,ekflag,awinFlag,parallelize){
  # require("sm")
  #simulate precipitation amounts and select dates
  #
  #get month for a given julian day
  lpyear = dat.d$year[min(which(leap_year(dat.d$year)))]
  aday <- ymd(paste(lpyear, 1, 1, sep="-")) #jan 1 of a leap year to have a 366-day year
  it1 = which(dat.d$date == aday)
  it2 = it1+366-1
  jdaymth <- dat.d$month[it1:it2] #calculations are based on a 366-day year
  #
  #get dates in window for each julian day 1-366
  Xdates=getDatesInWindow(syr,eyr,wwidth,leapflag=T)
  #
  if (ekflag){

    #get the bandwidth for each julian day
    bSJ <- vector()
    diwprcp <- vector()              #days in window with precipitation

    jday = 1
    for (jday in 1:366){

      diw <- na.omit(Xdates[jday,])  #dates in window for a given julian day
      idxlist <- vector()

      iday = 1
      for (iday in 1:length(diw)){
        idxlist[iday] = which(dat.d$date1==diw[iday])
      } #iday
      baprcp <- dat.d$prcp[idxlist]       #basin average precipitation
                                          #also includes 0 prcp amount
                                          #for days within the window

      pamt <- baprcp[which(baprcp>=0.01)] #precipitation amount vector
      wwidth.adapt = wwidth
      #
      while(length(pamt) < 2){

        wwidth.adapt = wwidth.adapt + 1

        #get dates in window for each julian day 1-366
        Xdates.adapt=getDatesInWindow(syr,eyr,wwidth.adapt,leapflag=T)
        diw <- na.omit(Xdates.adapt[jday,])  #dates in window for a given julian day
        idxlist <- vector()

        for (iday in 1:length(diw)){
          idxlist[iday]=which(dat.d$date1==diw[iday])
        } #iday
        baprcp <- dat.d$prcp[idxlist]

        pamt <- baprcp[which(baprcp>=0.01)] #precipitation amount vector

      }

      # if(awinFlag == T & wwidth.adapt != wwidth){
      #   cat(paste0("\n Window width too small on Julian day ", jday,", increased window to ", wwidth.adapt*2+1, " days\n"))
      # }

      diwprcp[jday] = length(pamt)
      logpamt <- log(pamt)                #log-transformed precipitation amount vector
      bSJ[jday] = hsj(logpamt)          #Sheather-Jones plug-in bandwidth

    } #jday
  } #ekflag
  #
  #simulate precipitation amount and select prcp date
  # Xpdate <- Xpamt <- matrix(NA, nrow=nsim*366, ncol=nrealz)

  if(parallelize == T){

    Xpdate <- Xpamt <- matrix(NA, nrow=nsim*366, ncol=1)

    # library(foreach)
    # library(doParallel)
    cl <- makePSOCKcluster(detectCores()-1)
    registerDoParallel(cl)

    # startTime <- Sys.time() #benchmark run time

    irealz = 1
    result <- foreach(irealz=1:nrealz, .packages='foreach', .export=c('repan', 'getDatesInWindow')) %dopar% {
    # for (irealz in 1:nrealz){

      message(paste0("-- Starting trace number ", irealz, " --"))
      xp = Xjday[,irealz]
      nxp = length(xp)

      ixp = 1
      # foreach(ixp=1:nxp, .combine = c) %do% {
      for (ixp in 1:nxp){

        jd = xp[ixp]

        if (is.na(jd)){
          Xpamt[ixp] = 0.0
        }else{

          diw <- na.omit(Xdates[jd,])  #dates in window for a given julian day
          idxlist <- vector()

          iday = 1
          for (iday in 1:length(diw)){
            idxlist[iday] = which(dat.d$date1 == diw[iday])
          } #iday

          baprcp <- dat.d$prcp[idxlist]       #e.g., basin average precipitation
          #also includes 0 prcp amount
          #days within the window

          np = length(which(baprcp>=0.01))      #number of prcp days in window
          pdate <- diw[which(baprcp>=0.01)]   #dates in window where prcp occurred
          pamt <- baprcp[which(baprcp>=0.01)] #precipitation amount vector

          wwidth.adapt = wwidth
          #
          while(np < 2){

            wwidth.adapt = wwidth.adapt + 1

            #get dates in window for each julian day 1-366
            Xdates.adapt = getDatesInWindow(syr,eyr,wwidth.adapt,leapflag=T)

            diw <- na.omit(Xdates.adapt[jd,])  #dates in window for a given julian day
            idxlist <- vector()

            iday = 1
            for (iday in 1:length(diw)){
              idxlist[iday]=which(dat.d$date1==diw[iday])
            } #iday

            baprcp <- dat.d$prcp[idxlist]       #e.g., basin average precipitation
            #also includes 0 prcp amount
            #days within the window

            np = length(which(baprcp>=0.01))      #number of prcp days in window
            pdate <- diw[which(baprcp>=0.01)]   #dates in window where prcp occurred
            pamt <- baprcp[which(baprcp>=0.01)] #precipitation amount vector

            message(paste0("\n Window width too small on Julian day ", jd,", increasing window to ", wwidth.adapt*2+1, " days"))

          }

          # if(awinFlag == T & wwidth.adapt != wwidth){
          #   print(paste0("\n Window width too small on Julian day ", jd,", increased window to ", wwidth.adapt*2+1, " days\n"))
          # }

          logpamt <- log(pamt)                #log-transformed prcp amount vector
          aindex = sample(1:np, 1)               #randomly pick a prcp day
          Xpdate[ixp] = pdate[aindex]
          ybar = logpamt[aindex]
          Xpamt[ixp] = exp(ybar)
          if (ekflag){
            rek = repan(1)                      #simulate a random number from the EKD
            Xpamt[ixp] = exp(ybar+rek*bSJ[jd])
          } #ekflag
        }
      } #ixp

      message("\n")

      list(Xpamt, Xpdate)

    } #irealz

    stopCluster(cl)

    # endTime = Sys.time()

    # timeP = difftime(endTime, startTime, units='mins')

    Xpamt = as.data.frame(do.call(cbind,lapply(result,function(x){x[[1]]})))
    Xpdate = as.data.frame(do.call(cbind,lapply(result,function(x){x[[2]]})))

    #
    #default
    olist=list("Xpamt"=Xpamt,"Xpdate"=Xpdate)
    if (ekflag) olist=list("Xpamt"=Xpamt,"Xpdate"=Xpdate,"bSJ"=bSJ)
    return(olist)

  }else{ #non-parallel loop

      Xpdate <- Xpamt <- matrix(NA, nrow=nsim*366, ncol=nrealz)

      # startTime <- Sys.time() #benchmark run time

      irealz = 1
      for (irealz in 1:nrealz){

      message(paste0("-- Starting trace number ", irealz, " --"))
      xp = Xjday[,irealz]
      nxp = length(xp)

      ixp = 1
      for (ixp in 1:nxp){

        jd = xp[ixp]

        if (is.na(jd)){
          Xpamt[ixp,irealz] = 0.0
        }else{

          diw <- na.omit(Xdates[jd,])  #dates in window for a given julian day
          idxlist <- vector()

          iday = 1
          for (iday in 1:length(diw)){
            idxlist[iday] = which(dat.d$date1 == diw[iday])
          } #iday

          baprcp <- dat.d$prcp[idxlist]       #e.g., basin average precipitation
                                              #also includes 0 prcp amount
                                              #days within the window

          np = length(which(baprcp>=0.01))      #number of prcp days in window
          pdate <- diw[which(baprcp>=0.01)]   #dates in window where prcp occurred
          pamt <- baprcp[which(baprcp>=0.01)] #precipitation amount vector

          wwidth.adapt = wwidth
          #
          while(np < 2){

            wwidth.adapt = wwidth.adapt + 1

            #get dates in window for each julian day 1-366
            Xdates.adapt = getDatesInWindow(syr,eyr,wwidth.adapt,leapflag=T)

            diw <- na.omit(Xdates.adapt[jd,])  #dates in window for a given julian day
            idxlist <- vector()

            iday = 1
            for (iday in 1:length(diw)){
              idxlist[iday]=which(dat.d$date1==diw[iday])
            } #iday

            baprcp <- dat.d$prcp[idxlist]       #e.g., basin average precipitation
            #also includes 0 prcp amount
            #days within the window

            np = length(which(baprcp>=0.01))      #number of prcp days in window
            pdate <- diw[which(baprcp>=0.01)]   #dates in window where prcp occurred
            pamt <- baprcp[which(baprcp>=0.01)] #precipitation amount vector

            message(paste0("\n Window width too small on Julian day ", jd,", increasing window to ", wwidth.adapt*2+1, " days"))

          }

          # if(awinFlag == T & wwidth.adapt != wwidth){
          #   print(paste0("\n Window width too small on Julian day ", jd,", increased window to ", wwidth.adapt*2+1, " days\n"))
          # }

          logpamt <- log(pamt)                #log-transformed prcp amount vector
          aindex = sample(1:np, 1)               #randomly pick a prcp day
          Xpdate[ixp,irealz] = pdate[aindex]
          ybar = logpamt[aindex]
          Xpamt[ixp,irealz] = exp(ybar)
          if (ekflag){
            rek = repan(1)                      #simulate a random number from the EKD
            Xpamt[ixp,irealz] = exp(ybar+rek*bSJ[jd])
          } #ekflag
        }
      } #ixp

      message("\n")
    } #irealz

    # endTime = Sys.time()

    # timeNP = difftime(endTime, startTime, units='mins')

    #
    #default
    olist=list("Xpamt"=Xpamt,"Xpdate"=Xpdate)
    if (ekflag) olist=list("Xpamt"=Xpamt,"Xpdate"=Xpdate,"bSJ"=bSJ)
    return(olist)

  } #end non-parallel loop

} #end function

