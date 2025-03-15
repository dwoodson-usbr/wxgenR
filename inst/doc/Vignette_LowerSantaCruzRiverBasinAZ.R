## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(wxgenR)
library(lubridate)
library(dplyr)
library(tidyr)
library(ggpubr)
library(moments)
library(seas)

data(LowerSantaCruzRiverBasinAZ)

head(LowerSantaCruzRiverBasinAZ)

## ----results = 'hide'---------------------------------------------------------
nsim = 5   #number of simulation years
nrealz = 2 #number of traces in ensemble

startTime <- Sys.time() #benchmark run time

z = wx(trainingData = LowerSantaCruzRiverBasinAZ,
       syr = 1970, eyr = 1974, smo = 10, emo = 9,
       nsim = nsim, nrealz = nrealz, aseed = 123,
       wwidth = c(7,1,3), unitSystem = "U.S. Customary", ekflag = TRUE,
       awinFlag = T, tempPerturb = TRUE, pcpOccFlag = FALSE,
       numbCores = 2)

endTime = Sys.time()


## -----------------------------------------------------------------------------
glimpse(z)

## -----------------------------------------------------------------------------
#parse variables from wx() output
dat.d = z$dat.d
simyr1 = z$simyr1
X = z$X
Xseas = z$Xseas
Xpdate = z$Xpdate
Xpamt = z$Xpamt
Xtemp = z$Xtemp

#write simulation output
#
it1 <- seq(1, length(X[,1]), 366)
it2 = it1+366-1

#initialize storage
sim.pcp = matrix(NA, nrow = nsim*366, ncol = nrealz+3)
sim.tmp = matrix(NA, nrow = nsim*366, ncol = nrealz+3)
sim.szn = matrix(NA, nrow = nsim*366, ncol = nrealz+3)

#loop through realization
irealz = 1
for (irealz in 1:nrealz){
  
  outmat <- vector()

  #loop through simulation years
  isim = 1
  for (isim in 1:nsim){
    leapflag = FALSE
    ayr = simyr1[isim, irealz]
    if (lubridate::leap_year(ayr)) leapflag = TRUE
    col1 = rep(isim, 366)        #column 1, simulation year
    d1 = ayr*10^4+01*10^2+01; d2 = ayr*10^4+12*10^2+31
    i1 = which(dat.d$date1 == d1)
    i2 = which(dat.d$date1 == d2)
    col2 = dat.d$date1[i1:i2]   #column 2, simulation date
    if (leapflag == FALSE) col2 = c(col2,NA)
    i1 = it1[isim]
    i2 = it2[isim]
    col3 = Xseas[i1:i2, irealz]  #column 3, simulation season
    col4 = X[i1:i2, irealz]      #column 4, precipitation occurrence
    col5 = Xpdate[i1:i2, irealz] #column 5, precipation resampling date
    col6 = Xpamt[i1:i2, irealz]  #column 6, resampled precipitation amount
    col7 = Xtemp[i1:i2, irealz]  #column7, simulated temperature

    #create time series of 'simulation day'
    sim.yr = rep(isim, length(col2))
    sim.month = month(ymd(col2))
    sim.day = day(ymd(col2))

    outmat = rbind(outmat, cbind(sim.yr, sim.month, sim.day, col6, col7, col3))
  } #isim

  colnames(outmat) = c("simulation year", "month", "day", "prcp", "temp", "season")
  
  if(irealz == 1){
    sim.pcp[,1:3] = outmat[,1:3]
    sim.tmp[,1:3] = outmat[,1:3]
    sim.szn[,1:3] = outmat[,1:3]
  }
  
    sim.pcp[,irealz+3] = outmat[,4]
    sim.tmp[,irealz+3] = outmat[,5]
    sim.szn[,irealz+3] = outmat[,6]

} #irealz


## -----------------------------------------------------------------------------

# df = sim.pcp

formatting = function(df){
  
  df = as.data.frame(df)
  
  colnames(df) = c("simulation year", "month", "day", paste0("Trace_", 1:nrealz))

  #remove 366 days for non-leap years
  df = drop_na(df, c(month, day))
  
  #assign simulation year to start at the same time as training data
  df$`simulation year` = df$`simulation year` + dat.d$year[1] - 1    
  
  #format date
  df$Date = ymd(paste(df$`simulation year`, df$month, df$day, sep = "-"))
  
  #remove years that aren't leap years
  # df = drop_na(df, Date)
  
  df = df %>%
    mutate(yday = as.numeric(yday(Date)),
           week = as.numeric(week(Date))) %>%
    relocate(c(Date,yday,week), .after = day) %>%
    pivot_longer(cols = starts_with("Trace_"), 
                 names_to = "variable", 
                 values_to = "value")
  return(df)
}

## -----------------------------------------------------------------------------
sim.pcp = formatting(sim.pcp)
sim.tmp = formatting(sim.tmp)
sim.szn = formatting(sim.szn)

## -----------------------------------------------------------------------------

colnames(dat.d)[11] = "yday"
obs.pcp = dat.d[,c(1:3,8:9,11,4)]
obs.tmp = dat.d[,c(1:3,8:9,11,5)]


## -----------------------------------------------------------------------------
#plot simulated daily data
# simDat = sim.tmp
# obsDat = obs.tmp
# Tag = "Temp"

dailyPlot = function(simDat, obsDat, Tag){
  
  simD = simDat %>%
    drop_na() %>%
    group_by(variable, yday) %>%
    summarise(
      mean = mean(value, na.rm = T),
      max = max(value, na.rm = T),
      sd = sd(value, na.rm = T),
      skew = skewness(value, na.rm = T)
      ) %>%
    ungroup()
    
  simDq <- simD %>%
    group_by(yday) %>%
    summarise(
      mean_q5 = quantile(mean, 0.05, na.rm = T),
      mean_med = median(mean, na.rm = T),
      mean_q95 = quantile(mean, 0.95, na.rm =T),
      max_q5 = quantile(max, 0.05, na.rm = T),
      max_med = median(max, na.rm = T),
      max_q95 = quantile(max, 0.95, na.rm = T),
      sd_q5 = quantile(sd, 0.05, na.rm = T),
      sd_med = median(sd),
      sd_q95 = quantile(sd, 0.95, na.rm = T),
      skew_q5 = quantile(skew, 0.05, na.rm = T),
      skew_med = median(skew, na.rm = T),
      skew_q95 = quantile(skew, 0.95, na.rm = T)
      ) %>%
      drop_na() %>%
    ungroup()
  
  if(Tag == "Temp"){
    obs <- obsDat %>%
      drop_na() %>%
      group_by(yday) %>%
      summarise(
        mean = mean(temp, na.rm = T),
        max = max(temp, na.rm = T),
        sd = sd(temp, na.rm = T),
        skew = skewness(temp, na.rm = T)
        ) %>%
      ungroup()
  } else if(Tag == "Precip"){
    obs <- obsDat %>%
      drop_na() %>%
      group_by(yday) %>%
      summarise(
        mean = mean(prcp, na.rm = T),
        max = max(prcp, na.rm = T),
        sd = sd(prcp, na.rm = T),
        skew = skewness(prcp, na.rm = T)
        ) %>%
      ungroup()
    }

  colnames(obs)[-1] = paste0("obs_", colnames(obs)[-1])
  
  df.comb = left_join(simDq, obs, by = "yday")
  
  #plotting --------------------------------
  lgdLoc = c(0.8, 0.9)
  
  if(Tag == "Temp"){
    yLabel = "Daily Temperature "
    units = "(°F)"
  } else if(Tag == "Precip"){
    yLabel = "Daily Precipitation "
    units = "(inches)"
  }
  
  trnAlpha = 0.65
  
  #daily mean
  p1 = ggplot(df.comb) +
    geom_ribbon(aes(x = yday, ymin = mean_q5, ymax = mean_q95, fill = "95% Confidence"), alpha = 0.25) +
    geom_line(aes(x = yday, y = mean_med, color = "Simulation Median"), size = 1, alpha = 0.8) +
    geom_line(aes(x = yday, y = obs_mean, color = "Training Data"), size = 0.3, alpha = trnAlpha, linetype = "solid") +
    geom_point(aes(x = yday, y = obs_mean, color = "Training Data"), size = 0.6, alpha = trnAlpha) +
    
    # Corrected scale_colour_manual()
    scale_colour_manual(values = c("Training Data" = "blue", "Simulation Median" = "red")) +
    
    # Ensure fill is mapped for the ribbon
    scale_fill_manual(values = c("95% Confidence" = "grey")) +
  
    theme_classic() +
    theme(axis.title = element_text(face = "bold"),
          panel.grid.major = element_line(),
          legend.title = element_blank(),
          legend.position = lgdLoc,
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.key = element_blank()) +
    xlab("Day of Year") + ylab(paste0("Mean ", yLabel, units))
  
  #daily SD
  p2 = ggplot(df.comb) +
    geom_ribbon(aes(x = yday, ymin = sd_q5, ymax = sd_q95, fill = "95% Confidence"), alpha = 0.25) +
    geom_line(aes(x = yday, y = sd_med, color = "Simulation Median"), size = 1, alpha = 0.8) +
    geom_line(aes(x = yday, y = obs_sd, color = "Training Data"), size = 0.3, alpha = trnAlpha, linetype = "solid") +
    geom_point(aes(x = yday, y = obs_sd, color = "Training Data"), size = 0.6, alpha = trnAlpha) +
  
    scale_colour_manual(values = c("Training Data" = "blue", "Simulation Median" = "red")) +
    scale_fill_manual(values = c("95% Confidence" = "grey")) +
  
    theme_classic() +
    theme(axis.title = element_text(face = "bold"),
          panel.grid.major = element_line(),
          legend.title = element_blank(),
          legend.position = lgdLoc,
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.key = element_blank()) +
    xlab("Day of Year") + ylab(paste0("Std. Deviation of ", yLabel, units))

  #daily skew
  p3 = ggplot(df.comb) +
    geom_ribbon(aes(x = yday, ymin = skew_q5, ymax = skew_q95, fill = "95% Confidence"), alpha = 0.25) +
    geom_line(aes(x = yday, y = skew_med, color = "Simulation Median"), size = 1, alpha = 0.8) +
    geom_line(aes(x = yday, y = obs_skew, color = "Training Data"), size = 0.3, alpha = trnAlpha, linetype = "solid") +
    geom_point(aes(x = yday, y = obs_skew, color = "Training Data"), size = 0.6, alpha = trnAlpha) +
  
    scale_colour_manual(values = c("Training Data" = "blue", "Simulation Median" = "red")) +
    scale_fill_manual(values = c("95% Confidence" = "grey")) +
  
    theme_classic() +
    theme(axis.title = element_text(face = "bold"),
          panel.grid.major = element_line(),
          legend.title = element_blank(),
          legend.position = lgdLoc,
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.key = element_blank()) +
    xlab("Day of Year") + ylab(paste0("Skew of ", yLabel, " (-)"))

    
  
  #daily Max
  p4 = ggplot(df.comb) +
    geom_ribbon(aes(x = yday, ymin = max_q5, ymax = max_q95, fill = "95% Confidence"), alpha = 0.25) +
    geom_line(aes(x = yday, y = max_med, color = "Simulation Median"), size = 1, alpha = 0.8) +
    geom_line(aes(x = yday, y = obs_max, color = "Training Data"), size = 0.3, alpha = trnAlpha, linetype = "solid") +
    geom_point(aes(x = yday, y = obs_max, color = "Training Data"), size = 0.6, alpha = trnAlpha) +
  
    scale_colour_manual(values = c("Training Data" = "blue", "Simulation Median" = "red")) +
    scale_fill_manual(values = c("95% Confidence" = "grey")) +
  
    theme_classic() +
    theme(axis.title = element_text(face = "bold"),
          panel.grid.major = element_line(),
          legend.title = element_blank(),
          legend.position = lgdLoc,
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.key = element_blank()) +
    xlab("Day of Year") + ylab(paste0("Maximum ", yLabel, units))
  
  p.comb = ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2, common.legend = TRUE, legend = "bottom")
  
  print(p.comb)
  
  # p.out = paste0(tempdir(), "/outputPlots/dailyStats_", Tag, ".png")
  # ggsave(filename = p.out, plot = p.comb, device = "png")

}


## ----fig.width=8, fig.height=8------------------------------------------------
dailyPlot(sim.pcp, obs.pcp, "Precip")

## ----fig.width=8, fig.height=8------------------------------------------------
dailyPlot(sim.tmp, obs.tmp, "Temp")

## -----------------------------------------------------------------------------
#plot simulated daily data
# simDat = sim.tmp
# obsDat = obs.tmp
# Tag = "Temp"

monthlyPlot = function(simDat, obsDat, Tag){
  
  if(Tag == "Temp"){
    
    simM = simDat %>%
      drop_na() %>%
      group_by(variable, month, `simulation year`) %>%
      summarise(
        mean = mean(value, na.rm = T),
        max = max(value, na.rm = T),
        sd = sd(value, na.rm = T),
        skew = skewness(value, na.rm = T)
        ) %>%
      ungroup()
    
    simMM <- simM %>%
      group_by(variable, month) %>%
      summarise(
        mean=mean(mean),
        max=mean(max),
        sd=sqrt(mean(sd^2)),
        skew=mean(skew, na.rm=T)
        ) %>%
      ungroup()
    
    obs <- obsDat %>%
      drop_na() %>%
      group_by(month, year) %>%
      summarise(
        mean = mean(temp, na.rm = T),
        max = max(temp, na.rm = T),
        sd = sd(temp, na.rm = T),
        skew = skewness(temp, na.rm = T)
        ) %>%
      ungroup()
    
    obsMM <- obs %>%
      group_by(month) %>%
      summarise(
        mean = mean(mean, na.rm = T),
        max = mean(max, na.rm = T),
        sd = sqrt(mean(sd^2)),
        skew = mean(skew, na.rm=T)
        ) %>%
      mutate(variable = "Observed") %>%
      relocate(variable) %>%
      ungroup()
    
  # colnames(obsMM)[-1] = paste0("obs_", colnames(obsMM)[-1])
  
  }else if(Tag == "Precip"){
    
    simM = simDat %>%
      drop_na() %>%
      group_by(variable, month, `simulation year`) %>%
      summarise(
        sum = sum(value, na.rm = T),
        max = max(value, na.rm = T),
        sd = sd(value, na.rm = T),
        skew = skewness(value, na.rm = T)
        ) %>%
      ungroup()  
    
      simMM <- simM %>%
        group_by(variable, month) %>%
        summarise(
          sum=mean(sum),
          max=mean(max),
          sd=sqrt(mean(sd^2)),
          skew=mean(skew, na.rm=T)
          ) %>%
        ungroup()
      
      obs <- obsDat %>%
        drop_na() %>%
        group_by(month, year) %>%
        summarise(
          sum = sum(prcp, na.rm = T),
          max = max(prcp, na.rm = T),
          sd = sd(prcp, na.rm = T),
          skew = skewness(prcp, na.rm = T)
          ) %>%
        ungroup()
    
      obsMM <- obs %>%
        group_by(month) %>%
        summarise(
          sum = mean(sum, na.rm = T),
          max = mean(max, na.rm = T),
          sd = sqrt(mean(sd^2)),
          skew = mean(skew, na.rm=T)
          ) %>%
        mutate(variable = "Observed") %>%
        relocate(variable) %>%
        ungroup()
    
  # colnames(obsMM)[-1] = paste0("obs_", colnames(obsMM)[-1])
  
  }
  
  df.comb = rbind(obsMM, simMM)
  
  #plotting --------------------------------
  if(Tag == "Temp"){
    
    p1 = ggplot(df.comb) + 
      geom_boxplot(data = subset(df.comb, variable != "Observed"), aes(x = month, y = mean, group = month)) +
      geom_line(data = subset(df.comb, variable == "Observed"), size = 0.5, aes(x = month, y = mean, color = "Observed")) +
      geom_point(data = subset(df.comb, variable == "Observed"), size = 1.5, aes(x = month, y = mean, color = "Observed")) +
      xlab("Month") + ylab("Temperature (°F)") +
      theme_classic() +
      theme(axis.title = element_text(face = "bold"),
            text=element_text(size=12),
            panel.grid.major = element_line(),
            legend.title=element_blank(),
            plot.title = element_text(size=10)
            ) +
      scale_x_continuous(breaks = 1:12) +
      ggtitle("Average Mean Monthly Temperature")
  
  }else if(Tag == "Precip"){
    
     p1 = ggplot(df.comb) + 
      geom_boxplot(data = subset(df.comb, variable != "Observed"), aes(x = month, y = sum, group = month)) +
      geom_line(data = subset(df.comb, variable == "Observed"), size = 0.5, aes(x = month, y = sum, color = "Observed")) +
      geom_point(data = subset(df.comb, variable == "Observed"), size = 1.5, aes(x = month, y = sum, color = "Observed")) +
      xlab("Month") + ylab("Precipitation (inches)") +
      theme_classic() +
      theme(axis.title = element_text(face = "bold"),
            text=element_text(size=12),
            panel.grid.major = element_line(),
            legend.title=element_blank(),
            plot.title = element_text(size=10)
            ) +
      scale_x_continuous(breaks = 1:12) +
      ggtitle("Average Total Monthly Precipitation")
    
  }
  
  if(Tag == "Temp"){
    yLabel = "Temperature "
    units = "(°F)"
  } else if(Tag == "Precip"){
    yLabel = "Precipitation "
    units = "(inches)"
  }
  
  #monthly SD
  p2 = ggplot(df.comb) + 
      geom_boxplot(data = subset(df.comb, variable != "Observed"), aes(x = month, y = sd, group = month)) +
      geom_line(data = subset(df.comb, variable == "Observed"), size = 0.5, aes(x = month, y = sd, color = "Observed")) +
      geom_point(data = subset(df.comb, variable == "Observed"), size = 1.5, aes(x = month, y = sd, color = "Observed")) +
      xlab("Month") + ylab(paste0("Standard Deviation ", units)) +
      theme_classic() +
      theme(axis.title = element_text(face = "bold"),
            text=element_text(size=12),
            panel.grid.major = element_line(),
            legend.title=element_blank(),
            plot.title = element_text(size=10)
            ) +
      scale_x_continuous(breaks = 1:12) +
      ggtitle(paste0("Average Standard Deviation in Monthly ", yLabel))
  
  
  #monthly Skew
  p3 = ggplot(df.comb) + 
      geom_boxplot(data = subset(df.comb, variable != "Observed"), aes(x = month, y = skew, group = month)) +
      geom_line(data = subset(df.comb, variable == "Observed"), size = 0.5, aes(x = month, y = skew, color = "Observed")) +
      geom_point(data = subset(df.comb, variable == "Observed"), size = 1.5, aes(x = month, y = skew, color = "Observed")) +
      xlab("Month") + ylab("Skew (-)") +
      theme_classic() +
      theme(axis.title = element_text(face = "bold"),
            text=element_text(size=12),
            panel.grid.major = element_line(),
            legend.title=element_blank(),
            plot.title = element_text(size=10)
            ) +
      scale_x_continuous(breaks = 1:12) +
      ggtitle(paste0("Average Skew in Monthly ", yLabel))
  
  
  #monthly max
  p4 = ggplot(df.comb) + 
      geom_boxplot(data = subset(df.comb, variable != "Observed"), aes(x = month, y = max, group = month)) +
      geom_line(data = subset(df.comb, variable == "Observed"), size = 0.5, aes(x = month, y = max, color = "Observed")) +
      geom_point(data = subset(df.comb, variable == "Observed"), size = 1.5, aes(x = month, y = max, color = "Observed")) +
      xlab("Month") + ylab(paste0("Maximum ", units)) +
      theme_classic() +
      theme(axis.title = element_text(face = "bold"),
            text=element_text(size=12),
            panel.grid.major = element_line(),
            legend.title=element_blank(),
            plot.title = element_text(size=10)
            ) +
      scale_x_continuous(breaks = 1:12) +
      ggtitle(paste0("Average Monthly Maximum ", yLabel))
  
  p.comb = ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2, common.legend = TRUE, legend = "bottom")
  
  print(p.comb)
  
  # p.out = paste0(tempdir(), "/outputPlots/monthlyStats_", Tag, ".png")
  # ggsave(filename = p.out, plot = p.comb, device = "png", height = 8, width = 8, units = "in")

}


## ----fig.width=8, fig.height=8------------------------------------------------
monthlyPlot(sim.pcp, obs.pcp, "Precip")

## ----fig.width=8, fig.height=8------------------------------------------------
monthlyPlot(sim.tmp, obs.tmp, "Temp")

## -----------------------------------------------------------------------------
#plot simulated daily data
simDat = sim.tmp
obsDat = obs.tmp
Tag = "Temp"

weeklyPlot = function(simDat, obsDat, Tag){
  
  if(Tag == "Temp"){
    
    simW = simDat %>%
      drop_na() %>%
      group_by(variable, week, `simulation year`) %>%
      summarise(
        mean = mean(value, na.rm = T),
        max = max(value, na.rm = T),
        sd = sd(value, na.rm = T),
        skew = skewness(value, na.rm = T)
        ) %>%
      ungroup()
    
    simWW <- simW %>%
      group_by(variable, week) %>%
      summarise(
        mean=mean(mean),
        max=mean(max),
        sd=sqrt(mean(sd^2)),
        skew=mean(skew, na.rm=T)
        ) %>%
      ungroup()
    
    obs <- obsDat %>%
      drop_na() %>%
      group_by(week, year) %>%
      summarise(
        mean = mean(temp, na.rm = T),
        max = max(temp, na.rm = T),
        sd = sd(temp, na.rm = T),
        skew = skewness(temp, na.rm = T)
        ) %>%
      ungroup()
    
    obsWW <- obs %>%
      group_by(week) %>%
      summarise(
        mean = mean(mean, na.rm = T),
        max = mean(max, na.rm = T),
        sd = sqrt(mean(sd^2)),
        skew = mean(skew, na.rm=T)
        ) %>%
      mutate(variable = "Observed") %>%
      relocate(variable) %>%
      ungroup()
    
  # colnames(obsMM)[-1] = paste0("obs_", colnames(obsMM)[-1])
  
  }else if(Tag == "Precip"){
    
    simW = simDat %>%
      drop_na() %>%
      group_by(variable, week, `simulation year`) %>%
      summarise(
        sum = sum(value, na.rm = T),
        max = max(value, na.rm = T),
        sd = sd(value, na.rm = T),
        skew = skewness(value, na.rm = T)
        ) %>%
      ungroup()  
    
      simWW <- simW %>%
        group_by(variable, week) %>%
        summarise(
          sum=mean(sum),
          max=mean(max),
          sd=sqrt(mean(sd^2)),
          skew=mean(skew, na.rm=T)
          ) %>%
        ungroup()
      
      obs <- obsDat %>%
        drop_na() %>%
        group_by(week, year) %>%
        summarise(
          sum = sum(prcp, na.rm = T),
          max = max(prcp, na.rm = T),
          sd = sd(prcp, na.rm = T),
          skew = skewness(prcp, na.rm = T)
          ) %>%
        ungroup()
    
      obsWW <- obs %>%
        group_by(week) %>%
        summarise(
          sum = mean(sum, na.rm = T),
          max = mean(max, na.rm = T),
          sd = sqrt(mean(sd^2)),
          skew = mean(skew, na.rm=T)
          ) %>%
        mutate(variable = "Observed") %>%
        relocate(variable) %>%
        ungroup()
    
  # colnames(obsMM)[-1] = paste0("obs_", colnames(obsMM)[-1])
  
  }
  
  df.comb = rbind(obsWW, simWW)
  
  #plotting --------------------------------
  if(Tag == "Temp"){
    
    p1 = ggplot(df.comb) + 
      geom_boxplot(data = subset(df.comb, variable != "Observed"), aes(x = week, y = mean, group = week)) +
      geom_line(data = subset(df.comb, variable == "Observed"), size = 0.5, aes(x = week, y = mean, color = "Observed")) +
      geom_point(data = subset(df.comb, variable == "Observed"), size = 1.5, aes(x = week, y = mean, color = "Observed")) +
      xlab("Week") + ylab("Temperature (°F)") +
      theme_classic() +
      theme(axis.title = element_text(face = "bold"),
            text=element_text(size=12),
            panel.grid.major = element_line(),
            legend.title=element_blank(),
            plot.title = element_text(size=10)
            ) +
      scale_x_continuous(breaks = seq(1,52,2)) +
      ggtitle("Average Mean Weekly Temperature")
  
  }else if(Tag == "Precip"){
    
     p1 = ggplot(df.comb) + 
      geom_boxplot(data = subset(df.comb, variable != "Observed"), aes(x = week, y = sum, group = week)) +
      geom_line(data = subset(df.comb, variable == "Observed"), size = 0.5, aes(x = week, y = sum, color = "Observed")) +
      geom_point(data = subset(df.comb, variable == "Observed"), size = 1.5, aes(x = week, y = sum, color = "Observed")) +
      xlab("Week") + ylab("Precipitation (inches)") +
      theme_classic() +
      theme(axis.title = element_text(face = "bold"),
            text=element_text(size=12),
            panel.grid.major = element_line(),
            legend.title=element_blank(),
            plot.title = element_text(size=10)
            ) +
      scale_x_continuous(breaks = seq(1,52,2)) +
      ggtitle("Average Total Weekly Precipitation")
    
  }
  
  if(Tag == "Temp"){
    yLabel = "Temperature "
    units = "(°F)"
  } else if(Tag == "Precip"){
    yLabel = "Precipitation "
    units = "(inches)"
  }
  
  #weekly SD
  p2 = ggplot(df.comb) + 
      geom_boxplot(data = subset(df.comb, variable != "Observed"), aes(x = week, y = sd, group = week)) +
      geom_line(data = subset(df.comb, variable == "Observed"), size = 0.5, aes(x = week, y = sd, color = "Observed")) +
      geom_point(data = subset(df.comb, variable == "Observed"), size = 1.5, aes(x = week, y = sd, color = "Observed")) +
      xlab("Week") + ylab(paste0("Standard Deviation ", units)) +
      theme_classic() +
      theme(axis.title = element_text(face = "bold"),
            text=element_text(size=12),
            panel.grid.major = element_line(),
            legend.title=element_blank(),
            plot.title = element_text(size=10)
            ) +
      scale_x_continuous(breaks = seq(1,52,2)) +
      ggtitle(paste0("Average Standard Deviation in Weekly ", yLabel))
  
  
  #weekly Skew
  p3 = ggplot(df.comb) + 
      geom_boxplot(data = subset(df.comb, variable != "Observed"), aes(x = week, y = skew, group = week)) +
      geom_line(data = subset(df.comb, variable == "Observed"), size = 0.5, aes(x = week, y = skew, color = "Observed")) +
      geom_point(data = subset(df.comb, variable == "Observed"), size = 1.5, aes(x = week, y = skew, color = "Observed")) +
      xlab("Week") + ylab("Skew (-)") +
      theme_classic() +
      theme(axis.title = element_text(face = "bold"),
            text=element_text(size=12),
            panel.grid.major = element_line(),
            legend.title=element_blank(),
            plot.title = element_text(size=10)
            ) +
      scale_x_continuous(breaks = seq(1,52,2)) +
      ggtitle(paste0("Average Skew in Weekly ", yLabel))
  
  
  #weekly max
  p4 = ggplot(df.comb) + 
      geom_boxplot(data = subset(df.comb, variable != "Observed"), aes(x = week, y = max, group = week)) +
      geom_line(data = subset(df.comb, variable == "Observed"), size = 0.5, aes(x = week, y = max, color = "Observed")) +
      geom_point(data = subset(df.comb, variable == "Observed"), size = 1.5, aes(x = week, y = max, color = "Observed")) +
      xlab("Week") + ylab(paste0("Maximum ", units)) +
      theme_classic() +
      theme(axis.title = element_text(face = "bold"),
            text=element_text(size=12),
            panel.grid.major = element_line(),
            legend.title=element_blank(),
            plot.title = element_text(size=10)
            ) +
      scale_x_continuous(breaks = seq(1,52,2)) +
      ggtitle(paste0("Average Weekly Maximum ", yLabel))
  
  p.comb = ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2, common.legend = TRUE, legend = "bottom")
  
  print(p.comb)
  
  # p.out = paste0(tempdir(), "/outputPlots/weeklyStats_", Tag, ".png")
  # ggsave(filename = p.out, plot = p.comb, device = "png", height = 8, width = 10, units = "in")

}


## ----fig.width=10, fig.height=8-----------------------------------------------
weeklyPlot(sim.pcp, obs.pcp, "Precip")

## ----fig.width=10, fig.height=8-----------------------------------------------
weeklyPlot(sim.tmp, obs.tmp, "Temp")

## -----------------------------------------------------------------------------
# dat.d = z$dat.d
# X = z$Xpamt
# syr = head(dat.d$year, 1)
# eyr = tail(dat.d$year, 1)
# wLO=0.05; wHI=0.95     #Set whisker percentile for boxplots

spellstats <- function(dat.d,X,syr,eyr,nrealz,nsim,wLO,wHI){
  #get these variables after running the driver_wx#.R code
  #dat.d <- z$dat.d
  #X <- z$X
  uyr=syr:eyr
  nyr=length(uyr)

  #Training Data
  Tdat=LowerSantaCruzRiverBasinAZ %>%
    mutate(occ = if_else(prcp>=0.01, 1, 0))
  
  #### END DATA PREPARATION BLOCK TO RUN STATS CODE BELOW ####
  
  #to get month sequence
  nobs=length(unique(Tdat$year))
  lpyear=dat.d$year[min(which(lubridate::leap_year(dat.d$year)))]
  aday <- ymd(paste(lpyear,1,1,sep="-")) #jan 1 of a leap year to have a 366-day year
  it1=which(dat.d$date==aday)
  it2=it1+366-1
  jdaymth <- dat.d$month[it1:it2]
  zz <- rep(jdaymth,nsim)
  yy <- rep(1:nsim,each=366)
  X1 <- cbind(yy,zz,X)
  
  #get spell length stats 3 statistics - mean, var and max
  Y <- array(NA,dim=c(nrealz,3,2,nsim,12)) #save dry and wet spell lengths by sim yr
  W <- array(NA,dim=c(3,2,nobs,12))        #save dry and wet spell lengths by obs yr
  Z <- array(NA,dim=c(3,2,nrealz,12))      #save average spell length by realz
  V <- array(NA,dim=c(3,2,12))             #save average spell length for obs
  
  fidx=seq(1,dim(X1)[1],by=366)
  
  for (irealz in 1:nrealz){
    for (isim in 1:nsim){
      i1=fidx[isim]
      i2=i1+366-1
      for (imth in 1:12){
        idxlist=i1 + which(X1[i1:i2,2]==imth) - 1
        s <- na.omit(X[idxlist,irealz])
        z.f <- spellLengths(s)
        if (length(z.f$`0`)>0){
          Y[irealz,1,1,isim,imth]=mean(z.f$`0`)
          Y[irealz,2,1,isim,imth]=var(z.f$`0`)
          Y[irealz,3,1,isim,imth]=max(z.f$`0`)
        }
        if (length(z.f$`1`)>0){
          Y[irealz,1,2,isim,imth]=mean(z.f$`1`)
          Y[irealz,2,2,isim,imth]=var(z.f$`1`)
          Y[irealz,3,2,isim,imth]=max(z.f$`1`)
        }
      } #imth
    }#isim
  } #irealz
  
  for (isim in 1:nobs){
    for (imth in 1:12){
      s <- dplyr::filter(Tdat, year==unique(Tdat$year)[isim], month==imth)$occ
      z.f <- spellLengths(s)
      if (length(z.f$`0`)>0){
        W[1,1,isim,imth]=mean(z.f$`0`)
        W[2,1,isim,imth]=var(z.f$`0`)
        W[3,1,isim,imth]=max(z.f$`0`)
      }
      if (length(z.f$`1`)>0){
        W[1,2,isim,imth]=mean(z.f$`1`)
        W[2,2,isim,imth]=var(z.f$`1`)
        W[3,2,isim,imth]=max(z.f$`1`)
      }
    } #imth
  }#isim
  
  for (imth in 1:12){
    Z[1,1,,imth]=apply(Y[,1,1,,imth],1,"mean",na.rm=T)
    Z[2,1,,imth]=apply(Y[,2,1,,imth],1,"mean",na.rm=T)
    Z[3,1,,imth]=apply(Y[,3,1,,imth],1,"mean",na.rm=T)
    Z[1,2,,imth]=apply(Y[,1,2,,imth],1,"mean",na.rm=T)
    Z[2,2,,imth]=apply(Y[,2,2,,imth],1,"mean",na.rm=T)
    Z[3,2,,imth]=apply(Y[,3,2,,imth],1,"mean",na.rm=T)
    V[1,1,imth]=mean(W[1,1,,imth],na.rm=T)
    V[2,1,imth]=mean(W[2,1,,imth],na.rm=T)
    V[3,1,imth]=mean(W[3,1,,imth],na.rm=T)
    V[1,2,imth]=mean(W[1,2,,imth],na.rm=T)
    V[2,2,imth]=mean(W[2,2,,imth],na.rm=T)
    V[3,2,imth]=mean(W[3,2,,imth],na.rm=T)
  } #imth
  
#Boxplots  
  d01 <- as.data.frame(Z[1,1,,])  #average mean dry spell length
  d02 <- as.data.frame(sqrt(Z[2,1,,]))  #average sd dry spell length
  d03 <- as.data.frame(Z[3,1,,])  #average max dry spell length
  
  d01T=V[1,1,]
  d02T=sqrt(V[2,1,])
  d03T=V[3,1,]

  RecRed = "red"
  RecBlue = "blue"
  
  # pdf(file=paste0(tempdir(), "/outputPlots/DryWetStats.pdf"), width=9, height=4)
  oldpar = par(mfrow=c(1,3), mar=c(2,2.5,2,1),
               oma=c(2,2,0,0), mgp=c(2,1,0), cex.axis=0.8)
  
  par(mfrow=c(1,3), mar=c(2,2.5,2,1), oma=c(2,2,0,0), mgp=c(2,1,0),cex.axis=0.8)
  #Mean
  #ylimit=range(c(d01,d01hist),na.rm=T)
  bb=boxplot(d01, plot=F, na.rm=T, names=1:12)
  ylimit=c(range(c(d01,d01T),na.rm=T))
  out=matrix(nrow=nsim, ncol=12)
  for(b in 1:12){
    x=d01[,b]
    quants=quantile(x, c(wLO,wHI),na.rm=T)
    bb$stats[c(1,5),b] = quants
    outs=which(x < quants[1] | x > quants[2])
    out[1:length(outs), b]= x[outs]
  }
  bxp(bb, ylim=ylimit,na.rm=T, outline=F,xlab="",ylab="")
  mtext("days", side=2, outer = T)
  title(main="average mean dry spell length")
  for(m in 1:12){points(rep(m, length(which(is.na(out[,m])==F))), out[!is.na(out[,m]),m])}
  points(1:12, d01T, pch=17, cex=1, col=RecRed)
  lines(1:12, d01T, col=RecRed)

  #Standard Deviation
  bb=boxplot(d02, plot=F, na.rm=T, names=1:12)
  ylimit=c(range(c(d02,d02T),na.rm=T))
  out=matrix(nrow=nsim, ncol=12)
  for(b in 1:12){
    x=d02[,b]
    quants=quantile(x, c(wLO,wHI),na.rm=T)
    bb$stats[c(1,5),b] = quants
    outs=which(x < quants[1] | x > quants[2])
    out[1:length(outs), b]= x[outs]
  }
  bxp(bb, ylim=ylimit,na.rm=T, outline=F,xlab="",ylab="")
  title(main="average sd dry spell length")
  mtext("month", side=1, outer = T, line=0.5)
  for(m in 1:12){points(rep(m, length(which(is.na(out[,m])==F))), out[!is.na(out[,m]),m])}
  points(1:12, d02T, pch=17, cex=1, col=RecRed)
  lines(1:12, d02T, col=RecRed)

  #Max Length
  bb=boxplot(d03, plot=F, na.rm=T, names=1:12)
  ylimit=c(range(c(d03,d03T),na.rm=T))
  out=matrix(nrow=nsim, ncol=12)
  for(b in 1:12){
    x=d03[,b]
    quants=quantile(x, c(wLO,wHI),na.rm=T)
    bb$stats[c(1,5),b] = quants
    outs=which(x < quants[1] | x > quants[2])
    out[1:length(outs), b]= x[outs]
  }
  bxp(bb, ylim=ylimit,xlab="",ylab="",na.rm=T, outline=F)
  title(main="average max dry spell length")
  for(m in 1:12){points(rep(m, length(which(is.na(out[,m])==F))), out[!is.na(out[,m]),m])}
  points(1:12, d03T, pch=17, cex=1, col=RecRed)
  lines(1:12, d03T, col=RecRed)
  
  
  dev.off()
  
  par(oldpar)
  ##########################################
}

## ----fig.width=8, fig.height=6------------------------------------------------
spellstats(dat.d = z$dat.d, X = z$Xpamt,
           syr = head(dat.d$year, 1),
           eyr = tail(dat.d$year, 1),
           nrealz = nrealz, nsim = nsim,
           wLO = 0.05, wHI = 0.95)

## -----------------------------------------------------------------------------

# setwd() to desired location for writeSim to save .csv files containing the simulated precipitation and temperature
# setwd(tempdir())
# 
# writeSim(wxOutput = z, nsim = nsim, nrealz = nrealz, debug = TRUE)



## -----------------------------------------------------------------------------
  #wxgenR weather generation run time:
  print(difftime(endTime, startTime, units='mins'))

