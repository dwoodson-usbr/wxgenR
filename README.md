# wxgenR: A Stochastic Weather Generator with Seasonality

A weather generator is a numerical tool that resamples an input timeseries many times, while preserving observed or projected characteristics of importance, such as the statistics of the transition between wet and dry days. The resulting large group, or ensemble, of likely rainfall and temperature timeseries represents a range of possible amounts, daily patterns, and seasonality. This weather generator is, to our knowledge, novel in that it includes seasons (up to 26) in training the simulation algorithm.

The goal of wxgenR is to provide users a tool that can simulate, with fidelity, an ensemble of precipitation and temperature based on training data that could include, for example, station based point measurements, grid cell values derived from models or remotely sensed data, or basin averages. The incorporation of seasonality as a covariate in the training algorithm allows users to examine the effects of shifts in seasonality due to climate warming (e.g., earlier snowmelt seasons or prolonged summer dry periods). wxgenR is an effective and robust scenario planning tool for a wide variety of potential futures.

## Features
- **Seasonality Support:** Incorporates up to 26 seasons, allowing for detailed examination of shifts in season start, end, and duration.
- **Flexible Input Data:** Handles single time series from sources such as station measurements, basin averages, grid cells, or model outputs.
- **Simulation Outputs:** Provides an ensemble of precipitation and temperature traces, preserving essential statistical characteristics of wet/dry day transitions.
- **Customizable Settings:** Includes configurable options like sampling window width, adaptive sampling, and temperature perturbation.

---

## Installation

### From CRAN
```R
install.packages("wxgenR")
```

### From GitHub
1. Clone the repository or download the code:
   ```bash
   git clone https://github.com/dwoodson-usbr/wxgenR.git
   ```
2. Install the package locally:
   ```R
   devtools::install_github("dwoodson-usbr/wxgenR")
   ```

---

## Getting Started

### Quickstart
1. **Load the packages**:
   ```R
   library(wxgenR)
   library(lubridate)
   library(dplyr)
   library(tidyr)
   library(reshape2)
   library(ggpubr)
   library(data.table)
   library(moments)
   library(seas)
   ```
2. **Prepare your data**: Input data must include columns for year, month, day, precipitation (`prcp`), temperature (`temp`), and season (`season`).

The package includes example data from the Blacksburg, VA NWS station for easy testing and learning. Access it with:
```R
data(BlacksburgVA)
```

3. **Run the weather generator**:

   ```R
   nsim = 20   #number of simulation years
   nrealz = 30 #number of traces in ensemble
      
   z = wx(trainingData = BlacksburgVA, syr = 2001, eyr = 2020,
          smo = 1, emo = 12,
          nsim = nsim, nrealz = nrealz, aseed = 123,
          wwidth = 7, unitSystem = "Metric", ekflag = TRUE,
          awinFlag = FALSE, tempPerturb = TRUE, pcpOccFlag = FALSE,
          numbCores = 2)
   
   # Access simulated precipitation and temperature
   simulatedPrecip = z$Xpamt
   simulatedTemp = z$Xtemp
   ```
   
4. **Analyze the outputs**: Use built-in or custom scripts to visualize and interpret simulated results.

---

## Detailed Documentation

### Input Data Requirements
- **Format**: Single time series in a dataframe or `.csv` file.
- **Columns**:
  - `year`, `month`, `day`: Calendar date information.
  - `prcp`: Daily precipitation values.
  - `temp`: Daily temperature values.
  - `season`: Numeric indices representing seasons (e.g., 1 = winter, 2 = spring).

### Key Parameters
- **`nsim`**: Length of simulation in years.
- **`nrealz`**: Number of ensemble traces to simulate.
- **`wwidth`**: Sampling window width for daily resampling.
- **`tempPerturb`**: Adds temperature variability using monthly residuals.
- **`ekflag`**: Enables precipitation pertubation via Epanechnikov Kernel.

---

## Visualization

### Example: Daily Precipitation and Temperature

First, use modified approach from writeSim function to post-process/format output.

```R
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


```
Format dataframes for simulated precip, temperature, and season.

```R
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
    reshape2::melt(id = 1:6)
    
  return(df)
}
```

Format training data
```R
sim.pcp = formatting(sim.pcp)
sim.tmp = formatting(sim.tmp)
sim.szn = formatting(sim.szn)

colnames(dat.d)[11] = "yday"
obs.pcp = dat.d[,c(1:3,8:9,11,4)]
obs.tmp = dat.d[,c(1:3,8:9,11,5)]
```

Now we plot the daily results.
```R
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
  geom_ribbon(aes(x = yday, ymin = mean_q5, ymax = mean_q95), alpha = 0.25) +
  geom_line(aes(x = yday, y = mean_med, color = "red"), size = 1, alpha = 0.8) +
  geom_line(aes(x = yday, y = obs_mean), size = 0.3, alpha = trnAlpha, linetype = "solid", color = "blue") +
  geom_point(aes(x = yday, y = obs_mean), size = 0.6, alpha = trnAlpha, color = "blue") +
  scale_colour_manual(values =c('blue'='blue','red'='red', 'grey' = 'grey'), labels = c('Training Data','Simulation Median', '95% Confidence')) +
  theme_classic() +
  theme(axis.title = element_text(face = "bold"),
        # text=element_text(size=14),
        panel.grid.major = element_line(),
        legend.title=element_blank(),
        legend.position = lgdLoc,
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank()) +
  xlab("Day of Year") + ylab(paste0("Mean ", yLabel, units))   
  
  #daily SD
  p2 = ggplot(df.comb) +
  geom_ribbon(aes(x = yday, ymin = sd_q5, ymax = sd_q95), alpha = 0.25) +
  geom_line(aes(x = yday, y = sd_med, color = "red"), size = 1, alpha = 0.8) +
  geom_line(aes(x = yday, y = obs_sd), size = 0.3, alpha = trnAlpha, linetype = "solid", color = "blue") +
  geom_point(aes(x = yday, y = obs_sd), size = 0.6, alpha = trnAlpha, color = "blue") +
  scale_colour_manual(values =c('blue'='blue','red'='red', 'grey' = 'grey'), labels = c('Training Data','Simulation Median', '95% Confidence')) +
  theme_classic() +
  theme(axis.title = element_text(face = "bold"),
        # text=element_text(size=14),
        panel.grid.major = element_line(),
        legend.title=element_blank(),
        legend.position = lgdLoc,
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank()) +
  xlab("Day of Year") + ylab(paste0("Std. Deviation of ", yLabel, units))   
  
  #daily skew
  p3 = ggplot(df.comb) +
  geom_ribbon(aes(x = yday, ymin = skew_q5, ymax = skew_q95), alpha = 0.25) +
  geom_line(aes(x = yday, y = skew_med, color = "red"), size = 1, alpha = 0.8) +
  geom_line(aes(x = yday, y = obs_skew), size = 0.3, alpha = trnAlpha, linetype = "solid", color = "blue") +
  geom_point(aes(x = yday, y = obs_skew), size = 0.6, alpha = trnAlpha, color = "blue") +
  scale_colour_manual(values =c('blue'='blue','red'='red', 'grey' = 'grey'), labels = c('Training Data','Simulation Median', '95% Confidence')) +
  theme_classic() +
  theme(axis.title = element_text(face = "bold"),
        # text=element_text(size=14),
        panel.grid.major = element_line(),
        legend.title=element_blank(),
        legend.position = lgdLoc,
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank()) +
  xlab("Day of Year") + ylab(paste0("Skew of ", yLabel, " (-)"))  
    
  #daily Max
  p4 = ggplot(df.comb) +
  geom_ribbon(aes(x = yday, ymin = max_q5, ymax = max_q95), alpha = 0.25) +
  geom_line(aes(x = yday, y = max_med, color = "red"), size = 1, alpha = 0.8) +
  geom_line(aes(x = yday, y = obs_max), size = 0.3, alpha = trnAlpha, linetype = "solid", color = "blue") +
  geom_point(aes(x = yday, y = obs_max), size = 0.6, alpha = trnAlpha, color = "blue") +
  scale_colour_manual(values =c('blue'='blue','red'='red', 'grey' = 'grey'), labels = c('Training Data','Simulation Median', '95% Confidence')) +
  theme_classic() +
  theme(axis.title = element_text(face = "bold"),
        # text=element_text(size=14),
        panel.grid.major = element_line(),
        legend.title=element_blank(),
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

```

Plot precip and temperature.
```R
dailyPlot(sim.pcp, obs.pcp, "Precip")
dailyPlot(sim.tmp, obs.tmp, "Temp")
```

We can also plot monthly stats:
```R

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
  
  p.out = paste0(tempdir(), "/outputPlots/monthlyStats_", Tag, ".png")
  # ggsave(filename = p.out, plot = p.comb, device = "png", height = 8, width = 8, units = "in")

}
```

```R
monthlyPlot(sim.pcp, obs.pcp, "Precip")
monthlyPlot(sim.tmp, obs.tmp, "Temp")
```

---

## References

For theoretical background and applications, refer to:

Bearup, L., Gangopadhyay, S., & Mikkelson, K. (2021). Hydroclimate Analysis Lower Santa Cruz River Basin Study (Technical Memorandum No ENV-2020-056). Bureau of Reclamation. <https://www.usbr.gov/lc/phoenix/programs/lscrbasin/LSCRBStudy.html>

Gangopadhyay, S., Bearup, L. A., Verdin, A., Pruitt, T., Halper, E., & Shamir, E. (2019, December 1). A collaborative stochastic weather generator for climate impacts assessment in the Lower Santa Cruz River Basin, Arizona. Fall Meeting 2019, American Geophysical Union. <https://ui.adsabs.harvard.edu/abs/2019AGUFMGC41G1267G>

Rajagopalan, B., Lall, U., and Tarboton, D. G.: Nonhomogeneous Markov Model for Daily Precipitation, Journal of Hydrologic Engineering, 1, 33–40, <https://doi.org/10.1061/(ASCE)1084-0699(1996)1:1(33)>, 1996.

Verdin, A., Rajagopalan, B., Kleiber, W., and Katz, R. W.: Coupled stochastic weather generation using spatial and generalized linear models, Stoch Environ Res Risk Assess, 29, 347–356, <https://doi.org/10.1007/s00477-014-0911-6>, 2015.

Verdin, A., Rajagopalan, B., Kleiber, W., Podestá, G., and Bert, F.: A conditional stochastic weather generator for seasonal to multi-decadal simulations, Journal of Hydrology, 556, 835–846, <https://doi.org/10.1016/j.jhydrol.2015.12.036>, 2018.

---

## Contributing

We welcome contributions! To contribute:
1. Fork the repository.
2. Create a new branch for your feature/fix.
3. Submit a pull request.

---

## Disclaimer

This information is preliminary and is subject to revision. It is being provided to meet the need for timely best science. The information is provided on the condition that neither the U.S. Bureau of Reclamation nor the U.S. Government may be held liable for any damages resulting from the authorized or unauthorized use of the information.

## License

This project is licensed under the CC0 1.0 Universal.

For more details, visit our [GitHub Repository](https://github.com/dwoodson-usbr/wxgenR).
