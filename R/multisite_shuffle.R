#' Multisite shuffling of wxgenR simulation results
#'
#' Runs a postprocessor on `wx` simulation results to shuffle
#' multiple stations' simulations such that wxgenR can be used in multisite applications.
#' Specifically, the `multisite_shuffle` function uses an approach developed by Iman and Conover (1982) and later applied by Clark et al. (2004)
#' to capture the rank correlation among observed station data and introduce it to those stations' simulations. The Iman and Conover approach is implemented
#' using the `cornode` function from the `mc2d` package.
#'
#' @param list.sta.wx A list containing multiple stations' observed and simulated data to be used in multisite shuffling.
#' Each list element should be named for the station's data it holds and should contain the dataframe output from the `wx` for that station.
#' @param numbCores Enable parallel computing for multisite shuffling,
#'  set number of cores to enable (must be a positive integer greater than or equal to 2).
#'   Turned off by default; if set to 0 or 1 it will run as single thread.
#'    Use function 'detectCores()' from 'parallel' package to show the number of available cores on your machine.
#' @param aseed Specify a seed for reproducibility.
#'
#' @return Returns a list containing results and metadata from the multisite shuffling in 'long' format for easy analysis and visualization.
#' \itemize{
#'   \item shuffledResultsOnly - Dataframe containing shuffled simulations for all traces and stations ('sim_prcp' and 'sim_temp').
#'   \item shuffledResultsAndObs - Dataframe containing shuffled simulations as well as corresponding observations/training data ('prcp' and 'temp' are the observed data).
#'   \item shuffledAndUnshuffled - Dataframe containing shuffled simulations, unshuffled simulations, observations. Unshuffled vs and shuffled are indicated by the 'Tag' variable.
#' }
#'
#' @examples
#' \donttest{
#' # Simulated example with two stations
#'
#' data(BOCO_sims)
#'
#' ms = multisite_shuffle(BOCO_sims, numbCores = 2, aseed = 123)
#'
#' print(ms)
#' }
#'
#' @references
#' Iman, Ronald & Conover, William. (1982). A Distribution-Free Approach to Inducing Rank Correlation Among Input Variates. Communications in Statistics-simulation and Computation - COMMUN STATIST-SIMULAT COMPUT. 11. 311-334. 10.1080/03610918208812265.
#'
#' Clark, M., Gangopadhyay, S., Hay, L., Rajagopalan, B., & Wilby, R. (2004). The Schaake Shuffle: A Method for Reconstructing Space–Time Variability in Forecasted Precipitation and Temperature Fields. Journal of Hydrometeorology, 5(1), 243-262. https://doi.org/10.1175/1525-7541(2004)005<0243:TSSAMF>2.0.CO;2
#'
#' R. Pouillot, M.-L. Delignette-Muller (2010), Evaluating variability and uncertainty in microbial quantitative risk assessment using two R packages. International Journal of Food Microbiology. 142(3):330-40
#'
#' @export
#'
#' @importFrom dplyr mutate filter bind_rows bind_cols arrange left_join group_by summarise select case_when
#' @importFrom tidyr pivot_longer pivot_wider separate starts_with
#' @importFrom lubridate ymd leap_year
#' @importFrom magrittr %>%
#' @importFrom mc2d cornode
#' @importFrom foreach foreach %dopar%
#' @importFrom doRNG registerDoRNG
#' @import parallel
#' @import doParallel
#'
#'

"multisite_shuffle" = function(list.sta.wx, numbCores = NULL, aseed = NULL){

  ### extract data list of stations' individual simulations
  stations = names(list.sta.wx)

  df.combo = NULL
  df.hist = NULL

  sta = stations[1]
  for(sta in stations){

    df = list.sta.wx[[sta]]

    data = df$dat.d
    psim = df$Xpamt
    tsim = df$Xtemp

    monthly_climatology = data %>%
      group_by(month) %>%
      summarise(monthly_temp = mean(temp, na.rm = TRUE), .groups = "drop")

    # Extract unique years from data
    years = unique(data$year)

    # Identify leap and non-leap years
    leap_years = years[leap_year(years)]
    non_leap_years = years[!leap_year(years)]

    # Create the extra NA row at the end of each non-leap year
    extra_na_rows = data.frame(
      year = non_leap_years,
      month = 12,
      day = 32,  # Use an artificial placeholder day
      prcp = NA,
      temp = NA,
      season = case_when(
        12 %in% c(12, 1, 2) ~ 1,
        12 %in% c(3, 4, 5) ~ 2,
        12 %in% c(6, 7, 8) ~ 3,
        12 %in% c(9, 10, 11) ~ 4
      )
    )

    # Merge with original data
    data_fixed = data[,c("year","month","day","prcp","temp","season")] %>%
      bind_rows(extra_na_rows) %>%
      arrange(year, month, day)

    # Step 1: Rename psim and tsim columns
    colnames(psim) = paste0("sim_prcp_", seq_len(ncol(psim)))
    colnames(tsim) = paste0("sim_temp_", seq_len(ncol(tsim)))

    # Bind precipitation simulations
    df.c.p = cbind(data_fixed, psim)

    # Bind temperature simulations
    df.c.t = cbind(data_fixed, tsim)

    # Convert precipitation simulations to long format
    df.prcp.long = df.c.p %>%
      pivot_longer(cols = starts_with("sim_prcp_"),
                   names_to = "simulation",
                   values_to = "sim_prcp") %>%
      mutate(simulation = gsub("sim_prcp_", "sim_", simulation))  # Rename for clarity

    # Convert temperature simulations to long format
    df.temp.long = df.c.t %>%
      pivot_longer(cols = starts_with("sim_temp_"),
                   names_to = "simulation",
                   values_to = "sim_temp") %>%
      mutate(simulation = gsub("sim_temp_", "sim_", simulation))  # Ensure consistency

    df.long = left_join(df.prcp.long, df.temp.long,
                        by = c("year", "month", "day", "prcp", "temp", "season", "simulation")) %>%
      left_join(monthly_climatology, by = "month") %>%
      mutate(temp = ifelse(is.na(temp) | is.nan(temp), monthly_temp, temp),
             sim_temp = ifelse(is.na(sim_temp) | is.nan(sim_temp), monthly_temp, sim_temp),
             prcp = ifelse(is.na(prcp) | is.nan(prcp), 0, prcp),  # Replace missing precipitation with 0
             sim_prcp = ifelse(is.na(sim_prcp) | is.nan(sim_prcp), 0, sim_prcp)  # Replace missing precipitation with 0
      ) %>%
      dplyr::select(-monthly_temp) %>%  # Remove helper column
      filter(!(month == 12 & day == 32)) %>% #remove dummy date for non-leap years
      mutate(date = ymd(paste(year,month,day, sep = "-")),
             station = sta
      )

    # Count missing values
    missing_temp = sum(is.na(df.long$temp) | is.nan(df.long$temp))
    missing_sim_temp = sum(is.na(df.long$sim_temp) | is.nan(df.long$sim_temp))
    missing_prcp = sum(is.na(df.long$prcp) | is.nan(df.long$prcp))
    missing_sim_prcp = sum(is.na(df.long$sim_prcp) | is.nan(df.long$sim_prcp))

    # Print warnings if any missing values are found
    if (missing_temp > 0 | missing_sim_temp > 0 | missing_prcp > 0 | missing_sim_prcp > 0) {
      warning("Missing values found:\n",
              if (missing_temp > 0) paste("temp:", missing_temp, "\n") else "",
              if (missing_sim_temp > 0) paste("sim_temp:", missing_sim_temp, "\n") else "",
              if (missing_prcp > 0) paste("prcp:", missing_prcp, "\n") else "",
              if (missing_sim_prcp > 0) paste("sim_prcp:", missing_sim_prcp, "\n") else "")
    } else {
      message(paste0(sta, " - All checks passed: No missing values in temp, sim_temp, prcp, or sim_prcp."))
    }

    df.combo = rbind(df.combo, df.long)

    # Extract relevant columns
    df.hist.sta = filter(df.long, simulation == "sim_1") %>%
      # dplyr::select(prcp)
      dplyr::select(prcp, temp)
    # df.hist.sta = data[, c("prcp", "temp")]

    # Rename columns with station name
    colnames(df.hist.sta) = paste0(colnames(df.hist.sta), ".", sta) #c(paste0("prcp.", sta), paste0("temp.", sta))

    # Bind to df.hist (first iteration assigns, subsequent iterations bind columns)
    if (is.null(df.hist)) {

      metaData = df.long %>%
        filter(simulation == "sim_1") %>%
        dplyr::select(date, year, month)

      df.hist = bind_cols(metaData, df.hist.sta)

    } else {
      df.hist = bind_cols(df.hist, df.hist.sta)
    }

  }

  df.combo$Tag = "Unshuffled"

  ### calculate historical correlation between stations -----------------------------------
  list.hist.cor = list()

  yr = 1992
  for(yr in years){

    list.hist.cor[[as.character(yr)]] = list()  # Initialize year-level list

    mo = 2
    for(mo in 1:12){

      df.hist.mo = df.hist %>%
        filter(year == yr & month == mo)

      hist.cor = cor(dplyr::select(df.hist.mo, -c("date", "year", "month")), method = "spearman", use = "complete.obs")

      list.hist.cor[[as.character(yr)]][[as.character(mo)]] = hist.cor  # Assign correlation matrix
      # list.hist.cor[[yr]][[mo]] = hist.cor
    }

  }

  ### begin shuffling procedure ------------------------------------------------------------
  if (is.null(numbCores) || !is.numeric(numbCores) || numbCores < 2) {
    numbCores = 1
  }
  ### single thread version
  if(numbCores == 1){
    #Pull out sim data for all stations by trace, then shuffle, repeat for next trace
    simz = unique(df.combo$simulation)
    df.shuff = NULL
    list.shuff.cor = list()

    sim = simz[1]
    for(sim in simz){
      print(sim)
      list.shuff.cor[[sim]] = list()  # Initialize year-level list

      df.sim = df.combo %>%
        filter(simulation == sim)

      #shuffle entire time series -----------------------------------------------------

      #now loop thru stations to reorder data
      df.combo.sim = NULL
      sta = stations[1]
      for(sta in stations){

        df.sta = df.sim %>%
          filter(station == sta) %>%
          dplyr::select(sim_prcp, sim_temp)

        # Rename columns with station name
        colnames(df.sta) = paste0(colnames(df.sta), ".", sta) #c(paste0("sim_prcp.", sta), paste0("sim_temp.", sta))

        df.combo.sim = bind_cols(df.combo.sim, df.sta)

      }

      metaData = df.sim %>%
        filter(station == sta) %>%
        dplyr::select(date, year, month)

      df.combo.sim = df.combo.sim %>%
        bind_cols(metaData)

      yr = years[1]
      for(yr in years){

        list.shuff.cor[[sim]][[as.character(yr)]] = list()  # Initialize year-level list

        mo = 1
        for(mo in 1:12){

          df.combo.sim.mo = df.combo.sim %>%
            filter(year == yr & month == mo) %>%
            dplyr::select(-c(date, year, month))
          # dplyr::select(starts_with("sim_prcp") | starts_with("sim_temp"))

          df.dates.mo = df.combo.sim %>%
            filter(year == yr & month == mo) %>%
            dplyr::select(date)

          # df.combo.sim.mo.prcp = df.combo.sim.mo %>%
          #   dplyr::select(starts_with("sim_prcp"))

          target.cor = NA
          while(anyNA(target.cor)) {
            target.cor.yr = sample(years, 1)
            target.cor = list.hist.cor[[as.character(target.cor.yr)]][[mo]]
          }

          invisible(capture.output({
            df.shuffle = cornode(as.matrix(df.combo.sim.mo),
                                 target = target.cor, result = TRUE, outrank = F, seed = aseed)
          }))

          invisible(capture.output({
            shuffle.cor = cor(cornode(as.matrix(df.combo.sim.mo),
                                      target = target.cor, result = TRUE, outrank = FALSE, seed = aseed), method = "spearman")
          }))

          df.shuffle = data.frame(date = df.dates.mo, df.shuffle)
          # df.combo.sim.shuff = df.combo.sim.mo

          # for (sta in colnames(df.combo.sim.mo)) {
          #   rank_order = df.shuffle[[sta]]
          #   temp_col = gsub("prcp", "temp", sta)
          #   df.combo.sim.shuff[[sta]] = df.combo.sim.mo[[sta]][rank_order]
          #   df.combo.sim.shuff[[temp_col]] = df.combo.sim.mo[[temp_col]][rank_order]
          # }

          # df.shuffle.long = df.combo.sim.shuff %>%
          df.shuffle.long = df.shuffle %>%
            # dplyr::select(-c("year", "month")) %>%
            pivot_longer(cols = -date, names_to = "variable", values_to = "value") %>%
            separate(variable, into = c("variable", "station"), sep = "\\.") %>%
            pivot_wider(names_from = "variable", values_from = "value") %>%
            mutate(simulation = sim, Tag = "Shuffled") %>%
            arrange(date)

          #store
          df.shuff = bind_rows(df.shuff, df.shuffle.long)
          list.shuff.cor[[sim]][[as.character(yr)]][[as.character(mo)]] = shuffle.cor

        }

      }

    }

  }else if(numbCores >= 2){

    if(numbCores > detectCores()){
      message(paste0("numbCores is set to more than the ", detectCores(), " available cores on your machine. Setting numbCores to ", detectCores()-1, " cores (leaving one core free)."))
      message(paste0("If you wish, interupt code to set a new value for numbCores that is between 2 and ", detectCores(), ", otherwise the simulation will continue."))
    }

    # Set up parallel backend
    cl = makeCluster(numbCores)
    registerDoParallel(cl)
    registerDoRNG(aseed)  # Set seed for reproducibility

    # Unique simulation values
    simz = unique(df.combo$simulation)

    # Parallelized loop
    results = foreach(sim = simz, .combine = function(x, y) {
      list(df.shuff = bind_rows(x$df.shuff, y$df.shuff),
           list.shuff.cor = c(x$list.shuff.cor, y$list.shuff.cor))
    }, .packages = c("dplyr", "tidyr", "mc2d")) %dorng% {

      print(sim)
      list.shuff.cor = list()  # Initialize year-level list
      df.shuff = NULL

      df.sim = df.combo %>%
        filter(simulation == sim)

      # Shuffle entire time series
      df.combo.sim = NULL
      for (sta in stations) {

        df.sta = df.sim %>%
          filter(station == sta) %>%
          dplyr::select(sim_prcp, sim_temp)

        colnames(df.sta) = paste0(colnames(df.sta), ".", sta)

        df.combo.sim = bind_cols(df.combo.sim, df.sta)
      }

      metaData = df.sim %>%
        filter(station == sta) %>%
        dplyr::select(date, year, month)

      df.combo.sim = df.combo.sim %>%
        bind_cols(metaData)

      yr = years[1]
      for (yr in years) {

        list.shuff.cor[[as.character(yr)]] = list()

        mo = 1
        for (mo in 1:12) {

          df.combo.sim.mo = df.combo.sim %>%
            filter(year == yr & month == mo) %>%
            dplyr::select(-c(date, year, month))
          # dplyr::select(starts_with("sim_prcp") | starts_with("sim_temp"))

          df.dates.mo = df.combo.sim %>%
            filter(year == yr & month == mo) %>%
            dplyr::select(date)

          # df.combo.sim.mo.prcp = df.combo.sim.mo %>%
          #   dplyr::select(starts_with("sim_prcp"))

          target.cor = NA
          while(anyNA(target.cor)) {
            target.cor.yr = sample(years, 1)
            target.cor = list.hist.cor[[as.character(target.cor.yr)]][[mo]]
          }

          invisible(capture.output({
            df.shuffle = cornode(as.matrix(df.combo.sim.mo),
                                 target = target.cor, result = TRUE, outrank = F, seed = aseed)
          }))

          invisible(capture.output({
            shuffle.cor = cor(cornode(as.matrix(df.combo.sim.mo),
                                      target = target.cor, result = TRUE, outrank = FALSE, seed = aseed), method = "spearman")
          }))

          df.shuffle = data.frame(date = df.dates.mo, df.shuffle)
          # df.combo.sim.shuff = df.combo.sim.mo

          # for (sta in colnames(df.combo.sim.mo)) {
          #   rank_order = df.shuffle[[sta]]
          #   temp_col = gsub("prcp", "temp", sta)
          #   df.combo.sim.shuff[[sta]] = df.combo.sim.mo[[sta]][rank_order]
          #   df.combo.sim.shuff[[temp_col]] = df.combo.sim.mo[[temp_col]][rank_order]
          # }

          # df.shuffle.long = df.combo.sim.shuff %>%
          df.shuffle.long = df.shuffle %>%
            # dplyr::select(-c("year", "month")) %>%
            pivot_longer(cols = -date, names_to = "variable", values_to = "value") %>%
            separate(variable, into = c("variable", "station"), sep = "\\.") %>%
            pivot_wider(names_from = "variable", values_from = "value") %>%
            mutate(simulation = sim, Tag = "Shuffled") %>%
            arrange(date)

          df.shuff = bind_rows(df.shuff, df.shuffle.long)
          list.shuff.cor[[as.character(yr)]][[as.character(mo)]] = shuffle.cor
        }
      }

      list(df.shuff = df.shuff, list.shuff.cor = list.shuff.cor)
    }

    # Stop parallel backend
    stopCluster(cl)

    # Extract results
    df.shuff = results$df.shuff
    list.shuff.cor = results$list.shuff.cor

  }

  df.shuff = data.frame(df.shuff)

  df.shuff.full = df.shuff %>%
    left_join(
      df.combo %>% dplyr::select(year, month, day, prcp, temp, season, date, station, simulation),
      by = c("date", "station", "simulation")
    ) %>%
    arrange(date, station)

  df.merge = bind_rows(df.combo, df.shuff.full)

  out = list(shuffledResultsOnly = df.shuff,
             shuffledResultsAndObs = df.shuff.full,
             shuffledAndUnshuffled = df.merge
             )

  return(out)

}
