# wxgenR: Weather Generator for Simulating Precipitation and Temperature

**wxgenR** is a novel weather generator designed to simulate ensembles of precipitation and temperature data for regions with distinct seasonality. By incorporating up to 26 seasons as training variables, it offers users a robust scenario-planning tool for analyzing potential impacts of phenomena like climate change on seasonality.

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
1. **Load the package**:
   ```R
   library(wxgenR)
   ```
2. **Prepare your data**: Input data must include columns for year, month, day, precipitation (`prcp`), temperature (`temp`), and season (`season`).
3. **Run the weather generator**:
   ```R
   result = wx(trainingData = myData, syr = 1990, eyr = 2020, nsim = 5, nrealz = 50)
   ```
4. **Analyze the outputs**: Use built-in or custom scripts to visualize and interpret simulated results.

### Sample Data
The package includes example data from the Blacksburg, VA NWS station for easy testing and learning. Access it with:
```R
data(BlacksburgVA)
```

---

## Usage Example

Here’s how to simulate 10 years of weather data with a 20-trace ensemble:
```R
library(wxgenR)
data(BlacksburgVA)

simulatedData = wx(
  trainingData = BlacksburgVA, 
  syr = 1991, eyr = 2020, 
  nsim = 10, nrealz = 20, 
  wwidth = 7, tempPerturb = TRUE, ekflag = TRUE
)

# Access simulated precipitation and temperature
simulatedPrecip = simulatedData$Xpamt
simulatedTemp = simulatedData$Xtemp
```

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
- **`wwidth`**: Sampling window width for daily resampling (default: 7).
- **`tempPerturb`**: Adds temperature variability using monthly residuals (default: FALSE).

---

## Visualization

### Example: Daily Precipitation and Temperature
```R
library(ggplot2)

# Plotting daily precipitation
ggplot(simulatedPrecip, aes(x = Date, y = Trace_1)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Simulated Daily Precipitation", y = "Precipitation (inches)")

# Plotting daily temperature
ggplot(simulatedTemp, aes(x = Date, y = Trace_1)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Simulated Daily Temperature", y = "Temperature (°F)")
```

---

## Performance

Using 5 years of training data on a standard laptop:
- **Simulation length**: 5 years
- **Ensemble size**: 10 traces
- **Run time**: ~4 seconds with parallel computing enabled (`numbCores = 11`).

---

## References

For theoretical background and applications, refer to:

1. Rajagopalan, B., et al. (1996). *Nonhomogeneous Markov Model for Daily Precipitation*.  
2. Verdin, A., et al. (2015, 2018). *Coupled stochastic weather generation*.

---

## Contributing

We welcome contributions! To contribute:
1. Fork the repository.
2. Create a new branch for your feature/fix.
3. Submit a pull request.

---

## License

This project is licensed under the MIT License.

For more details, visit our [GitHub Repository](https://github.com/dwoodson-usbr/wxgenR).
