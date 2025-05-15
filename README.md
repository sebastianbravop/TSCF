# TSCF: A statistic-based tool for smoothing temperature sensitivity effects on matric potential dielectric readings

The Temperature Sensitivity Correction Function (TSCF) is introduced as straightforward method to localise, quantify, and smooth the daily-scale periodic noise. This function is based on a dynamic moving average approach using a combination of two centred simple moving average (SMA) with different window size (k). First, a fast SMA (mp_f) with a window size (k1), in this example equal to a hourly period, is used to decrease data noise normally produced in this type of high resolution observations. This magnitude allows for keeping short reactions times produced by fast wetting fronts. The size of k1 can also be set to a different resolution depending on the data resolution and noise intensity (e.g., 3 hours). Secondly, a slow SMA (mp_s) with a window equal to a daily period (k2) is used to smooth the daily oscillation produced by temperature sensitivity effects. The resulting corrected time series is a dynamic combination of fast and slow smoothing routines.

### How does it work?

To produce an algorithm able to recognise and implement a proper combination of these time series, we subtracted the mp_f from the mp_s to obtain a convergence-divergence noise signal (CDNS). This application resembles the Moving Average Convergence-Divergence (MACD) oscillator proposed by Appel (2005) for economic sciences, with the simplification of using the simple unweighted mean (Savitzky & Golay, 1964) instead of the exponential moving average. The CDNS function allows for localising and quantifying the intensity of the noise signal dynamically.

We avoid smoothing important wetting processes during the dry season. When a wetting front reaches the matric potential sensor under lower water content (more negative potentials), the observed potential will drop towards zero and will change the noise periodic pattern, producing an easily-noticeable peak in the CDNS spectrum. The capacity to recognise these peaks is best when using non-linear magnitudes of matric potential (e.g., kPa). This enables us to set the two conditions used for combining both mp_f and mp_s time series into a corrected series of data over an iterative process: i) the matric potential threshold (Tau, τ) at which the temperature sensitivity effects become significant (or visible), and ii) the change in the noise pattern when a wetting process is perceived.

The successful dynamic combination of mp_f and mp_s ensures that fast wetting fronts are not smoothed and therefore described in high resolution.

### Authors:

[Sebastián Bravo Peña (1,2) ![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0009-0005-3970-6187)\
[Meindert Commelin (1) ![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0001-7460-1915)\
[Ole Wendroth (3)](https://scholar.google.com/citations?user=cf-NkNEAAAAJ&hl=en)

1: Soil Physics and Land Management group; Wageningen University, Wageningen, The Netherlands\
2: Instituto de Ingeniería Agraria y Suelos, Universidad Austral de Chile, Valdivia, Chile\
3: Department of Plant and Soil Sciences, University of Kentucky, Lexington, United States\
Contact: [sebastian.bravopena\@wur.nl](mailto:sebastian.bravopena@wur.nl){.email}

## Performing calculations

```{r, include = TRUE, echo = TRUE, eval = FALSE}
# This approach first requires the computation of the noise signal. For which 
# we proposed the Convergence-Divergence Noise Signal (CDNS) function. 

#' @title Convergence-Divergence Noise Signal (CDNS): rollmean.cdns() 
#' @description Computes centered fast and slow simple moving averages (SMA) to 
#' calculate the noise signal. 
#' @param data Data frame with columns order as follows:
#' date = 1st col, date-time
#' num  = 2nd col, number of observation 
#' mp   = 3th col, absolute value of pressure head
#' @param k1 odd-number window size for fast SMA
#' @param k2 odd-number window size for slow SMA

# Creating a temporary data frame with the new indexes (extremes discarded)
# the idx will start from 1 + (k-1)/2, and end nrow(data) - (k-1)/2
# so the first value of the new df can calculate the MA from data rows  
# Two data frames will run in parallel
rollmean.cdns <- function(data, k1, k2) {
  
  # FAST simple moving average (SMA)
  idx_f <- seq(from = 1+(k1-1)/2, to = (nrow(data)-(k1-1)/2), by = 1)
  df_f  <- data.frame(idx_f) # creates the first column in the df_f
  
  # Calculating fast SMA average for head
  for (i in (k1-1)/2:ncol(data)) {
    df_f[[paste0(names(data)[i], "_f")]] <- NA
    
    for (j in idx_f) {
      if (any(!is.na(data[(j-(k1-1)/2):(j+(k1-1)/2), i]))) {
        df_f[j-(k1-1)/2, paste0(names(data)[i], "_f")] <- mean((data[(j-(k1-1)/2):(j+(k1-1)/2), i]), na.rm = TRUE)
      } else {
        df_f[j-(k1-1)/2, paste0(names(data)[i], "_f")] <- data[(j-(k1-1)/2), i]
      }
    }
  }
  # SLOW simple moving average (SMA)
  idx_s <- seq(from = 1+(k2-1)/2, to = (nrow(data)-(k2-1)/2), by = 1)
  df_s  <- data.frame(idx_s) # this creates the first column in the df_s
  
  # Calculating moving average for head
  for (i in 3) {
    df_s[[paste0(names(data)[i], "_s")]] <- NA
    
    for (j in idx_s) {
      if (any(!is.na(data[(j-(k2-1)/2):(j+(k2-1)/2), i]))) {
        df_s[j-(k2-1)/2, paste0(names(data)[i], "_s")] <- mean((data[(j-(k2)/2):(j+(k2-1)/2), i]), na.rm = TRUE)
      } else {
        df_s[j-(k2-1)/2, paste0(names(data)[i], "_s")] <- data[(j-(k2-1)/2), i]
      }
    }
  }
  
  N <- nrow(data)
  data$wc   <- c(data$wc[1:floor(k1/2)], df_f$wc_f, data$wc[(N-floor(k1/2)+1):N])
  data$temp <- c(data$temp[1:floor(k1/2)], df_f$temp_f, data$temp[(N-floor(k1/2)+1):N])
  data$mp_f <- c(data$mp[1:floor(k1/2)], df_f$mp_f, data$mp[(N-floor(k1/2)+1):N])
  data$mp_s <- c(data$mp[1:floor(k2/2)], df_s$mp_s, data$mp[(N-floor(k2/2)+1):N])
  data$macd <- data$mp_s-data$mp_f
  
  return(data)
}

#' @title Temperature Sensitivity Correction Function (TSCF): tscf.function()
#' @param data Data frame returned by the CDNS function.
#' @param tau Matric potential threshold for temperature sensitivity
#' @param signal Threshold to discriminate noise from wetting front signal change

tscf.function <- function(data, tau, k1, k2, signal = 0.3) {
  # Initialize variables
  tau <- tau # Set the temperature sensitivity threshold
  df$jump <- 0 # State variable to jump from mp_s to mp_f (0: no jump, 1 = jump)
  df$sens <- 0 # State variable to define sensitivity range (0: h <= 200, 1: h > 200)
  df$mp_corrected <- NaN # Corrected time series combining MACD
  jumped <- FALSE
  
  # Loop through the data
  for (i in k2:(nrow(df) - k2 + 1)) {
    
    # Check for NA/NaN values
    if (all(!is.na(df$macd[(i - k2 + 1):(i)])) &&
        all(!is.na(df$macd[(i + 1):(i + k2 - 1)])) &&
        all(!is.na(df$mp_f[i])) &&
        all(!is.na(df$mp_s[i])))  {
      
      # We set the jump backward from mp_s to mp_f
      if ((max(df$macd[(i - k2 + 1):(i)], na.rm = TRUE)) < 
          (signal * max(abs(df$macd[(i + 1):(i + k2 - 1)]), na.rm = TRUE)) &&
          
          (-10 < df$macd[i]) & (df$macd[i] < 10) &&
          
          (0.99 < (df$mp_f[i]/df$mp_s[i]) & (df$mp_f[i]/df$mp_s[i]) < 1.01) &&
          
          sd(df$macd[(i - k2 + 1):(i)], na.rm = TRUE) < 
          sd(df$macd[(i + 1):(i + k2 - 1)], na.rm = TRUE) &&
          
          any(diff(df$mp_s[i:(i + k2 - 1)]) <= 0))  {
        
        # Set jump to 1 if conditions are met
        df$jump[i] <- 1
        jumped <- 1
      } else {
        df$jump[i] <- 0 
      }  
      
      # We set conditions for the sensitivity range
      if (df$mp_f[i] > tau) {
        df$sens[i] <- 1 
      } else {
        df$jump[i] <- 0
        jumped <- 0
        df$sens[i] <- 0 
      }
      # Check if jumped should be deactivated
      if (jumped == 1) {
        df$jump[i] <- 1
        if (i + k1 <= nrow(df)) {
          # Deactivate jumped if mp_s starts increasing
          if (df$mp_s[i + k1] > df$mp_s[i + k1 - 1]) {
            jumped <- 0
          }
        }
      }
      
    } else {
      df$jump[i] <- NA
      df$sens[i] <- NA
      df$mp_corrected[i] <- NA
    }
  }

  for (i in 1:nrow(df))  {
    # Combining MACD into a corrected time series
    if ((df$sens[i] == 1 && df$jump[i] == 1) || df$sens[i] == 0) {
      df$mp_corrected[i] <- df$mp_f[i]
    } else {
      df$mp_corrected[i] <- df$mp_s[i]
    }
  }
  # remove unnecessary columns
  df$jump <- NULL
  df$sens <- NULL
  
  return(df)
}
```
