rm(list = ls())

# ---- Day-scale Temperature Sensitivity Correction ----
#' @author Sebastián Bravo Peña, Meindert Comellin, & Ole Wendroth
 
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
  
  # Calculating fast SMA average for head, theta and temp
  for (i in 3:ncol(data)) {
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
  df_s  <- data.frame(idx_s) # this create the first column in the df_s
  
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
  # remove uncessary columns
  df$jump <- NULL
  df$sens <- NULL
  
  return(df)
}

# ---- Example run ----
df <- read.csv(file = "Example_data.csv", header = TRUE, sep = ";") 

# Parameters
k1 <- 7
k2 <- 145
tau <- 200

# Run both functions
df <- rollmean.cdns(df, k1, k2)
df <- tscf.function(df, k1, k2, tau)

summary(df) # print summary

# Visualizing
par(mfrow = c(1,1), mar = c(4, 4, 1, 1))
plot(df$mp_f, type = 'l', lwd = 1, log = 'y', pch = ".", ylim = c(10,5000),
     ylab = 'Matric potential , |kPa|', col = '#EEB949', xlab = 'Time , dd-mm-yyyy', xaxt ='n')
lines(df$mp_s, col = 'red', lwd = 1, lty = 1)
lines(df$mp_corrected, col = 'blue', lwd = 1, lty = "dotted")
abline(h = tau, lty = 3, col = 'grey')
legend('topleft', legend = c(expression("h"[f]), expression("h"[s]), expression("h"[corrected])), 
       lty = c(1,1,3), col = c('#EEB949', 'red', 'blue'), 
       lwd = c(1,1,1), cex = 0.8, bty = 'n')
axis(side = 1, at = c(1, 13249, 26359,  39319, 52561), 
     lab = expression("01-06-2020", "01-09-2020", "01-12-2020", "01-03-2021", "01-06-2021"), cex.axis = 1)
text(y = 250, x = 2000, expression(paste(tau, ' = -200 kPa')), col = 'darkgrey', cex = 0.8)

# Zoomed-in 1
plot(df$mp_f, type = 'l', lwd = 1, log = 'y', lty = 1, xlim = c(26500,28500),ylim = c(500,3000),
     ylab = 'Matric potential , |kPa|', col = '#EEB949', xlab = 'Time index')
lines(df$mp_s, col = 'red', lwd = 0.6, lty = 1)
lines(df$mp_corrected, col = 'blue', lwd = 1, lty = 5)
abline(h = tau, lty = 3, col = 'grey')

# Zoomed-in 2
plot(df$mp_f, type = 'l', lwd = 1, log = 'y', lty = 1 ,xlim = c(27500,29500), ylim = c(10,500),
     ylab = 'Matric potential , |kPa|', col = '#EEB949', xlab = 'Time index')
lines(df$mp_s, col = 'red', lwd = 0.6, lty = 1)
lines(df$mp_corrected, col = 'blue', lwd = 1, lty = 5)
abline(h = tau, lty = 3, col = 'grey')
