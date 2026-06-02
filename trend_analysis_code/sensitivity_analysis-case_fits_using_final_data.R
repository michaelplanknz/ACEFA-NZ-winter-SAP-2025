# Script for running analyses on the final dataset for COVID-19 cases 
# fitting up to the dates which data was available for in real-time during the season

##############################################################################################################
# Create directory for saving model fits for sensitivity analyses
dir.create(paste('fitted_stan_models/nz_sens/', sep=""))
date_column <- "notification_date"

###############################################################################################################
# Loading required packages
library(ggplot2)
library(rstan)
library(patchwork)
library(idpalette)

# Loading required functions
source('R/ps_single_analysis_scripts.R')
source('R/extra_functions_in_script.R')

################################################################################
# Set dates to consider for model fitting

# The final dataset
final_date <- as.Date("2025-10-23")

# The maximum dates used in real-time analyses during the season
cov_dates <- c(as.Date("2025-06-01"),
               as.Date("2025-06-08"),
               as.Date("2025-06-15"),
               as.Date("2025-06-29"),
               as.Date("2025-07-20"),
               as.Date("2025-07-27"),
               as.Date("2025-07-31"),
               as.Date("2025-08-07"),
               as.Date("2025-08-13"),
               as.Date("2025-08-21"),
               as.Date("2025-08-27"),
               as.Date("2025-09-04"),
               as.Date("2025-09-10"),
               as.Date("2025-09-18"),
               as.Date("2025-09-25"))

###############################################################################################################
## Loading and processing case data (this data is not publicly available)
# Load the final case data 
df_cov <- read.csv(paste("processed-data/2025/", final_date"/SARSCOV2-PCR-only-case-count-", final_date, ".csv", sep=""))
df_cov <- df_cov[df_cov$test_type=="PCR",]

# Set limits on dates to consider
max_date <- final_date
min_date <- max_date - 365*3
df_cov <- df_cov[df_cov$notification_date<=max_date & df_cov$notification_date>min_date,]

## Some processing in case case data isn't in order
df_cov[,date_column] <- as.Date(df_cov[,date_column])
df_cov$time_index <- as.numeric(df_cov[,date_column]) - min(as.numeric(df_cov[,date_column]))+1
df_cov <- df_cov[order(df_cov$time_index),]
#####################################################################################################################
# Set some stan settings
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = 4)

# Loading Stan model
ps_single_mod <- stan_model('stan/ps_single_final.stan')

######################################################################################################################
## Fit models to COVID-19 case data
################################################################################

for(i in 1:length(cov_dates)){
  
  origin_date <- cov_dates[i]
  print(Sys.time())
  print(origin_date)
  
  df <- df_cov[df_cov$location =="NZ" & df_cov$notification_date<=origin_date,]
  
  knots <- get_knots(df$time_index, days_per_knot = 5, spline_degree = 3)
  
  mod_data <- list(num_data = nrow(df),
                   num_knots = length(knots),
                   knots = knots,
                   spline_degree=3,
                   Y = df$cases,
                   X = df$time_index,
                   week_effect = 7,
                   DOW = (df$time_index %% 7)+1) 
  
  mod_fit <- sampling(ps_single_mod,
                      iter= 2500,
                      warmup = 500,
                      chains=4,
                      data = mod_data)

  saveRDS(mod_fit, paste('fitted_stan_models/nz_sens/',origin_date,'-cov_case_sens.rds', sep=""))
  print(Sys.time())
}

######################################################################################################################
## Produce csv of modelled outputs
################################################################################

# Produce COVID outputs

b_cov <- 0.27
n_cov <- 0.89

gammaDist <- function(b, n, a){
  val <- (b**n) * (a**(n-1)) * exp(-b*a) / gamma(n)
  val[val==Inf] <- 0
  val
} 

cov_inc <- data.frame()
cov_inc_dow <- data.frame()
cov_gr <- data.frame()
cov_Rt <- data.frame()


for(i in 1:length(cov_dates)){
  # Loop through dates considered
  origin_date <- cov_dates[i]
  print(origin_date)
  
  df <- df_cov[df_cov$location =="NZ" & df_cov$notification_date<=origin_date,]

  # Reload the stan model for the time series
  mod_fit <- readRDS(paste('fitted_stan_models/nz_sens/',origin_date,'-cov_case_sens.rds', sep=""))

  # Get all outputs (function in separate R script)
  outputs <- get_all_outputs(df, mod_fit, location="NZ", gamma_dist = gammaDist, b=b_cov, n=n_cov, tau_max = 21, pathogen ="SARS-CoV-2")
  
  outputs[[1]]$origin_date <- origin_date
  outputs[[2]]$origin_date <- origin_date
  outputs[[3]]$origin_date <- origin_date
  outputs[[4]]$origin_date <- origin_date
  
  cov_inc <- rbind(cov_inc, outputs[[1]])
  cov_inc_dow <- rbind(cov_inc_dow, outputs[[2]])
  cov_gr <- rbind(cov_gr, outputs[[3]])
  cov_Rt <- rbind(cov_Rt, outputs[[4]])
  
}

# Save all outputs as csv files
write.csv(cov_inc, paste('smoothed_estimates/nz-cov_inc_sens.csv', sep=""))
write.csv(cov_inc_dow, paste('smoothed_estimates/nz-cov_inc_dow_sens.csv', sep=""))
write.csv(cov_gr, paste('smoothed_estimates/nz-cov_gr_sens.csv', sep=""))
write.csv(cov_Rt, paste('smoothed_estimates/nz-cov_Rt_sens.csv', sep=""))

