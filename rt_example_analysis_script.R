# Example script that produced past and current trend analysis in real-time
# during the season

##############################################################################################################
# Define origin_date (used for labelling and reading in data)
origin_date <- as.Date("2025-10-23")
dir.create(paste('figure/', origin_date, sep=""))
dir.create(paste('fitted_stan_models/', origin_date, sep=""))
first_date <- origin_date - 180 # For plotting
date_column <- "notification_date"

###############################################################################################################
# Loading required packages
library(ggplot2)
library(rstan)
library(patchwork)
#library(idpalette)

# Loading required functions
#source('R/ps_analysis_scripts.R')
source('R/ps_single_analysis_scripts.R')
source('R/extra_functions_in_script.R')

###############################################################################################################
## Loading and processing case data (this data is not publicly available)
# Load the case data 
df_cov <- read.csv(paste("data/SARSCOV2-PCR-only-case-count-", origin_date, ".csv", sep=""))
df_cov <- df_cov[df_cov$test_type=="PCR",]

# Set limits on dates to consider (Note some data does not have three years of data
max_date <- origin_date
min_date <- max_date - 365*3
df_cov <- df_cov[df_cov$notification_date<=max_date & df_cov$notification_date>min_date,]

## Some processing in case case data isn't in order
df_cov[,date_column] <- as.Date(df_cov[,date_column])
df_cov$time_index <- as.numeric(df_cov[,date_column]) - min(as.numeric(df_cov[,date_column]))+1
df_cov <- df_cov[order(df_cov$time_index),]

#####################################################################################################################
## Loading and processing hospitalisation data (this data is not publicly available)
# Load the hospitalisation data
df_hosp <- read.csv(paste("data/NZ-hospitalisations-count-", origin_date, ".csv", sep=""))

# Set limits on dates to consider (Note some data does not have three years of data
max_date <- origin_date
min_date <- origin_date - 3*365
df_hosp <- df_hosp[df_hosp$admission_date<=max_date & df_hosp$admission_date>min_date,]

## Some processing in case hospitalisation data isn't in order
df_hosp$admission_date <- as.Date(df_hosp$admission_date)
df_hosp$time_index <- as.numeric(df_hosp$admission_date) - min(as.numeric(df_hosp$admission_date))+1
df_hosp <- df_hosp[order(df_hosp$time_index),]

#####################################################################################################################
# Setting some stan settings
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = 4)

# Loading Stan models
ps_single_mod <- stan_model('stan/ps_single_final.stan')
ps_single_mod_priors <- stan_model('stan/ps_single_final_priors.stan')


#############################################################################################################################################
## Fitting to COVID-19 case data
for(i in "NZ"){
  print(Sys.time())
  print("SARS-CoV-2")
  print(i)
  
  # Ensure NZ data only (some processed datasets used by ACEFA have the eight Australian jurisidictions in them)
  df <- df_cov[df_cov$location ==i,]
  
  # Set the boundaries for the series of equally spaced knots
  knots <- get_knots(df$time_index, days_per_knot = 5, spline_degree = 3)
  
  # Format data into format required for the stan model
  mod_data <- list(num_data = nrow(df),
                   num_knots = length(knots),
                   knots = knots,
                   spline_degree=3,
                   Y = df$cases,
                   X = df$time_index,
                   week_effect = 7,
                   DOW = (df$time_index %% 7)+1) 
  
  # Fit the stan model
  mod_fit <- sampling(ps_single_mod,
                      iter= 2500,
                      warmup = 500,
                      chains=4,
                      data = mod_data)
  
  # Save the stan model into a directory 
  saveRDS(mod_fit, paste('fitted_stan_models/', origin_date, '/',origin_date,'-', i,'-cov_fit.rds', sep=""))
  print(Sys.time())
}

#############################################################################################################################################
## Fitting to COVID-19 hospitalisation data

for(i in c("SARSCOV2")){
  
  print(Sys.time())
  print("SARS-CoV-2 Hospitalisations")
  print(i)
  
  df <- df_hosp[df_hosp$pathogen ==i,]
  
  # Figure out maximum date (sometimes differences in maximum dates between pathogens in the hospitalisation data)
  max_datei <- max(df$admission_date)
  df <- df[df$admission_date<= max_datei,]
  
  # Set the boundaries for the series of equally spaced knots
  knots <- get_knots(df$time_index, days_per_knot = 5, spline_degree = 3)
  
  # Format data into format required for the stan model
  mod_data <- list(num_data = nrow(df),
                    num_knots = length(knots),
                    knots = knots,
                    spline_degree=3,
                    Y = df$hospitalisations, 
                    X = df$time_index,
                    week_effect = 1,
                    DOW = (df$time_index %% 1)+1,
                    phi_mean = 0,
                    phi_sd = 50,
                    tau_mean = 0,
                    tau_sd = 50) 
  
  # Fit the stan model
  mod_fit <- sampling(ps_single_mod_priors,
                      iter= 2500,
                      warmup = 500,
                      chains=4,
                      data = mod_data)
  
  # Save the stan model into a directory 
  saveRDS(mod_fit, paste('fitted_stan_models/', origin_date, '/',origin_date,'-NZ-', i,'-hosp_fit.rds', sep=""))
  print(Sys.time())
  
}

#############################################################################################################################################
## Fitting to influenza and RSV hospitalisation data

# Extract posterior of model fit to COVID-19 hospitalisation data
post <- rstan::extract(mod_fit)

# Get the mean and standard deviation of phi parameter (used to define prior for models being fit)
phi_mn <- mean(post$phi)
phi_sd <- sd(post$phi)

#Get the mean and standard deviation of tau parameter (used to define prior for models being fit)
tau_mn <- mean(post$tau)
tau_sd <- sd(post$tau)

for(i in c("RSV","flu")){
  
  print(Sys.time())
  print("SARS-CoV-2 Hospitalisations")
  print(i)
  
  df <- df_hosp[df_hosp$pathogen ==i,]
  
  # Figure out maximum date (sometimes differences in maximum dates between pathogens in the hospitalisation data)
  max_datei <- max(df$admission_date)
  df <- df[df$admission_date<= max_datei,]
  # Set minimum date (influenza and RSV hospitalisation data is only available continously from 1 January 2025)
  df <- df[df$admission_date>as.Date("2025-01-01"),]
  
  # Set the boundaries for the series of equally spaced knots
  knots <- get_knots(df$time_index, days_per_knot = 5, spline_degree = 3)
  
  # Format data into format required for the stan model
  mod_data <- list(num_data = nrow(df),
                   num_knots = length(knots),
                   knots = knots,
                   spline_degree=3,
                   Y = df$hospitalisations, 
                   X = df$time_index,
                   week_effect = 1,
                   DOW = (df$time_index %% 1)+1,
                   phi_mean = phi_mn,
                   phi_sd = phi_sd,
                   tau_mean = tau_mn,
                   tau_sd = tau_sd) 
  
  # Fit the stan model
  mod_fit <- sampling(ps_single_mod_priors,
                      iter= 2500,
                      warmup = 500,
                      chains=4,
                      data = mod_data)
  
  # Save the stan model into a directory 
  saveRDS(mod_fit, paste('fitted_stan_models/', origin_date, '/',origin_date,'-NZ-', i,'-hosp_fit.rds', sep=""))
  print(Sys.time())
  
}

################################################################################
## Produce figures for inclusion in reports and make csv of modelled outputs
################################################################################

################################################################################
#SARS-CoV-2 case time series figures

# Parameters for generation interval (for reproduction number estimates in reports only)
b_cov <- 0.27
n_cov <- 0.89

gammaDist <- function(b, n, a){
  val <- (b**n) * (a**(n-1)) * exp(-b*a) / gamma(n)
  val[val==Inf] <- 0
  val
} 

# Define empty data.frames
cov_inc <- data.frame()
cov_inc_dow <- data.frame()
cov_gr <- data.frame()
cov_Rt <- data.frame()

# Loop through jurisdictions (just NZ in this case) and run function producing modelled outputs and figure
for(i in "NZ"){
  print(i)
  df <- df_cov[df_cov$location ==i,]
  
  # reload the stan model for the time series
  mod_fit <- readRDS(paste('fitted_stan_models/', origin_date, '/',origin_date,"-", i,'-cov_fit.rds', sep=""))
  
  # Get all outputs (function in separate R script)
  outputs <- get_all_outputs(df, mod_fit, location=i, gamma_dist = gammaDist, b=b_cov, n=n_cov, tau_max = 21, pathogen ="SARS-CoV-2")
  
  cov_inc <- rbind(cov_inc, outputs[[1]])
  cov_inc_dow <- rbind(cov_inc_dow, outputs[[2]])
  cov_gr <- rbind(cov_gr, outputs[[3]])
  cov_Rt <- rbind(cov_Rt, outputs[[4]])
}

######################################################################################################################################################
# Hospitalisation outputs

# Define empty data.frames
hosp_inc <- data.frame()
hosp_gr <- data.frame()

# Loop through the pathogens (for NZ hospitalisation data)
for(i in c("SARSCOV2","RSV","flu")){
  print(i)
  
  # reload the stan model for the time series
  mod_fit <- readRDS(paste('fitted_stan_models/', origin_date, '/',origin_date,'-NZ-', i,'-hosp_fit.rds', sep=""))
  
  # Limit data to correct time period
  df <- df_hosp[df_hosp$pathogen ==i,]
  if(i == "SARSCOV2"){
    max_datei <- max(df$admission_date)
    
  } else{
    max_datei <- max(df$admission_date)
    df <- df[df$admission_date>as.Date("2025-01-01"),]
  }
  
  df <- df[df$admission_date<= max_datei,]
  
  # Create dummy variables to work with the function (that was written for case data)
  df$cases <- df$hospitalisations
  df$notification_date <- df$admission_date
  
  # Get all outputs (function in separate R script)
  outputs <- get_all_outputs_hosp(df, mod_fit, location="NZ", gamma_dist = gammaDist, b=b_cov, n=n_cov, tau_max = 21, pathogen =i, dow="No")
  
  hosp_inc <- rbind(hosp_inc, outputs[[1]])
  hosp_gr <- rbind(hosp_gr, outputs[[3]])
  
}

################################################################################
## Format and save modelled outputs and data for future comparison plots

# Relevel some factors
hosp_inc$pathogen <- as.factor(hosp_inc$pathogen)
df_hosp$pathogen <- as.factor(df_hosp$pathogen)
levels(hosp_inc$pathogen) <- c("Influenza", "RSV", "SARS-CoV-2")
levels(df_hosp$pathogen) <- c("Influenza", "RSV", "SARS-CoV-2")

# Format the hospitalisation data to match what was used to produce modelled outputs
df1 <- df_hosp[df_hosp$pathogen=="Influenza",]
max_datei <- max(df1$admission_date)
df1 <- df1[df1$admission_date<= max_datei,]

df2 <- df_hosp[df_hosp$pathogen=="RSV",]
max_datei <- max(df2$admission_date)
df2 <- df2[df2$admission_date<= max_datei,]

df3 <- df_hosp[df_hosp$pathogen=="SARS-CoV-2",]
max_datei <- max(df3$admission_date)
df3 <- df3[df3$admission_date<= max_datei,]

df_hosp <- rbind(df1, df2, df3)

# Save formatted data used in models
write.csv(df_hosp,paste("smoothed_estimates/nz_df_hosp",origin_date,".csv", sep=""), row.names=FALSE)
write.csv(df_cov,paste("smoothed_estimates/nz_df_cov",origin_date,".csv", sep=""), row.names=FALSE)

# Save modelled outputs for NZ cases
write.csv(cov_inc,paste("smoothed_estimates/nz_cov_inc",origin_date,".csv", sep=""), row.names=FALSE)
write.csv(cov_inc_dow,paste("smoothed_estimates/nz_cov_inc_dow",origin_date,".csv", sep=""), row.names=FALSE)
write.csv(cov_gr,paste("smoothed_estimates/nz_cov_gr",origin_date,".csv", sep=""), row.names=FALSE)

# Save modelled outputs for NZ hospitalisations
write.csv(hosp_inc,paste("smoothed_estimates/nz_hosp_inc",origin_date,".csv", sep=""), row.names=FALSE)
write.csv(hosp_gr,paste("smoothed_estimates/nz_hosp_gr",origin_date,".csv", sep=""), row.names=FALSE)

