
# Script for running sensitivity analyses on influena and RSV hospitalisation data 
# fitting assuming different prior for parameter tau
# Produced supplementary figure 1 at end of the script

##############################################################################################################
# Define origin_date (using final dataset)
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
library(idpalette)

# Loading required functions
source('R/ps_single_analysis_scripts.R')
source('R/extra_functions_in_script.R')



#####################################################################################################################
# Load the hospitalisation data
df_hosp <- read.csv(paste("processed-data/2025/", origin_date"/NZ-hospitalisations-count-", origin_date, ".csv", sep=""))

# Set limits on dates to consider 
max_date <- origin_date
min_date <- origin_date - 3*365

# Some processing in case hospitalisation data isn't in order
df_hosp <- df_hosp[df_hosp$admission_date<=max_date & df_hosp$admission_date>min_date,]
df_hosp$admission_date <- as.Date(df_hosp$admission_date)
df_hosp$time_index <- as.numeric(df_hosp$admission_date) - min(as.numeric(df_hosp$admission_date))+1
df_hosp <- df_hosp[order(df_hosp$time_index),]


#####################################################################################################################
# Set some stan settings
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = 4)

# Loading Stan model
ps_single_mod_priors <- stan_model('stan/ps_single_final_priors.stan')


####################################################################################################
# Read in model fit to final SARS-CoV-2 hospitalisation data (output from main analysis script)
mod_fit <- readRDS(paste('fitted_stan_models/', origin_date, '/',origin_date,'-NZ-', "SARSCOV2",'-hosp_fit.rds', sep=""))

post <- rstan::extract(mod_fit)

phi_mn <- mean(post$phi)
phi_sd <- sd(post$phi)

tau_mn <- mean(post$tau)
tau_sd <- sd(post$tau)

######################################################################################################################
## Fit models to hospitalisation data
################################################################################

##############################################################################################
# Fit models assuming tau_mn is two times the original

for(i in c("RSV","flu")){
  
  print(Sys.time())
  print("SARS-CoV-2 Hospitalisations")
  print(i)
  
  df <- df_hosp[df_hosp$pathogen ==i,]
  max_datei <- max(df$admission_date)
  df <- df[df$admission_date<= max_datei,]
  df <- df[df$admission_date>as.Date("2025-01-01"),]
  
  knots <- get_knots(df$time_index, days_per_knot = 5, spline_degree = 3)
  
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
                   tau_mean = tau_mn*2,
                   tau_sd = tau_sd) 
  
  
  mod_fit <- sampling(ps_single_mod_priors,
                      iter= 2500,
                      warmup = 500,
                      chains=4,
                      data = mod_data)
  
  saveRDS(mod_fit, paste('fitted_stan_models/nz_sens/',origin_date,'-NZ-', i,'-hosp_fit_sens2.rds', sep=""))
  print(Sys.time())
  
}

##############################################################################################
# Fit models assuming tau_mn is 0.5 times the original

for(i in c("RSV","flu")){
  
  print(Sys.time())
  print("SARS-CoV-2 Hospitalisations")
  print(i)
  
  df <- df_hosp[df_hosp$pathogen ==i,]
  max_datei <- max(df$admission_date)
  df <- df[df$admission_date<= max_datei,]
  df <- df[df$admission_date>as.Date("2025-01-01"),]
  
  knots <- get_knots(df$time_index, days_per_knot = 5, spline_degree = 3)
  
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
                   tau_mean = tau_mn*0.5,
                   tau_sd = tau_sd) 
  
  
  mod_fit <- sampling(ps_single_mod_priors,
                      iter= 2500,
                      warmup = 500,
                      chains=4,
                      data = mod_data)
  
  saveRDS(mod_fit, paste('fitted_stan_models/nz_sens/',origin_date,'-NZ-', i,'-hosp_fit_sens05.rds', sep=""))
  print(Sys.time())
  
}


######################################################################################################################
## Produce modelled outputs for plotting
################################################################################

hosp_inc <- data.frame()
hosp_inc_dow <- data.frame()
hosp_gr <- data.frame()

# Loop through RSV and influenza for models with tau_mn = 2*original
for(i in c("RSV","flu")){
  print(i)
  
  mod_fit <- readRDS(paste('fitted_stan_models/nz_sens/',origin_date,'-NZ-', i,'-hosp_fit_sens2.rds', sep=""))
  
  df <- df_hosp[df_hosp$pathogen ==i,]

  max_datei <- max(df$admission_date)
  df <- df[df$admission_date>as.Date("2025-01-01"),]
  
  df <- df[df$admission_date<= max_datei,]
  df$cases <- df$hospitalisations
  df$notification_date <- df$admission_date
  
  outputs <- get_all_outputs_hosp(df, mod_fit, location="NZ", gamma_dist = gammaDist, b=b_cov, n=n_cov, tau_max = 21, pathogen =i, dow="No")
  
  outputs[[1]]$label <- 2
  outputs[[3]]$label <- 2
  hosp_inc <- rbind(hosp_inc, outputs[[1]])
  hosp_gr <- rbind(hosp_gr, outputs[[3]])
  
}

# Loop through RSV and influenza for models with tau_mn = 0.5*original
for(i in c("RSV","flu")){
  print(i)
  
  mod_fit <- readRDS(paste('fitted_stan_models/nz_sens/',origin_date,'-NZ-', i,'-hosp_fit_sens05.rds', sep=""))
  
  df <- df_hosp[df_hosp$pathogen ==i,]
  
  max_datei <- max(df$admission_date)
  df <- df[df$admission_date>as.Date("2025-01-01"),]
  
  df <- df[df$admission_date<= max_datei,]
  df$cases <- df$hospitalisations
  df$notification_date <- df$admission_date
  
  outputs <- get_all_outputs_hosp(df, mod_fit, location="NZ", gamma_dist = gammaDist, b=b_cov, n=n_cov, tau_max = 21, pathogen =i, dow="No")
  
  outputs[[1]]$label <- 0.5
  outputs[[3]]$label <- 0.5
  hosp_inc <- rbind(hosp_inc, outputs[[1]])
  #hosp_inc_dow <- rbind(cov_inc_dow, outputs[[2]])
  hosp_gr <- rbind(hosp_gr, outputs[[3]])
  
}

# Loop through RSV and influenza for models with tau_mn = original
for(i in c("RSV","flu")){
  print(i)
  
  mod_fit <- readRDS(paste('fitted_stan_models/', origin_date, '/',origin_date,'-NZ-', i,'-hosp_fit.rds', sep=""))
  
  df <- df_hosp[df_hosp$pathogen ==i,]
  
  max_datei <- max(df$admission_date)
  df <- df[df$admission_date>as.Date("2025-01-01"),]

  
  df <- df[df$admission_date<= max_datei,]
  df$cases <- df$hospitalisations
  df$notification_date <- df$admission_date
  
  outputs <- get_all_outputs_hosp(df, mod_fit, location="NZ", gamma_dist = gammaDist, b=b_cov, n=n_cov, tau_max = 21, pathogen =i, dow="No")
  
  outputs[[1]]$label <- 1
  outputs[[3]]$label <- 1
  
  hosp_inc <- rbind(hosp_inc, outputs[[1]])
  hosp_gr <- rbind(hosp_gr, outputs[[3]])
  
}

######################################################################################################################
## Produce supplementary figure 1
################################################################################
cols3 <- RColorBrewer::brewer.pal(3,"Dark2")

# Relabelling some variables
hosp_inc$label <- as.factor(hosp_inc$label)
relevel(hosp_inc$label, ref="1")
levels(hosp_inc$label) <- c("Original","0.5 times original","2 times original" )

# Relabelling some variables
hosp_gr$label <- as.factor(hosp_gr$label)
relevel(hosp_gr$label, ref="1")
levels(hosp_gr$label) <- c("Original","0.5 times original","2 times original" )

# Plotting panels
plt1a<- ggplot(hosp_inc[hosp_inc$pathogen=="flu" &hosp_inc$time>as.Date("2025-03-01"),])+
  geom_line(aes(x=time, y=y, col=factor(label), group=label))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=factor(label), group=label), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=factor(label), group=label), alpha=0.2)+
  geom_point(data = df_hosp[df_hosp$pathogen=="flu" &df_hosp$admission_date>as.Date("2025-03-01"),], aes(x=admission_date, y=hospitalisations ), col="black",shape=16, size=0.8)+
  geom_line(data = df_hosp[df_hosp$pathogen=="flu" &df_hosp$admission_date>as.Date("2025-03-01"),], aes(x=admission_date, y=hospitalisations ), col="black",linewidth=0.2)+
  scale_color_manual("",values=cols3)+
  scale_fill_manual("",values=cols3)+
  theme_bw(base_size = 14)+
  ylab("Hospitalisations" )+
  xlab("Date")+
  coord_cartesian()+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b")+
  theme(strip.background = element_rect(fill="white"),
        legend.position = "none")


plt1b <- ggplot(hosp_gr[hosp_gr$pathogen=="flu" &hosp_gr$time>as.Date("2025-03-01"),])+
  geom_line(aes(x=time, y=y, col=factor(label), group=label))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=factor(label), group=label), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=factor(label), group=label), alpha=0.2)+
  scale_color_manual("",values=cols3)+
  scale_fill_manual("",values=cols3)+
  theme_bw(base_size = 14)+
  geom_hline(yintercept = 0, linetype="dashed")+
  ylab("Hospitalisations" )+
  xlab("Date")+
  coord_cartesian()+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b")+
  theme(strip.background = element_rect(fill="white"),
        legend.position = "none")

plt1c <- ggplot(hosp_inc[hosp_inc$pathogen=="RSV" &hosp_inc$time>as.Date("2025-03-01"),])+
  geom_line(aes(x=time, y=y, col=factor(label), group=label))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=factor(label), group=label), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=factor(label), group=label), alpha=0.2)+
  geom_point(data = df_hosp[df_hosp$pathogen=="RSV" &df_hosp$admission_date>as.Date("2025-03-01"),], aes(x=admission_date, y=hospitalisations ), col="black",shape=16, size=0.8)+
  geom_line(data = df_hosp[df_hosp$pathogen=="RSV" &df_hosp$admission_date>as.Date("2025-03-01"),], aes(x=admission_date, y=hospitalisations ), col="black",linewidth=0.2)+
  scale_color_manual("",values=cols3)+
  scale_fill_manual("",values=cols3)+
  theme_bw(base_size = 14)+
  ylab("Hospitalisations" )+
  xlab("Date")+
  coord_cartesian()+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b")+
  theme(strip.background = element_rect(fill="white"),
        legend.position = "none")


plt1d <- ggplot(hosp_gr[hosp_gr$pathogen=="RSV" &hosp_gr$time>as.Date("2025-03-01"),])+
  geom_line(aes(x=time, y=y, col=factor(label), group=label))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=factor(label), group=label), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=factor(label), group=label), alpha=0.2)+
  scale_color_manual("Mean tau in\nprior distribution",values=cols3)+
  scale_fill_manual("Mean tau in\nprior distribution",values=cols3)+
  geom_hline(yintercept = 0, linetype="dashed")+
  theme_bw(base_size = 14)+
  ylab("Hospitalisations" )+
  xlab("Date")+
  coord_cartesian()+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b")+
  theme(strip.background = element_rect(fill="white"),
        legend.position = "bottom")



plt1a <- plt1a+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.tag.position = c(0.0,0.98))+
  labs(tag="A")+
  annotate("label",label="Influenza", y=Inf, x = as.Date("2025-03-15"), fill= "white", color="black", vjust=1.2, size=5, hjust=0.0)


plt1b <- plt1b+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())

plt1c <- plt1c+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.tag.position = c(0.0,0.98))+
  labs(tag="B")+
  annotate("label",label="RSV", y=Inf, x = as.Date("2025-03-15"), fill= "white", color="black", vjust=1.2, size=5, hjust=0.0)




plt1a+plt1b+plt1c+plt1d+plot_layout(nrow=4, heights=c(1,0.5,1,0.5))

ggsave(paste('figure/', "paper_NZ",'/', 'revisions_prior_sens', '.png', sep=""), width=8, height=10)
