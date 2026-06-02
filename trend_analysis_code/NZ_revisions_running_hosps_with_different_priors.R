
setwd("C:/Users/EALESO/R Projects/Winter 2025 reporting")

##############################################################################################################
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
source('R/ps_analysis_scripts.R')
source('R/ps_single_analysis_scripts.R')
source('R/format_data.R')
source('R/extra_functions_in_script.R')

###############################################################################################################
# Load the case data 

df_cov <- read.csv(paste("data/SARSCOV2-PCR-only-case-count-", origin_date, ".csv", sep=""))
df_cov <- df_cov[df_cov$test_type=="PCR",]

# Set limits on dates to consider
max_date <- origin_date
min_date <- max_date - 365*3
df_cov <- df_cov[df_cov$notification_date<=max_date & df_cov$notification_date>min_date,]


df_cov[,date_column] <- as.Date(df_cov[,date_column])

## Will the data have to be ordered at all? 
df_cov$time_index <- as.numeric(df_cov[,date_column]) - min(as.numeric(df_cov[,date_column]))+1


df_cov <- df_cov[order(df_cov$time_index),]

#####################################################################################################################
# Load the hospitalisation data
df_hosp <- read.csv(paste("data/NZ-hospitalisations-count-", origin_date, ".csv", sep=""))
#df_hospSub <- read.csv(paste("data/NZ-hospitalisations-flu-type-count-", origin_date, ".csv", sep="")) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

max_date <- origin_date
min_date <- origin_date - 3*365

df_hosp <- df_hosp[df_hosp$admission_date<=max_date & df_hosp$admission_date>min_date,]
df_hosp$admission_date <- as.Date(df_hosp$admission_date)
df_hosp$time_index <- as.numeric(df_hosp$admission_date) - min(as.numeric(df_hosp$admission_date))+1
df_hosp <- df_hosp[order(df_hosp$time_index),]


#####################################################################################################################
# Set some stan settings
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = 4)

# Loading Stan models
ps_single_mod <- stan_model('stan/ps_single_final.stan')
ps_single_mod_priors <- stan_model('stan/ps_single_final_priors.stan')
ps_inf_mod <- stan_model('stan/ps_influenza_finalV2.stan')


####################################################################################################
mod_fit <- readRDS(paste('fitted_stan_models/', origin_date, '/',origin_date,'-NZ-', "SARSCOV2",'-hosp_fit.rds', sep=""))

post <- rstan::extract(mod_fit)

phi_mn <- mean(post$phi)
phi_sd <- sd(post$phi)

tau_mn <- mean(post$tau)
tau_sd <- sd(post$tau)

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
                      iter= 2000,
                      warmup = 500,
                      chains=4,
                      data = mod_data)
  
  saveRDS(mod_fit, paste('fitted_stan_models/nz_sens/',origin_date,'-NZ-', i,'-hosp_fit_sens2.rds', sep=""))
  print(Sys.time())
  
}


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
                      iter= 2000,
                      warmup = 500,
                      chains=4,
                      data = mod_data)
  
  saveRDS(mod_fit, paste('fitted_stan_models/nz_sens/',origin_date,'-NZ-', i,'-hosp_fit_sens05.rds', sep=""))
  print(Sys.time())
  
}


################################################################################
# Preparing material for figure

hosp_inc <- data.frame()
hosp_inc_dow <- data.frame()
hosp_gr <- data.frame()

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
  #hosp_inc_dow <- rbind(cov_inc_dow, outputs[[2]])
  hosp_gr <- rbind(hosp_gr, outputs[[3]])
  
}

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


for(i in c("RSV","flu")){
  print(i)
  
  mod_fit <- readRDS(paste('fitted_stan_models/', origin_date, '/',origin_date,'-NZ-', i,'-hosp_fit.rds', sep=""))
  
  df <- df_hosp[df_hosp$pathogen ==i,]
  if(i == "SARSCOV2"){
    max_datei <- max(df$admission_date)
    
  } else{
    max_datei <- max(df$admission_date)
    df <- df[df$admission_date>as.Date("2025-01-01"),]
  }
  
  df <- df[df$admission_date<= max_datei,]
  df$cases <- df$hospitalisations
  df$notification_date <- df$admission_date
  
  outputs <- get_all_outputs_hosp(df, mod_fit, location="NZ", gamma_dist = gammaDist, b=b_cov, n=n_cov, tau_max = 21, pathogen =i, dow="No")
  
  outputs[[1]]$label <- 1
  outputs[[3]]$label <- 1
  
  hosp_inc <- rbind(hosp_inc, outputs[[1]])
  #hosp_inc_dow <- rbind(cov_inc_dow, outputs[[2]])
  hosp_gr <- rbind(hosp_gr, outputs[[3]])
  
}


###################################################################################
# Plot figure

cols3 <- RColorBrewer::brewer.pal(3,"Dark2")

hosp_inc$label <- as.factor(hosp_inc$label)
relevel(hosp_inc$label, ref="1")
levels(hosp_inc$label) <- c("Original","0.5 times original","2 times original" )

hosp_gr$label <- as.factor(hosp_gr$label)
relevel(hosp_gr$label, ref="1")
levels(hosp_gr$label) <- c("Original","0.5 times original","2 times original" )

plt1a<- ggplot(hosp_inc[hosp_inc$pathogen=="flu" &hosp_inc$time>as.Date("2025-03-01"),])+
  geom_line(aes(x=time, y=y, col=factor(label), group=label))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=factor(label), group=label), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=factor(label), group=label), alpha=0.2)+
  #geom_line(data = new_data_hosp[new_data_hosp$pathogen=="Influenza" &new_data_hosp$mask>0 & new_data_hosp$time>as.Date("2025-03-01"),], aes(x=time, y=hospitalisations, group=time), linewidth=0.2)+
  #geom_point(data = data_hosp[data_hosp$pathogen=="Influenza" &data_hosp$mask>0 & data_hosp$time>as.Date("2025-03-01"),], aes(x=time, y=hospitalisations),shape=21, size=0.8, fill="white")+
  #geom_point(data = data_hosp_final[data_hosp_final$pathogen=="Influenza" &data_hosp_final$time>as.Date("2025-03-01"),], aes(x=time, y=hospitalisations),shape=21, size=0.8, fill="black")+
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
  #geom_line(data = new_data_hosp[new_data_hosp$pathogen=="Influenza" &new_data_hosp$mask>0 & new_data_hosp$time>as.Date("2025-03-01"),], aes(x=time, y=hospitalisations, group=time), linewidth=0.2)+
  #geom_point(data = data_hosp[data_hosp$pathogen=="Influenza" &data_hosp$mask>0 & data_hosp$time>as.Date("2025-03-01"),], aes(x=time, y=hospitalisations),shape=21, size=0.8, fill="white")+
  #geom_point(data = data_hosp_final[data_hosp_final$pathogen=="RSV" &data_hosp_final$time>as.Date("2025-03-01"),], aes(x=time, y=hospitalisations),shape=21, size=0.8, fill="black")+
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
  #geom_line(data = new_data_hosp[new_data_hosp$pathogen=="RSV" &new_data_hosp$mask>0 & new_data_hosp$time>as.Date("2025-03-01"),], aes(x=time, y=hospitalisations, group=time), linewidth=0.2)+
  #geom_point(data = data_hosp[data_hosp$pathogen=="RSV" &data_hosp$mask>0 & data_hosp$time>as.Date("2025-03-01"),], aes(x=time, y=hospitalisations),shape=21, size=0.8, fill="white")+
  #geom_point(data = data_hosp_final[data_hosp_final$pathogen=="RSV" &data_hosp_final$time>as.Date("2025-03-01"),], aes(x=time, y=hospitalisations),shape=21, size=0.8, fill="black")+
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
  #geom_line(data = new_data_hosp[new_data_hosp$pathogen=="RSV" &new_data_hosp$mask>0 & new_data_hosp$time>as.Date("2025-03-01"),], aes(x=time, y=hospitalisations, group=time), linewidth=0.2)+
  #geom_point(data = data_hosp[data_hosp$pathogen=="RSV" &data_hosp$mask>0 & data_hosp$time>as.Date("2025-03-01"),], aes(x=time, y=hospitalisations),shape=21, size=0.8, fill="white")+
  #geom_point(data = data_hosp_final[data_hosp_final$pathogen=="RSV" &data_hosp_final$time>as.Date("2025-03-01"),], aes(x=time, y=hospitalisations),shape=21, size=0.8, fill="black")+
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


################################################################################
#

mod_fit <- readRDS(paste('fitted_stan_models/', origin_date, '/',origin_date,'-NZ-', "SARSCOV2",'-hosp_fit.rds', sep=""))

post <- rstan::extract(mod_fit)
phi_mn <- mean(post$phi)
phi_sd <- sd(post$phi)
tau_mn <- mean(post$tau)
tau_sd <- sd(post$tau)

rsv1 <- readRDS(paste('fitted_stan_models/', origin_date, '/',origin_date,'-NZ-', "RSV",'-hosp_fit.rds', sep=""))
rsv2 <- readRDS(paste('fitted_stan_models/nz_sens/',origin_date,'-NZ-', "RSV",'-hosp_fit_sens2.rds', sep=""))
rsv05 <- readRDS(paste('fitted_stan_models/nz_sens/',origin_date,'-NZ-', "RSV",'-hosp_fit_sens05.rds', sep=""))

inf1 <- readRDS(paste('fitted_stan_models/', origin_date, '/',origin_date,'-NZ-', "flu",'-hosp_fit.rds', sep=""))
inf2 <- readRDS(paste('fitted_stan_models/nz_sens/',origin_date,'-NZ-', "flu",'-hosp_fit_sens2.rds', sep=""))
inf05 <- readRDS(paste('fitted_stan_models/nz_sens/',origin_date,'-NZ-', "flu",'-hosp_fit_sens05.rds', sep=""))

post_rsv1 <- rstan::extract(rsv1)
post_rsv2 <- rstan::extract(rsv2)
post_rsv05 <- rstan::extract(rsv05)

post_inf1 <- rstan::extract(inf1)
post_inf2 <- rstan::extract(inf2)
post_inf05 <- rstan::extract(inf05)


ggplot(mapping=aes(x=post_rsv2$tau))+
  geom_histogram()

ggplot(mapping=aes(x=post_inf05$phi))+
  geom_histogram()


# Three posterior samples (replace with your vectors)
post1 <- rnorm(5000, 2.0, 0.4)
post2 <- rnorm(5000, 2.3, 0.3)
post3 <- rnorm(5000, 1.8, 0.5)

# Prior parameters (replace with your values)
prior_mean <- 2.0
prior_sd   <- 0.8


# Combine posteriors into one data frame ----------------------------------

library(ggplot2)

df <- data.frame(
  value = c(post_rsv1$tau, post_rsv05$tau, post_rsv2$tau),
  posterior = factor(
    rep(c("Original", "0.5 times original", "2 times original"),
        times = c(length(post_rsv1$tau), length(post_rsv05$tau), length(post_rsv2$tau)))
  )
)


# Plot overlapping histograms + prior density -----------------------------

ggplot(df, aes(x = value, fill = posterior, colour = posterior)) +
  geom_histogram(
    aes(y=after_stat(count)),
    binwidth = 0.0025,
    alpha = 0.3,
    position = "identity"
  ) +
  
  # prior normal density overlay
  stat_function(
    fun = function(x) dnorm(x, tau_mn, tau_sd) *6000* 0.0025,
    colour = cols3[3],
    linewidth = 1.5,
    linetype = "dashed"
  ) +
  stat_function(
    fun = function(x) dnorm(x, tau_mn*2, tau_sd) *6000* 0.0025,
    colour = cols3[2],
    linewidth = 1.5,
    linetype = "dashed"
  ) +
  stat_function(
    fun = function(x) dnorm(x, tau_mn*0.5, tau_sd) *6000* 0.0025,
    colour = cols3[1],
    linewidth = 1.5,
    linetype = "dashed"
  ) +
  
  labs(
    x = "Tau",
    y = "Density"
  ) +
  
  scale_color_manual(values=c(cols3[1],cols3[2],cols3[3]))+
  scale_fill_manual(values=c(cols3[1],cols3[2],cols3[3]))+

  theme_bw(base_size = 14)
