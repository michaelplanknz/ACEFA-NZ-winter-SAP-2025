# This script produces figures 1--4 from the manuscript.
# As data can not be shared publicly all lines involving the datasets have been 
# commented out (so the code can be run). Accordingly, the figures produced only
# do not include data and only display the modelled outputs.

##########################################################################################
library(ggplot2)
library(patchwork)
##############################################################################################################
final_date <- as.Date("2025-10-23")

unq_dates <- c(as.Date("2025-06-05"),
               as.Date("2025-06-12"),
               as.Date("2025-06-26"),
               as.Date("2025-07-03"),
               as.Date("2025-07-24"),
               as.Date("2025-07-31"),
               as.Date("2025-08-07"),
               as.Date("2025-08-14"),
               as.Date("2025-08-21"),
               as.Date("2025-08-28"),
               as.Date("2025-09-04"),
               as.Date("2025-09-11"),
               as.Date("2025-09-18"),
               as.Date("2025-09-25"),
               as.Date("2025-10-02"),
               final_date)

unq_dates <- as.Date(unq_dates)

#######################################################################################
# Loading in the modelled estimates

df_gr <- data.frame()
df_inc <- data.frame()
df_dow <- data.frame()

hosp_gr <- data.frame()
hosp_inc <- data.frame()

#data_hosp <- data.frame()
#data_cov <- data.frame()


for(i in 1:length(unq_dates) ){
  print(unq_dates[i])
  
  #temp_df <- read.csv(paste("smoothed_estimates/nz_df_cov",unq_dates[i],".csv", sep=""))
  #temp_df$origin <- unq_dates[i]
  #data_cov <- rbind(data_cov, temp_df)
  
  #temp_df <- read.csv(paste("smoothed_estimates/nz_df_hosp",unq_dates[i],".csv", sep=""))
  #temp_df$origin <- unq_dates[i]
  #data_hosp <-rbind(data_hosp, temp_df)
  
  temp_df <- read.csv(paste("smoothed_estimates/nz_cov_inc",unq_dates[i],".csv", sep=""))
  temp_df$origin <- unq_dates[i]
  df_inc <-rbind(df_inc, temp_df)
  
  temp_df <- read.csv(paste("smoothed_estimates/nz_cov_inc_dow",unq_dates[i],".csv", sep=""))
  temp_df$origin <- unq_dates[i]
  df_dow <-rbind(df_dow, temp_df)
  
  temp_df <- read.csv(paste("smoothed_estimates/nz_cov_gr",unq_dates[i],".csv", sep=""))
  temp_df$origin <- unq_dates[i]
  df_gr <-rbind(df_gr, temp_df)
  
  temp_df <- read.csv(paste("smoothed_estimates/nz_hosp_inc",unq_dates[i],".csv", sep=""))
  temp_df$origin <- unq_dates[i]
  hosp_inc <-rbind(hosp_inc, temp_df)
  
  temp_df <- read.csv(paste("smoothed_estimates/nz_hosp_gr",unq_dates[i],".csv", sep=""))
  temp_df$origin <- unq_dates[i]
  hosp_gr <-rbind(hosp_gr, temp_df)
  
}

df_inc$time <- as.Date(df_inc$time)
df_inc$origin <- as.Date(df_inc$origin)

df_gr$time <- as.Date(df_gr$time)
df_gr$origin <- as.Date(df_gr$origin)

df_dow$time <- as.Date(df_dow$time)
df_dow$origin <- as.Date(df_dow$origin)

#data_cov$time <- as.Date(data_cov$notification_date)
#data_cov$origin <- as.Date(data_cov$origin)
#data_cov <- data_cov[data_cov$location == "NZ",]

hosp_inc[hosp_inc$pathogen=="flu",]$pathogen <- "Influenza"
hosp_inc[hosp_inc$pathogen=="SARSCOV2",]$pathogen <- "SARS-CoV-2"
hosp_inc$pathogen <- as.factor(hosp_inc$pathogen)

hosp_gr[hosp_gr$pathogen=="flu",]$pathogen <- "Influenza"
hosp_gr[hosp_gr$pathogen=="SARSCOV2",]$pathogen <- "SARS-CoV-2"
hosp_gr$pathogen <- as.factor(hosp_gr$pathogen)


#data_hosp$pathogen <- as.factor(data_hosp$pathogen)
#data_hosp$admission_date <- as.Date(data_hosp$admission_date)
#data_hosp$time <- data_hosp$admission_date

########################################################################################
# Overall epidemic trends (Figures 1&2)
########################################################################################

################################################################################
## Figure 1

df_final <- df_inc[df_inc$origin == final_date,]
df_gr_final <- df_gr[df_gr$origin == final_date,]
df_dow_final <- df_dow[df_dow$origin == final_date,]

#data_cov_final <- data_cov[data_cov$origin == final_date,]

p1 <- ggplot(df_dow_final[df_dow_final$time>as.Date("2025-01-01"),])+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  #geom_point(data=data_cov_final[data_cov_final$notification_date>as.Date("2025-01-01"),], aes(x=time, y=cases), size=0.5)+
  theme_bw(base_size = 14)+
  ylab("SARS-CoV-2 cases" )+
  xlab("Date")+
  coord_cartesian(xlim=c(as.Date("2025-01-11"), as.Date("2025-10-01")))+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b")+
  theme(strip.background = element_rect(fill="white"),
        legend.position = "bottom")

p2 <- ggplot(df_final[df_final$time>as.Date("2025-01-01"),])+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  theme_bw(base_size = 14)+
  ylab("Modelled cases" )+
  xlab("Date")+
  coord_cartesian(xlim=c(as.Date("2025-01-11"), as.Date("2025-10-01")))+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b")+
  theme(strip.background = element_rect(fill="white"),
        legend.position = "bottom")

p3 <- ggplot(df_gr_final[df_gr_final$time>as.Date("2025-01-01"),])+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_hline(yintercept = 0, linetype="dashed")+
  ylab("Growth rate" )+
  xlab("Date")+
  coord_cartesian(xlim=c(as.Date("2025-01-11"), as.Date("2025-10-01")))+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b")+
  theme(strip.background = element_rect(fill="white"),
        legend.position = "bottom")

p1 <- p1 +
  labs(tag="A")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

p2 <- p2 +
  labs(tag="B")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

p3 <- p3 +
  labs(tag="C")


p1 + p2 + p3 +plot_layout(nrow=3)
ggsave(paste('figures','/', 'All_final_draft', '.png', sep=""), width=8, height=8)


################################################################################
# Figure 2 

#data_hosp$admission_date <- as.Date(data_hosp$admission_date)
hosp_inc$time <- as.Date(hosp_inc$time)
hosp_gr$time <- as.Date(hosp_gr$time)


hosp_inc_final <- hosp_inc[hosp_inc$origin == final_date,]
hosp_gr_final <- hosp_gr[hosp_gr$origin == final_date,]

#data_hosp_final <- data_hosp[data_hosp$origin == final_date,]
hosp_inc_final$time <- as.Date(hosp_inc_final$time)


h1a <- ggplot(hosp_inc_final[hosp_inc_final$time>as.Date("2025-01-01") &hosp_inc_final$pathogen=="SARS-CoV-2",])+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  #geom_point(data=data_hosp_final[data_hosp_final$pathogen=="SARS-CoV-2" &data_hosp_final$admission_date>as.Date("2025-01-01"),], aes(x=admission_date, y=hospitalisations), size=0.5)+
  #geom_line(data=data_hosp_final[data_hosp_final$pathogen=="SARS-CoV-2" &data_hosp_final$admission_date>as.Date("2025-01-01"),], aes(x=admission_date, y=hospitalisations), linewidth=0.2)+
  theme_bw(base_size = 14)+
  ylab("Hospitalisations" )+
  xlab("Date")+
  coord_cartesian(xlim=c(as.Date("2025-01-11"), as.Date("2025-10-01")))+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b")+
  theme(strip.background = element_rect(fill="white"),
        legend.position = "bottom")

h2a <- ggplot(hosp_gr_final[hosp_gr_final$time>as.Date("2025-01-01") &hosp_gr_final$pathogen=="SARS-CoV-2",])+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_y_continuous(breaks=c(-0.04,0,0.04))+
  ylab("Growth rate" )+
  xlab("Date")+
  coord_cartesian(xlim=c(as.Date("2025-01-11"), as.Date("2025-10-01")))+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b")+
  theme(strip.background = element_rect(fill="white"),
        legend.position = "bottom")

h1b <- ggplot(hosp_inc_final[hosp_inc_final$time>as.Date("2025-01-01") &hosp_inc_final$pathogen=="Influenza",])+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  #geom_point(data=data_hosp_final[data_hosp_final$pathogen=="Influenza" &data_hosp_final$admission_date>as.Date("2025-01-01"),], aes(x=admission_date, y=hospitalisations), size=0.5)+
  #geom_line(data=data_hosp_final[data_hosp_final$pathogen=="Influenza" &data_hosp_final$admission_date>as.Date("2025-01-01"),], aes(x=admission_date, y=hospitalisations), linewidth=0.2)+
  theme_bw(base_size = 14)+
  ylab("Hospitalisations" )+
  xlab("Date")+
  coord_cartesian(xlim=c(as.Date("2025-01-11"), as.Date("2025-10-01")))+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b")+
  theme(strip.background = element_rect(fill="white"),
        legend.position = "bottom")

h2b <- ggplot(hosp_gr_final[hosp_gr_final$time>as.Date("2025-01-01") &hosp_gr_final$pathogen=="Influenza",])+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_y_continuous(breaks=c(-0.04,0,0.04))+
  ylab("Growth rate" )+
  xlab("Date")+
  coord_cartesian(xlim=c(as.Date("2025-01-11"), as.Date("2025-10-01")))+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b")+
  theme(strip.background = element_rect(fill="white"),
        legend.position = "bottom")

h1c <- ggplot(hosp_inc_final[hosp_inc_final$time>as.Date("2025-01-01") &hosp_inc_final$pathogen=="RSV",])+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  #geom_point(data=data_hosp_final[data_hosp_final$pathogen=="RSV" &data_hosp_final$admission_date>as.Date("2025-01-01"),], aes(x=admission_date, y=hospitalisations), size=0.5)+
  #geom_line(data=data_hosp_final[data_hosp_final$pathogen=="RSV" &data_hosp_final$admission_date>as.Date("2025-01-01"),], aes(x=admission_date, y=hospitalisations), linewidth=0.2)+
  theme_bw(base_size = 14)+
  ylab("Hospitalisations" )+
  xlab("Date")+
  coord_cartesian(xlim=c(as.Date("2025-01-11"), as.Date("2025-10-01")))+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b")+
  theme(strip.background = element_rect(fill="white"),
        legend.position = "bottom")

h2c <- ggplot(hosp_gr_final[hosp_gr_final$time>as.Date("2025-01-01") &hosp_gr_final$pathogen=="RSV",])+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_hline(yintercept = 0, linetype="dashed")+
  ylab("Growth rate" )+
  xlab("Date")+
  scale_y_continuous(breaks=c(-0.04,0,0.04))+
  coord_cartesian(xlim=c(as.Date("2025-01-11"), as.Date("2025-10-01")))+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b")+
  theme(strip.background = element_rect(fill="white"),
        legend.position = "bottom")

h1a <- h1a+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.tag.position = c(0.0,0.98))+
  labs(tag="A")+
  annotate("label",label="SARS-Cov-2", y=Inf, x = as.Date("2025-03-15"), fill= "white", color="black", vjust=1.2, size=5, hjust=0.0)


h2a <- h2a+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())

h1b <- h1b+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.tag.position = c(0.0,0.98))+
  labs(tag="B")+
  annotate("label",label="Influenza", y=Inf, x = as.Date("2025-03-15"), fill= "white", color="black", vjust=1.2, size=5, hjust=0.0)


h2b <- h2b+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())

h1c <- h1c+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.tag.position = c(0.0,0.98))+
  labs(tag="C")+
  annotate("label",label="RSV", y=Inf, x = as.Date("2025-03-15"), fill= "white", color="black", vjust=1.2, size=5, hjust=0.0)


h1a+h2a+h1b+h2b+h1c+h2c+plot_layout(nrow=6, heights=c(1,0.5,1,0.5,1,0.5))

ggsave(paste('figure','/', 'Hosp_final_draft', '.png', sep=""), width=8, height=10)

##################################################################################################################
# Real-time analysis figures (Figures 3&4)
################################################################################

# Save original data (for plotting purposes)
#data_df_og <- data_cov
#data_hosp_og <- data_hosp

################################################################################
## Some formatting to produce variable mask 
## (to avoid plotting modelled outputs over each other)


initial_date <- unq_dates[1]
df <- df_inc[df_inc$time>df_inc$origin-35,]
df_gr <- df_gr[df_gr$time>df_gr$origin-35,]
df_dow <- df_dow[df_dow$time>df_dow$origin-35,]
#data_df <- data_cov[data_cov$time>data_cov$origin-35,]

hosp_gr <- hosp_gr[hosp_gr$time>hosp_gr$origin-35,]
hosp_inc <- hosp_inc[hosp_inc$time>hosp_inc$origin-35,]
#data_hosp <- data_hosp[data_hosp$time>data_hosp$origin-35,]


df$mask <- 0
df_gr$mask <- 0
df_dow$mask <- 0
#data_df$mask <- 0

hosp_inc$mask <- 0
hosp_gr$mask <- 0
#data_hosp$mask <- 0


for(i in 1:length(unq_dates) ){
  print(i)
  
  if(i < length(unq_dates)){
    min_date <- max(df[df$origin == rev(unq_dates)[i +1],]$time)
    max_date <- max(df[df$origin == rev(unq_dates)[i],]$time)
    
    df[df$origin == rev(unq_dates)[i] &df$time>=min_date &df$time>=(max_date-35),]$mask <- i 
    df_gr[df_gr$origin == rev(unq_dates)[i] &df_gr$time>=min_date &df_gr$time>=(max_date-35),]$mask <- i 
    df_dow[df_dow$origin == rev(unq_dates)[i] &df_dow$time>=min_date &df_dow$time>=(max_date-35),]$mask <- i 
    #data_df[data_df$origin == rev(unq_dates)[i] &data_df$time>min_date &data_df$time>=(max_date-35),]$mask <- i 
    
  } else{

    df[df$origin == rev(unq_dates)[i],]$mask <- i 
    df_gr[df_gr$origin == rev(unq_dates)[i],]$mask <- i 
    df_dow[df_dow$origin == rev(unq_dates)[i] ,]$mask <- i 
    #data_df[data_df$origin == rev(unq_dates)[i],]$mask <- i 
    
  }
  
 
  
}


for(i in 1:(length(unq_dates)-2) ){
  print(i)
  for(j in c("SARS-CoV-2","RSV","Influenza")){
    print(j)
    
    if(i<(length(unq_dates)-2)){
      min_date <- max(hosp_inc[hosp_inc$pathogen==j & hosp_inc$origin == rev(unq_dates)[i +1],]$time)
      max_date <- max(hosp_inc[hosp_inc$pathogen==j & hosp_inc$origin == rev(unq_dates)[i],]$time)
      
      hosp_inc[hosp_inc$pathogen==j & hosp_inc$origin == rev(unq_dates)[i] &hosp_inc$time>=min_date &hosp_inc$time>=(max_date-35),]$mask <- i 
      hosp_gr[hosp_gr$pathogen==j & hosp_gr$origin == rev(unq_dates)[i] &hosp_gr$time>=min_date &hosp_gr$time>=(max_date-35),]$mask <- i 
      #data_hosp[data_hosp$pathogen==j & data_hosp$origin == rev(unq_dates)[i] &data_hosp$time>min_date &data_hosp$time>=(max_date-35),]$mask <- i 
      
      
    } else{
      hosp_inc[hosp_inc$pathogen==j & hosp_inc$origin == rev(unq_dates)[i],]$mask <- i 
      hosp_gr[hosp_gr$pathogen==j & hosp_gr$origin == rev(unq_dates)[i],]$mask <- i 
      #data_hosp[data_hosp$pathogen==j & data_hosp$origin == rev(unq_dates)[i],]$mask <- i 
      
    }
  }
  

}


##########################################################
# Figure 3
##########################################################

cols1 <- RColorBrewer::brewer.pal(8,"Dark2")
cols16 <- c(cols1, cols1)

rt2 <- ggplot(df[df$mask>0 &df$time>as.Date("2025-03-01") & df$origin!=final_date,])+
  geom_line(aes(x=time, y=y, col=factor(origin), group=origin))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=factor(origin), group=origin), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=factor(origin), group=origin), alpha=0.2)+
  geom_ribbon(data = df_final[df_final$time>as.Date("2025-03-01"),], aes(x=time, y=y, ymin=lb_95, ymax=ub_95),alpha=0,color="black", linetype="dashed")+
  #geom_point(data = data_df_og[data_df_og$time>as.Date("2025-03-01") & data_df_og$origin %in% unq_dates[1:length(unq_dates)],], aes(x=time, y=cases, col=as.factor(origin), group=time ), shape=16, size=0.8)+
  #geom_line(data = data_df_og[data_df_og$time>as.Date("2025-03-01")& data_df_og$origin %in% unq_dates[1:length(unq_dates)],], aes(x=time, y=cases, col=as.factor(origin), group=time ))+
  #geom_point(data = data_cov_final[data_cov_final$pathogen=="SARS-CoV-2"&data_cov_final$time>as.Date("2025-03-01"),], aes(x=time, y=cases ), col="black",shape=16, size=0.8)+
  #geom_point(data = data_cov_final[data_cov_final$time>as.Date("2025-03-01"),], aes(x=time, y=cases),shape=16, size=0.8, fill="black")+
  scale_color_manual("",values=cols16)+
  scale_fill_manual("",values=cols16)+
  theme_bw(base_size = 14)+
  ylab("SARS-CoV-2 cases" )+
  xlab("Date")+
  coord_cartesian(xlim=c(initial_date-14, final_date-22 ))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b")+
  theme(strip.background = element_rect(fill="white"),
        legend.position = "none")

rt1<-ggplot(df_dow[df_dow$mask>0 &df_dow$time>as.Date("2025-03-01") &df_dow$origin!=final_date,])+
  geom_line(aes(x=time, y=y, col=factor(origin), group=origin))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=factor(origin), group=origin), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=factor(origin), group=origin), alpha=0.2)+
  geom_ribbon(data = df_dow_final[df_dow_final$time>as.Date("2025-03-01"),], aes(x=time, y=y, ymin=lb_95, ymax=ub_95),alpha=0,color="black", linetype="dashed")+
  #geom_point(data = data_df_og[data_df_og$time>as.Date("2025-03-01") & data_df_og$origin %in% unq_dates[1:length(unq_dates)],], aes(x=time, y=cases, col=as.factor(origin), group=time ), shape=16, size=0.8)+
  #geom_line(data = data_df_og[data_df_og$time>as.Date("2025-03-01")& data_df_og$origin %in% unq_dates[1:length(unq_dates)],], aes(x=time, y=cases, col=as.factor(origin), group=time ))+
  #geom_point(data = data_cov_final[data_cov_final$pathogen=="SARS-CoV-2"&data_cov_final$time>as.Date("2025-03-01"),], aes(x=time, y=cases ), col="black",shape=16, size=0.8)+
  #geom_point(data = data_cov_final[data_cov_final$time>as.Date("2025-03-01"),], aes(x=time, y=cases),shape=16, size=0.8, fill="black")+
  scale_color_manual("",values=cols16)+
  scale_fill_manual("",values=cols16)+
  theme_bw(base_size = 14)+
  ylab("SARS-CoV-2 cases" )+
  xlab("Date")+
  coord_cartesian(xlim=c(initial_date-14, final_date-22 ))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b")+
  theme(strip.background = element_rect(fill="white"),
        legend.position = "none")

rt3<-ggplot(df_gr[df_gr$mask>0 &df_gr$time>as.Date("2025-03-01") &df_gr$origin!=final_date,])+
  geom_line(aes(x=time, y=y, col=factor(origin), group=origin))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=factor(origin), group=origin), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=factor(origin), group=origin), alpha=0.2)+
  geom_ribbon(data = df_gr_final[df_gr_final$time>as.Date("2025-03-01"),], aes(x=time, y=y, ymin=lb_95, ymax=ub_95),alpha=0,color="black", linetype="dashed")+
  scale_color_manual("",values=cols16)+
  scale_fill_manual("",values=cols16)+
  geom_hline(yintercept = 0, linetype="dotted")+
  theme_bw(base_size = 14)+
  ylab("Growth rate" )+
  xlab("Date")+
  coord_cartesian(xlim=c(initial_date-14, final_date-22 ))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b")+
  theme(strip.background = element_rect(fill="white"),
        legend.position = "none")

rt3 <- rt3+
  labs(tag="C")
rt1 <- rt1+
  theme(axis.text.x = element_blank(),
                 axis.title.x = element_blank())+
  labs(tag = "A")

rt2 <- rt2+
  theme(axis.text.x = element_blank(),
                 axis.title.x = element_blank())+
  labs(tag = "B")
rt1 + rt2 + rt3 + plot_layout(nrow=3)
ggsave(paste('figure','/', 'real_time_casesALT', '.png', sep=""), width=8, height=8)


################################################################################
## Figure 4

rtH1a <- ggplot(hosp_inc[hosp_inc$pathogen=="SARS-CoV-2"& hosp_inc$mask>0 &hosp_inc$time>as.Date("2025-03-01")&hosp_inc$origin!=final_date,])+
  geom_line(aes(x=time, y=y, col=factor(origin), group=origin))+
  facet_wrap(.~pathogen, scale="free_y", axes = "all")+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=factor(origin), group=origin), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=factor(origin), group=origin), alpha=0.2)+
  geom_ribbon(data = hosp_inc_final[hosp_inc_final$pathogen=="SARS-CoV-2"& hosp_inc_final$time>as.Date("2025-03-01"),], aes(x=time, y=y, ymin=lb_95, ymax=ub_95),alpha=0,color="black", linetype="dashed")+
  #geom_line(data = data_hosp_og[data_hosp_og$pathogen=="SARS-CoV-2" &data_hosp_og$time>as.Date("2025-03-01")& data_hosp_og$origin %in% rev(unq_dates)[1:(length(unq_dates)-1) ],], aes(x=time, y=hospitalisations, col=as.factor(origin), group=time ), linewidth=0.2)+
  #geom_point(data = data_hosp_og[data_hosp_og$pathogen=="SARS-CoV-2" &data_hosp_og$time>as.Date("2025-03-01")& data_hosp_og$origin %in% rev(unq_dates)[1:(length(unq_dates)-2) ],], aes(x=time, y=hospitalisations, col=as.factor(origin), group=time ),shape=16, size=0.8)+
  #geom_point(data = data_hosp_final[data_hosp_final$pathogen=="SARS-CoV-2" &data_hosp_final$time>as.Date("2025-03-01"),], aes(x=time, y=hospitalisations),shape=16, size=0.8, fill="black")+
  scale_color_manual("",values=cols16)+
  scale_fill_manual("",values=cols16)+
  theme_bw(base_size = 14)+
  ylab("Hospitalisations" )+
  xlab("Date")+
  coord_cartesian(xlim=c(initial_date-14, final_date-22 ))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b")+
  theme(strip.background = element_rect(fill="white"),
        legend.position = "none")

rtH1b <- ggplot(hosp_inc[hosp_inc$pathogen=="Influenza"& hosp_inc$mask>0 &hosp_inc$time>as.Date("2025-03-01")&hosp_inc$origin!=final_date,])+
  geom_line(aes(x=time, y=y, col=factor(origin), group=origin))+
  facet_wrap(.~pathogen, scale="free_y", axes = "all")+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=factor(origin), group=origin), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=factor(origin), group=origin), alpha=0.2)+
  geom_ribbon(data = hosp_inc_final[hosp_inc_final$pathogen=="Influenza"& hosp_inc_final$time>as.Date("2025-03-01"),], aes(x=time, y=y, ymin=lb_95, ymax=ub_95),alpha=0,color="black", linetype="dashed")+
  #geom_line(data = data_hosp_og[data_hosp_og$pathogen=="Influenza" &data_hosp_og$time>as.Date("2025-03-01")& data_hosp_og$origin %in% rev(unq_dates)[1:(length(unq_dates)-2)],], aes(x=time, y=hospitalisations, col=as.factor(origin), group=time ))+
  #geom_point(data = data_hosp_og[data_hosp_og$pathogen=="Influenza" &data_hosp_og$time>as.Date("2025-03-01")& data_hosp_og$origin %in% rev(unq_dates)[1:(length(unq_dates)-2)],], aes(x=time, y=hospitalisations, col=as.factor(origin), group=time ),shape=16, size=0.8)+
  #geom_point(data = data_hosp_final[data_hosp_final$pathogen=="Influenza" &data_hosp_final$time>as.Date("2025-03-01"),], aes(x=time, y=hospitalisations ), col="black",shape=16, size=0.8)+
  #geom_point(data = data_hosp_final[data_hosp_final$pathogen=="Influenza" &data_hosp_final$time>as.Date("2025-03-01"),], aes(x=time, y=hospitalisations),shape=16, size=0.8, fill="black")+
  scale_color_manual("",values=cols16)+
  scale_fill_manual("",values=cols16)+
  theme_bw(base_size = 14)+
  ylab("Hospitalisations" )+
  xlab("Date")+
  scale_y_continuous(breaks=c(0,2,4,6,8,10,12))+
  coord_cartesian(xlim=c(initial_date-14, final_date-22 ))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b")+
  theme(strip.background = element_rect(fill="white"),
        legend.position = "none")


rtH1c <- ggplot(hosp_inc[hosp_inc$pathogen=="RSV"& hosp_inc$mask>0 &hosp_inc$time>as.Date("2025-03-01")&hosp_inc$origin!=final_date,])+
  geom_line(aes(x=time, y=y, col=factor(origin), group=origin))+
  facet_wrap(.~pathogen, scale="free_y", axes = "all")+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=factor(origin), group=origin), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=factor(origin), group=origin), alpha=0.2)+
  geom_ribbon(data = hosp_inc_final[hosp_inc_final$pathogen=="RSV"& hosp_inc_final$time>as.Date("2025-03-01"),], aes(x=time, y=y, ymin=lb_95, ymax=ub_95),alpha=0,color="black", linetype="dashed")+
  #geom_line(data = data_hosp_og[data_hosp_og$pathogen=="RSV" &data_hosp_og$time>as.Date("2025-03-01")& data_hosp_og$origin %in% rev(unq_dates)[1:(length(unq_dates)-2)],], aes(x=time, y=hospitalisations, col=as.factor(origin), group=time ))+
  #geom_point(data = data_hosp_og[data_hosp_og$pathogen=="RSV" &data_hosp_og$time>as.Date("2025-03-01")& data_hosp_og$origin %in% rev(unq_dates)[1:(length(unq_dates)-2)],], aes(x=time, y=hospitalisations, col=as.factor(origin), group=time ),shape=16, size=0.8)+
  #geom_point(data = data_hosp_final[data_hosp_final$pathogen=="RSV" &data_hosp_final$time>as.Date("2025-03-01"),], aes(x=time, y=hospitalisations ), col="black",shape=16, size=0.8)+
  #geom_point(data = data_hosp_final[data_hosp_final$pathogen=="RSV" &data_hosp_final$time>as.Date("2025-03-01"),], aes(x=time, y=hospitalisations),shape=16, size=0.8, fill="black")+
  scale_color_manual("",values=cols16)+
  scale_fill_manual("",values=cols16)+
  theme_bw(base_size = 14)+
  ylab("Hospitalisations" )+
  xlab("Date")+
  coord_cartesian(xlim=c(initial_date-14, final_date-22 ))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b")+
  theme(strip.background = element_rect(fill="white"),
        legend.position = "none")

rtH2a <- ggplot(hosp_gr[hosp_gr$pathogen=="SARS-CoV-2" &hosp_gr$mask>0 &hosp_gr$time>as.Date("2025-03-01")&hosp_gr$origin!=final_date,])+
  geom_line(aes(x=time, y=y, col=factor(origin), group=origin))+
  facet_wrap(.~pathogen, axes = "all")+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=factor(origin), group=origin), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=factor(origin), group=origin), alpha=0.2)+
  geom_ribbon(data = hosp_gr_final[hosp_gr_final$pathogen=="SARS-CoV-2" &hosp_gr_final$time>as.Date("2025-03-01"),], aes(x=time, y=y, ymin=lb_95, ymax=ub_95),alpha=0,color="black", linetype="dashed")+
  scale_color_manual("",values=cols16)+
  scale_fill_manual("",values=cols16)+
  geom_hline(yintercept = 0, linetype="dotted")+
  theme_bw(base_size = 14)+
  ylab("Growth rate" )+
  xlab("Date")+
  scale_y_continuous(breaks=c(-0.06,-0.04,-0.02,0,0.02,0.04,0.06))+
  coord_cartesian(xlim=c(initial_date-14, final_date-22 ))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b")+
  theme(strip.background = element_rect(fill="white"),
        legend.position = "none")

rtH2b <- ggplot(hosp_gr[hosp_gr$pathogen=="Influenza" &hosp_gr$mask>0 &hosp_gr$time>as.Date("2025-03-01")&hosp_gr$origin!=final_date,])+
  geom_line(aes(x=time, y=y, col=factor(origin), group=origin))+
  facet_wrap(.~pathogen, axes = "all")+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=factor(origin), group=origin), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=factor(origin), group=origin), alpha=0.2)+
  geom_ribbon(data = hosp_gr_final[hosp_gr_final$pathogen=="Influenza" &hosp_gr_final$time>as.Date("2025-03-01"),], aes(x=time, y=y, ymin=lb_95, ymax=ub_95),alpha=0,color="black", linetype="dashed")+
  scale_color_manual("",values=cols16)+
  scale_fill_manual("",values=cols16)+
  geom_hline(yintercept = 0, linetype="dotted")+
  theme_bw(base_size = 14)+
  ylab("Growth rate" )+
  xlab("Date")+
  scale_y_continuous(breaks=c(-0.06,-0.04,-0.02,0,0.02,0.04,0.06))+
  coord_cartesian(xlim=c(initial_date-14, final_date-22 ))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b")+
  theme(strip.background = element_rect(fill="white"),
        legend.position = "none")

rtH2c <- ggplot(hosp_gr[hosp_gr$pathogen=="RSV" &hosp_gr$mask>0 &hosp_gr$time>as.Date("2025-03-01") &hosp_gr$origin!=final_date,])+
  geom_line(aes(x=time, y=y, col=factor(origin), group=origin))+
  facet_wrap(.~pathogen, axes = "all")+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=factor(origin), group=origin), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=factor(origin), group=origin), alpha=0.2)+
  geom_ribbon(data = hosp_gr_final[hosp_gr_final$pathogen=="RSV" &hosp_gr_final$time>as.Date("2025-03-01"),], aes(x=time, y=y, ymin=lb_95, ymax=ub_95),alpha=0,color="black", linetype="dashed")+
  scale_color_manual("",values=cols16)+
  scale_fill_manual("",values=cols16)+
  geom_hline(yintercept = 0, linetype="dotted")+
  theme_bw(base_size = 14)+
  ylab("Growth rate" )+
  xlab("Date")+
  scale_y_continuous(breaks=c(-0.06,-0.04,-0.02,0,0.02,0.04,0.06))+
  coord_cartesian(xlim=c(initial_date-14, final_date-22 ))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b")+
  theme(strip.background = element_rect(fill="white"),
        legend.position = "none")


rtH1a <- rtH1a+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  labs(tag="A")+
  annotate("label",label="SARS-Cov-2", y=Inf, x = as.Date("2025-09-01"), fill= "white", color="black", vjust=1.2, size=5, hjust=0.0)


rtH1b <- rtH1b+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  labs(tag="B")+
  annotate("label",label="Influenza", y=Inf, x = as.Date("2025-09-01"), fill= "white", color="black", vjust=1.2, size=5, hjust=0.0)


rtH1c <- rtH1c+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  labs(tag="C")+
  annotate("label",label="RSV", y=Inf, x = as.Date("2025-09-01"), fill= "white", color="black", vjust=1.2, size=5, hjust=0.0)


rtH2a <- rtH2a+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

rtH2b <- rtH2b+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

rtH2c <- rtH2c+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())


rtH1a+rtH2a+rtH1b + rtH2b +rtH1c + rtH2c+ plot_layout(nrow=6, heights=c(2,1,2,1,2,1))
ggsave(paste('figure','/', 'real_time_hosps_alt', '.png', sep=""), width=8, height=12)


