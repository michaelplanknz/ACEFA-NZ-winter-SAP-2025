
get_vals <- function(mod_gr, mod_Rt, pathogen, date, state){
  
  gr <- mod_gr[mod_gr$time == date,]
  Rt <- mod_Rt[mod_Rt$time == date,]
  data.frame(gr = gr$y,
             gr_lb = gr$lb_95,
             gr_ub = gr$ub_95,
             dt = log(2)/gr$y,
             dt_lb = log(2)/gr$lb_95,
             dt_ub = log(2)/gr$ub_95,
             p_gr = gr$prop,
             Rt = Rt$y,
             Rt_lb = Rt$lb_95,
             Rt_ub = Rt$ub_95,
             p_Rt = Rt$prop,
             state = state,
             date = date,
             pathogen = pathogen)
  
}

get_all_outputs <- function(df, mod_fit, location, gamma_dist, b, n, tau_max=14, pathogen, date_column="notification_date", dow="Yes"){
  
  # Modelled cases with DOW effect
  if(dow=="Yes"){
    mod_inc_dow <- ps_single_incidence_dow(mod_fit,
                                           week_effect = 7,
                                           DOW = (df$time_index %% 7)+1,
                                           X=df$time_index, num_days=nrow(df), time_labels = df[,date_column])
  } else{
    mod_inc_dow <- data.frame()
  }
  
  
  # Modelled cases
  mod_inc <- ps_single_incidence(mod_fit, df$time_index, num_days=nrow(df), time_labels = df[,date_column])
  
  # Growth rate
  mod_gr <- ps_single_growth_rate(mod_fit, df$time_index, num_days=nrow(df), time_labels = df[,date_column])
  
  # Rt
  mod_Rt <- ps_single_Rt(mod_fit, df$time_index, num_days=nrow(df), time_labels = df[,date_column],
                         tau_max=tau_max,
                         gi_dist = function(x) gammaDist(b=b, n=n, x))
  
  
  ymax1 <- mod_inc[mod_inc$time> (max(mod_inc$time)-180),]$ub_95
  df_tmp <- df[df$notification_date> (max(df$notification_date)-180),]
  ymax2 <- max(df_tmp[df_tmp$cases < max(df_tmp$cases),]$cases)
  
  ymax <- max(ymax1, ymax2) 
  
  plt1 <- ggplot(mod_inc)+
    geom_line(aes(x=time, y=y))+
    geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
    geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
    theme_bw(base_size = 14)+
    geom_point(data=df, aes(x=notification_date, y=cases), size=0.8)+
    geom_line(data=df, aes(x=notification_date, y=cases), linewidth=0.2)+
    ylab(paste(pathogen," cases", sep='') )+
    xlab("Date")+
    coord_cartesian(xlim=c(first_date, origin_date), ylim=c(0,ymax))+
    scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")+
    geom_label(data=data.frame(), aes(x = first_date, y = Inf, label = paste(pathogen," cases", sep='')),hjust=-0,vjust=1.2, fill = "white", size=5)+
    theme( axis.text.x=element_blank(),
           axis.title.x = element_blank())
  
  
  plt2 <- ggplot(mod_gr)+
    geom_line(aes(x=time, y=y))+
    geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
    geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
    theme_bw(base_size = 14)+
    geom_hline(yintercept = 0, linetype="dashed")+
    #geom_hline(yintercept = log(2)/50, linetype="dotted", color='darkgreen')+
    #geom_hline(yintercept = -log(2)/50, linetype="dotted", color='darkgreen')+
    annotate("rect", xmin = as.Date("2020-01-01"), xmax = as.Date("2045-06-01"), ymin = -log(2)/50, ymax = log(2)/50, color="red3",
             fill = NA,alpha = 0.1, linetype = "dotted")+
    scale_y_continuous(
      "Growth rate", 
      sec.axis = sec_axis(~., breaks=c(log(2)/7, log(2)/14, log(2)/21,log(2)/50,0,-log(2)/50,-log(2)/21, -log(2)/14, -log(2)/7), labels = c( "7", "14","21","50", expression(infinity/-infinity),"-50","-21", "-14", "-7"), name = "Doubling(+) / Halving(-) time (days)")
    )+
    xlab("Date")+
    coord_cartesian(ylim=c(-0.1,0.1), xlim=c(first_date, origin_date))+
    scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")+
    geom_label(data=data.frame(), aes(x = first_date, y = Inf, label ="Growth rate"),hjust=-0,vjust=1.2, fill = "white", size=5)+
    theme( axis.text.x=element_blank(),
           axis.title.x = element_blank())
  
  plt3 <- ggplot(mod_Rt)+
    geom_line(aes(x=time, y=y))+
    geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
    geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
    theme_bw(base_size = 14)+
    geom_hline(yintercept = 1, linetype="dashed")+
    ylab("Effective reproduction number")+
    xlab("Date")+
    coord_cartesian(xlim=c(first_date, origin_date))+
    scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")+
    geom_label(data=data.frame(), aes(x = first_date, y = Inf, label = paste(expression(R[t])) ),parse=T,hjust=-0,vjust=1.2, fill = "white", size=5)+
    theme(axis.title.x = element_blank())
  
  
  
  plt1 / plt2 / plt3 
  ggsave(paste('figure/', origin_date,'/', location, pathogen, '.png', sep=""), width=10, height=10)
  
  mod_inc$location <- location
  mod_inc_dow$location <- location
  mod_gr$location <- location
  mod_Rt$location <- location
  
  return(list(mod_inc, mod_inc_dow, mod_gr, mod_Rt))
}

get_all_outputs_hosp <- function(df, mod_fit, location, gamma_dist, b, n, tau_max=14, pathogen, date_column="admission_date", dow="Yes"){
  
  # Modelled cases with DOW effect
  if(dow=="Yes"){
    mod_inc_dow <- ps_single_incidence_dow(mod_fit,
                                           week_effect = 7,
                                           DOW = (df$time_index %% 7)+1,
                                           X=df$time_index, num_days=nrow(df), time_labels = df[,date_column])
  } else{
    mod_inc_dow <- data.frame()
  }
  
  
  # Modelled cases
  mod_inc <- ps_single_incidence(mod_fit, df$time_index, num_days=nrow(df), time_labels = df[,date_column])
  
  # Growth rate
  mod_gr <- ps_single_growth_rate(mod_fit, df$time_index, num_days=nrow(df), time_labels = df[,date_column])
  
  
  
  plt1 <- ggplot(mod_inc)+
    geom_line(aes(x=time, y=y))+
    geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
    geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
    theme_bw(base_size = 14)+
    geom_point(data=df, aes(x=notification_date, y=cases), size=0.8)+
    geom_line(data=df, aes(x=notification_date, y=cases), linewidth=0.2)+
    ylab(paste(pathogen," hospitalisations", sep='') )+
    xlab("Date")+
    coord_cartesian(xlim=c(first_date, origin_date), ylim=c(0,max(mod_inc[mod_inc$time> (max(mod_inc$time)-180),]$ub_95, df[df$notification_date> (max(df$notification_date)-180),]$cases )))+
    scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")+
    geom_label(data=data.frame(), aes(x = first_date, y = Inf, label = paste(pathogen," hospitalisations", sep='')),hjust=-0,vjust=1.2, fill = "white", size=5)+
    theme( axis.text.x=element_blank(),
           axis.title.x = element_blank())
  
  
  plt2 <- ggplot(mod_gr)+
    geom_line(aes(x=time, y=y))+
    geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
    geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
    theme_bw(base_size = 14)+
    geom_hline(yintercept = 0, linetype="dashed")+
    #geom_hline(yintercept = log(2)/50, linetype="dotted", color='darkgreen')+
    #geom_hline(yintercept = -log(2)/50, linetype="dotted", color='darkgreen')+
    annotate("rect", xmin = as.Date("2020-01-01"), xmax = as.Date("2045-06-01"), ymin = -log(2)/50, ymax = log(2)/50, color="red3",
             fill = NA,alpha = 0.1, linetype = "dotted")+
    scale_y_continuous(
      "Growth rate", 
      sec.axis = sec_axis(~., breaks=c(log(2)/7, log(2)/14, log(2)/21,log(2)/50,0,-log(2)/50,-log(2)/21, -log(2)/14, -log(2)/7), labels = c( "7", "14","21","50", expression(infinity/-infinity),"-50","-21", "-14", "-7"), name = "Doubling(+) / Halving(-) time (days)")
    )+
    xlab("Date")+
    coord_cartesian(ylim=c(-0.1,0.1), xlim=c(first_date, origin_date))+
    scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")+
    geom_label(data=data.frame(), aes(x = first_date, y = Inf, label ="Growth rate"),hjust=-0,vjust=1.2, fill = "white", size=5)+
    theme( axis.text.x=element_blank(),
           axis.title.x = element_blank())
  
  plt1 / plt2
  ggsave(paste('figure/', origin_date,'/', location, pathogen,'hosp', '.png', sep=""), width=10, height=10)
  
  mod_inc$location <- location
  #mod_inc_dow$location <- location
  mod_gr$location <- location
  
  mod_inc$pathogen <- pathogen
  mod_gr$pathogen <- pathogen

  
  return(list(mod_inc, mod_inc_dow, mod_gr))
}


get_all_outputs_inf <- function(df, mod_fit, location, gamma_dist, b, n, tau_max=14, pathogen, num_path=3, pathogen_names=c("Influenza A H3N2", "Influenza A H1N1", "Influenza B")){
  
  
  mod_inc <- ps_incidence(mod_fit, df$time_index, num_days=nrow(df), time_labels = df[,date_column],
                          num_path = num_path,
                          pathogen_names = pathogen_names)
  
  mod_inc_dow <- ps_incidence_dow(mod_fit, X=df$time_index, num_days=nrow(df), time_labels = df[,date_column],
                                  week_effect = 7,
                                  DOW = (df$time_index %% 7)+1,
                                  num_path = num_path,
                                  pathogen_names = pathogen_names)
  
  
  mod_gr <- ps_growth_rate(mod_fit, df$time_index, num_days=nrow(df), time_labels = df[,date_column],
                           num_path = num_path,
                           pathogen_names = pathogen_names)
  
  if(num_path==3){
    mod_prop <- ps_proportion(mod_fit, df$time_index, num_days=nrow(df), time_labels = df[,date_column],
                              num_path = num_path,
                              comb_num=list(c(1,2), c(3), c(1), c(2), c(1)),
                              comb_den=list(c(1,2,3), c(1,2,3), c(1,2,3), c(1,2,3), c(1,2)),
                              comb_names=c("Influenza A", "Influenza B", "Influenza A H3N2", "Influenza A H1N1", "H3N2 vs B"))
  } else if(num_path==2){
    mod_prop <- ps_proportion(mod_fit, df$time_index, num_days=nrow(df), time_labels = df[,date_column],
                              num_path = num_path,
                              comb_num=list(c(1), c(2)),
                              comb_den=list(c(1,2), c(1,2)),
                              comb_names=c("Influenza A", "Influenza B"))
  }
  
  
  mod_Rt <- ps_Rt(mod_fit, df$time_index, num_days=nrow(df), time_labels = df[,date_column],
                  num_path = num_path,
                  tau_max=tau_max,
                  gi_dist = function(x) gammaDist(b=b, n=n, x),
                  pathogen_names = pathogen_names)
  
  
  
  inf0 <- ggplot(mod_inc_dow[mod_inc_dow$pathogen=="Total",])+
    geom_line(aes(x=time, y=y))+
    geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
    geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
    theme_bw(base_size = 14)+
    geom_point(data=df, aes(x=notification_date, y=cases), size=0.8)+
    geom_line(data=df, aes(x=notification_date, y=cases), linewidth=0.2)+
    ylab("Influenza cases")+
    xlab("Date")+
    coord_cartesian(xlim=c(first_date, origin_date))+
    scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")+
    theme( axis.text.x=element_blank(),
           axis.title.x = element_blank())
  
  inf1 <- ggplot(mod_inc[mod_inc$pathogen=="Total",])+
    geom_line(aes(x=time, y=y))+
    geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
    geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
    theme_bw(base_size = 14)+
    geom_point(data=df, aes(x=notification_date, y=cases), size=0.8)+
    geom_line(data=df, aes(x=notification_date, y=cases), linewidth=0.2)+
    ylab("Influenza cases")+
    xlab("Date")+
    coord_cartesian(xlim=c(first_date, origin_date), ylim=c(0,max(mod_inc[mod_inc$time> (max(mod_inc$time)-180),]$ub_95, df[df$notification_date> (max(df$notification_date)-180),]$cases )))+
    scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")+
    geom_label(data=data.frame(), aes(x = first_date, y = Inf, label ="Influenza cases"),hjust=-0,vjust=1.2, fill = "white", size=5)+
    theme( axis.text.x=element_blank(),
           axis.title.x = element_blank())
  
  inf2 <- ggplot(mod_gr[mod_gr$pathogen=="Total",])+
    geom_line(aes(x=time, y=y))+
    geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
    geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
    theme_bw(base_size = 14)+
    geom_hline(yintercept = 0, linetype="dashed")+
    annotate("rect", xmin = as.Date("2020-01-01"), xmax = as.Date("2045-06-01"), ymin = -log(2)/50, ymax = log(2)/50, color="red3",
             fill = NA,alpha = 0.1, linetype = "dotted")+
    scale_y_continuous(
      "Growth rate", 
      sec.axis = sec_axis(~., breaks=c(log(2)/7, log(2)/14, log(2)/21,log(2)/50,0,-log(2)/50,-log(2)/21, -log(2)/14, -log(2)/7), labels = c( "7", "14","21","50", expression(infinity/-infinity),"-50","-21", "-14", "-7"), name = "Doubling(+) / Halving(-) time (days)")
    )+
    xlab("Date")+
    coord_cartesian(ylim=c(-0.1,0.1), xlim=c(first_date, origin_date))+
    geom_label(data=data.frame(), aes(x = first_date, y = Inf, label ="Growth rate"),hjust=-0,vjust=1.2, fill = "white", size=5)+
    scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")
  
  
  inf3 <- ggplot(mod_Rt[mod_Rt$pathogen=="Total",])+
    geom_line(aes(x=time, y=y))+
    geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
    geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
    theme_bw(base_size = 14)+
    geom_hline(yintercept = 1, linetype="dashed")+
    ylab("Effective reproduction number")+
    xlab("Date")+
    coord_cartesian(xlim=c(first_date, origin_date))+
    scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")+
    geom_label(data=data.frame(), aes(x = first_date, y = Inf, label = paste(expression(R[t])) ),parse=T,hjust=-0,vjust=1.2, fill = "white", size=5)+
    theme(axis.title.x = element_blank())
  
  
  infSub0 <- ggplot(mod_inc_dow[mod_inc_dow$pathogen!="Total",])+
    geom_line(aes(x=time, y=y, color=pathogen))+
    geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=pathogen), alpha=0.2)+
    geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=pathogen), alpha=0.2)+
    theme_bw(base_size = 14)+
    scale_fill_brewer(palette = "Dark2")+
    scale_color_brewer(palette = "Dark2")+
    theme(legend.position = "none")+
    ylab("Modelled influenza cases")+
    xlab("Date")+
    coord_cartesian(xlim=c(first_date, origin_date))+
    scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")+
    theme( axis.text.x=element_blank(),
           axis.title.x = element_blank())
  
  infSub1 <- ggplot(mod_inc[mod_inc$pathogen!="Total",])+
    geom_line(aes(x=time, y=y, color=pathogen))+
    geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=pathogen), alpha=0.2)+
    geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=pathogen), alpha=0.2)+
    theme_bw(base_size = 14)+
    scale_fill_brewer(palette = "Dark2")+
    scale_color_brewer(palette = "Dark2")+
    theme(legend.position = "none")+
    ylab("Modelled influenza cases")+
    xlab("Date")+
    coord_cartesian(xlim=c(first_date, origin_date), ylim=c(0,max(mod_inc[mod_inc$time> (max(mod_inc$time)-180),]$ub_95, df[df$notification_date> (max(df$notification_date)-180),]$cases )))+
    scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")+
    geom_label(data=data.frame(), aes(x = first_date, y = Inf, label ="Influenza cases"),hjust=-0,vjust=1.2, fill = "white", size=5)+
    theme( axis.text.x=element_blank(),
           axis.title.x = element_blank())
  
  infSub2 <- ggplot(mod_gr[mod_gr$pathogen!="Total",])+
    geom_line(aes(x=time, y=y, color=pathogen))+
    geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=pathogen), alpha=0.2)+
    geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=pathogen), alpha=0.2)+
    theme_bw(base_size = 14)+
    scale_fill_brewer(palette = "Dark2")+
    scale_color_brewer(palette = "Dark2")+
    geom_hline(yintercept = 0, linetype="dashed")+
    scale_y_continuous(
      "Growth rate", 
      sec.axis = sec_axis(~., breaks=c(log(2)/7, log(2)/14, log(2)/21,0, -log(2)/21, -log(2)/14, -log(2)/7), labels = c( "7", "14","21", expression(infinity/-infinity),"-21", "-14", "-7"), name = "Doubling(+) / Halving(-) time (days)")
    )+
    xlab("Date")+
    coord_cartesian(ylim=c(-0.12,0.12), xlim=c(first_date, origin_date))+
    scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")+
    geom_label(data=data.frame(), aes(x = first_date, y = Inf, label ="Growth rate"),hjust=-0,vjust=1.2, fill = "white", size=5)+
    theme(legend.position = "none",
          axis.text.x=element_blank(),
          axis.title.x = element_blank())
  
  infSub3 <- ggplot(mod_Rt[mod_Rt$pathogen!="Total",])+
    geom_line(aes(x=time, y=y, color=pathogen))+
    geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=pathogen), alpha=0.2)+
    geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=pathogen), alpha=0.2)+
    theme_bw(base_size = 14)+
    geom_hline(yintercept = 1, linetype="dashed")+
    scale_fill_brewer("Sub-type", palette = "Dark2")+
    scale_color_brewer("Sub-type", palette = "Dark2")+
    theme(legend.position = "bottom",
          axis.title.x = element_blank())+
    ylab("Effective reproduction number")+
    xlab("Date")+
    coord_cartesian(xlim=c(first_date, origin_date))+
    geom_label(data=data.frame(), aes(x = first_date, y = Inf, label = paste(expression(R[t])) ),parse=T,hjust=-0,vjust=1.2, fill = "white", size=5)+
    scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")
  
  
  inf1 / inf2 / inf3
  ggsave(paste('figure/', origin_date,'/', location,"-Influenza", '.png', sep=""), width=10, height=10)
  
  infSub1 / infSub2 / infSub3
  ggsave(paste('figure/', origin_date,'/', location,"-InfluenzaSub", '.png', sep=""), width=10, height=10)
  
  mod_inc$location <- location
  mod_inc_dow$location <- location
  mod_gr$location <- location
  mod_Rt$location <- location
  mod_prop$location <- location
  
  return(list(mod_inc, mod_inc_dow, mod_gr, mod_Rt, mod_prop))
}

slightly_formatted_table <- function(df){
  
  df_tab <- data.frame()
  for(i in 1:nrow(df)){
    
    gr=paste(round.cust(df$gr[i],digits=2), " (", round.cust(df$gr_lb[i],digits = 2 ),", ",round.cust(df$gr_ub[i], digits=2), ")",sep="")
    
    
    format_dt <- function(dt){
      dt = round(dt, digits=0)
      #if(dt>100){
      #  dt = ">100"
      #} else if(dt< -100){
      #  dt = "<-100"
      #} else{
      #  dt = round(dt, digits=0)
      #}
    }
    d_time=paste(format_dt(df$dt[i]), " (", format_dt(df$dt_lb[i]),", ",format_dt(df$dt_ub[i]), ")",sep="")
    
    format_Pr <- function(p){
      p = round(p, digits=2)
      if(p>0.99){
        p = ">0.99"
      } else if(p<0.01){
        p = "<0.01"
      } else{
        p = p
      }
    }
    prob_r =format_Pr(df$p_gr[i])
    
    Rt=paste(round(df$Rt[i],digits=2), " (", round(df$Rt_lb[i],digits = 2 ),", ",round(df$Rt_ub[i], digits=2), ")",sep="")
    
    
    prob_Rt =format_Pr(df$p_Rt[i])
    
    if(length(df$prop)==0){
      
      row <- data.frame(Pathogen = df$pathogen[i],
                        Jurisdiction = df$state[i],
                        MaxDate = df$date[i],
                        gr=gr,
                        d_time=d_time,
                        prob_r = prob_r,
                        Rt = Rt,
                        prob_Rt = prob_Rt)
      
    } else{
      prop = format_Pr(df$prop[i])
      row <- data.frame(Pathogen = df$pathogen[i],
                        Jurisdiction = df$state[i],
                        MaxDate = df$date[i],
                        gr=gr,
                        d_time=d_time,
                        prob_r = prob_r,
                        prop = prop,
                        Rt = Rt,
                        prob_Rt = prob_Rt)
      
    }
    
    
    
    df_tab <- rbind(df_tab, row)
    
    
  }
  df_tab
  
  
}
