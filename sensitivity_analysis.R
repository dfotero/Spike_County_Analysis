create_plots_sd <- function(geographical_data)
{
  
  lag       <- 30
  threshold <- c(1,1.25,1.5,1.75,2,2.25,2.5,2.75,3)
  influence <- 0.3
  
  y <- geographical_data$Total_Overdose_Deaths
  
  theResult<-0
  
  for(i in 1:length(threshold))
  {
    result <- ThresholdingAlgo(y,lag,threshold[i],influence)
    df <- data.frame(date = geographical_data$DOD_4_FD, 
                   y = y, 
                   avgFilter = result$avgFilter, 
                   upperThreshold = result$avgFilter + threshold[i] * result$stdFilter, 
                   signals = result$signals,
                   YOD = geographical_data$YOD)
    df$year <- year(df$date)
    df_modified <- df
    df_modified$modified_y <- ifelse(df$signals == 0, 0, df$y)
    df_modified$y_minus_avgFilter_or_zero <- ifelse(df$signals == 0, 0, df$y - df$avgFilter)
    df_modified$threshold<-threshold[i]
    
    if(i==1)
    {
      theResult<-df_modified
    }
    else
    {
      theResult<-rbind(theResult,df_modified)
    }  
  }
  
  for (i in 1:length(threshold)) {
    variable_name <- paste0("plot", i)
    thedf<-theResult[theResult$year==2022 & theResult$threshold==threshold[i],]
    long_df <- thedf %>%
      gather(key = "type", value = "value", -date, -signals, -year, -y_minus_avgFilter_or_zero)
    plot <- ggplot(long_df, aes(x = date, y = y_minus_avgFilter_or_zero)) +
      geom_step(color = "red") +
      scale_x_date(date_breaks = "1 month", date_labels = "%b") +
      scale_y_continuous(limits = c(-10, 10)) +
      labs(x = "", y = "") +
      theme_minimal() +
      ggtitle(threshold[i])
    assign(variable_name,plot)
  }
  combined_plot<-plot1 / plot2/ plot3/ plot4/ plot5 / plot6/ plot7/ plot8/ plot9 
  #combined_plot<-plot1 / plot2/ plot3/ plot4
  
  return(combined_plot)
 
}

ThresholdingAlgo2 <- function(y,lag,threshold,influence) {
  signals <- rep(0,length(y))
  filteredY <- y[1:lag]
  avgFilter <- NULL
  stdFilter <- NULL
  spikes<-rep(0,length(y))
  dips<-rep(0,length(y))
  avgFilter[lag] <- mean(y[1:lag], na.rm=TRUE)
  stdFilter[lag] <- sd(y[1:lag], na.rm=TRUE)
  for (i in (lag+1):length(y)){
    if (abs(y[i]-avgFilter[i-1]) > threshold*stdFilter[i-1]) {
      if (y[i] > avgFilter[i-1]) {
        signals[i] <- 1
        spikes[i]<-y[i]
        dips[i]<-0
        } 
      else {
        signals[i] <- -1
        dips[i]<-y[i]
        spikes[i]<-0
      }
      filteredY[i] <- influence*y[i]+(1-influence)*filteredY[i-1]
    } else {
      signals[i] <- 0
      spikes[i]<-0
      dips[i]<-0
      filteredY[i] <- y[i]
    }
    avgFilter[i] <- mean(filteredY[(i-lag+1):i], na.rm=TRUE)
    stdFilter[i] <- sd(filteredY[(i-lag+1):i], na.rm=TRUE)
  }
  return(list("signals"=signals,"avgFilter"=avgFilter,"stdFilter"=stdFilter,"spikes"=spikes,"dips"=dips))
}

sens_analysis <- function(geographical_data)
{
  
  lag       <- c(7,30,90)
  threshold<-2
  
  #lag<-c(7,15,30,60,90,120,150,180)
  #threshold <- c(2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3)
  influence<-0.3
  #influence<-c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95)
  
  y <- geographical_data$Total_Overdose_Deaths
  
  theResult<-0
  
  for(i in 1:length(lag))
  {
    for(j in 1:length(threshold))
    {
      for( k in 1:length(influence))
      {
        result <- ThresholdingAlgo2(y,lag[i],threshold[j],influence[k])
        df <- data.frame(date = geographical_data$DOD_4_FD, 
                         y = y, 
                         avgFilter = result$avgFilter, 
                         upperThreshold = result$avgFilter + threshold[j] * result$stdFilter, 
                         signals = result$signals,
                         spikes=result$spikes,
                         dips=result$dips,
                         YOD = geographical_data$YOD)
        df$year <- year(df$date)
        df_modified <- df
        df_modified$modified_y <- ifelse(df$signals == 0, 0, df$y)
        df_modified$y_minus_avgFilter_or_zero <- ifelse(df$signals == 0, 0, df$y - df$avgFilter)
        df_modified$lag<-lag[i]
        df_modified$threshold<-threshold[j]
        df_modified$influence<-influence[k]
        df_modified$percSpike<-ifelse(df_modified$signals==1,df_modified$y_minus_avgFilter_or_zero/df_modified$avgFilter,0)
        df_modified<-df_modified[df_modified$year>=2017 & df_modified$year<=2023,]
        #df_modified<-df_modified %>%
        #group_by(year,lag,threshold,influence) %>%
        #group_by(lag,threshold,influence) %>%
        #summarise(Count = sum(y_minus_avgFilter_or_zero > 0),Average_spike=mean(spikes[spikes > 0]),Average_dip=mean(dips[dips>0]),Total_spikes=sum(spikes),Total_deaths=sum(y),Average_per_inc_spike=mean(percSpike[percSpike>0]),Average_diff=mean(y_minus_avgFilter_or_zero[y_minus_avgFilter_or_zero>0]))
        #df_modified$perc_sike_death<-df_modified$Total_spikes/df_modified$Total_deaths
        if(i==1 & j==1 & k==1)
        {
          theResult<-df_modified
        }
        else
        {
          print(paste(i,j,k))
          theResult<-rbind(theResult,df_modified)
        }   
      }
    }

  }
  
   #combined_plot<-plot1 / plot2/ plot3/ plot4
  #countS<-theResult %>% group_by(year,lag,threshold,influence) %>% summarise(Count = sum(y_minus_avgFilter_or_zero != 0))
  
  #plot<-ggplot(countS, aes(x = threshold, y = Count, group = year, colour = year)) + geom_line() + theme_minimal() + labs(title = "Number of spikes and dips varying threshold value",
         #x = "Threshold",
         #y = "Number of Spikes and Dips ",
         #colour = "Year")
  
  
  #return(countS)
  return(theResult)
  
}

addDrugsMA<-function(drug_categories,sens)
{
  #theDrugs<-drug_categories %>% group_by(DOD_4_FD,County_name)%>%summarise(ps_only=sum(ps_only),coc_only=sum(coc_only),stim_only=sum(stim_only),opioid_and_stim_only=sum(opioid_and_stim_only),stim_involved=sum(stim_involved),opioid_and_stim_involved=sum(opioid_and_stim_involved))
  theDrugs<-drug_categories %>% group_by(DOD_4_FD)%>%summarise(ps_only=sum(ps_only),stim_only=sum(stim_only),opioid_and_stim_only=sum(opioid_and_stim_only),stim_involved=sum(stim_involved),opioid_and_stim_involved=sum(opioid_and_stim_involved),fentanyl_involved=sum(fentanyl),heroin_involved=sum(heroin),rx_opioids_involved=sum(rx_opioids))
  
  theData<- sens %>% left_join(theDrugs,by=c("date"="DOD_4_FD"))
  
  theData$ps_only<-ifelse(is.na(theData$ps_only),0,theData$ps_only)
  theData$stim_only<-ifelse(is.na(theData$stim_only),0,theData$stim_only)
  theData$opioid_and_stim_only<-ifelse(is.na(theData$opioid_and_stim_only),0,theData$opioid_and_stim_only)
  theData$stim_involved<-ifelse(is.na(theData$stim_involved),0,theData$stim_involved)
  theData$opioid_and_stim_involved<-ifelse(is.na(theData$opioid_and_stim_involved),0,theData$opioid_and_stim_involved)
  theData$fentanyl_involved<-ifelse(is.na(theData$fentanyl_involved),0,theData$fentanyl_involved)
  theData$heroin_involved<-ifelse(is.na(theData$heroin_involved),0,theData$heroin_involved)
  theData$rx_opioids_involved<-ifelse(is.na(theData$rx_opioids_involved),0,theData$rx_opioids_involved)
  
  theData$non_stim<-theData$y-theData$stim_involved
  theData$stim_plus<-theData$stim_involved-theData$stim_only
  theData$non_fen<-theData$y-theData$fentanyl_involved
  theData$non_her<-theData$y-theData$heroin_involved
  theData$non_rx<-theData$y-theData$rx_opioids_involved
  
  return(theData)
}

sens_analysis_which_counties <- function(geographical_data,counties)
{
  
  #lag       <- 30
  #threshold<-2.5
  
  lag<-c(7,15,30,60,90,120,150,180)
  threshold <- c(2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3)
  influence<-0.3
  #influence<-c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95)
  y <- geographical_data$Total_Overdose_Deaths
  
  theResult<-0
  
  for(i in 1:length(lag))
  {
    for(j in 1:length(threshold))
    {
      for( k in 1:length(influence))
      {
        result <- ThresholdingAlgo2(y,lag[i],threshold[j],influence[k])
        df <- data.frame(date = geographical_data$DOD_4_FD, 
                         y = y, 
                         avgFilter = result$avgFilter, 
                         upperThreshold = result$avgFilter + threshold[j] * result$stdFilter, 
                         signals = result$signals,
                         spikes=result$spikes,
                         dips=result$dips,
                         YOD = geographical_data$YOD)
        df$year <- year(df$date)
        df_modified <- df
        df_modified$modified_y <- ifelse(df$signals == 0, 0, df$y)
        df_modified$y_minus_avgFilter_or_zero <- ifelse(df$signals == 0, 0, df$y - df$avgFilter)
        df_modified$lag<-lag[i]
        df_modified$threshold<-threshold[j]
        df_modified$influence<-influence[k]
        df_modified<-df_modified[df_modified$signals==1,]
        df_modified<- df_modified %>% inner_join(counties,by=c("date"="DOD_4_FD"))
        df_modified<-df_modified%>% group_by(year,County_name,lag,threshold,influence) %>% summarise(Tot_deaths = sum(Total_Overdose_Deaths))
        
        if(i==1 & j==1 & k==1)
        {
          theResult<-df_modified
        }
        else
        {
          print(paste(i,j,k))
          theResult<-rbind(theResult,df_modified)
        }   
      }
    }
    
  }
  
  #combined_plot<-plot1 / plot2/ plot3/ plot4
  #countS<-theResult %>% group_by(year,lag,threshold,influence) %>% summarise(Count = sum(y_minus_avgFilter_or_zero != 0))
  
  #plot<-ggplot(countS, aes(x = threshold, y = Count, group = year, colour = year)) + geom_line() + theme_minimal() + labs(title = "Number of spikes and dips varying threshold value",
  #x = "Threshold",
  #y = "Number of Spikes and Dips ",
  #colour = "Year")
  
  
  #return(countS)
  return(theResult)
  
}

percStatTestMA<-function(spikes,lag,threshold,year)
{
  theData<-spikes[spikes$lag==lag & spikes$threshold==threshold & spikes$year==year,]
  
  theData$perc_stim<-theData$stim_involved/theData$y
  theData$perc_fentanyl<-theData$fentanyl_involved/theData$y
  theData$perc_heroin<-theData$heroin_involved/theData$y
  theData$perc_rx<-theData$rx_opioids_involved/theData$y
  
  theData$signals<-ifelse(theData$signals>0,1,0)
  
  dataset<-theData
  #theRes<-data.frame("County_name"="","num_spike"=0,"num_non_spike"=0,"num_spike_g0"=0,"num_non_spike_g0"=0,"mean_stim_spike"=0,"mean_stim_non_spike"=0,"sd_stim_spike"=0,"sd_stim_non_spike"=0,"Ttest_stim"=0,"mean_fen_spike"=0,"mean_fen_non_spike"=0,"sd_fen_spike"=0,"sd_fen_non_spike"=0,"Ttest_fen"=0,"mean_her_spike"=0,"mean_her_non_spike"=0,"sd_her_spike"=0,"sd_her_non_spike"=0,"Ttest_her"=0)
  theRes<-data.frame("Drug"="","min"=0,"max"=0,"pvalue"=0)
  ind<-1
  dataSpike<-dataset[dataset$signals==1,]
  datanonSpike<-dataset[dataset$signals==0,]
    
  numDataSpike<-nrow(dataSpike)
  numDataNonSpike<-nrow(datanonSpike)
    
  dataset<-dataset[!is.nan(dataset$perc_stim),]
    
  dataSpike<-dataset[dataset$signals==1,]
  datanonSpike<-dataset[dataset$signals==0,]
    
  numDataSpikenon0<-nrow(dataSpike)
  numDataNonSpikenon0<-nrow(datanonSpike)

  theTestStim<-t.test(dataSpike$perc_stim, datanonSpike$perc_stim, var.equal = FALSE) 
  theTestFen<-t.test(dataSpike$perc_fen, datanonSpike$perc_fen, var.equal = FALSE) 
  theTestHer<-t.test(dataSpike$perc_her, datanonSpike$perc_her, var.equal = FALSE) 
  theTestRx<-t.test(dataSpike$perc_rx, datanonSpike$perc_rx, var.equal = FALSE) 
        #theRes[ind,]<-c(theCounties[i],"Stimulants",theTestStim$conf.int[1],theTestStim$conf.int[2],theTestStim$p.value)
  theRes[ind,]<-c("Stimulants",theTestStim$estimate[1],theTestStim$estimate[2],theTestStim$p.value)
  ind<-ind+1
  theRes[ind,]<-c("Fentanyl",theTestFen$estimate[1],theTestFen$estimate[2],theTestFen$p.value)
  ind<-ind+1
  theRes[ind,]<-c("Heroin",theTestHer$estimate[1],theTestHer$estimate[2],theTestHer$p.value)
  ind<-ind+1
  theRes[ind,]<-c("Rx_Opioids",theTestRx$estimate[1],theTestRx$estimate[2],theTestRx$p.value)
  ind<-ind+1

  
  return(theRes)
}

graphPercDrugs<-function(spikes,lag,threshold)
{
  theData<-spikes[spikes$lag==lag & spikes$threshold==threshold & spikes$year>=2017 & spikes$year<=2022,]
  theData$stim_spike<-ifelse(theData$signals==1,theData$stim_involved,0)
  theData$fentanyl_spike<-ifelse(theData$signals==1,theData$fentanyl_involved,0)
  theData$heroin_spike<-ifelse(theData$signals==1,theData$heroin_involved,0)
  theData$rx_spike<-ifelse(theData$signals==1,theData$rx_opioids_involved,0)
  
  theData$stim_non_spike<-ifelse(theData$signals==1,0,theData$stim_involved)
  theData$fentanyl_non_spike<-ifelse(theData$signals==1,0,theData$fentanyl_involved)
  theData$heroin_non_spike<-ifelse(theData$signals==1,0,theData$heroin_involved)
  theData$rx_non_spike<-ifelse(theData$signals==1,0,theData$rx_opioids_involved)
  
  theData$non_spike_deaths<-theData$y-theData$spikes
  
  theData<-theData%>%group_by(year)%>%summarise(total_deaths=sum(y),deaths_spikes=sum(spikes),deaths_non_spikes=sum(non_spike_deaths),total_stim_involved=sum(stim_involved),total_fentanyl_involved=sum(fentanyl_involved),total_heroin_involved=sum(heroin_involved),total_rx_involved=sum(rx_opioids_involved),stim_spikes=sum(stim_spike),fentanyl_spikes=sum(fentanyl_spike),heroin_spikes=sum(heroin_spike),rx_spikes=sum(rx_spike),stim_non_spikes=sum(stim_non_spike),fentanyl_non_spikes=sum(fentanyl_non_spike),heroin_non_spikes=sum(heroin_non_spike),rx_non_spikes=sum(rx_non_spike))
  
  #theData$perc_deaths_spikes<-theData$deaths_spikes/theData$total_deaths
  
  theData$perc_stim_spike<-theData$stim_spikes/theData$deaths_spikes
  theData$perc_fentanyl_spike<-theData$fentanyl_spikes/theData$deaths_spikes
  theData$perc_heroin_spike<-theData$heroin_spikes/theData$deaths_spikes
  theData$perc_rx_spike<-theData$rx_spikes/theData$deaths_spikes
  
  
  theData$perc_stim_non_spike<-theData$stim_non_spikes/theData$deaths_non_spikes
  theData$perc_fentanyl_non_spike<-theData$fentanyl_non_spikes/theData$deaths_non_spikes
  theData$perc_heroin_non_spike<-theData$heroin_non_spikes/theData$deaths_non_spikes
  theData$perc_rx_non_spike<-theData$rx_non_spikes/theData$deaths_non_spikes
  
  
  df_long <- theData %>% pivot_longer(
    #cols = c(perc_stim,perc_fentanyl,perc_stim_spike,perc_fentanyl_spike), # Only specify columns to be transformed
    #cols = c(stim_only,stim_plus,non_stim), # Only specify columns to be transformed
    cols = c(perc_stim_spike,perc_fentanyl_spike,perc_heroin_spike,perc_rx_spike,perc_stim_non_spike,perc_fentanyl_non_spike,perc_heroin_non_spike,perc_rx_non_spike),
    names_to = "categories",
    values_to = "deaths"
  )
  
  
  
  df_long$categories <- factor(df_long$categories , levels = c("perc_stim_spike","perc_fentanyl_spike","perc_heroin_spike","perc_rx_spike","perc_stim_non_spike","perc_fentanyl_non_spike","perc_heroin_non_spike","perc_rx_non_spike"))
  
  
  ggplot(df_long, aes(x=year, y=deaths, color=categories)) +
    geom_line(linewidth=1.5) +
    labs(title="Overdose deaths and the percentage of substances involved for the state of Massachusetts",
         x="Year",
         y="Percentages",
         color="Category") +
    scale_x_continuous(breaks = df_long$year) +
    scale_color_brewer(palette = "Paired")+
    theme_minimal(base_size = 20)
}


plotSurface<-function(data)
{
  library(plotly)
  library(htmlwidgets)
  library(akima)
  
  data<-data[data$lag==30 & data$year==2022,]
  theData<-with(data, interp(x = threshold, y = influence, z = Count))
  
  plot<-plot_ly() %>%
    add_surface(x = theData$x, y = theData$y, z = theData$z) %>%
    layout(title = "Number of spikes and dips for 2022 with a lag of 30",
           scene = list(xaxis = list(title = 'Threshold'),
                        yaxis = list(title = 'Influence'),
                        zaxis = list(title = 'Total number of spikes and dips')))
  
  return(plot)
}

plotDaysSpikes<-function(data)
{
  #theData<-data[data$threshold==2 & data$lag==30 & data$year==2022,]
  theData<-data
  theData$spikes<-ifelse(theData$spikes>0,theData$spikes,NA)
  ggplot(theData, aes(x = date)) +
    geom_line(aes(y = y, color = 'Deaths')) +
    geom_point(aes(y = spikes, color = 'Spikes'), size = 2) +
    geom_line(aes(y = avgFilter, color = 'Moving Average (1 months)')) +
    theme_minimal() +
    scale_color_manual(name = 'Labels',
                       labels = c('Deaths', 'Moving Average (15 days)', 'Spikes'),
                       values = c("gray", "lightseagreen", "red")
    ) +
    labs(
      title = paste("Time Series Daily Overdose Death Data with Adjusted Moving Average of 15 days and Spikes,", year_of_interest),
      x = "Day",
      y = "Raw Counts of Overdose Deaths"
    ) + 
    ylim(0,6)
}