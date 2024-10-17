library(randtests)
library(lme4)
library(data.table)
library(dplyr)
library(emmeans)
library(ggplot2)
library(grid)
library(forcats)
library(cowplot)
library(tableone)
library(ggpubr)
library(sf)
library(tigris)
library(ggpattern)


county_analysis <- function(geographical_data)
{
  
  
  #lag<-c(7,15,30,60,90,120,150,180)
  #threshold <- c(1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3)
  lag<-c(7,30,60,90)
  threshold<-c(2,2.5)
  influence<-0.3
  
  theCounties<-count(geographical_data,"County_name"=geographical_data$County_name)
  theCounties<-na.omit(theCounties)
  theCounties<-theCounties$County_name
  #influence<-c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95)
  geographical_data<-na.omit(geographical_data)
  
  theResult<-0
  
  for(j in 1:length(lag))
    {
     
      for(k in 1:length(threshold))
      {
        for(l in 1:length(influence))
        {
          if(j==2)
          {
            print("error")
          }
          for(i in 1:length(theCounties))
          {
      
          y <- geographical_data[geographical_data$County_name==theCounties[i],]
       
          theDeaths<-y$Total_Overdose_Deaths
          result <- ThresholdingAlgo2(theDeaths,lag[j],threshold[k],influence[l])
          df <- data.frame(date = y$DOD_4_FD, 
                           y = theDeaths, 
                           avgFilter = result$avgFilter, 
                           upperThreshold = result$avgFilter + threshold * result$stdFilter, 
                           signals = result$signals,
                           spikes=result$spikes,
                           dips=result$dips,
                           YOD = y$YOD)
          df$year <- year(df$date)
          df_modified <- df
          df_modified$modified_y <- ifelse(df$signals == 0, 0, df$y)
          df_modified$y_minus_avgFilter_or_zero <- ifelse(df$signals == 0, 0, df$y - df$avgFilter)
          df_modified$County_name<-theCounties[i]
          df_modified$lag<-lag[j]
          df_modified$threshold<-threshold[k]
          df_modified$influence<-influence[l]
          df_modified$percSpike<-ifelse(df_modified$signals==1,df_modified$y_minus_avgFilter_or_zero/df_modified$avgFilter,0)
          #df_modified<-df_modified %>%
          #group_by(year,County_name,lag,threshold,influence) %>%
          #summarise(Count = sum(y_minus_avgFilter_or_zero > 0),Average_spike=mean(spikes[spikes > 0]),Average_dip=mean(dips[dips>0]),Total_spikes=sum(spikes),Total_deaths=sum(y),Average_per_inc_spike=mean(percSpike[percSpike>0]))
          #df_modified$perc_sike_death<-df_modified$Total_spikes/df_modified$Total_deaths

        
          if(i==1 & j==1 & k==1 & l==1)
          {
            theResult<-df_modified
          }
          else
          {
            
            theResult<-rbind(theResult,df_modified)
          }
          
        }
      }
    }
    print(j)
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


dtw<-function(x,y)
{
  theM<-matrix(0,length(x),length(y))
  theM[1,1]<-abs(x[1]-y[1])
  theB<-array(0,dim=c(length(x),length(y),2))
  
  for(i in 1:length(x))
  {
    for(j in 1:length(y))
    {
      if(i>1 & j>1)
      {
        theM[i,j]<-abs(x[i]-y[j])+min(theM[i-1,j],theM[i,j-1],theM[i-1,j-1])
        
        if(theM[i-1,j-1]<=min(theM[i-1,j],theM[i,j-1]))
        {
          theB[i,j,1]<-i-1
          theB[i,j,2]<-j-1
        }
        else if(theM[i-1,j]<=min(theM[i-1,j-1],theM[i,j-1]))
        {
          theB[i,j,1]<-i-1
          theB[i,j,2]<-j
        }
        else
        {
          theB[i,j,1]<-i
          theB[i,j,2]<-j-1
        }
        
      }
      else if(i>1 & j==1)
      {
        theM[i,j]<-abs(x[i]-y[j])+theM[i-1,j]
        theB[i,j,1]<-i-1
        theB[i,j,2]<-j
      }
      else if(i==1 & j>1)
      {
        theM[i,j]<-abs(x[i]-y[j])+theM[i,j-1]
        theB[i,j,1]<-i
        theB[i,j,2]<-j-1
      }
    }
  }
  
  return(theM[length(x),length(y)])
}

estimateDistance<-function(spikes,year)
{
  theCounties<-count(spikes,County_name)
  theCounties<-na.omit(theCounties)
  theCounties<-theCounties$County_name
  
  theDistances <- matrix(0, nrow = length(theCounties), ncol = length(theCounties))
  
  rownames(theDistances) <- theCounties
  colnames(theDistances) <- theCounties
  
  for(i in 1:length(theCounties))
  {
    for(j in i: length(theCounties))
    {
      county1<-spikes[spikes$YOD==year & spikes$County_name==theCounties[i],]
      county2<-spikes[spikes$YOD==year & spikes$County_name==theCounties[j],]
      
      theDistances[i,j]<-dtw(county1$y,county2$y)
      #theDistances[i,j]<-dtw(county1$spikes,county2$spikes)
      theDistances[j,i]<-theDistances[i,j]
    }
  }

  return(theDistances)  
}

estimateDistanceYear<-function(spikes)
{
  theCounties<-count(spikes,County_name)
  theCounties<-na.omit(theCounties)
  theCounties<-theCounties$County_name
  theYears<-c(2018,2019,2020,2021,2022)
  
  theDistances <- matrix(0, nrow = length(theCounties), ncol = length(theYears))
  
  rownames(theDistances) <- theCounties
  colnames(theDistances) <- theYears
  
  for(i in 1:length(theCounties))
  {
    for(j in 1: (length(theYears)))
    {
      county1<-spikes[spikes$YOD==theYears[j]-1 & spikes$County_name==theCounties[i],]
      county2<-spikes[spikes$YOD==theYears[j] & spikes$County_name==theCounties[i],]
      
      theDistances[i,j]<-dtw(county1$spikes,county2$spikes)
    }
  }
  
  return(theDistances)  
}


#Stimulants

addDrugs<-function(drug_categories,county_sens)
{
  #theDrugs<-drug_categories %>% group_by(DOD_4_FD,County_name)%>%summarise(ps_only=sum(ps_only),coc_only=sum(coc_only),stim_only=sum(stim_only),opioid_and_stim_only=sum(opioid_and_stim_only),stim_involved=sum(stim_involved),opioid_and_stim_involved=sum(opioid_and_stim_involved))
  theDrugs<-drug_categories %>% group_by(DOD_4_FD,County_name)%>%summarise(ps_only=sum(ps_only),stim_only=sum(stim_only),opioid_and_stim_only=sum(opioid_and_stim_only),stim_involved=sum(stim_involved),opioid_and_stim_involved=sum(opioid_and_stim_involved),fentanyl_involved=sum(fentanyl),heroin_involved=sum(heroin),rx_opioids_involved=sum(rx_opioids))
  
  theData<- county_sens %>% left_join(theDrugs,by=c("date"="DOD_4_FD","County_name"="County_name"))
  
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

calculateStatisticsDrugs<-function(county_sens_drugs)
{
  
  county_sens_drugs$stim_only<-county_sens_drugs$stim_only*ifelse(county_sens_drugs$signals==1,1,0)
  county_sens_drugs$stim_plus<-county_sens_drugs$stim_plus*ifelse(county_sens_drugs$signals==1,1,0)
  county_sens_drugs$non_stim<-county_sens_drugs$non_stim*ifelse(county_sens_drugs$signals==1,1,0)
  
  df_modified<-county_sens_drugs %>%
    group_by(year,County_name,lag,threshold,influence) %>%
    summarise(total_deaths = sum(spikes),stim_only=sum(stim_only),stim_plus=sum(stim_plus),non_stim=sum(non_stim))
  
  df_modified$per_stim<-ifelse(df_modified$total_deaths>0,(df_modified$stim_only+df_modified$stim_plus)/df_modified$total_deaths,0)
  df_modified$per_stim_only<-ifelse(df_modified$total_deaths>0,df_modified$stim_only/df_modified$total_deaths,0)
  df_modified$per_stim_plus<-ifelse(df_modified$total_deaths>0,df_modified$stim_plus/df_modified$total_deaths,0)
  df_modified$per_non_stim<-ifelse(df_modified$total_deaths>0,df_modified$non_stim/df_modified$total_deaths,0)
  
  
  
  return(df_modified)
}


#RANDOM TESTS

randRunTest<-function(spikes,lag,threshold,year)
{
  theData<-spikes[spikes$lag==lag & spikes$threshold==threshold & spikes$year>=year &spikes$year<=2022,]
  #theData<-spikes[spikes$lag==lag & spikes$threshold==threshold & spikes$year==year,]
  
  theCounties<-count(spikes,County_name)
  theCounties<-na.omit(theCounties)
  theCounties<-theCounties$County_name
  
  #theData$perc_deaths_spikes<-theData$deaths_spikes/theData$total_deaths

  thePValues<-data.frame("County_name"="","p_value"=0)
  
  for(i in 1:length(theCounties))
  {
    dataset<-theData[theData$County_name==theCounties[i],]
    dataset<-dataset$signals
    dataset<-ifelse(dataset>0,1,0)
    
    if(sum(dataset)>0)
    {
      result_runs <- runs.test(dataset,threshold = 0.5)
      thePValues[i,]<-c(theCounties[i],result_runs$p.value)
      
    }
    else
    {
      thePValues[i,]<-c(theCounties[i],-1)
    }
  }
  
  return(thePValues)
  
}

bartlesTest<-function(spikes,lag,threshold,year)
{
  theData<-spikes[spikes$lag==lag & spikes$threshold==threshold & spikes$year>=year &spikes$year<=2022,]
  #theData<-spikes[spikes$lag==lag & spikes$threshold==threshold & spikes$year==year,]
  
  theCounties<-count(spikes,County_name)
  theCounties<-na.omit(theCounties)
  theCounties<-theCounties$County_name
  
  #theData$perc_deaths_spikes<-theData$deaths_spikes/theData$total_deaths
  
  thePValues<-data.frame("County_name"="","p_value"=0)
  
  for(i in 1:length(theCounties))
  {
    dataset<-theData[theData$County_name==theCounties[i],]
    dataset<-dataset$signals
    dataset<-ifelse(dataset>0,1,0)
    
    if(sum(dataset)>0)
    {
      result_runs <- bartels.rank.test(dataset)
      thePValues[i,]<-c(theCounties[i],result_runs$p.value)
      
    }
    else
    {
      thePValues[i,]<-c(theCounties[i],-1)
    }
  }
  
  return(thePValues)
  
}

autoTest<-function(spikes,lag,threshold,year)
{
  theData<-spikes[spikes$lag==lag & spikes$threshold==threshold & spikes$year==year,]
  
  theCounties<-count(spikes,County_name)
  theCounties<-na.omit(theCounties)
  theCounties<-theCounties$County_name
  
  #theData$perc_deaths_spikes<-theData$deaths_spikes/theData$total_deaths
  
  thePValues<-data.frame("County_name"="","max_corr"=0)
  
  for(i in 1:length(theCounties))
  {
    dataset<-theData[theData$County_name==theCounties[i],]
    dataset<-dataset$signals
    dataset<-ifelse(dataset>0,1,0)
    
    if(sum(dataset)>0)
    {
      ts_data <- ts(dataset) 
      acf_result <- acf(ts_data, lag.max=30, plot = FALSE)
      thePValues[i,]<-c(theCounties[i],max(abs(acf_result$acf[-1])))
      
    }
    else
    {
      thePValues[i,]<-c(theCounties[i],-1)
    }
  }
  
  return(thePValues)
  
}

poissonMonthTest<-function(spikes,lag,threshold,year)
{
  theData<-spikes[spikes$lag==lag & spikes$threshold==threshold & spikes$year>=year &spikes$year<=2022,]
  
  theCounties<-count(spikes,County_name)
  theCounties<-na.omit(theCounties)
  theCounties<-theCounties$County_name
  
  theData$month<-month(theData$date)
  #theData$perc_deaths_spikes<-theData$deaths_spikes/theData$total_deaths
  
  thePValues<-data.frame("County_name"="","p_value"=0,"lambda"=0)
  
  for(i in 1:length(theCounties))
  {
    dataset<-theData[theData$County_name==theCounties[i],]
    dataset$signals<-ifelse(dataset$signals>0,1,0)
    
    dataset<-dataset%>%group_by(year,month)%>%summarise(numSpikes=sum(signals))
    
    if(sum(dataset$numSpikes)>0)
    {
      lambda_estimate <- mean(dataset$numSpikes)
      poisson_test <- sum((dataset$numSpikes - lambda_estimate)^2 / lambda_estimate)
      
      # Compare this statistic to a chi-squared distribution
      p_value <- 1-pchisq(poisson_test, df = length(dataset$numSpikes) - 1, lower.tail = FALSE)
      
      thePValues[i,]<-c(theCounties[i],p_value,lambda_estimate)
      
    }
    else
    {
      thePValues[i,]<-c(theCounties[i],-1,0)
    }
  }
  
  return(thePValues)
  
}


#T-Test for percentages

percStatTest<-function(spikes,lag,threshold,year)
{
  theData<-spikes[spikes$lag==lag & spikes$threshold==threshold & spikes$year==year,]
  
  theCounties<-count(spikes,County_name)
  theCounties<-na.omit(theCounties)
  theCounties<-theCounties$County_name
  
  theData$perc_stim<-theData$stim_involved/theData$y
  theData$perc_fentanyl<-theData$fentanyl_involved/theData$y
  theData$perc_heroin<-theData$heroin_involved/theData$y
  theData$perc_rx<-theData$rx_opioids_involved/theData$y
  
  theData$signals<-ifelse(theData$signals>0,1,0)
 

  #theRes<-data.frame("County_name"="","num_spike"=0,"num_non_spike"=0,"num_spike_g0"=0,"num_non_spike_g0"=0,"mean_stim_spike"=0,"mean_stim_non_spike"=0,"sd_stim_spike"=0,"sd_stim_non_spike"=0,"Ttest_stim"=0,"mean_fen_spike"=0,"mean_fen_non_spike"=0,"sd_fen_spike"=0,"sd_fen_non_spike"=0,"Ttest_fen"=0,"mean_her_spike"=0,"mean_her_non_spike"=0,"sd_her_spike"=0,"sd_her_non_spike"=0,"Ttest_her"=0)
  theRes<-data.frame("County_name"="","Drug"="","min"=0,"max"=0,"pvalue"=0)
  ind<-1
  for(i in 1:length(theCounties))
  {
    
    dataset<-theData[theData$County_name==theCounties[i],]
    
    dataSpike<-dataset[dataset$signals==1,]
    datanonSpike<-dataset[dataset$signals==0,]
    
    numDataSpike<-nrow(dataSpike)
    numDataNonSpike<-nrow(datanonSpike)
    
    dataset<-dataset[!is.nan(dataset$perc_stim),]
    
    dataSpike<-dataset[dataset$signals==1,]
    datanonSpike<-dataset[dataset$signals==0,]
    
    numDataSpikenon0<-nrow(dataSpike)
    numDataNonSpikenon0<-nrow(datanonSpike)
    
    
    
    if(nrow(dataset)>0)
    {

      if(nrow(dataSpike)>0 & nrow(datanonSpike>0))
      {
        theTestStim<-t.test(dataSpike$perc_stim, datanonSpike$perc_stim, var.equal = FALSE) 
        theTestFen<-t.test(dataSpike$perc_fen, datanonSpike$perc_fen, var.equal = FALSE) 
        theTestHer<-t.test(dataSpike$perc_her, datanonSpike$perc_her, var.equal = FALSE) 
        theTestRx<-t.test(dataSpike$perc_rx, datanonSpike$perc_rx, var.equal = FALSE) 
        #theRes[ind,]<-c(theCounties[i],"Stimulants",theTestStim$conf.int[1],theTestStim$conf.int[2],theTestStim$p.value)
        theRes[ind,]<-c(theCounties[i],"Stimulants",theTestStim$estimate[1],theTestStim$estimate[2],theTestStim$p.value)
        ind<-ind+1
        theRes[ind,]<-c(theCounties[i],"Fentanyl",theTestFen$estimate[1],theTestFen$estimate[2],theTestFen$p.value)
        ind<-ind+1
        theRes[ind,]<-c(theCounties[i],"Heroin",theTestHer$estimate[1],theTestHer$estimate[2],theTestHer$p.value)
        ind<-ind+1
        theRes[ind,]<-c(theCounties[i],"Rx_Opioids",theTestRx$estimate[1],theTestRx$estimate[2],theTestRx$p.value)
        ind<-ind+1
        #theRes[i,]<-c(theCounties[i],numDataSpike,numDataNonSpike,numDataSpikenon0,numDataNonSpikenon0,mean(dataSpike$perc_stim),mean(datanonSpike$perc_stim),sd(dataSpike$perc_stim),sd(datanonSpike$perc_stim),theTestStim$p.value,mean(dataSpike$perc_fen),mean(datanonSpike$perc_fen),sd(dataSpike$perc_fen),sd(datanonSpike$perc_fen),theTestFen$p.value,mean(dataSpike$perc_her),mean(datanonSpike$perc_her),sd(dataSpike$perc_her),sd(datanonSpike$perc_her),theTestHer$p.value)
      }
      else
      {
        #theRes[i,]<-c(theCounties[i],numDataSpike,numDataNonSpike,numDataSpikenon0,numDataNonSpikenon0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1)
        
      }
      
    }
    else
    {
      #theRes[i,]<-c(theCounties[i],numDataSpike,numDataNonSpike,numDataSpikenon0,numDataNonSpikenon0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1)
    }
  }
  return(theRes)
}

logTest<-function(drug_categories,county_sens,census,threshold,lag,year,county)
{
  theDrugs<-drug_categories %>%select(DOD_4_FD,YOD,County_name,AGE1_CALC,SEX,RACEGROUP,cocaine,psychostimulants,fentanyl,heroin,rx_opioids)
  theSens<-county_sens[county_sens$threshold==threshold & county_sens$lag==lag,]
  theSens<-theSens%>%select(date,y,County_name,signals)
  theSens$signals<-ifelse(theSens$signals>0,1,0)
  
  theData<- theDrugs %>% left_join(theSens,by=c("DOD_4_FD"="date","County_name"="County_name"))
  #theData<-theData[theData$County_name!="Barnstable" & theData$County_name!="Dukes" & theData$County_name!="Nantucket",] #Comment for State analysis
  theData<-theData[theData$County_name %in% c("Bristol","Essex","Hampden","Middlesex","Norfolk","Plymouth","Suffolk","Worcester"),] #Comment for State analysis
  #theData<-theData[theData$County_name %in% c("Middlesex","Suffolk","Worcester"),] #Comment for State analysis
  
  
  theData<-theData %>% left_join(census,by=c("County_name"="County_name"))
  theData$RACEGROUP<-ifelse(theData$RACEGROUP=="Black NH","Black",ifelse(theData$RACEGROUP=="White NH","White",ifelse(theData$RACEGROUP=="Hispanic","Hispanic","Other")))
  theData<-theData[theData$SEX=="Male" | theData$SEX=="Female",]
  #theData$Median.Age<-theData$Median.Age/100
  theData$AGE1_CALC<-theData$AGE1_CALC/100
  theData$Median.Age<-theData$Median.Age/100
  theData$popSize<-theData$popSize/1000000
  theData$Income<-theData$Income/100000
  theData$White<-theData$White/100
  if(county=="none")
  {
    if(year>0)
    {
      theData<-theData[theData$YOD==year,]
      model <- glmer(signals ~ cocaine+SEX+RACEGROUP+AGE1_CALC+(1+cocaine|County_name),data=theData,family = binomial(link = "logit"))
      
      #model <- glmer(signals ~ cocaine+(1+cocaine|County_name),data=theData,family = binomial(link = "logit"))
      
      #model <- glmer(signals ~ cocaine+(1+cocaine|County_name),data=theData,family = binomial(link = "logit"))
      
      #model <- glmer(signals ~ cocaine+psychostimulants+fentanyl+heroin+rx_opioids+Income+Median.Age+White+SEX+RACEGROUP+AGE1_CALC+(1+cocaine+psychostimulants+fentanyl+heroin+rx_opioids|County_name),data=theData,family = binomial(link = "logit"))
      
      #model <- glmer(signals ~ cocaine+fentanyl+heroin+rx_opioids+Male+Rurality+White+SEX+Median.Age+RACEGROUP+AGE1_CALC+(1+cocaine+fentanyl+heroin+rx_opioids|County_name),data=theData,family = binomial(link = "logit"))
  
    }
    else
    { 
      theData<-theData[theData$YOD>=2020 & theData$YOD<=2023,]
      
      #model <- glmer(signals ~ cocaine+psychostimulants+fentanyl+heroin+rx_opioids+popSize+Income+Median.Age+White+SEX+RACEGROUP+AGE1_CALC+YOD+(1+cocaine+psychostimulants+fentanyl+heroin+rx_opioids|County_name),data=theData,family = binomial(link = "logit"))
      model <- glmer(signals ~ rx_opioids+popSize+SEX+RACEGROUP+AGE1_CALC+(0+rx_opioids|County_name),data=theData,family = binomial(link = "logit"))
      
      
    }
      #model <- glmer(signals ~ cocaine+psychostimulants+fentanyl+heroin+rx_opioids+Male+Rurality+White+SEX+RACEGROUP+(1+cocaine+psychostimulants+fentanyl+heroin+rx_opioids|County_name),data=theData,family = binomial(link = "logit"))
    
    #model <- glmer(signals ~ cocaine+fentanyl+heroin+rx_opioids+Male+Rurality+White+SEX+Median.Age+RACEGROUP+(1+cocaine+fentanyl+heroin+rx_opioids|County_name),data=theData,family = binomial(link = "logit"))
    
    #model <- glmer(signals ~ cocaine+fentanyl+heroin+rx_opioids+Male+Rurality+White+SEX+RACEGROUP+(1+cocaine+fentanyl+heroin+rx_opioids|County_name),data=theData,family = binomial(link = "logit"))
    
    
  }
  else
  {
    theData<-theData[theData$County_name==county,]
    #model <- glmer(signals ~cocaine+psychostimulants+fentanyl+heroin+rx_opioids+(1+cocaine+psychostimulants+fentanyl+heroin+rx_opioids|SEX)+(1+cocaine+psychostimulants+fentanyl+heroin+rx_opioids|RACEGROUP),data=theData,family = binomial(link = "logit"), control = glmerControl(optimizer = "nloptwrap"))
    model <- glm(signals ~cocaine+psychostimulants+fentanyl+heroin+rx_opioids,data=theData, family = "binomial")
    
  }
  return(model)
  
}

logTestState<-function(drug_categories,state_sens,census,threshold,lag,year)
{
  theDrugs<-drug_categories %>%select(DOD_4_FD,YOD,AGE1_CALC,SEX,RACEGROUP,County_name,cocaine,psychostimulants,fentanyl,heroin,rx_opioids)
  theSens<-state_sens[state_sens$threshold==threshold & state_sens$lag==lag,]
  theSens<-theSens%>%select(date,signals)
  theSens$signals<-ifelse(theSens$signals>0,1,0)
  
  theData<- theDrugs %>% left_join(theSens,by=c("DOD_4_FD"="date"))
  theData$RACEGROUP<-ifelse(theData$RACEGROUP=="Black NH","Black",ifelse(theData$RACEGROUP=="White NH","White",ifelse(theData$RACEGROUP=="Hispanic","Hispanic","Other")))
  theData$SEX<-as.factor(theData$SEX)
  theData$RACEGROUP<-as.factor(theData$RACEGROUP)
  theData<-theData %>% left_join(census,by=c("County_name"="County_name"))
  theData$AGE1_CALC<-theData$AGE1_CALC/100
  theData$Median.Age<-theData$Median.Age/100
  theData<-theData[theData$YOD>=2020 & theData$YOD<=2023,]
  theData<-theData[theData$SEX=="Male" | theData$SEX=="Female",]
  theData$popSize<-theData$popSize/1000000
  theData$Income<-theData$Income/100000
  theData$RACEGROUP <- relevel(theData$RACEGROUP, ref = "White")
  
  if(year>0)
  {
    theData<-theData[theData$YOD==year,]
    model <- glm(signals ~cocaine+psychostimulants+fentanyl+heroin+rx_opioids+popSize+Income+Median.Age+White+SEX+RACEGROUP+AGE1_CALC,data=theData, family = "binomial")
  }
  else
    {
      model <- glm(signals ~cocaine+psychostimulants+fentanyl+heroin+rx_opioids+popSize+Income+Median.Age+White+SEX+RACEGROUP+AGE1_CALC+YOD,data=theData, family = "binomial")
      #model <- glmer(signals ~ cocaine+psychostimulants+fentanyl+heroin+rx_opioids+popSize+Income+Median.Age+White+SEX+RACEGROUP+AGE1_CALC+YOD+(1|County_name),data=theData,family = binomial(link = "logit"))
      
  }
  return(model)
  
}

estimateCI<-function(model,year)
{
  #cDrugs<-c("cocaine","psychostimulants","fentanyl","heroin","rx_opioids","AGE1_CALC","SEXM","RACEGROUPHispanic","RACEGROUPOther","RACEGROUPWhite")
  cDrugs<-c("cocaine","psychostimulants","fentanyl","heroin","rx_opioids")
  #cDrugs<-c("cocaine","fentanyl","heroin","rx_opioids")
  
  fixed_effects<-fixef(model)
  random_effects<-ranef(model)$County_name
  theRes<-0
  for(i in 1:length(cDrugs))
  {
    theRS<-random_effects[,i+1]
    adjustedRS<-theRS+fixed_effects[i+1]

    se_random_slopes <- attr(ranef(model, condVar = TRUE)$County_name, "postVar")[i+1, i+1, ]
    se_fixed_effects<-sqrt(vcov(model)[i+1, i+1])
    the_se<-sqrt(se_fixed_effects^2 + se_random_slopes^2)
    
    ci_lower <- adjustedRS - 1.96 *the_se
    ci_upper <- adjustedRS + 1.96 *the_se
    random_slopes_df <- data.frame(group = rownames(random_effects),
                                   random_slope = adjustedRS,
                                   ci_lower = ci_lower,
                                   ci_upper = ci_upper,
                                   drug=cDrugs[i])
    if(i==1)
    {
      theRes<-random_slopes_df
    }
    else
    {
      theRes<-rbind(theRes,random_slopes_df)
    }

  }
  # 
  # random_effects<-ranef(model)$'County_name:SEX'
  # for(i in 1:length(cDrugs))
  # {
  #   theRS<-random_effects[,i+1]
  #   adjustedRS<-theRS+fixed_effects[i+1]
  #   adjustedRS<-theRS
  #   
  #   se_random_slopes <- attr(ranef(model, condVar = TRUE)$'County_name:SEX', "postVar")[i+1, i+1, ]
  #   ci_lower <- adjustedRS - 1.96 * sqrt(se_random_slopes)
  #   ci_upper <- adjustedRS + 1.96 * sqrt(se_random_slopes)
  #   random_slopes_df <- data.frame(group = rownames(random_effects),
  #                                  random_slope = theRS,
  #                                  ci_lower = ci_lower,
  #                                  ci_upper = ci_upper,
  #                                  drug=cDrugs[i])
  # 
  #     theRes<-rbind(theRes,random_slopes_df)
  # 
  # }
  # 
  # random_effects<-ranef(model)$'County_name:RACEGROUP'
  # for(i in 1:length(cDrugs))
  # {
  #   theRS<-random_effects[,i+1]
  #   adjustedRS<-theRS+fixed_effects[i+1]
  #   adjustedRS<-theRS
  #   
  #   se_random_slopes <- attr(ranef(model, condVar = TRUE)$'County_name:RACEGROUP', "postVar")[i+1, i+1, ]
  #   ci_lower <- adjustedRS - 1.96 * sqrt(se_random_slopes)
  #   ci_upper <- adjustedRS + 1.96 * sqrt(se_random_slopes)
  #   random_slopes_df <- data.frame(group = rownames(random_effects),
  #                                  random_slope = theRS,
  #                                  ci_lower = ci_lower,
  #                                  ci_upper = ci_upper,
  #                                  drug=cDrugs[i])
  #   
  #   theRes<-rbind(theRes,random_slopes_df)
  #   
  # }
  
  theRes$year<-year
  return(theRes)
  
}

estimateCI_County<-function(model)
{
  #cDrugs<-c("cocaine","psychostimulants","fentanyl","heroin","rx_opioids","AGE1_CALC","SEXM","RACEGROUPHispanic","RACEGROUPOther","RACEGROUPWhite")
  cDrugs<-c("cocaine","psychostimulants","fentanyl","heroin","rx_opioids")
  
  fixed_effects<-fixef(model)
  #fixed_effects<-fixef(model)-fixef(model)
  random_effects<-ranef(model)$SEX
  theRes<-0
  for(i in 1:length(cDrugs))
  {
    theRS<-random_effects[,i+1]
    adjustedRS<-theRS+fixed_effects[i+1]
    #adjustedRS<-theRS
    
    se_random_slopes <- attr(ranef(model, condVar = TRUE)$SEX, "postVar")[i+1, i+1, ]
    ci_lower <- adjustedRS - 1.81 * sqrt(se_random_slopes)
    ci_upper <- adjustedRS + 1.81 * sqrt(se_random_slopes)
    random_slopes_df <- data.frame(group = rownames(random_effects),
                                   random_slope = theRS,
                                   ci_lower = ci_lower,
                                   ci_upper = ci_upper,
                                   drug=cDrugs[i])
    if(i==1)
    {
      theRes<-random_slopes_df
    }
    else
    {
      theRes<-rbind(theRes,random_slopes_df)
    }
    
  }
  
  random_effects<-ranef(model)$'RACEGROUP'
  for(i in 1:length(cDrugs))
  {
    theRS<-random_effects[,i+1]
    adjustedRS<-theRS+fixed_effects[i+1]
    #adjustedRS<-theRS
    se_random_slopes <- attr(ranef(model, condVar = TRUE)$'RACEGROUP', "postVar")[i+1, i+1, ]
    ci_lower <- adjustedRS - 1.96 * sqrt(se_random_slopes)
    ci_upper <- adjustedRS + 1.96 * sqrt(se_random_slopes)
    random_slopes_df <- data.frame(group = rownames(random_effects),
                                   random_slope = theRS,
                                   ci_lower = ci_lower,
                                   ci_upper = ci_upper,
                                   drug=cDrugs[i])
    
    theRes<-rbind(theRes,random_slopes_df)
    
  }
  return(theRes)
  
}

estimateCI_State<-function(model,year)
{

  cDrugs<-c("cocaine","psychostimulants","fentanyl","heroin","rx_opioids")
  
  theConfInt<-confint(model)
  theConfInt<-data.frame(theConfInt)
  theConfInt$drugs <- rownames(theConfInt)
  theConfInt$estimates<-model$coefficients
  
  theRes<-theConfInt[theConfInt$drugs %in% cDrugs,]
  names(theRes)[1]<-"ci_lower"
  names(theRes)[2]<-"ci_upper"
  theRes$year<-year
  return(theRes)
  
}

estimateCIOneDrug<-function(model,year,drug)
{
  #cDrugs<-c("cocaine","psychostimulants","fentanyl","heroin","rx_opioids","AGE1_CALC","SEXM","RACEGROUPHispanic","RACEGROUPOther","RACEGROUPWhite")

  #cDrugs<-c("cocaine","fentanyl","heroin","rx_opioids")
  
  fixed_effects<-fixef(model)
  random_effects<-ranef(model)$County_name
  theRes<-0
  i<-1
  theRS<-random_effects[,i+1]
  adjustedRS<-theRS+fixed_effects[i+1]
    
  se_random_slopes <- attr(ranef(model, condVar = TRUE)$County_name, "postVar")[i+1, i+1, ]
  se_fixed_effects<-sqrt(vcov(model)[i+1, i+1])
  the_se<-sqrt(se_fixed_effects^2 + se_random_slopes^2)
    
  ci_lower <- adjustedRS - 1.96 *the_se
  ci_upper <- adjustedRS + 1.96 *the_se
  random_slopes_df <- data.frame(group = rownames(random_effects),
                                   random_slope = adjustedRS,
                                   ci_lower = ci_lower,
                                   ci_upper = ci_upper,
                                   drug=drug)
  theRes<-random_slopes_df
  theRes$year<-year
  return(theRes)
  
}

estimateCIOneDrugNoInt<-function(model,year,drug)
{
  #cDrugs<-c("cocaine","psychostimulants","fentanyl","heroin","rx_opioids","AGE1_CALC","SEXM","RACEGROUPHispanic","RACEGROUPOther","RACEGROUPWhite")
  
  #cDrugs<-c("cocaine","fentanyl","heroin","rx_opioids")
  
  fixed_effects<-fixef(model)
  random_effects<-ranef(model)$County_name
  theRes<-0
  i<-0
  theRS<-random_effects[,i+1]
  adjustedRS<-theRS+fixed_effects[i+1]
  
  se_random_slopes <- attr(ranef(model, condVar = TRUE)$County_name, "postVar")[i+1, i+1, ]
  se_fixed_effects<-sqrt(vcov(model)[i+1, i+1])
  the_se<-sqrt(se_fixed_effects^2 + se_random_slopes^2)
  
  ci_lower <- adjustedRS - 1.96 *the_se
  ci_upper <- adjustedRS + 1.96 *the_se
  random_slopes_df <- data.frame(group = rownames(random_effects),
                                 random_slope = adjustedRS,
                                 ci_lower = ci_lower,
                                 ci_upper = ci_upper,
                                 drug=drug)
  theRes<-random_slopes_df
  theRes$year<-year
  return(theRes)
  
}


logCI<-function(drug_categories,county_sens,census,threshold,lag,year)
{
  theDrugs<-drug_categories %>%select(DOD_4_FD,YOD,County_name,AGE1_CALC,SEX,RACEGROUP,cocaine,psychostimulants,fentanyl,heroin,rx_opioids)
  theSens<-county_sens[county_sens$threshold==threshold & county_sens$lag==lag,]
  theSens<-theSens%>%select(date,y,County_name,signals)
  theSens$signals<-ifelse(theSens$signals>0,1,0)
  
  theData<- theDrugs %>% left_join(theSens,by=c("DOD_4_FD"="date","County_name"="County_name"))
  theData<-theData[theData$County_name %in% c("Bristol","Essex","Hampden","Middlesex","Norfolk","Plymouth","Suffolk","Worcester"),] #Comment for State analysis
 
  
  theData<-theData %>% left_join(census,by=c("County_name"="County_name"))
  theData$RACEGROUP<-ifelse(theData$RACEGROUP=="Black NH","Black",ifelse(theData$RACEGROUP=="White NH","White",ifelse(theData$RACEGROUP=="Hispanic","Hispanic","Other")))
  theData<-theData[theData$SEX=="Male" | theData$SEX=="Female",]
  theData$AGE1_CALC<-theData$AGE1_CALC/100
  theData$Median.Age<-theData$Median.Age/100
  theData$popSize<-theData$popSize/1000000
  theData$Income<-theData$Income/100000
  theData$White<-theData$White/100
  theData$RACEGROUP<-as.factor(theData$RACEGROUP)
  theData$RACEGROUP <- relevel(theData$RACEGROUP, ref = "White")
  
  if(year>0)
  {
    theData<-theData[theData$YOD==year,]
    modelC <- glmer(signals ~ cocaine+Income+SEX+RACEGROUP+AGE1_CALC+(1+cocaine|County_name),data=theData,family = binomial(link = "logit"))
    ciC<-estimateCIOneDrug(modelC,year,"cocaine")
    
    modelF <- glmer(signals ~ fentanyl+popSize+Income+Median.Age+White+SEX+RACEGROUP+AGE1_CALC+(0+fentanyl|County_name),data=theData,family = binomial(link = "logit"))
    ciF<-estimateCIOneDrugNoInt(modelF,year,"fentanyl")
    
    modelH <- glmer(signals ~ heroin+Income+Median.Age+White+SEX+RACEGROUP+AGE1_CALC+(1+heroin|County_name),data=theData,family = binomial(link = "logit"))
    ciH<-estimateCIOneDrug(modelH,year,"heroin")
    
    modelP <- glmer(signals ~ psychostimulants+Income+White+SEX+RACEGROUP+AGE1_CALC+(0+psychostimulants|County_name),data=theData,family = binomial(link = "logit"))
    ciP<-estimateCIOneDrugNoInt(modelP,year,"psychostimulants")
    
    modelR <- glmer(signals ~ rx_opioids+popSize+SEX+RACEGROUP+AGE1_CALC+(0+rx_opioids|County_name),data=theData,family = binomial(link = "logit"))
    
    ciR<-estimateCIOneDrugNoInt(modelR,year,"rx_opioids")

    theCIs<-rbind(ciC,ciF,ciH,ciP,ciR)
  }
  else
  { 
    theData<-theData[theData$YOD>=2020 & theData$YOD<=2023,]
      
    modelC <- glmer(signals ~ cocaine+Income+SEX+RACEGROUP+AGE1_CALC+(1+cocaine|County_name),data=theData,family = binomial(link = "logit"))
    ciC<-estimateCIOneDrug(modelC,year,"cocaine")
    
    modelF <- glmer(signals ~ fentanyl+popSize+Income+Median.Age+White+SEX+RACEGROUP+AGE1_CALC+(1+fentanyl|County_name),data=theData,family = binomial(link = "logit"))
    ciF<-estimateCIOneDrug(modelF,year,"fentanyl")
    
    modelH <- glmer(signals ~ heroin+Income+Median.Age+White+SEX+RACEGROUP+AGE1_CALC+(1+heroin|County_name),data=theData,family = binomial(link = "logit"))
    ciH<-estimateCIOneDrug(modelH,year,"heroin")
    
    modelP <- glmer(signals ~ psychostimulants+Income+White+SEX+RACEGROUP+AGE1_CALC+(1+psychostimulants|County_name),data=theData,family = binomial(link = "logit"))
    ciP<-estimateCIOneDrug(modelP,year,"psychostimulants")
    
    modelR <- glmer(signals ~ rx_opioids+SEX+RACEGROUP+AGE1_CALC+(0+rx_opioids|County_name),data=theData,family = binomial(link = "logit"))
    
    ciR<-estimateCIOneDrugNoInt(modelR,year,"prescription opioids")
    
    theCIs<-rbind(ciC,ciF,ciH,ciP,ciR)
  }

  return(theCIs)
}

logDrug<-function(drug_categories,county_sens,census,threshold,lag,drug)
{
  theDrugs<-drug_categories %>%select(DOD_4_FD,YOD,County_name,AGE1_CALC,SEX,RACEGROUP,cocaine,psychostimulants,fentanyl,heroin,rx_opioids)
  theSens<-county_sens[county_sens$threshold==threshold & county_sens$lag==lag,]
  theSens<-theSens%>%select(date,y,County_name,signals)
  theSens$signals<-ifelse(theSens$signals>0,1,0)
  
  theData<- theDrugs %>% left_join(theSens,by=c("DOD_4_FD"="date","County_name"="County_name"))
  theData<-theData[theData$County_name %in% c("Bristol","Essex","Hampden","Middlesex","Norfolk","Plymouth","Suffolk","Worcester"),] #Comment for State analysis
  
  
  theData<-theData %>% left_join(census,by=c("County_name"="County_name"))
  theData$RACEGROUP<-ifelse(theData$RACEGROUP=="Black NH","Black",ifelse(theData$RACEGROUP=="White NH","White",ifelse(theData$RACEGROUP=="Hispanic","Hispanic","Other")))
  theData$RACEGROUP<-as.factor(theData$RACEGROUP)
  theData$RACEGROUP <- relevel(theData$RACEGROUP, ref = "White")
  theData<-theData[theData$SEX=="Male" | theData$SEX=="Female",]
  theData$AGE1_CALC<-theData$AGE1_CALC/100
  theData$Median.Age<-theData$Median.Age/100
  theData$popSize<-theData$popSize/1000000
  theData$Income<-theData$Income/100000
  theData$White<-theData$White/100

 
  theData<-theData[theData$YOD>=2020 & theData$YOD<=2023,]
  
  if(drug=="cocaine")
  {
    model <- glmer(signals ~ cocaine+popSize+Income+Median.Age+White+SEX+RACEGROUP+AGE1_CALC+(1+cocaine|County_name),data=theData,family = binomial(link = "logit"))
   
  }
  else if(drug=="fentanyl")
  {
    model <- glmer(signals ~ fentanyl+popSize+Income+Median.Age+White+SEX+RACEGROUP+AGE1_CALC+(1+fentanyl|County_name),data=theData,family = binomial(link = "logit"))
  }
  else if(drug=="heroin")
  {
    model <- glmer(signals ~ heroin+popSize+Income+Median.Age+White+SEX+RACEGROUP+AGE1_CALC+(1+heroin|County_name),data=theData,family = binomial(link = "logit"))
  }
  else if(drug=="psychostimulants")
  {
    model <- glmer(signals ~ psychostimulants+popSize+Income+Median.Age+White+SEX+RACEGROUP+AGE1_CALC+(1+psychostimulants|County_name),data=theData,family = binomial(link = "logit"))
  }
  else
  {
    #model <- glmer(signals ~ rx_opioids+popSize+Income+Median.Age+White+SEX+RACEGROUP+AGE1_CALC+(1+rx_opioids|County_name),data=theData,family = binomial(link = "logit"))
    model <- glmer(signals ~ rx_opioids+(0+rx_opioids|County_name),data=theData,family = binomial(link = "logit"))
    
  }
  return(model)
}
#PLOTS#
graphCounties<-function(spikes,year)
{
  theCounties<-count(spikes,County_name)
  theCounties<-na.omit(theCounties)
  theCounties<-theCounties$County_name
  
  theData<-spikes[spikes$year==year & spikes$threshold==2 & spikes$lag==30,]
  
  for (i in 1:length(theCounties)) {
    variable_name <- paste0("plot", i)
    thedf<-theData[theData$County_name==theCounties[i],]
    long_df <- thedf %>%
      gather(key = "type", value = "value", -date,-y, -signals, -year, -y_minus_avgFilter_or_zero,-spikes)
    plot <- ggplot(long_df, aes(x = date, y = spikes)) +
      geom_step(color = "red") +
      scale_x_date(date_breaks = "1 month", date_labels = "%b") +
      scale_y_continuous(limits = c(-10, 10)) +
      labs(x = "", y = "") +
      theme_minimal() +
      ggtitle(theCounties[i])
    assign(variable_name,plot)
  }
  #combined_plot<-plot1 / plot2/ plot3/ plot4/ plot5 / plot6/ plot7/ plot8/ plot9 /plot10/ plot11/ plot12/ plot13/ plot14  
  combined_plot<-plot11 / plot12
  
  return(combined_plot)
}

graphCountiesStim<-function(spikes,lag,threshold,year)
{
  theData<-spikes[spikes$lag==lag & spikes$threshold==threshold & spikes$year==year,]
  
  theCounties<-count(spikes,County_name)
  theCounties<-na.omit(theCounties)
  theCounties<-theCounties$County_name
  
  upperlim<-max(theData$y)+1
  for (i in 1:length(theCounties)) 
  {
    variable_name <- paste0("plot", i)
    thedf<-theData[theData$County_name==theCounties[i],]
    
    df_long <- thedf %>% 
      pivot_longer(
        cols = c(stim_involved,non_stim), # Only specify columns to be transformed
        names_to = "drugs",
        values_to = "deaths"
      )
    
    df_long$drugs <- factor(df_long$drugs , levels = c("non_stim", "stim_involved"))
    
    df_long$deaths<-df_long$deaths*ifelse(df_long$signals==1,1,0)
    
    plot<-ggplot(data = df_long, aes(x = date, y = deaths, fill = drugs)) +
      geom_bar(stat = "identity") +
      labs(x = "Date", y = "Number of deaths in spikes", title = theCounties[i]) +
      scale_x_date(date_breaks = "1 month", date_labels = "%b") +
      scale_y_continuous(limits = c(0,upperlim )) +
      scale_fill_manual(values=c("non_stim"="red","stim_involved"="blue")) +
      labs(x = "", y = "") +
      theme_minimal()
    
    assign(variable_name,plot)
    
  }
  
  combined_plot<-plot1 / plot2/ plot3/ plot4/ plot5 / plot6/ plot7/ plot8/ plot9 /plot10/ plot11/ plot12/ plot13/ plot14  
  return(combined_plot)
}

graphCountiesFen<-function(spikes,lag,threshold,year)
{
  theData<-spikes[spikes$lag==lag & spikes$threshold==threshold & spikes$year==year,]
  
  theCounties<-count(spikes,County_name)
  theCounties<-na.omit(theCounties)
  theCounties<-theCounties$County_name
  
  upperlim<-max(theData$y)+1
  for (i in 1:length(theCounties)) 
  {
    variable_name <- paste0("plot", i)
    thedf<-theData[theData$County_name==theCounties[i],]
    
    df_long <- thedf %>% 
      pivot_longer(
        cols = c(fentanyl_involved,non_fen), # Only specify columns to be transformed
        names_to = "drugs",
        values_to = "deaths"
      )
    
    df_long$drugs <- factor(df_long$drugs , levels = c("non_fen", "fentanyl_involved"))
    
    df_long$deaths<-df_long$deaths*ifelse(df_long$signals==1,1,0)
    
    plot<-ggplot(data = df_long, aes(x = date, y = deaths, fill = drugs)) +
      geom_bar(stat = "identity") +
      labs(x = "Date", y = "Number of deaths in spikes", title = theCounties[i]) +
      scale_x_date(date_breaks = "1 month", date_labels = "%b") +
      scale_y_continuous(limits = c(0,upperlim )) +
      scale_fill_manual(values=c("non_fen"="red","fentanyl_involved"="blue")) +
      labs(x = "", y = "") +
      theme_minimal()
    
    assign(variable_name,plot)
    
  }
  
  combined_plot<-plot1 / plot2/ plot3/ plot4/ plot5 / plot6/ plot7/ plot8/ plot9 /plot10/ plot11/ plot12/ plot13/ plot14  
  return(combined_plot)
}


graphCountiesDrugs_stats_years<-function(spikes,lag,threshold)
{
  theData<-spikes[spikes$lag==lag & spikes$threshold==threshold & spikes$year>=2017 & spikes$year<=2022,]
  
  theCounties<-theData%>%group_by(County_name)%>%summarize(lag=sum(lag))
  theCounties<-na.omit(theCounties)
  theCounties<-theCounties$County_name
  
  #upperlimit<-max(theData$total_deaths)
  upperlimit<-1
  for (i in 1:length(theCounties)) 
  {
    variable_name <- paste0("plot", i)
    thedf<-theData[theData$County_name==theCounties[i],]
    
    df_long <- thedf %>% 
      pivot_longer(
        cols = c(per_stim_only,per_stim_plus,per_non_stim), # Only specify columns to be transformed
        #cols = c(stim_only,stim_plus,non_stim), # Only specify columns to be transformed
        names_to = "drugs",
        values_to = "percentages"
      )
    
    df_long$drugs <- factor(df_long$drugs , levels = c("per_non_stim", "per_stim_plus", "per_stim_only"))
    
    plot<-ggplot(data = df_long, aes(x = year, y = percentages, fill = drugs)) +
      geom_bar(stat = "identity") +
      labs(x = "Date", y = "Percentages of deaths in spikes", title = theCounties[i]) +
      scale_y_continuous(limits = c(0,upperlimit)) +
      scale_fill_manual(values=c("per_non_stim"="red",  "per_stim_plus"="violet","per_stim_only"="blue")) +
      labs(x = "", y = "") +
      theme_minimal()
    
    assign(variable_name,plot)
    
  }
  
  combined_plot<-plot1 / plot2/ plot3/ plot4/ plot5 / plot6/ plot7/ plot8/ plot9 /plot10/ plot11/ plot12/ plot13/ plot14  
  return(combined_plot)
}

graphTotalDrugs<-function(spikes,lag,threshold)
{
  theData<-spikes[spikes$lag==lag & spikes$threshold==threshold & spikes$year>=2012 & spikes$year<=2022,]
  theData$stim_spike<-ifelse(theData$signals==1,theData$stim_involved,0)
  theData$fentanyl_spike<-ifelse(theData$signals==1,theData$fentanyl_involved,0)
  
  theData<-theData%>%group_by(year)%>%summarise(total_deaths=sum(y),deaths_spikes=sum(spikes),total_stim_involved=sum(stim_involved),total_fentanyl_involved=sum(fentanyl_involved),stim_spikes=sum(stim_spike),fentanyl_spikes=sum(fentanyl_spike))
  
  df_long <- theData %>% pivot_longer(
      cols = c(total_deaths,deaths_spikes,total_stim_involved,total_fentanyl_involved,stim_spikes,fentanyl_spikes), # Only specify columns to be transformed
      #cols = c(stim_only,stim_plus,non_stim), # Only specify columns to be transformed
      names_to = "categories",
      values_to = "deaths"
    )
  
  
  
  df_long$categories <- factor(df_long$categories , levels = c("total_deaths","total_stim_involved","total_fentanyl_involved","deaths_spikes","stim_spikes","fentanyl_spikes"))
  
  
  ggplot(df_long, aes(x=year, y=deaths, color=categories)) +
    geom_line(linewidth=1.5) +
    labs(title="Comparison of total deaths and deaths in spikes including stimulants and fentanyl",
         x="Year",
         y="Deaths",
         color="Category") +
    scale_x_continuous(breaks = df_long$year) +
    scale_color_brewer(palette = "Paired")+
    theme_minimal()
}

graphPercDrugs<-function(spikes,lag,threshold,County)
{
  theData<-spikes[spikes$lag==lag & spikes$threshold==threshold & spikes$year>=2012 & spikes$year<=2022,]
  theData<-theData[theData$County_name==County,]
  theData$stim_spike<-ifelse(theData$signals==1,theData$stim_involved,0)
  theData$fentanyl_spike<-ifelse(theData$signals==1,theData$fentanyl_involved,0)

  theData$stim_non_spike<-ifelse(theData$signals==1,0,theData$stim_involved)
  theData$fentanyl_non_spike<-ifelse(theData$signals==1,0,theData$fentanyl_involved)
  
  theData$non_spike_deaths<-theData$y-theData$spikes
    
  theData<-theData%>%group_by(year)%>%summarise(total_deaths=sum(y),deaths_spikes=sum(spikes),deaths_non_spikes=sum(non_spike_deaths),total_stim_involved=sum(stim_involved),total_fentanyl_involved=sum(fentanyl_involved),stim_spikes=sum(stim_spike),fentanyl_spikes=sum(fentanyl_spike),stim_non_spikes=sum(stim_non_spike),fentanyl_non_spikes=sum(fentanyl_non_spike))
  
  #theData$perc_deaths_spikes<-theData$deaths_spikes/theData$total_deaths
  
  theData$perc_stim<-theData$total_stim_involved/theData$total_deaths
  theData$perc_fentanyl<-theData$total_fentanyl_involved/theData$total_deaths
  
  theData$perc_stim_spike<-theData$stim_spikes/theData$deaths_spikes
  theData$perc_fentanyl_spike<-theData$fentanyl_spikes/theData$deaths_spikes
  
  theData$perc_stim_non_spike<-theData$stim_non_spikes/theData$deaths_non_spikes
  theData$perc_fentanyl_non_spike<-theData$fentanyl_non_spikes/theData$deaths_non_spikes
  
  
  df_long <- theData %>% pivot_longer(
    #cols = c(perc_stim,perc_fentanyl,perc_stim_spike,perc_fentanyl_spike), # Only specify columns to be transformed
    #cols = c(stim_only,stim_plus,non_stim), # Only specify columns to be transformed
    cols = c(perc_stim_spike,perc_fentanyl_spike,perc_stim_non_spike,perc_fentanyl_non_spike),
    names_to = "categories",
    values_to = "deaths"
  )
  
  
  
  df_long$categories <- factor(df_long$categories , levels = c("perc_stim_spike","perc_fentanyl_spike","perc_stim_non_spike","perc_fentanyl_non_spike"))
  
  
  ggplot(df_long, aes(x=year, y=deaths, color=categories)) +
    geom_line(linewidth=1.5) +
    labs(title=County,
         x="Year",
         y="Percentages",
         color="Category") +
    scale_x_continuous(breaks = df_long$year) +
    scale_color_brewer(palette = "Paired")+
    theme_minimal()
}

graphDendrogram<-function(distances,year)
{
  dist_matrix<-as.dist(distances)
  plotClusters<-hclust(dist_matrix,method="complete")
  plot(plotClusters,xlab="Counties",main=paste("Cluster Dendrogram",year))
}
graphHeatMap<-function(distances,year)
{
  row_clusters <- hclust(dist(distances))
  col_clusters <- hclust(dist(t(distances)))
  m<-melt(distances[row_clusters$order,col_clusters$order])
  ggplot(m, aes(Var1, Var2, fill=value)) + geom_tile()+theme_minimal()+ggtitle(paste("Heat Map ",year))+xlab("Counties")+ylab("Counties")
}

graphBarCountyState<-function(drug_categories,county_sens,state_sens,lag,threshold,year,county)
{
  theDrugs<-drug_categories %>%select(DOD_4_FD,YOD,County_name,AGE1_CALC,SEX,RACEGROUP,cocaine,psychostimulants,fentanyl,heroin,rx_opioids)
  theSens<-county_sens[county_sens$threshold==threshold & county_sens$lag==lag,]
  theSens<-theSens%>%select(date,County_name,signals)
  theSens$signals<-ifelse(theSens$signals>0,1,0)
  
  theData<- theDrugs %>% left_join(theSens,by=c("DOD_4_FD"="date","County_name"="County_name"))
  theData<-theData[theData$YOD==year &theData$County_name==county,]
  
  theCounty<-theData%>%group_by(signals)%>%summarise(Deaths=n(),cocaine=sum(cocaine),psychostimulants=sum(psychostimulants),fentanyl=sum(fentanyl),heroin=sum(heroin),rx_opioids=sum(rx_opioids))
  
  theSensState<-state_sens[state_sens$threshold==threshold & state_sens$lag==lag,]
  theSensState<-theSensState%>%select(date,signals)
  theSensState$signals<-ifelse(theSensState$signals>0,1,0)
  
  theDataState<- theDrugs %>% left_join(theSensState,by=c("DOD_4_FD"="date"))
  theDataState<-theDataState[theDataState$YOD==year,]
  
  theState<-theDataState%>%group_by(signals)%>%summarise(Deaths=n(),cocaine=sum(cocaine),psychostimulants=sum(psychostimulants),fentanyl=sum(fentanyl),heroin=sum(heroin),rx_opioids=sum(rx_opioids))
  
  theCounty$cocaine<-theCounty$cocaine/theCounty$Deaths
  theCounty$psychostimulants<-theCounty$psychostimulants/theCounty$Deaths
  theCounty$fentanyl<-theCounty$fentanyl/theCounty$Deaths
  theCounty$heroin<-theCounty$heroin/theCounty$Deaths
  theCounty$rx_opioids<-theCounty$rx_opioids/theCounty$Deaths
  
  theState$cocaine<-theState$cocaine/theState$Deaths
  theState$psychostimulants<-theState$psychostimulants/theState$Deaths
  theState$fentanyl<-theState$fentanyl/theState$Deaths
  theState$heroin<-theState$heroin/theState$Deaths
  theState$rx_opioids<-theState$rx_opioids/theState$Deaths
  
  theCounty<-theCounty%>%select(signals,cocaine,psychostimulants,fentanyl,heroin,rx_opioids)
  theState<-theState%>%select(signals,cocaine,psychostimulants,fentanyl,heroin,rx_opioids)
  
  theCounty<-theCounty %>% pivot_longer(cols = -signals, names_to = "drug", values_to = "percentage")
  theState<-theState %>% pivot_longer(cols = -signals, names_to = "drug", values_to = "percentage")
  
  theCounty$Region<-county
  theState$Region<-"Massachusetts"
  
  plot_data<-rbind(theCounty,theState)
  
  ggplot(plot_data, aes(x = drug, y = percentage, fill = Region)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ signals, labeller = labeller(signals = c(`0` = "Typical Day", `1` = "Spike"))) +
    labs(title = paste("Percentage of Deaths Involving Each Drug for",year),
         x = "Drug",
         y = "Percentage of Deaths",
         fill = "Region") +
    theme_minimal()
  
}

graphForestPlot<-function(CIs,theDrug)
{
  CIs$ci_lower<-exp(CIs$ci_lower)
  CIs$ci_upper<-exp(CIs$ci_upper)
  CIs$random_slope<-exp(CIs$random_slope)
  
  CIs<-CIs[CIs$drug==theDrug,]
  
  CIs <- CIs %>%
    mutate(color = case_when(
      ci_lower > 1 ~ "blue",
      ci_upper < 1 ~ "purple",
      TRUE ~ "grey"
    ))
  theCIs <- CIs %>%
    mutate(year_legend = factor(year),
      year = fct_rev(factor(year)),
           group = fct_rev(factor(group)))
  
  theCIs$Estimate_CI <- sprintf("%.2f (%.2f, %.2f)", theCIs$random_slope, theCIs$ci_lower, theCIs$ci_upper)
  
  pairwise_plot<-ggplot(theCIs, aes(y = year, x = random_slope , color = color, shape = year_legend, group = year)) +
    geom_line(size = 1) + 
    geom_vline(xintercept = 1, colour = "gray50", linetype = "dashed") +
    scale_x_continuous(
      breaks = seq(0, 5, by = 0.2),
      trans = "log10"
    ) +
    geom_point(size = 4) + 
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0.2, size = 1) +
    facet_wrap(~fct_rev(group), scales = "free_y", ncol=1, strip.position = "left") + 
    labs(title = "More on typical days                                                         More on spike days", x = "Adjusted odds ratio", y = "Counties") +
    theme_minimal()+
    theme(
      legend.position = 'bottom',
      legend.title=element_blank(),
      panel.spacing = unit(1, "lines"),
      axis.text.y = element_blank(),  
      strip.text = element_text(size = 14),  
      strip.text.y.left = element_text(angle = 0, hjust = 0, vjust = 1),
      plot.title = element_text(size = 20, hjust = 0.38, vjust = 2) 
    ) +
    theme(axis.title.x = element_text(face = "bold", margin = margin(t=10)), axis.title.y = element_text(face = "bold", margin=margin(r=10))) +
    guides(color = "none",
      shape = guide_legend(title = "Year"),
           strip.placement = "outside")+
    #scale_color_manual(values = c("grey"))
  scale_color_manual(values = c("blue", "grey", "purple"))
  
  return(pairwise_plot)
}

graphForestPlot_State<-function(CIs)
{
  CIs$ci_lower<-exp(CIs$ci_lower)
  CIs$ci_upper<-exp(CIs$ci_upper)
  CIs$estimates<-exp(CIs$estimates)
  CIs$drugs<-ifelse(CIs$drugs=="rx_opioids","rx opioids",CIs$drugs)
 CIs <- CIs %>%
    mutate(color = case_when(
      ci_lower > 1 ~ "black",
      ci_upper < 1 ~ "black",
      TRUE ~ "grey"
    ))
  theCIs <- CIs %>%
    mutate(year_legend = factor(year),
           year = fct_rev(factor(year)),
           group = fct_rev(factor(drugs)))
  
  theCIs$Estimate_CI <- sprintf("%.2f (%.2f, %.2f)", theCIs$estimates, theCIs$ci_lower, theCIs$ci_upper)
  
  pairwise_plot<-ggplot(theCIs, aes(y = year, x = estimates , color = color, shape = year_legend, group = year)) +
    geom_line(size = 1) + 
    geom_vline(xintercept = 1, colour = "gray50", linetype = "dashed") +
    scale_x_continuous(
      breaks = seq(-0.5, 3.5, by = 0.5),
      limits=c(-0.5,3.5)
      #trans = "log10"
    ) +
    geom_point(size = 4) + 
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0.2, size = 1) +
    facet_wrap(~fct_rev(group), scales = "free_y", ncol=1, strip.position = "left") + 
    labs(title = "Massachusetts \n \n                                   More on non-spike days              More on spike days", x = "Adjusted odds ratio", y = "Drug") +
    theme_minimal()+
    theme(
      legend.position = 'bottom',
      legend.title=element_blank(),
      legend.text = element_text(size = 16),
      panel.spacing = unit(1, "lines"),
      axis.text.y = element_blank(),  
      strip.text = element_text(size = 14),  
      strip.text.y.left = element_text(angle = 0, hjust = 0, vjust = 1),
      plot.title = element_text(size = 20, hjust = 0, vjust = 2) ,
      plot.title.position = "plot"
    ) +
    theme(axis.title.x = element_text(face = "bold", size = 16,margin = margin(t=10)), axis.title.y = element_text(face = "bold", size = 16, margin=margin(r=10))) +
    guides(color = "none",
           shape = guide_legend(title = "Year"),
           strip.placement = "outside")+
    scale_color_manual(values = c("grey"))
    #scale_color_manual(values = c("black", "grey", "black"))
  
  #annotate_figure(x, top = text_grob("Odds ratio confidence intervals for Massachussets", face = "bold", size = 20, hjust = 1.2))
  return(pairwise_plot)
}

graphForest_all<-function(CIs)
{
  CIs$ci_lower<-exp(CIs$ci_lower)
  CIs$ci_upper<-exp(CIs$ci_upper)
  CIs$random_slope<-exp(CIs$random_slope)
  
  #-----------COCAINE---------------
  
  theCIs<- CIs[CIs$drug =="cocaine",] 
  
  theCIs <- theCIs %>%
    mutate(year_legend = factor(year),
           year = fct_rev(factor(year)),
           group = fct_rev(factor(group)))
  
  theCIs <- theCIs %>%
    mutate(color = case_when(
      ci_lower > 1 ~ "black",
      ci_upper < 1 ~ "black",
      TRUE ~ "grey"
    ))
  
  theCIs$Estimate_CI <- sprintf("%.2f (%.2f, %.2f)", theCIs$random_slope, theCIs$ci_lower, theCIs$ci_upper)
  
  pCocaine<-ggplot(theCIs, aes(y = year, x = random_slope , color = color, shape = year_legend, group = year)) +
    geom_line(size = 1) + 
    geom_vline(xintercept = 1, colour = "gray50", linetype = "dashed") +
    scale_x_continuous(
      breaks = seq(-0.5, 4, by = 0.5),
      limits=c(-0.5,3.5),
    trans = "log10"
    ) +
    geom_point(size = 4) + 
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0.2, size = 1) +
    facet_wrap(~fct_rev(group), scales = "free_y", ncol=1, strip.position = "left") + 
    labs(title = "A) Cocaine \n  \n          Non-Spike                        Spike",x =  element_blank(), y =  element_blank()) +
    theme_minimal()+
    theme(
      legend.position = 'bottom',
      legend.title=element_blank(),
      legend.text = element_text(size = 14),
      panel.spacing = unit(1, "lines"),
      axis.text.y = element_blank(),  
      strip.text = element_text(size = 14),
      #strip.text.y.left = element_text(angle = 0, hjust = 0, vjust = 1),
      plot.title = element_text(size = 16)  
    ) +
    theme(axis.title.x = element_text(face = "bold", margin = margin(t=10)), axis.title.y = element_text(face = "bold", margin=margin(r=10))) +
    guides(color = "none",
           shape = guide_legend(title = "Year"),
           strip.placement = "outside")+
    #scale_color_manual(values = c("grey"))
    scale_color_manual(values = c("black", "grey", "black"))
    
  
  #-----------HEROIN---------------
  
  theCIs<- CIs[CIs$drug =="heroin",] 
  
  theCIs <- theCIs %>%
    mutate(year_legend = factor(year),
           year = fct_rev(factor(year)),
           group = fct_rev(factor(group)))
  
  theCIs <- theCIs %>%
    mutate(color = case_when(
      ci_lower > 1 ~ "black",
      ci_upper < 1 ~ "black",
      TRUE ~ "grey"
    ))
  
  theCIs$Estimate_CI <- sprintf("%.2f (%.2f, %.2f)", theCIs$random_slope, theCIs$ci_lower, theCIs$ci_upper)
  
  pHeroin<-ggplot(theCIs, aes(y = year, x = random_slope , color = color, shape = year_legend, group = year)) +
    geom_line(size = 1) + 
    geom_vline(xintercept = 1, colour = "gray50", linetype = "dashed") +
    scale_x_continuous(
      breaks = seq(-0.5, 8, by = 0.5),
      limits=c(-0.5,8),
    trans = "log10"
    ) +
    geom_point(size = 4) + 
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0.2, size = 1) +
    facet_wrap(~fct_rev(group), scales = "free_y", ncol=1, strip.position = "left") + 
    labs(title = "C) Heroin \n \n        Non-Spike                     Spike",x =  element_blank(), y =  element_blank()) +
    theme_minimal()+
    theme(
      legend.position = 'bottom',
      legend.text = element_text(size = 14),
      legend.title=element_blank(),
      panel.spacing = unit(1, "lines"),
      axis.text.y = element_blank(),  
      strip.text = element_text(size = 14),  
      #strip.text.y.left = element_text(angle = 0, hjust = 0, vjust = 1),
      plot.title = element_text(size = 16) 
    ) +
    theme(axis.title.x = element_text(face = "bold", margin = margin(t=10)), axis.title.y = element_text(face = "bold", margin=margin(r=10))) +
    guides(color = "none",
           shape = guide_legend(title = "Year"),
           strip.placement = "outside")+
    scale_color_manual(values = c("black", "grey", "black"))
  
  
  #-----------RX OPIOIDS---------------
    
  theCIs<- CIs[CIs$drug =="rx_opioids",] 
  
  theCIs <- theCIs %>%
    mutate(year_legend = factor(year),
           year = fct_rev(factor(year)),
           group = fct_rev(factor(group)))
  
  theCIs <- theCIs %>%
    mutate(color = case_when(
      ci_lower > 1 ~ "black",
      ci_upper < 1 ~ "black",
      TRUE ~ "grey"
    ))
  
  theCIs$Estimate_CI <- sprintf("%.2f (%.2f, %.2f)", theCIs$random_slope, theCIs$ci_lower, theCIs$ci_upper)
  
  pRx<-ggplot(theCIs, aes(y = year, x = random_slope , color = color, shape = year_legend, group = year)) +
    geom_line(size = 1) + 
    geom_vline(xintercept = 1, colour = "gray50", linetype = "dashed") +
    scale_x_continuous(
      breaks = seq(-0.5, 3.5, by = 0.5),
      limits=c(-0.5,3.5),
    trans = "log10"
    ) +
    geom_point(size = 4) + 
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0.2, size = 1) +
    facet_wrap(~fct_rev(group), scales = "free_y", ncol=1, strip.position = "left") + 
    labs(title = "D) RX Opioids \n  \n        Non-Spike                          Spike",x = element_blank(), y =  element_blank()) +
    theme_minimal()+
    theme(
      legend.position = 'bottom',
      legend.text = element_text(size = 14),
      legend.title=element_blank(),
      panel.spacing = unit(1, "lines"),
      axis.text.y = element_blank(),  
      strip.text = element_text(size = 14),  
      #strip.text.y.left = element_text(angle = 0, hjust = 0, vjust = 1),
      plot.title = element_text(size = 16) 
    ) +
    theme(axis.title.x = element_text(face = "bold", margin = margin(t=10)), axis.title.y = element_text(face = "bold", margin=margin(r=10))) +
    guides(color = "none",
           shape = guide_legend(title = "Year"),
           strip.placement = "outside")+
  scale_color_manual(values = c("black", "grey", "black"))
  
  
  #-----------FENTANYL---------------
  
  theCIs<- CIs[CIs$drug =="fentanyl",] 
  
  theCIs <- theCIs %>%
    mutate(year_legend = factor(year),
           year = fct_rev(factor(year)),
           group = fct_rev(factor(group)))
  
  theCIs <- theCIs %>%
    mutate(color = case_when(
      ci_lower > 1 ~ "black",
      ci_upper < 1 ~ "black",
      TRUE ~ "grey"
    ))
  
  theCIs$Estimate_CI <- sprintf("%.2f (%.2f, %.2f)", theCIs$random_slope, theCIs$ci_lower, theCIs$ci_upper)
  
  pFentanyl<-ggplot(theCIs, aes(y = year, x = random_slope , color = color, shape = year_legend, group = year)) +
    geom_line(size = 1) + 
    geom_vline(xintercept = 1, colour = "gray50", linetype = "dashed") +
    scale_x_continuous(
      breaks = seq(-0.5, 3.5, by = 0.5),
      limits=c(-0.5,3.5),
      trans = "log10"
    ) +
    geom_point(size = 4) + 
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0.2, size = 1) +
    facet_wrap(~fct_rev(group), scales = "free_y", ncol=1, strip.position = "left") + 
    labs(title = "B) Fentanyl \n  \n          Non-Spike                        Spike",x =  element_blank(), y =  element_blank()) +
    theme_minimal()+
    theme(
      legend.position = 'bottom',
      legend.title=element_blank(),
      legend.text = element_text(size = 14),
      panel.spacing = unit(1, "lines"),
      axis.text.y = element_blank(),  
      strip.text = element_text(size = 14),
      #strip.text.y.left = element_text(angle = 0, hjust = 0, vjust = 1),
      plot.title = element_text(size = 16)  
    ) +
    theme(axis.title.x = element_text(face = "bold", margin = margin(t=10)), axis.title.y = element_text(face = "bold", margin=margin(r=10))) +
    guides(color = "none",
           shape = guide_legend(title = "Year"),
           strip.placement = "outside")+
    #scale_color_manual(values = c("grey"))
    scale_color_manual(values = c("black", "grey", "black"))
  
  
  #-----------PSYCHOSTIMULANTS---------------
  
  theCIs<- CIs[CIs$drug =="psychostimulants",] 
  
  theCIs <- theCIs %>%
    mutate(year_legend = factor(year),
           year = fct_rev(factor(year)),
           group = fct_rev(factor(group)))
  
  theCIs <- theCIs %>%
    mutate(color = case_when(
      ci_lower > 1 ~ "black",
      ci_upper < 1 ~ "black",
      TRUE ~ "grey"
    ))
  
  theCIs$Estimate_CI <- sprintf("%.2f (%.2f, %.2f)", theCIs$random_slope, theCIs$ci_lower, theCIs$ci_upper)
  
  pPsycho<-ggplot(theCIs, aes(y = year, x = random_slope , color = color, shape = year_legend, group = year)) +
    geom_line(size = 1) + 
    geom_vline(xintercept = 1, colour = "gray50", linetype = "dashed") +
    scale_x_continuous(
      breaks = seq(-0.5, 4, by = 0.5),
      limits=c(-0.5,4.5),
      trans = "log10"
    ) +
    geom_point(size = 4) + 
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0.2, size = 1) +
    facet_wrap(~fct_rev(group), scales = "free_y", ncol=1, strip.position = "left") + 
    labs(title = "E) Psychostimulants \n  \n          Non-Spike                        Spike",x =  element_blank(), y =  element_blank()) +
    theme_minimal()+
    theme(
      legend.position = 'bottom',
      legend.title=element_blank(),
      legend.text = element_text(size = 14),
      panel.spacing = unit(1, "lines"),
      axis.text.y = element_blank(),  
      strip.text = element_text(size = 14),
      #strip.text.y.left = element_text(angle = 0, hjust = 0, vjust = 1),
      plot.title = element_text(size = 16)  
    ) +
    theme(axis.title.x = element_text(face = "bold", margin = margin(t=10)), axis.title.y = element_text(face = "bold", margin=margin(r=10))) +
    guides(color = "none",
           shape = guide_legend(title = "Year"),
           strip.placement = "outside")+
    #scale_color_manual(values = c("grey"))
    scale_color_manual(values = c("black", "grey", "black"))
  
  
    
  #thePlot<-ggarrange(pCocaine,pHeroin,pRx, nrow = 1,common.legend = TRUE,legend="bottom")
  thePlot<-ggarrange(pCocaine,pFentanyl,pHeroin, nrow = 1,common.legend = TRUE,legend="bottom")
  #thePlot<-ggarrange(pRx,pPsycho, nrow = 1,common.legend = TRUE,legend="bottom")
  
    final_plot <- annotate_figure(thePlot,
                                  left = text_grob("Counties", rot = 90, vjust = 1,size = 16),
                                  bottom = text_grob("Adjusted odds ratio", vjust = -3,size = 16))
    
  return(final_plot)
}


graphForestPlot_OneYear<-function(CIs)
{
  CIs$ci_lower<-exp(CIs$ci_lower)
  CIs$ci_upper<-exp(CIs$ci_upper)
  CIs$estimates<-exp(CIs$random_slope)
  
  CIs <- CIs %>%
    mutate(color = case_when(
      ci_lower > 1 ~ "black",
      ci_upper < 1 ~ "black",
      TRUE ~ "grey"
    ))
  theCIs <- CIs %>%
    mutate(Drug_legend = factor(drug),
           drug = fct_rev(factor(drug)),
           group = fct_rev(factor(group)))
  
  theCIs$Estimate_CI <- sprintf("%.2f (%.2f, %.2f)", theCIs$estimates, theCIs$ci_lower, theCIs$ci_upper)
  
  pairwise_plot<-ggplot(theCIs, aes(y = drug, x = estimates , color = color, shape = Drug_legend, group = drug)) +
    geom_line(size = 1) + 
    geom_vline(xintercept = 1, colour = "gray50", linetype = "dashed") +
    scale_x_continuous(
      breaks = seq(0.5, 2, by = 0.5),
      limits=c(0.5,2.5),
      trans = "log10"
    ) +
    geom_point(size = 4) + 
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0.2, size = 1) +
    facet_wrap(~fct_rev(group), scales = "free_y", ncol=1, strip.position = "left") + 
    labs(title = "Results by County \n \n                                   More on non-spike days              More on spike days", x = "Adjusted odds ratio", y = "County") +
    theme_minimal()+
    theme(
      legend.position = 'bottom',
      legend.title=element_blank(),
      legend.text = element_text(size = 16),
      panel.spacing = unit(1, "lines"),
      axis.text.y = element_blank(),  
      strip.text = element_text(size = 14),  
      strip.text.y.left = element_text(angle = 0, hjust = 0, vjust = 1),
      plot.title = element_text(size = 20, hjust = 0, vjust = 2) ,
      plot.title.position = "plot"
    ) +
    theme(axis.title.x = element_text(face = "bold", size = 16,margin = margin(t=10)), axis.title.y = element_text(face = "bold", size = 16, margin=margin(r=10))) +
    guides(color = "none",
           shape = guide_legend(title = "Drug"),
           strip.placement = "outside")+
    #scale_color_manual(values = c("grey"))
    scale_color_manual(values = c("black", "grey", "black"))
  
  #annotate_figure(x, top = text_grob("Odds ratio confidence intervals for Massachussets", face = "bold", size = 20, hjust = 1.2))
  return(pairwise_plot)
}

graphTable<-function(CIs,theDrug)
{
  theCIs$ci_lower<-exp(theCIs$ci_lower)
  theCIs$ci_upper<-exp(theCIs$ci_upper)
  theCIs$random_slope<-exp(theCIs$random_slope)
  
  theCIs<-theCIs[theCIs$drug==theDrug,]
  
  theCIs <- theCIs %>%
    mutate(color = case_when(
      ci_lower > 1 ~ "blue",
      ci_upper < 1 ~ "lightblue",
      TRUE ~ "grey"
    ))
  theCIs <- theCIs %>%
    mutate(year = fct_rev(factor(year)),
           group = fct_rev(factor(group)))
  
  theCIs$Estimate_CI <- sprintf("%.2f (%.2f, %.2f)", theCIs$random_slope, theCIs$ci_lower, theCIs$ci_upper)
  
  dat_table <- theCIs %>%
    select(year, group, Estimate_CI) %>%
    tidyr::pivot_longer(c(Estimate_CI), names_to = "stat") %>%
    mutate(stat = factor(stat, levels = c("Estimate_CI")))
  
  pairwise_table <- ggplot(dat_table, aes(stat, year, label = value)) +
    geom_text(size = 4) +
    scale_x_discrete(position = "top", labels = c("aOR (95% CI)")) +
    facet_wrap(~fct_rev(group), scales = "free_y", ncol=1, strip.position = "left", labeller = labeller(group = c(Bristol = "", Essex= "", Hampden = "", Middlesex = '', Suffolk = "", Worcester = ''))) +
    labs(y = NULL, x = NULL) +
    theme_classic() +
    theme(
      panel.spacing = unit(1, "lines"),
      strip.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 12, face ="bold"),
      strip.text.x = element_text(face="bold", size = 12),
      strip.text.y.left = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_text(face = "bold", size = 11),
      plot.title = element_text(size = 11)
    ) + ggtitle("") 
  return(pairwise_table)
}

tableOne<-function(drug_categories,county_sens,lag,threshold)
{
  
  theDrugs<-drug_categories %>%select(DOD_4_FD,YOD,County_name,AGE1_CALC,SEX,RACEGROUP,cocaine,psychostimulants,fentanyl,heroin,rx_opioids)
  theSens<-county_sens[county_sens$threshold==threshold & county_sens$lag==lag,]
  theSens<-theSens%>%select(date,County_name,signals)
  theSens$signals<-ifelse(theSens$signals>0,1,0)
  
  theData<- theDrugs %>% left_join(theSens,by=c("DOD_4_FD"="date","County_name"="County_name"))
  #theData<-theData[theData$County_name!="Barnstable" & theData$County_name!="Dukes" & theData$County_name!="Nantucket",]
  #theData<-theData[theData$County_name %in% c("Bristol","Essex","Hampden","Middlesex","Suffolk","Worcester"),]
  theData<-theData[theData$YOD>=2020 & theData$YOD<=2023,]
  theData$RACEGROUP<-ifelse(theData$RACEGROUP=="Black NH","Black",ifelse(theData$RACEGROUP=="White NH","White",ifelse(theData$RACEGROUP=="Hispanic","Hispanic","Other")))
  
  
  listVars<-c("AGE1_CALC")
  theData$SEX<-as.factor(theData$SEX)
  theData$RACEGROUP<-as.factor(theData$RACEGROUP)
  theData$County_name<-as.factor(theData$County_name)
  catVars<-c("cocaine","psychostimulants","fentanyl","heroin","rx_opioids","SEX","RACEGROUP","County_name")
  
  #table1 <- CreateTableOne(vars = c(listVars, catVars), data = theData, factorVars = catVars,strata = c("signals"))
  table1 <- CreateTableOne(vars = c(listVars, catVars), data = theData, factorVars = catVars)
  
  return(table1)
}

graphMap<-function()
{
  counties <- counties(state = "MA", cb = TRUE, class = "sf")
  
  # Example data frame
  county_data <- data.frame(
    NAME = c("Barnstable", "Berkshire", "Bristol", "Dukes", "Essex", "Franklin", "Hampden", "Hampshire", "Middlesex", "Nantucket", "Norfolk", "Plymouth", "Suffolk", "Worcester"),
    value = c(143.2,190.7,194.4,97.1,134.7,154.9,215.9,81.9,92.6,90.9,82.5,130.7,217.6,150.9)  # Replace this with your actual datat
  )
  
  county_data$value<-county_data$value/4
  
  # Merge the data with the shapefile
  counties <- counties %>%
    left_join(county_data, by = "NAME")
  
  # Plot the map with labels
  ggplot(data = counties) +
    geom_sf(aes(fill = value)) +
    scale_fill_gradient(low = "white", high = "grey60", name = "Mortality rate") +
    geom_sf_text(aes(label = paste(NAME, round(value, 1), sep = "\n")), size = 5, color = "black",fontface = "bold") +
    theme_minimal() +
    labs(title = "Massachusetts counties with overdose mortality rate (per 100,000 people)",size=10) +
    theme(legend.position = "bottom",
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(size = 20, face = "bold"),  
          legend.title = element_text(size = 16))
}

graphMapDrugs<-function()
{
  
  counties <- counties(state = "MA", cb = TRUE, class = "sf")
  
  county_data <- data.frame(
    NAME = c("Barnstable", "Berkshire", "Bristol", "Dukes", "Essex", "Franklin", "Hampden", "Hampshire", "Middlesex", "Nantucket", "Norfolk", "Plymouth", "Suffolk", "Worcester"),
    value = c(0, 0, 1,0, 1, 0, 1, 0,1,0, 1,1, 1, 1),
    additional_info = c("", "", "Cocaine", "", "Fentanyl", "", "Cocaine", "", "", "", "Psychostimulants", "Fentanyl \n Heroin", "Psychostimulants", ""),  # Replace with actual additional data
    group = c("Not included", "Not included", "Included", "Not included", "Included", "Not included", "Included","Not included", "Included", "Not included", "Included", "Included", "Included", "Included") 
  )
  
  # Divide value by 4
  county_data$value <- county_data$value / 4
  
  # Merge the data with the shapefile
  counties <- counties %>%
    left_join(county_data, by = "NAME")
  
  group_colors <- c("Not included" = "white", "Included" = "grey")
  
  
  # Plot the map with labels
  ggplot(data = counties) +
    geom_sf(aes(fill = as.factor(group))) +  # Fill based on the group
    scale_fill_manual(values = group_colors, name = "Group") +  # Use the custom color palette
    geom_sf_text(aes(label = paste(NAME)), size = 5, color = "black", fontface = "bold") +
    theme_minimal() +
    labs(title = "Massachusetts counties with overdose mortality rate (per 100,000 people)", size = 10) +
    theme(legend.position = "bottom",
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(size = 20, face = "bold"),
          legend.title = element_text(size = 16))
  
}