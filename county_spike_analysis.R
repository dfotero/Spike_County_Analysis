# ========================================
# Script Name: county_spike_analysis.R
# Author: Daniel Otero-Leon
# Date: 2024-11-20
# Description: This script analize the drugs involve in spikes for the state of Massachusetts and each of its counties
# ========================================

# =============================
# Section 1: Libraries
# =============================

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


# =============================
# Section 2: Data preparation
# =============================

#Total overdose deaths per day and county
daily_deaths_county <- read.csv("daily_deaths_county.csv")

#Total overdose deaths per day
daily_deaths_state <- read.csv("daily_deaths_state.csv")

#Each overdose death with the following information: demographics, date, drugs involved
deaths <- read.csv("deaths.csv")

#census information for each county in MA
census<-read.csv("census.csv")

# =============================
# Section 3: Spike identification
# =============================

#' Algorithm to identify spikes from a time series 
#'
#' This function detects anomalies (spikes and dips) in a time series based on a thresholding approach.
#' It identifies data points that deviate significantly from a rolling mean
#' 
#' @param y A numeric vector representing the time series data.
#' @param lag An integer specifying the lag window size for the rolling mean and standard deviation.
#' @param threshold A numeric value specifying the number of standard deviations from the mean to consider a point as an anomaly.
#' @param influence A numeric value (between 0 and 1) that determines the influence of anomalous points on the rolling statistics.
#' @return A list containing the following elements:
#'   \item{signals}{A numeric vector indicating anomalies. Values are 1 for positive anomalies (spikes), -1 for negative anomalies (dips), and 0 for normal points.}
#'   \item{avgFilter}{A numeric vector representing the rolling mean of the filtered time series.}
#'   \item{stdFilter}{A numeric vector representing the rolling standard deviation of the filtered time series.}
#'   \item{spikes}{A numeric vector where spikes (positive anomalies) are recorded. Non-spike points are 0.}
#'   \item{dips}{A numeric vector where dips (negative anomalies) are recorded. Non-dip points are 0.}
#' @export

ThresholdingAlgo <- function(y,lag,threshold,influence) {
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

#' Analyze State-Level Overdose Death Trends
#'
#' This function performs a multi-parameter analysis of daily overdose deaths in a given state. 
#' It applies the `ThresholdingAlgo` function to detect anomalies (spikes and dips) based on 
#' different combinations of lag, threshold, and influence parameters. The results are organized
#' into a data frame for further analysis.
#'
#' @param daily_deaths_state A data frame containing state-level overdose death data. 
#'   Must include a column `Total_Overdose_Deaths` with daily overdose death counts.
#' @return A data frame with the following columns:
#'   \item{date}{The date corresponding to each observation.}
#'   \item{y}{The original daily overdose death counts.}
#'   \item{avgFilter}{The rolling average filter calculated using `ThresholdingAlgo`.}
#'   \item{upperThreshold}{The upper threshold for anomaly detection (average + threshold * std).}
#'   \item{signals}{Anomaly signals: 1 for spikes, -1 for dips, 0 for normal.}
#'   \item{spikes}{The value of the spikes (positive anomalies).}
#'   \item{dips}{The value of the dips (negative anomalies).}
#'   \item{year}{The year extracted from the `date` column.}
#'   \item{modified_y}{The death count only for anomalies, 0 otherwise.}
#'   \item{y_minus_avgFilter_or_zero}{Difference between actual deaths and the rolling average for anomalies, 0 otherwise.}
#'   \item{lag}{The lag parameter used for the calculation.}
#'   \item{threshold}{The threshold parameter used for the calculation.}
#'   \item{influence}{The influence parameter used for the calculation.}
#'   \item{percSpike}{Percentage increase of spikes over the rolling average.}
#'
#'
#' @export
stateAnalysis <- function(daily_deaths_state) {
  lag <- c(7, 30, 90)
  threshold <- 2
  influence <- 0.3
  
  y <- daily_deaths_state$Total_Overdose_Deaths
  theResult <- 0
  
  for (i in 1:length(lag)) {
    for (j in 1:length(threshold)) {
      for (k in 1:length(influence)) {
        result <- ThresholdingAlgo(y, lag[i], threshold[j], influence[k])
        df <- data.frame(
          date = daily_deaths_state$DOD_4_FD,
          y = y,
          avgFilter = result$avgFilter,
          upperThreshold = result$avgFilter + threshold[j] * result$stdFilter,
          signals = result$signals,
          spikes = result$spikes,
          dips = result$dips,
          YOD = geographical_data$YOD
        )
        df$year <- year(df$date)
        df_modified <- df
        df_modified$modified_y <- ifelse(df$signals == 0, 0, df$y)
        df_modified$y_minus_avgFilter_or_zero <- ifelse(df$signals == 0, 0, df$y - df$avgFilter)
        df_modified$lag <- lag[i]
        df_modified$threshold <- threshold[j]
        df_modified$influence <- influence[k]
        df_modified$percSpike <- ifelse(df_modified$signals == 1, df_modified$y_minus_avgFilter_or_zero / df_modified$avgFilter, 0)
        df_modified <- df_modified[df_modified$year >= 2017 & df_modified$year <= 2023, ]
        if (i == 1 && j == 1 && k == 1) {
          theResult <- df_modified
        } else {
          print(paste(i, j, k))
          theResult <- rbind(theResult, df_modified)
        }
      }
    }
  }
  return(theResult)
}

#' Analyze County-Level Overdose Death Trends
#'
#' This function performs anomaly detection and analysis on county-level overdose death data
#' using the `ThresholdingAlgo` function. It iterates over multiple counties and parameter
#' combinations to identify spikes and dips in overdose death trends and aggregates the results
#' into a comprehensive data frame.
#'
#' @param daily_deaths_county A data frame containing county-level overdose death data. 
#'   Must include columns `County_name` (county names), `Total_Overdose_Deaths` (daily overdose death counts), 
#'   and `DOD_4_FD` (date of death).
#' @return A data frame with the following columns:
#'   \item{date}{The date corresponding to each observation.}
#'   \item{y}{The original daily overdose death counts for each county.}
#'   \item{avgFilter}{The rolling average filter calculated using `ThresholdingAlgo`.}
#'   \item{upperThreshold}{The upper threshold for anomaly detection (average + threshold * std).}
#'   \item{signals}{Anomaly signals: 1 for spikes, -1 for dips, 0 for normal.}
#'   \item{spikes}{The value of the spikes (positive anomalies).}
#'   \item{dips}{The value of the dips (negative anomalies).}
#'   \item{YOD}{Additional geographical information from the data.}
#'   \item{year}{The year extracted from the `date` column.}
#'   \item{modified_y}{The death count only for anomalies, 0 otherwise.}
#'   \item{y_minus_avgFilter_or_zero}{Difference between actual deaths and the rolling average for anomalies, 0 otherwise.}
#'   \item{County_name}{The name of the county corresponding to the data.}
#'   \item{lag}{The lag parameter used for the calculation.}
#'   \item{threshold}{The threshold parameter used for the calculation.}
#'   \item{influence}{The influence parameter used for the calculation.}
#'   \item{percSpike}{Percentage increase of spikes over the rolling average.}
#'
#' @export
countyAnalysis <- function(daily_deaths_county)
{

  lag<-c(7,30,60,90)
  threshold<-c(2,2.5)
  influence<-0.3
  
  theCounties<-count(daily_deaths_county,"County_name"=daily_deaths_county$County_name)
  theCounties<-na.omit(theCounties)
  theCounties<-theCounties$County_name
  daily_deaths_county<-na.omit(daily_deaths_county)
  
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
          
          y <- daily_deaths_county[daily_deaths_county$County_name==theCounties[i],]
          
          theDeaths<-y$Total_Overdose_Deaths
          result <- ThresholdingAlgo(theDeaths,lag[j],threshold[k],influence[l])
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
  return(theResult)
}


state_spikes<-stateAnalysis(daily_deaths_state)
county_spikes<-countyAnalysis(daily_deaths_county)

# ===============================
# Section 4: Logistic Regressions
# ===============================


#' Logistic Regression for State-Level Drug Sensitivity Analysis
#'
#' This function performs a logistic regression analysis to assess the relationship 
#' between drug use, demographic, and socioeconomic factors, and the occurrence of 
#' spikes in a given state. It integrates drug-related data, 
#' census data, and spikes in the state.
#'
#' @param drug_categories A data frame containing drug-related information
#' @param state_sens A data frame containing spikes for the state.
#' @param census A data frame containing census data with demographic and socioeconomic information.
#' @param threshold Numeric. The threshold used in the spike detection algorithm.
#' @param lag Numeric. The lag value used in the spike detection algorithm
#' @param year Numeric. A specific year for analysis. If set to 0, all years from 2020 to 2023 are included in the analysis.
#' @return A fitted logistic regression model with predictors related to drug categories, 
#' demographic factors, and socioeconomic variables. The dependent variable is the presence of a spike.
#'
#' @export

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

  }
  return(model)
  
}

#' Logistic Mixed-Effects Analysis for County-Level Overdose Data
#'
#' This function performs a logistic mixed-effects analysis of county-level overdose data 
#' to evaluate the relationship between drug categories, demographic and socioeconomic factors, 
#' and the presence of spikes. It models drug-specific effects across selected counties 
#' using mixed-effects models with random intercepts and slopes.
#'
#' @param drug_categories A data frame containing drug-related information
#' @param county_sens A data frame containing spikes for each county.
#' @param census A data frame containing census data with demographic and socioeconomic information.
#' @param threshold Numeric. The threshold used in the spike detection algorithm.
#' @param lag Numeric. The lag value used in the spike detection algorithm
#' @param year Numeric. A specific year for analysis. If set to 0, all years from 2020 to 2023 are included in the analysis.
#' @return A data frame containing the confidence intervals for each drug category, modeled separately. 
#' Columns in the output include:
#'   \itemize{
#'     \item \code{Drug}: The drug category being modeled.
#'     \item \code{Lower_CI}: Lower bound of the confidence interval.
#'     \item \code{Estimate}: Model estimate for the drug effect.
#'     \item \code{Upper_CI}: Upper bound of the confidence interval.
#'   }
#' @export

logTestCounties<-function(drug_categories,county_sens,census,threshold,lag,year)
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


log_model_state<-logTestState(deaths,state_spikes,census,threshold,lag,year)
confidence_interval_counties<-logTestCounties(deaths,county_sspikes,census,threshold,lag,year)

# ===============================
# Section 5: Graphics
# ===============================

#' Generate Forest Plots for Drug-Specific Adjusted Odds Ratios
#'
#' This function creates forest plots for adjusted odds ratios (AORs) of multiple drug categories 
#' (e.g., cocaine, heroin, fentanyl, prescription opioids, and psychostimulants). The plots visualize 
#' the random slope estimates and confidence intervals for each drug category across different years and groups.
#'
#' @param CIs A data frame containing confidence intervals and model estimates for drug categories. 
#'   The data frame must include the following columns:
#' @return A \code{ggplot} object containing annotated forest plots for selected drug categories.
#'   Each plot visualizes the adjusted odds ratios for a drug across years and groups.
#' @export

graphForest<-function(CIs)
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
