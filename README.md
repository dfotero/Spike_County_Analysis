#### **County Spike Analysis Script for Overdose Trends in Massachusetts**

**Author**: Daniel Otero-Leon  
**Date**: 2024-11-20  

This script, `county_spike_analysis.R`, is a tool designed to analyze spikes in overdose deaths for the state of Massachusetts and its counties. It integrates spike detection in time series, logistic regression modeling, and visualization techniques to uncover insights into the drugs involved in overdose spikes.

---

### **Key Features**

1. **Spike Detection**  
   - Utilizes a custom thresholding algorithm to identify significant deviations (spikes) in overdose trends across the state and counties.

2. **Data Preparation**  
   - Processes multiple datasets, including:
     - Daily overdose deaths at the county and state levels.
     - Demographic and socioeconomic census data.
     - Overdose death records with detailed drug information.

3. **State and County-Level Analysis**  
   - Applies logistic regression models to evaluate the relationships between demographic factors, socioeconomic variables, and drug involvement in overdose spikes.
   - Implements mixed-effects models for county-level analysis to account for variations across counties.

4. **Interactive Visualizations**  
   - Generates forest plots to visualize adjusted odds ratios (AORs) and confidence intervals for various drug categories (e.g., cocaine, heroin, fentanyl, psychostimulants, and prescription opioids).

---

### **Main Functions**

- **`ThresholdingAlgo`**  
  - Detects anomalies in a time series based on deviations from rolling averages.

- **`stateAnalysis` and `countyAnalysis`**  
  - Analyze overdose death trends at the state and county levels using spike detection algorithms.

- **`logTestState` and `logTestCounties`**  
  - Perform logistic regression analyses to evaluate the impact of drug involvement, demographics, and socioeconomic factors on overdose spikes.

- **`graphForest`**  
  - Creates annotated forest plots to visualize AORs for different drugs across years and groups.

---

### **Required Libraries**

This script relies on several R packages, including:
- **Data Manipulation**: `dplyr`, `data.table`, `sf`
- **Visualization**: `ggplot2`, `ggpubr`, `cowplot`
- **Modeling**: `lme4`, `randtests`, `emmeans`
