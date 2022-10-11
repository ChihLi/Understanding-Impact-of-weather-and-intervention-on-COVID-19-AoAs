# Reproducibility for Sung (2022) Estimating Functional Parameters for Understanding the impact of Weather and Government Interventions on COVID‐19 Outbreak
This folder consists of the data and R code for the manuscript [Estimating Functional Parameters for Understanding the impact of Weather and Government Interventions on COVID‐19 Outbreak](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-16/issue-4/Estimating-functional-parameters-for-understanding-the-impact-of-weather-and/10.1214/22-AOAS1601.short) by Sung (2022). 

* The code folder reproduces the results in the Sections 4 and 5 of the manuscript. 
  * The required R packages include `plgp`, `ggplot2`, `gbm`, `gridExtra`, `PerformanceAnalytics`, and `ggpubr`
  * `main.R` is the master file, which reproduces all the figures of the manuscript.
  * Note: reproduces all the results may take more than 24 hours, depending on the computer resourses.
* The data folder contains the required data:
  * `climate_city_test.csv`: weather and government interventions in the test data
  * `climate_city.csv`: weather and government interventions in the training data
  * `us-counties.csv`: US county level covid data, including daily cases and deaths
  * `city_population.csv`: city population
  * `county.txt`: link city to county
