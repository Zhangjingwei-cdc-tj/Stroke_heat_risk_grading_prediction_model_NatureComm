# Author: Jingwei Zhang
# Date: July 28，2025
# Purpose: Statistical Model for the Development and 
#          Validation of a Stroke Heat Risk Grading Prediction Model

# NOTE: this code is a guide for transparency and 
#       reproducibility and is not able to be run

# Key packages used
library(dlnm)       # For distributed lag nonlinear models (DLNM)
library(splines)    # For spline functions
library(survival)   # For survival analysis (used in case-crossover)


# 1.Construction of the Stroke Heat Risk Grading Prediction Model  
# NOTE: time-series analyses were conducted using quasi-Poisson regression combined with DLNM, 
#       separately within age- and gender-specific subgroups, based on stroke mortality, 
#       meteorological, and air pollution data during the summer months from 2013 to 2018.
#      


# Model temperature using a natural cubic spline with three internal knots
perc <- c(10,75,90)
argvar <- list(fun="ns", df=3, knots=quantile(
                                              #Tmeanc is the daily mean temperature
                                              #data1 is county-level time-series data 
                                              #on stroke mortality from 2013 to 2018
                                              data1$Tmeanc, perc/100,na.rm=T))

# Define lag structure on log scale
arglag <- list(knots=logknots(3, nk=2, df=4))

# Create crossbasis function for temperature
cb <- crossbasis(data1$Tmeanc, lag=3, argvar=argvar, arglag=arglag)

# Fit the quasi-Poisson model with adjustment variables
# NOTE: Fit the quasi-Poisson model with adjustment variables 
#       for one selected county as an example
model <- glm(
  # indicating the county-level daily stroke mortality 
  str ~ 
    cb + 
    
    # adjust for daily mean relative humidity with 3 df
    ns(Rhu_mean, df=3) + 
    
    # adjust for daily mean wind speed with 3 df
    ns(Wind_mean, df=3) +   
    
    # adjust for daily average concentrations of PM2.5 and O3
    PM2.5 + M8H_O3 +   
    
    # adjust for day of week
    as.factor(dow) +     
    
    # adjust for long-term trend (7 df per year for 6 years)(2013-2018)
    ns(time, df=7*6),      
  #data1 is county-level time-series data on stroke mortality from 2013 to 2018
  data=data1,
  family = quasipoisson(link = "log"),
  control = glm.control(epsilon = 10E-8, maxit = 5000))

# 2.Validation of the Stroke Heat Risk Grading Prediction Model
# NOTE: As a validation analysis, case-crossover analyses were conducted based on stroke mortality, 
#       meteorological, and air pollution data during the summer months from 2019 to 2022.

# Fit conditional logistic regression
model <- clogit(
  case ~ 
    #an ordinal variable with “1” for identified "moderate risk" level, 
    #“2” for identified "high risk" level,
    #“3” for identified "extremely high risk" level,
    #and “0” for identified "low risk" level
    risk_level + 
    
    # adjust for daily mean relative humidity with 3 df
    ns(Rhu_mean, df = 3) + 
    
    # adjust for daily mean wind speed with 3 df
    ns(Wind_mean, df = 3) + 
    
    # adjust for daily average concentrations of PM2.5 and O3
    PM2.5 + M8H_O3 + 
    
    # Match on case-control sets
    strata(ID),  
  #data2 is individual-level stroke mortality data from 2019 to 2022
  data = data2,
  #breslow method is a maximum likelihood approach for matched case-control studies 
  method = "breslow"
)

# Calculate OR for each risk level
or_results <- as.data.frame(summary(model)$coefficients) %>%
  mutate(
    OR = exp(coef),
    OR_low = exp(coef - 1.96 * `se(coef)`),
    OR_high = exp(coef + 1.96 * `se(coef)`)
  )
