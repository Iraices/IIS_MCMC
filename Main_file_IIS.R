setwd("")

install.packages("extrafontdb")
library(extrafont)
font_import()
loadfonts()

library(rstan)
library(readxl)
library(dplyr)
library(LaplacesDemon)
library(optimx)
library(lattice)
library(latticeExtra)
library(broom)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')

## 
source('Prior_Predictive_Check.R')

source('generate_sets_of_priors.R')

source('Functions_For_IIS.R')

source('Simulation_IIS.R')

source('Test_function.R')


