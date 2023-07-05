# R Package for batch installation
packages <- c("deSolve","bbmle","stats","ggplot2","broom","tidyverse", "mgcv")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(packages)
# R Package for batch calls
lapply(packages, function(x){library(x, character.only = T)}) 

