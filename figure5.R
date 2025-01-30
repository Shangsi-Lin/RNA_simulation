library(tidyverse)
library(dplyr)
library(readxl)
library(writexl)

ppm = function(observed, theo){
  if(abs((observed - theo) / theo * 10^6) > 10) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

dictionary = read_xlsx("Data/dictionary.xlsx")
dictionary = dictionary[-c(5, 9),]