library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)
library(readxl)
library(writexl)
library(formattable)
library(plotly)
library(ggrepel) #to prevent labels from covering each other
library(viridisLite) #provide the turbo 256 color scale

dictionary = read_xlsx("Data/dictionary.xlsx")

build_theo = function(known_sequence, reference, ladder_5){
  return_df = data.frame(matrix(ncol = 2))
  colnames(return_df) = c("base_name", "theoretical_mass")
  if(ladder_5 == TRUE){
    return_df[1,1] = "5'"
    return_df[1,2] = 18.015
    #return_df[1,2] = 97.9769
  } else {
    return_df[1,1] = "3'"
    return_df[1,2] = -61.95579
  }
  i = 1
  while(i <= nchar(known_sequence)) {
    next_base = substr(known_sequence, i, i)
    if(next_base == "("){
      next_base = substr(known_sequence, i, i+3)
    }
    for(j in 1:nrow(reference)){
      if(next_base == reference[j,1]){
        return_df[nrow(return_df) + 1,2] = return_df[nrow(return_df),2] + reference[j,2]
        return_df[nrow(return_df),1] = next_base
        break
      }
    }
    if(nchar(next_base) == 4){
      i = i+3
    } else{
      i = i+1
    }
  }
  return_df = return_df %>% 
    mutate(n_position = row_number() - 1)
  return(return_df)
}

prophet = function(df, theo_df){
  return_df = data.frame(matrix(ncol = 6))
  colnames(return_df) = c("base_name", "theoretical_mass", "n_position", "monoisotopic_mass", "apex_rt", "sum_intensity")
  return_df = return_df %>% 
    drop_na()
  for(i in 1 : nrow(theo_df)) {
    for(j in 1 : nrow(df)){
      if(ppm(df[j,1], theo_df[i,2])){
        temp_row = data.frame(matrix(ncol = 6))
        colnames(temp_row) = c("base_name", "theoretical_mass", "n_position", "monoisotopic_mass", "apex_rt", "sum_intensity")
        temp_row[1,1] = theo_df[i,1]
        temp_row[1,2] = theo_df[i,2]
        temp_row[1,3] = theo_df[i,3]
        temp_row[1,4] = df[j,1]
        temp_row[1,5] = df[j,2]
        temp_row[1,6] = df[j,3]
        return_df = rbind(return_df, temp_row)
        break
      }
    }
  }
  return(return_df)
}

ppm = function(observed, theo){
  if(abs((observed - theo) / theo * 10^6) > 10) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}