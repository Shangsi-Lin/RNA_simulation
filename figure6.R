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

load_df = function(file_location){
  df = read_xlsx(file_location) %>% 
    janitor::clean_names() %>% 
    drop_na() %>% 
    select(monoisotopic_mass, apex_rt, sum_intensity, relative_abundance, number_of_detected_intervals) %>% 
    arrange(desc(monoisotopic_mass)) %>%
    mutate(sum_intensity = as.numeric(sum_intensity), relative_abundance = as.numeric(relative_abundance), number_of_detected_intervals = as.numeric(number_of_detected_intervals))
}

build_theo = function(known_sequence, reference, ladder_5){
  return_df = data.frame(matrix(ncol = 2))
  colnames(return_df) = c("base_name", "theoretical_mass")
  if(ladder_5 == TRUE){
    return_df[1,1] = "5'"
    return_df[1,2] = 18.015
  } else {
    return_df[1,1] = "3'"
    return_df[1,2] = -61.95579
  }
  i = 1
  while(i <= nchar(known_sequence)) {
    next_base = substr(known_sequence, i, i)
    for(j in 1:nrow(reference)){
      if(next_base == reference[j,1]){
        return_df[nrow(return_df) + 1,2] = return_df[nrow(return_df),2] + reference[j,2]
        return_df[nrow(return_df),1] = next_base
        break
      }
    }
    i = i+1
  }
  return_df = return_df %>% 
    mutate(n_position = row_number() - 1)
  if(ladder_5){
    return_df[nrow(return_df), 2] = return_df[nrow(return_df), 2] - 79.97107
  }
  return(return_df)
}

prophet = function(df, theo_df){
  return_df = data.frame(matrix(ncol = 7))
  colnames(return_df) = c("base_name", "theoretical_mass", "n_position", "monoisotopic_mass", "sum_intensity", "apex_rt", "relative_abundance")
  return_df = return_df %>% 
    drop_na()
  for(i in 1 : nrow(theo_df)) {
    for(j in 1 : nrow(df)){
      if(ppm(df[j,1], theo_df[i,2])){
        temp_row = data.frame(matrix(ncol = 7))
        colnames(temp_row) = c("base_name", "theoretical_mass", "n_position", "monoisotopic_mass", "sum_intensity", "apex_rt", "relative_abundance")
        temp_row[1,1] = theo_df[i,1]
        temp_row[1,2] = theo_df[i,2]
        temp_row[1,3] = theo_df[i,3]
        temp_row[1,4] = df[j,1]
        temp_row[1,5] = df[j,3]
        temp_row[1,6] = df[j,2]
        temp_row[1,7] = df[j,4]
        return_df = rbind(return_df, temp_row)
        break
      }
    }
  }
  return(return_df)
}