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

# get the maximum mass in the dictionary for future processing.
mass_bound = 360

ppm = function(observed, theo){
  if(abs((observed - theo) / theo * 10^6) > 0.001) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

#filters the remaining data frame during blind sequencing to increse efficiency
filter_desc = function(df, mass_bound, version){
  if(version == TRUE){
    return(filter(df, monoisotopic_mass >= (df$monoisotopic_mass[1] - mass_bound)))
  } else {
    return(filter(df, monoisotopic_mass <= (df$monoisotopic_mass[1] + mass_bound)))
  }
}

#Match mass delta with the dictionary to do base calling, going descending
matcher_desc = function(match_df, match_dict, nth_attempt) {
  found_match = 0
  return_df = data.frame(matrix(ncol = 5)) # create the df to be returned
  colnames(return_df) = c("base_name", "monoisotopic_mass", "sum_intensity", "apex_rt", "n_iteration")
  for(i in 2:nrow(match_df)) {
    for(j in 1:nrow(match_dict)) {
      #if(ppm((match_df$monoisotopic_mass[1] - match_df$monoisotopic_mass[i]), match_dict$mass[j])){
      if(ppm((match_df$monoisotopic_mass[1] - match_dict$mass[j]), match_df$monoisotopic_mass[i])){
        found_match = found_match + 1 #numbers of matches found
        return_df[found_match,1] = dictionary$name[j]
        return_df[found_match,2] = match_df$monoisotopic_mass[i]
        return_df[found_match,3] = match_df$sum_intensity[i]
        return_df[found_match,4] = match_df$apex_rt[i]
        return_df[found_match,5] = nth_attempt
      }
    }
  }
  return_df = return_df %>%
    filter(sum_intensity == max(sum_intensity))
  return(return_df)
}

#Match mass delta with the dictionary to do base calling, going ascending
matcher_asce = function(match_df, match_dict, nth_attempt) {
  found_match = 0
  return_df = data.frame(matrix(ncol = 5)) # create the df to be returned
  colnames(return_df) = c("base_name", "monoisotopic_mass", "sum_intensity", "apex_rt", "n_iteration")
  for(i in 1:nrow(match_df - 1)) {
    for(j in 1:nrow(match_dict)) {
      #if(ppm((match_df$monoisotopic_mass[i] - match_df$monoisotopic_mass[nrow(match_df)]), match_dict$mass[j])){
      if(ppm((match_df$monoisotopic_mass[i] - match_dict$mass[j]), match_df$monoisotopic_mass[nrow(match_df)])){
        found_match = found_match + 1 #numbers of matches found
        return_df[found_match,1] = dictionary$name[j]
        return_df[found_match,2] = match_df$monoisotopic_mass[i]
        return_df[found_match,3] = match_df$sum_intensity[i]
        return_df[found_match,4] = match_df$apex_rt[i]
        return_df[found_match,5] = nth_attempt
      }
    }
  }
  return_df = return_df %>%
    filter(sum_intensity == max(sum_intensity))
  return(return_df)
}

#loop_down_match function
loop_down_match = function(df, mass_bound, dictionary, return_df, nth_attempt, begin){
  filtered_df = filter_desc(df, mass_bound, version = TRUE) #try to find the first match
  if(nrow(filtered_df) > 1){
    matched_row = matcher_desc(filtered_df, dictionary, nth_attempt)
    if(!is.na(matched_row[1,1])){
      if(begin == FALSE){
        temp_row = df[1,] %>% 
          mutate(n_iteration = nth_attempt, base_name = "High") %>% 
          select(base_name, monoisotopic_mass, sum_intensity, apex_rt, n_iteration) 
        return_df = rbind(return_df, temp_row) # include the beginning mass point
      }
      begin = TRUE
      return_df = rbind(return_df, matched_row) #add the found row for return
      temp_df = df %>% 
        filter(monoisotopic_mass <= matched_row[1,2])
      if(nrow(temp_df) > 1){
        return_df = loop_down_match(temp_df, mass_bound, dictionary, return_df, nth_attempt, begin)
      }
    }
  }
  return(return_df) 
}

#loop_up_match function
loop_up_match = function(df, mass_bound, dictionary, return_df, nth_attempt, begin){
  filtered_df = filter_desc(df, mass_bound, version = FALSE) #try to find the first match
  if(nrow(filtered_df) > 1){
    matched_row = matcher_asce(filtered_df, dictionary, nth_attempt)
    if(!is.na(matched_row[1,1])){ 
      if(begin == FALSE){
        temp_row = df[nrow(df),] %>% 
          mutate(n_iteration = nth_attempt, base_name = "High") %>% 
          select(base_name, monoisotopic_mass, sum_intensity, apex_rt, n_iteration) 
        return_df = rbind(return_df, temp_row) # include the beginning mass point
      }
      begin = TRUE
      return_df = rbind(return_df, matched_row) #add the found row for return
      temp_df = df %>% 
        filter(monoisotopic_mass >= matched_row[1,2])
      if(nrow(temp_df) > 1){
        return_df = loop_up_match(temp_df, mass_bound, dictionary, return_df, nth_attempt, begin)
      }
    }
  }
  return(return_df) 
} 