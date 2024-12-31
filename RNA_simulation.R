library(dplyr)
library(tidyverse)
library(readxl)
library(writexl)
library(ggrepel)
library(plotly)

load_dictionary = function(row_chosen){ # the dictionary file name has to be dictionary.xlsx
  print("Loading dictionary.xlsx")
  read_xlsx("dictionary.xlsx", n_max = row_chosen)
}

dictionary = load_dictionary(11)

simulate_RNA = function(length){ #This function generates a completely random RNA sequence with user defined length
  RNA = ""
  for(i in 1:length){
    new_base = dictionary[sample(1 : nrow(dictionary), 1), 1]
    RNA = paste0(RNA, new_base)
  }
  print(RNA)
  return(RNA)
}

perfect_ladder = function(RNA_sequence, RNA_number, intensity, standard_dev){ #This function produces perfect ladder dataset, both 3' and 5' for 
                                                     #the assigned RNA sequence with the assigned RNA number, intensity, and standard deviation for intensity
  perfect_3ladder = data.frame(
    base_name = character(),  
    monoisotopic_mass = numeric(),
    apex_rt = numeric(),
    sum_intensity = numeric(),
    number = numeric()
  )
  perfect_5ladder = data.frame(
    base_name = character(),  
    monoisotopic_mass = numeric(),
    apex_rt = numeric(),
    sum_intensity = numeric(),
    number = numeric()
  )
  
  reversed_RNA = paste(rev(strsplit(RNA_sequence, NULL)[[1]]), collapse = "") # reverse the input 5' RNA
  for(i in 1:nchar(reversed_RNA)){ #building the dataset containing the perfect 3' ladder
    perfect_3ladder[i,1] = substr(reversed_RNA, i, i)
    for(j in 1:nrow(dictionary)){
      if(substr(reversed_RNA, i, i) == dictionary[j, 1]){
        mass = dictionary[j, 2]
      }
    }
    if(i == 1){
      perfect_3ladder[i, 2] = mass - 61.95579
    } else {
      perfect_3ladder[i, 2] = perfect_3ladder[i-1, 2] + mass
    }
    perfect_3ladder[i, 3] = runif(1, min = i - 0.5, max = i + 0.5) # uniformly assign retention time
    perfect_3ladder[i, 4] = rnorm(1, mean = intensity, sd = standard_dev) # normally assign sum intensity
    perfect_3ladder[i, 5] = RNA_number
  }
  
  for(i in 1:(nchar(RNA_sequence) - 1)){ #building the dataset containing the perfect 5' ladder
    perfect_5ladder[i,1] = substr(RNA_sequence, i, i)
    for(j in 1:nrow(dictionary)){
      if(substr(RNA_sequence, i, i) == dictionary[j, 1]){
        mass = dictionary[j, 2]
      }
    }
    if(i == 1){
      perfect_5ladder[i, 2] = mass + 18.015
    } else {
      perfect_5ladder[i, 2] = perfect_5ladder[i-1, 2] + mass
    }
    perfect_5ladder[i, 3] = runif(1, min = i - 0.5, max = i + 0.5) # uniformly assign retention time
    perfect_5ladder[i, 4] = rnorm(1, mean = intensity, sd = standard_dev) # normally assign sum intensity
    perfect_5ladder[i, 5] = RNA_number
  }
  return(rbind(perfect_5ladder, perfect_3ladder))
}

# get the maximum mass in the dictionary for future processing.
mass_bound = max(dictionary$mass + 1)

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

ppm = function(observed, theo){
  if(abs((observed - theo) / theo * 10^6) > 1) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

blind_sequencing = function(df){
  n_iteration = 0
  output_df = data.frame(matrix(ncol = 5))
  colnames(output_df) = c("base_name", "monoisotopic_mass", "sum_intensity", "apex_rt", "n_iteration")
  output_df = output_df %>% 
    drop_na()
  exhaustive_level = 0
  sequenced_df = df %>% 
    select(monoisotopic_mass, apex_rt, sum_intensity, number) %>% 
    arrange(desc(monoisotopic_mass))
  while(nrow(sequenced_df) > 0){
    n_iteration = n_iteration + 1
    highest_intensity_mass = sequenced_df %>%
      filter(sum_intensity == max(sum_intensity)) %>% 
      filter(monoisotopic_mass == max(monoisotopic_mass)) #in case if there are two data points with same level of intensity
    if(highest_intensity_mass$sum_intensity < exhaustive_level) { #exhaustion reached
      break
    }
    loop_down_df = filter(sequenced_df, monoisotopic_mass <= highest_intensity_mass$monoisotopic_mass)
    loop_up_df = filter(sequenced_df, monoisotopic_mass >= highest_intensity_mass$monoisotopic_mass)
    temp_df = data.frame(matrix(ncol = 5))
    colnames(temp_df) = c("base_name", "monoisotopic_mass", "sum_intensity", "apex_rt", "n_iteration")
    temp_df = temp_df %>% 
      drop_na()
    temp_df = loop_up_match(loop_up_df, mass_bound, dictionary, temp_df, n_iteration, begin = FALSE) %>% 
      arrange(desc(monoisotopic_mass))
    output_df = rbind(output_df, temp_df)
    output_df = loop_down_match(loop_down_df, mass_bound, dictionary, output_df, n_iteration, begin = FALSE)
    suppressMessages({ #suppress the messages
      sequenced_df = anti_join(sequenced_df, output_df, by = "monoisotopic_mass")
      sequenced_df = anti_join(sequenced_df, highest_intensity_mass)
    })
  }
  
  output_df = output_df %>% #make sure every row is unique, remove the replication of two high mass points in one iteration
    distinct()
  
  return(output_df)
}
