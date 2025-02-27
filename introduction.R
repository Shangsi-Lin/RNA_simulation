library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)
library(readxl)
library(writexl)
library(formattable)
library(plotly)
library(ggrepel) #avoid label stacking when using ggplot
library(viridisLite) #turbo 256 color scale
library(patchwork) #for table presentation

# build the theoretical ladder based on ladder type and sample type
build_theo = function(known_sequence, reference, ladder_5, sample_type){
  return_df = data.frame(matrix(ncol = 2))
  colnames(return_df) = c("base_name", "theoretical_mass")
  if(sample_type == "natural_RNA"){
    if(ladder_5 == TRUE){
      return_df[1,1] = "5'"
      return_df[1,2] = 97.9769
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
    } else {
      return_df[nrow(return_df), 2] = return_df[nrow(return_df), 2] + 79.97107
    }
    return(return_df)
  } else if(sample_type == "synthetic_RNA"){
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
  } else if(sample_type == "peptide"){
    return_df[1,1] = "Water"
    return_df[1,2] = 18.015
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
    return(return_df)
  }
}


set_up = function(tRNA_name, reference_sequence_5, file_location, sample_type, intact_limit, interval_filter){
  assign("sample_name", tRNA_name, envir = .GlobalEnv)
  assign("reference_sequence_5", reference_sequence_5, envir = .GlobalEnv)
  assign("reference_sequence_3", paste(rev(strsplit(reference_sequence_5, NULL)[[1]]), collapse = ""), envir = .GlobalEnv)
  df = read_xlsx(file_location) %>% 
    janitor::clean_names() %>% 
    drop_na() %>% 
    select(monoisotopic_mass, apex_rt, sum_intensity, relative_abundance, number_of_detected_intervals) %>% 
    arrange(desc(monoisotopic_mass)) %>%
    mutate(sum_intensity = as.numeric(sum_intensity), relative_abundance = as.numeric(relative_abundance), number_of_detected_intervals = as.numeric(number_of_detected_intervals)) %>% 
    filter(number_of_detected_intervals >= interval_filter) %>% 
    select(-number_of_detected_intervals)
  assign("df", df, envir = .GlobalEnv)
  assign("intact_df", df %>% filter(monoisotopic_mass >= intact_limit), envir = .GlobalEnv)
  print(paste0("Starting the analysis on ", tRNA_name, " with file location in ", file_location, ", the sample is classified as ", sample_type, " with a intact mass higher than ", intact_limit))
  assign("theo_5", build_theo(reference_sequence_5, dictionary, TRUE, sample_type), envir = .GlobalEnv)
  assign("theo_3", build_theo(reference_sequence_3, dictionary, FALSE, sample_type), envir = .GlobalEnv)
  assign("file_location", file_location, envir = .GlobalEnv)
  assign("intact_weight", theo_5[nrow(theo_5), 2], envir = .GlobalEnv)
}

#do mass matching
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

#this function is to do homology search using your dataset against the theo_mass of your homology mass
homology = function(data_frame, theo_mass, description){
  #create return df
  return_df = data.frame(
    observed_mass = numeric(),
    theoretical_mass = numeric(),
    sum_intensity = numeric(),      
    relative_abundance = numeric(),
    apex_rt = numeric(),
    ppm = numeric(),
    identity = character()
  )
  #do matching according to the given theo_mass
  for(i in 1:nrow(data_frame)){
    if(ppm(theo_mass, data_frame[i, 1])){
      return_df[1,1] = data_frame[i,1]
      return_df[1,2] = theo_mass
      return_df[1,3] = data_frame[i,2]
      return_df[1,4] = data_frame[i,3]
      return_df[1,5] = data_frame[i,4]
      return_df[1,6] = round(abs((theo_mass-data_frame[i,1])/theo_mass*1000000), 3)
      return_df[1,7] = description
      return(return_df)
    }
  }
  return(return_df)
}

homology_search = function(sample_name, intact_mass, intact_datafile, mass_upper_bound, mass_lower_bound){
  assign("H2O", 18.01528, envir = .GlobalEnv)
  assign("K", 39.0983, envir = .GlobalEnv)
  assign("Na", 22.989769, envir = .GlobalEnv)
  assign("H", 1.00784, envir = .GlobalEnv)
  assign("m", 14.026, envir = .GlobalEnv)
  assign("C", 305.0413, envir = .GlobalEnv)
  assign("A", 329.0525, envir = .GlobalEnv)
  print(paste("The table below shows the result of homology search for", sample_name))
  df = read_xlsx(intact_datafile) %>% 
    janitor::clean_names() %>% 
    drop_na() %>% 
    select(monoisotopic_mass, sum_intensity, relative_abundance, apex_rt) %>% 
    mutate(monoisotopic_mass = as.numeric(monoisotopic_mass), sum_intensity = as.numeric(sum_intensity), relative_abundance = as.numeric(relative_abundance), apex_rt = as.numeric(apex_rt)) %>% 
    filter(monoisotopic_mass >= mass_lower_bound & monoisotopic_mass <= mass_upper_bound)
  
  if(nrow(df) == 0){
    print(paste("The function has been stopped because no mass between the lower bound and upper bound has been found in this dataset."))
    return(invisible(NULL))
  }
  
  #search for all potential isoform matches
  original = homology(df, intact_mass, paste0(sample_name, " intact")) #match original fmet
  methyl = homology(df, intact_mass + m, "+m") #match methyl adduct
  dimethyl = homology(df, intact_mass + 2*m, "+2m") #match dimethyl adduct
  trimethyl = homology(df, intact_mass + 3*m, "+3m") #match trimethyl adduct
  demethyl = homology(df, intact_mass - m, "-m") #match demethyl
  double_demethyl = homology(df, intact_mass - m*2, "-2m") #match two demethyl
  triple_demethyl = homology(df, intact_mass - m*3, "-3m") #match three demethyl
  hydride = homology(df, intact_mass + H2O, "+H2O") #match hydride form
  dehydride = homology(df, intact_mass - H2O, "-H2O") #match three dehydride 
  sodium = homology(df, intact_mass + Na - H, "+Na") #match sodium adduct - hydrogen
  potassium = homology(df, intact_mass + K - H, "+K") #match potassium adduct - hydrogen
  potassium_sodium = homology(df, intact_mass + K + Na - H - H, "+NaK") #match sodium&potassium adduct
  deA = homology(df, intact_mass - A, "-A") #the CC isoform
  deCA = homology(df, intact_mass - C - A, "-CA") #the C isoform
  deCCA = homology(df, intact_mass - C - C - A, "-CCA") #the without CCA isoform
  
  #combine all the isoforms found
  isoforms_df = rbind(original, methyl, dimethyl, trimethyl, demethyl, double_demethyl, triple_demethyl, hydride, dehydride, sodium, potassium, potassium_sodium, deA, deCA, deCCA)
  
  df = df %>% 
    mutate(identity = "unknown", theoretical_mass =  NA_real_, ppm =  NA_real_) %>% 
    rename(observed_mass = monoisotopic_mass) %>% 
    select(observed_mass, theoretical_mass, sum_intensity, relative_abundance, apex_rt, ppm, identity)
  
  if(nrow(isoforms_df) != 0){ #when there are matches between theoretical and observed mass
    for(i in 1:nrow(isoforms_df)){
      for(j in 1:nrow(df)){
        if(isoforms_df[i,1] == df[j, 1]){
          df[j,] = isoforms_df[i,]
          break
        }
      }
    }
    
    #add the difference between the observed mass with the intact mass of interest
    df = df %>% 
      mutate(difference_to_intact_mass = round((observed_mass - intact_mass), 3))
    
    formattable(df, align =c("l","c","c","c","c", "c", "c", "c", "r"), 
                list(sum_intensity= color_tile("#DeF7E9", "#71CA97"),
                     relative_abundance= color_tile("#DeF7E9", "#71CA97"),
                     ppm = color_tile("#ff7f7f", "red")))
    
    df = df %>% 
      mutate(identified = ifelse(identity == "unknown", "no", "yes")) %>% 
      mutate(identity = ifelse(identity == "unknown", "", identity))
    ggplot(df, aes(x = observed_mass, y = sum_intensity, color = identified)) +
      geom_point() +
      theme_bw() +
      labs(
        title = paste("Identities Assigned After Homology Search"),
        x = "Observed Mass(Da)",
        y = "Sum Intensity"
      ) +
      theme(plot.title = element_text(
        size = 10,        
        face = "bold",    
        hjust = 0.5       # Horizontal justification (0 = left, 0.5 = center, 1 = right)
      )) +
      geom_text_repel(
        aes(label = identity),             
        vjust = -1,                         
        size = 3,                           
        color = "black")
  } else { #when there are no match
    print("No matching to theoretical mass was found for this dataset")
    df = df %>% 
      select(observed_mass, sum_intensity, relative_abundance, apex_rt, identity) %>% 
      mutate(difference_to_intact_mass = round((observed_mass - intact_mass), 3))
    formattable(df, align =c("l","c","c","c","c", "c", "r"), 
                list(sum_intensity= color_tile("#DeF7E9", "#71CA97"),
                     relative_abundance= color_tile("#DeF7E9", "#71CA97")))
  }
}

build_theo_adduct = function(theo_df, adduct_mass){ #build the theoretical ladder for an adduct ladder
  return_df = theo_df
  for(i in 1:nrow(return_df)){
    return_df[i, 2] = return_df[i, 2] + adduct_mass
  }
  return(return_df)
}


theo_creater = function(df, theo_df, adduct_mass, adduct_name, length, ladder_5){
  adduct_theo_df = build_theo_adduct(theo_df, adduct_mass) 
  return_df = prophet(df, adduct_theo_df) %>% 
    mutate(ladder_type = adduct_name, log_intensity = log10(sum_intensity))
  if(ladder_5){
    return_df = mutate(return_df, position_5 = n_position)
  } else {
    return_df = mutate(return_df, position_5 = length - n_position + 1)
  }
  return(return_df)
}

plot_branching = function(df_name){
  ggplot(df_name, aes(x = monoisotopic_mass, y = apex_rt, color = log_intensity, shape = ladder_type)) +
    geom_point(size = 3) +
    geom_line() +
    scale_color_gradientn(colors = turbo(256)) +
    labs(x = "Monoisotopic Mass(Da)", y = "Retention Time(min)") +
    geom_text_repel(aes(label = base_name), color = "black", vjust = -1, size = 5) +
    theme(
      axis.title = element_text(size = 20), 
      axis.text = element_text(size = 20),   
      plot.title = element_text(size = 16)   
    ) +
    theme_bw() +
    ggtitle("")
}