library(dplyr)
library(tidyverse)
library(readxl)
library(writexl)

load_dictionary = function(row_chosen){ # the dictionary file name has to be dictionary.xlsx
  print("Loading dictionary.xlsx")
  read_xlsx("dictionary.xlsx", n_max = row_chosen)
}

dictionary = load_dictionary(4)

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