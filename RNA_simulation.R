library(dplyr)
library(tidyverse)
library(readxl)
library(writexl)
library(ggrepel)
library(plotly)

  #Load the dictionary, containing the most common 9 nucleotides(AUCG, their methylated form, and D)
dictionary = read_xlsx("Data/dictionary.xlsx") %>% 
  janitor::clean_names()

  #This function generates a completely random RNA sequence with user defined length variable
simulate_RNA = function(length){
  RNA = ""
  for(i in 1:length){
    new_base = dictionary[sample(1 : nrow(dictionary), 1), 1]
    RNA = paste0(RNA, new_base)
  }
  print(RNA)
  return(RNA)
}

  #This function generates the dataset containing the perfect 3' and 5' ladders 
#data based on the input sequence(5') of the corresponding RNA(name_of_RNA), the returned 
#dataset will contain corresponding base name at the fragment end, position, sequence
#of the fragment, ladder type(3' or 5'), and name of the RNA.
generate_ladder = function(input_suquence, name_of_RNA){
  return_df = tibble(
    fragment_sequence = character(),
    base_name = character(),
    position = integer(),
    ladder_type = character(),
    RNA_name = character()
  )
  for(i in 1:(nchar(input_suquence) - 1)){
    return_df[i, 1] = substring(input_suquence, 1, i)
    return_df[i, 2] = substring((return_df[i,1]), nchar(return_df[i,1]), nchar(return_df[i,1]))
    return_df[i, 3] = nchar(return_df[i,1])
    return_df[i, 4] = "5'"
    return_df[i, 5] = name_of_RNA
  }
  temp_n = nrow(return_df)
  reversed_test_RNA_sequence = paste(rev(strsplit(input_suquence, NULL)[[1]]), collapse = "") # reverse the input 5' RNA
  for(i in 1:nchar(reversed_test_RNA_sequence)){
    j = i + temp_n
    return_df[j, 1] = substring(reversed_test_RNA_sequence, 1, i)
    return_df[j, 2] = substring((return_df[j,1]), nchar(return_df[j,1]), nchar(return_df[j,1]))
    return_df[j, 3] = nchar(return_df[j,1])
    return_df[j, 4] = "3'"
    return_df[j, 5] = name_of_RNA
  }
  return(return_df)
}


# This function Generate mass column based on the given data frame(created by generate_ladder())
generate_mass = function(ladder_df){
  ladder_df = mutate(ladder_df, monoisotopic_mass = 0)
  first_line_3 = FALSE
  for(i in 1:nrow(ladder_df)){
    for(j in 1:nrow(dictionary)){
      if(dictionary[j, 1] == ladder_df[i, 2]){
        mass = dictionary[j, 2]
      }
    }
    if(i == 1){
      ladder_df[i, 6] = mass + 18.015
    } else if(!first_line_3 & ladder_df[i, 4] == "3'"){ #reached the first line of 3' 
      first_line_3 = TRUE
      ladder_df[i, 6] = mass - 61.95579
    } else {
      ladder_df[i, 6] = ladder_df[i-1, 6] + mass
    }
  }
  return(ladder_df)
}

# This function generates the perfect fragments data frame for an RNA of user's choice
perfect_fragments = function(input_sequence, name_of_RNA) {
  temp_df = generate_ladder("AAACCGUUACCAUUACUGAG", "RNA_A")
  temp_df_final = generate_mass(temp_df)
  return(temp_df_final)
}

### THE CODE ABOVE IS TO SIMULATE RNA LADDER MASSES
### THE CODE BELOW IS TO TRAIN MODEL TO PREDICT RETERNTION TIME BASED ON MASS

train_df = rbind(read_xlsx("Data/250115/Phe_1FA_10pmol_T01_250114_LSS.xlsx") %>% 
  drop_na() %>% 
  janitor::clean_names(),
  
  read_xlsx("Data/250115/Phe_1FA_15pmol_T01_250114_LSS.xlsx") %>% 
    drop_na() %>% 
    janitor::clean_names(),

  read_xlsx("Data/250115/Phe_1FA_20pmol_T01_250114_LSS.xlsx") %>% 
    drop_na() %>% 
    janitor::clean_names(),

  read_xlsx("Data/250115/Phe_50FA_15pmol_T01_250114_LSS.xlsx") %>% 
    drop_na() %>% 
    janitor::clean_names())

