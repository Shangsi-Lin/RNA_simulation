library(tidyverse)
library(dplyr)

blind_seq = function(df, dictionary){
  mass_bound = max(dictionary$mass + 1)
  n_iteration = 0
  output_df = data.frame(matrix(ncol = 5))
  colnames(output_df) = c("base_name", "monoisotopic_mass", "sum_intensity", "apex_rt", "n_iteration")
  output_df = output_df %>% 
    drop_na()
  sequenced_df = df
  exhaustive_level = quantile(df$sum_intensity, 0.001)
  while(nrow(sequenced_df) > 1){
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

output_ladders_found = function(output_df, length, df_location){
  # filter the output data frame from the nested algorithm based on the length user specified
  iterations_filtered = output_df %>% 
    group_by(n_iteration) %>% 
    count() %>% 
    filter(n > length) %>% 
    ungroup() %>% 
    mutate(ladder_number = paste0("ladder", row_number()))
  
  temp = output_df %>% 
    filter(n_iteration %in% iterations_filtered$n_iteration) %>% 
    mutate(ladder_number = "") 
  
  for(i in 1:nrow(temp)){
    for(j in 1:nrow(iterations_filtered)){
      if(temp[i,5] == iterations_filtered[j,1]){
        temp[i,6] = iterations_filtered[j,3]
      }
    }
  }
  
  order_temp = temp %>%
    group_by(n_iteration) %>% 
    arrange(monoisotopic_mass) %>% 
    ungroup() %>% 
    arrange(n_iteration)
    
  
  write_xlsx(order_temp, "Result/blind_sequencing_result_point_version.xlsx")
  
  temp = temp %>% 
    mutate(current_sequence = temp$base_name[1]) 
  for(i in 2:(nrow(temp))){
    if(temp[i, 5] == temp[i-1, 5]){
      temp[i, 7] = paste0(temp[i-1, 7], temp[i,1])
    } else {
      temp[i, 7] = temp[i, 1]
    }
  }
  
  for(i in 2:nrow(temp)){
    if(temp[i, 5] == temp[i-1, 5]){
      temp[i-1, 7] = ""
    }
  }
  
  temp = temp %>% filter(current_sequence != "") %>% 
    rename(sequence = current_sequence, ending_mass = monoisotopic_mass) %>% 
    mutate(sequence = gsub("High", "", sequence)) %>% 
    select(sequence, n_iteration, ladder_number, ending_mass)
  
  for(i in 1:nrow(temp)){
    temp[i,1] = paste(rev(strsplit(temp$sequence[i], NULL)[[1]]), collapse = "")
  }
  
  write_xlsx(temp, "Result/blind_sequencing_result_line_version.xlsx")
}

plot_ladders = function(input_df_location){
  temp = read_xlsx(input_df_location)
  ggplot(temp, aes(x = monoisotopic_mass, y = apex_rt, color = ladder_number)) +
    geom_line(aes(group = ladder_number),color = "gray", size = 1) +
    geom_point(
      alpha=0.9,
      size=3) +
    theme_bw() +
    labs(
      title = paste("All ladder generated for", file_location),
      x = "Monoisotopic Mass(Da)",
      y = "Rentention Time(min)"
    ) +
    theme(plot.title = element_text(
      size = 14,        
      face = "bold",    
      hjust = 0.5       # Horizontal justification (0 = left, 0.5 = center, 1 = right)
    )
    ) +
    geom_text_repel(
      aes(label = base_name),  
      size = 3,               
      color = "black",        
      max.overlaps = 70       
    )
}
