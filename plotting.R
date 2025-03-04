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

#This function plots the depth of the sequencing result, the x-axis will be the position of nucleotides, the y-axis will be
#the number of depth sequenced at that particular position.
#first variable: the df you want to draw based on, this should be in the syntax of output by function theo_creater.
#second variable: the positions that you don't want to be included in this depth plot.(Ex. the intact position)
plot_depth = function(data_frame, positions_to_eliminate){
  temp = data_frame %>%
    add_count(position_5, name = "count_appearances") %>% 
    distinct(position_5, .keep_all = TRUE) %>% 
    select(position_5, count_appearances) %>% 
    arrange(position_5) %>% 
    filter(!position_5 %in% positions_to_eliminate)
  
  temp = temp %>%
    rename(depth = count_appearances)
  
  # Draw the bar plot
  ggplot(temp, aes(x = position_5, y = depth)) +
    geom_bar(stat = "identity", fill = "green") + 
    labs(
      x = "Position",
      y = "Depth",
      title = "RNA-B Sequencing Depth at Each Position"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    )+
    scale_y_continuous(
      breaks = c(2, 4, 6, 8, 10, 12, 14)
    ) +
    scale_x_continuous(
      breaks = c(1, 5, 10, 15, 20),
      limits = c(1, NA)
    )
}

#This function is for plotting the homology search plot, user can later use lines to connect the dots in the same group to
#show their relationship.
#First variable: the data frame that this function is going to draw upon, this data frame needs to include columns named 
#monoisotopic_mass, sum_intensity, group(which mass group, Ex.M1), and label(Ex.M1IF1).
#Second variable: the colors that you want the groups to be.(Ex.c("M1" = "red", "M2" = "blue", "M3" = "green"))
plot_homology = function(data_frame, group_color){
  ggplot(data_frame, aes(x = as.numeric(monoisotopic_mass), y = as.numeric(sum_intensity), color = group)) +
    geom_point(size = 6) +  # Plot points with size 4
    geom_text(aes(label = label), vjust = -1, size = 4, color = "black") +  # Add labels above points
    labs(
      x = "Monoisotopic Mass",
      y = "Sum Intensity",
      color = "Group"
    ) +
    scale_color_manual(values = group_color) +  # Define colors for groups
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 12)
    )
}