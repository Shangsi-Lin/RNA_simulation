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
library(viridis)

plot_mixture = function(df_location){
  temp = read_xlsx(df_location) %>% 
    mutate(log_intensity = log10(sum_intensity))
  plot_ly(temp, 
          x = ~monoisotopic_mass, 
          y = ~apex_rt, 
          z = ~log_intensity, 
          color = ~sgRNA, 
          type = "scatter3d", 
          mode = "markers", 
          marker = list(size = 5)) %>%
    layout(scene = list(
      xaxis = list(title = "Monoisotopic Mass (Da)"),
      yaxis = list(title = "Time (min)"),
      zaxis = list(title = "Log Intensity")
    ))
}