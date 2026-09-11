library(ggthemes)
library(ggplot2)

# Set the theme of the figures
theme_set(theme_tufte(base_size = 7, base_family = "ArialMT"))
theme_update(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    plot.title = element_text(size = 8),
    axis.title = element_text(size = 8), 
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    strip.text = element_text(size = 8), 
    axis.line = element_line(size = 0.5),
    axis.ticks = element_line(size = 0.5),
    axis.ticks.length = unit(.1, "cm"), 
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
)

# palette
icgc <- c(
 "black"="#1c1f3f",
 "grey"="#d7d7de",
 "grey2"="#73788b",
 "blue"="#3979bb",
 "red"="#df2a55", 
 "red2"="#f5becc",
 "blue2"="#7ca7d2",
 "blue3"="#284d80",
 "white"="white"  
)

nsubclones_colours <- c(
    "0" = "white",
    "1" = "#FFEBEB", 
    "2" = "#FFC2C5", 
    "3" = "#FF8A97", 
    "4" = "#DF2A55"
)

gleason_colours <- c("#799FCB", "#AFC7D0", "#FEC9C9", "#F9665E")
names(gleason_colours) <- c("1", "2", "3", "4+")

trajectory_colours <- c(
  "Ordering 1" = "#8DA0CB", 
  "Ordering 2" = "#AAF0C9", 
  "Ordering 3" = "#C04667"
)

trajectory_colors <- trajectory_colours
