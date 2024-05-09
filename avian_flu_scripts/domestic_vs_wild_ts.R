## 17/04/2024
## Compare the timeseries for the cases in wild birds to the cases in domestic species. 

rm(list = ls())
library(tidyverse)
library(ggplot2)
library(terra)

domestic_ts <- read.csv("data/flu_data/prepped_data/time_series_domestic_hpai.csv")
#wild_ts <- read.csv("data/flu_data/prepped_data/weekly_counts_hpai_wild_europe.csv")
wild_ts <- read.csv("data/flu_data/prepped_data/weekly_counts_hpai_wild_europe_april2024.csv")

#rename the count column
domestic_ts <- rename(domestic_ts, Domestic = n)
wild_ts <- rename(wild_ts, Wild = n)

#merge with the wild bird data
dom_wild_ts <- left_join(wild_ts, domestic_ts)

dom_wild_long <- pivot_longer(dom_wild_ts, cols = c("Domestic", "Wild"))

dom_wild_long$weekdate <- as.Date(dom_wild_long$weekdate, format = "%Y-%m-%d")

# plot the time series
ggplot(data=dom_wild_long, aes(x= weekdate, y=value, colour=name)) +
  geom_line(linewidth = 0.6) 

## It is hard to see from this plot whether there is a consistent lag because of the length of the x-axis. 
## Look at in in subsections.

par(mar = c(0,0,0,0))
ggplot(data=dom_wild_long[which(dom_wild_long$year < 2008),], aes(x= weekdate, y=value, colour=name)) +
  geom_line(linewidth = 0.6) + 
  labs(y = "Number of weekly cases", x = "Year", size = 18) +
  theme(axis.text.x = element_text(angle=90, margin = margin(t = 0.1, r = 0.2, b = 0.2, l = 0.3, unit = "cm"), 
                                   face = "bold", size = 10, vjust = 0.5),
        axis.text.y = element_text(face = "bold", size = 10),
        plot.margin = margin(1,1,1,1, "cm"), 
        axis.title =  element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "bottom", 
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey"),
        panel.grid.minor = element_line(size = 0.5, linetype = 'dashed',
                                        colour = "grey90"))
#ggsave("plots/domestic_v_wild_pre2008.png")


ggplot(data=dom_wild_long[which(dom_wild_long$year > 2015 & dom_wild_long$year < 2019),], aes(x= weekdate, y=value, colour=name)) +
  geom_line(linewidth = 0.7) + 
  labs(y = "Number of weekly cases", x = "Year", size = 18) +
  theme(axis.text.x = element_text(angle=90, margin = margin(t = 0.1, r = 0.2, b = 0.2, l = 0.3, unit = "cm"), 
                                   face = "bold", size = 10, vjust = 0.5),
        axis.text.y = element_text(face = "bold", size = 10),
        plot.margin = margin(1,1,1,1, "cm"), 
        axis.title =  element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey"),
        panel.grid.minor = element_line(size = 0.5, linetype = 'dashed',
                                        colour = "grey90"))
#ggsave("plots/domestic_v_wild_2015_2019.png")


## And finally for the last part of the data
ggplot(data=dom_wild_long[which(dom_wild_long$year > 2019),], aes(x= weekdate, y=value, colour=name)) +
  geom_line(linewidth = 0.7) + 
  labs(y = "Number of weekly cases", x = "Year", size = 18) +
  theme(axis.text.x = element_text(angle=90, margin = margin(t = 0.1, r = 0.2, b = 0.2, l = 0.3, unit = "cm"), 
                                   face = "bold", size = 10, vjust = 0.5),
        axis.text.y = element_text(face = "bold", size = 10),
        plot.margin = margin(1,1,1,1, "cm"), 
        axis.title =  element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey"),
        panel.grid.minor = element_line(size = 0.5, linetype = 'dashed',
                                        colour = "grey90"))
ggsave("plots/domestic_v_wild_post2020.png")


