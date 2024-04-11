### Can compare the timeseries for the cases in wild birds to the cases in domestic species. 

domestic_ts <- read.csv("data/flu_data/prepped_data/time_series_domestic_hpai.csv")

#rename the count column
domestic_ts <- rename(domestic_ts, Domestic = n) %>%
  dplyr::select(!date)

#merge with the wild bird data
dom_wild_ts <- left_join(all_dates_subtype_data, domestic_ts)

dom_wild_ts <- rename(dom_wild_ts, Wild = n)

dom_wild_long <- pivot_longer(dom_wild_ts, cols = c("Domestic", "Wild"))

ggplot(data=dom_wild_long,
       aes(x=week_seq, y=value, colour=name)) +
  geom_line() ## output not shown, it's equivalent to the below graph (with a tiny difference in the legend title)


ggplot(dom_wild_ts, aes(week_seq)) + 
  geom_line(aes(y = n, colour = "n")) + 
  geom_line(aes(y = dom_n, colour = "dom_n"))

















trial_df <- hpai

class(trial_df$observation.date)
class(trial_df$week)

trial_df$week_date <- as.Date(trial_df$week)

## we need to aggregate by the time period we want to use

weekly_counts <- trial_df %>%
  count(week) ## has missing weeks as week is a character variable

library(dplyr)
data.frame(date=seq(as.Date("2005-01-01"), as.Date("2023-07-27"), by="day")) %>% 
  mutate(week_num=isoweek(date),year=year(date)) %>%
  group_by(year,week_num) %>% 
  summarise(weekdate=min(date)) -> week_calendar


dates <- merge(trial_df,week_calendar)
# And after you can plot with

count_dates <- dates %>% 
  count(weekdate)

weekly_counts <- left_join(week_calendar, count_dates)
ggplot(data=weekly_counts, aes(x=weekdate, y=n)) +
  geom_bar(stat="identity")


ggplot(dates, aes(x = weekdate)) +
  geom_bar() #------------------------+
  scale_x_date(date_breaks = "1 week", date_labels = "%d-%b")+
  theme(axis.text.x = element_text(angle = 90))
# however, this doesn't seem to work