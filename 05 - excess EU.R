
######################################################
##########          EXCESS ON COUNTS        ##########
######################################################
#####-------------- Directory and packages --------------#####

rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load data
load("fig/forecast/forecasts.EU.D.RData")

library(tidyverse)
require(ggplot2)
library(splines)
library(zoo)
library(data.table)
library(readxl)
library(foreign)
library(magic)
library(knitr)
library(kableExtra)
library(gt)
library(magrittr)
library(gridExtra)
library(thematic)


#####-------------- Plot excess death --------------#####

figEU80 <- db %>%
  filter(year>="2020-02-01", year<="2022-06-01") %>%
  filter(model=="Smooth trend model" & 
           country %in% 
           c("Austria","Bulgaria","Croatia","Czechia","Denmark","Estonia",
             "Finland","France","Germany","Greece","Hungary","Iceland",
             "Italy","Luxembourg","Netherlands","Portugal","Romania",
             "Slovenia","Spain","Sweden","Switzerland") |
           model=="Modulation model" & 
           country %in% c("Ireland","Lithuania","Norway","Poland")) %>%
  rename("observed"="deaths", "observed.fit"="deaths.fit") %>%
  mutate(expected = exp(eta), lwr = exp(eta.lo), upr = exp(eta.up)) %>%
  # Excess death by month
  mutate(excess = observed-expected,
         excess.lwr = observed-upr,
         excess.upr = observed-lwr) %>%
  select(country,year,observed,observed.fit,expected,lwr,upr) %>%
  
  ggplot(aes(x = year)) +
  # Differences
  stat_difference(aes(ymin=expected, ymax=observed), size=0.5, alpha=0.3) +
  # facet_grid(vars(country), scales = "free") +
  facet_wrap(vars(country), ncol = 3, scales = "free") +
  # Prediction intervals
  geom_smooth(aes(x=year, y=expected, ymin=lwr, ymax=upr), 
              stat='identity', linetype="solid", color="white", 
              size = 0.5, fill="white", alpha = 1) +
  geom_smooth(aes(x=year, y=expected, ymin=lwr, ymax=upr), 
              stat='identity', linetype="dashed", fill="lightgray", 
              size = 0.5, color="black", alpha = 1) +
  # Observed
  geom_line(aes(y = observed), size=0.5) +
  scale_x_date(date_labels="%Y", expand=c(0.05,0.05)) +
  scale_fill_manual(values = c(colorspace::lighten("#C32E5A"),
                               colorspace::lighten("#3D85F7"),"grey60"),
                    labels = c("excess", "deficit", "same")) +
  labs(title = "", x = "", y = "Death count") +
  guides(color = guide_legend(order = 1), fill = guide_legend(order = 2)) +
  theme_minimal() +
  theme(
    # Axis
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 13, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 13),
    # Facet
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    # Legend
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    legend.position="none", # "bottom"
    legend.key.width=unit(2.5, "cm"),
    # Title
    plot.title = element_text(hjust = 0.5, size = 15))


ggsave("fig/excess/FigEU 80.png", plot=figEU80, width=12, height=17, bg="white")


#####-------------- Compute excess death in periods --------------#####

excess.EU <- db %>%
  filter(year>="2020-02-01", year<="2022-06-01") %>%
  filter(model=="Smooth trend model" & 
         country %in% c("Austria","Bulgaria","Croatia","Czechia","Denmark",
                        "Estonia","Finland","France","Germany","Greece",
                        "Hungary","Iceland","Italy","Luxembourg",
                        "Netherlands","Portugal","Romania","Slovenia",
                        "Spain","Sweden","Switzerland") |
           model=="Modulation model" & 
         country %in% c("Ireland","Lithuania","Norway","Poland")) %>%
  # Transform back to original scale
  rename("observed"="deaths") %>%
  mutate(expected = exp(eta),
         lwr = exp(eta.lo),
         upr = exp(eta.up)) %>%
  # Excess death by month
  mutate(excess = observed-expected,
         excess.lwr = observed-upr,
         excess.upr = observed-lwr) %>%
  select(country,year,excess,excess.lwr,excess.upr) %>%
  # Excess death by pandemic periods
  group_by(country) %>%
  summarize(
    # For the table
    excess_2020 = paste0(
      comma(round(sum(excess[year>="2020-02-01" & year<="2020-06-01"], na.rm=TRUE)), accuracy=1), " (",
      comma(round(sum(excess.lwr[year>="2020-02-01" & year<="2020-06-01"], na.rm=TRUE)), accuracy=1), " ; ",
      comma(round(sum(excess.upr[year>="2020-02-01" & year<="2020-06-01"], na.rm=TRUE)), accuracy=1), ")"),
    excess_2021 = paste0(
      comma(round(sum(excess[year>="2020-07-01" & year<="2021-06-01"], na.rm=TRUE)), accuracy=1), " (",
      comma(round(sum(excess.lwr[year>="2020-07-01" & year<="2021-06-01"], na.rm=TRUE)), accuracy=1), " ; ",
      comma(round(sum(excess.upr[year>="2020-07-01" & year<="2021-06-01"], na.rm=TRUE)), accuracy=1), ")"),
    excess_2122 = paste0(
      comma(round(sum(excess[year>="2021-07-01" & year<="2022-06-01"], na.rm=TRUE)), accuracy=1), " (",
      comma(round(sum(excess.lwr[year>="2021-07-01" & year<="2022-06-01"], na.rm=TRUE)), accuracy=1), " ; ",
      comma(round(sum(excess.upr[year>="2021-07-01" & year<="2022-06-01"], na.rm=TRUE)), accuracy=1), ")"),
    # For the plot
    ex_2020 = sum(excess[year>="2020-02-01" & year<="2020-06-01"], na.rm=TRUE),
    ex_2020_lw = sum(excess.lwr[year>="2020-02-01" & year<="2020-06-01"], na.rm=TRUE),
    ex_2020_up = sum(excess.upr[year>="2020-02-01" & year<="2020-06-01"], na.rm=TRUE),
    ex_2021 = sum(excess[year>="2020-07-01" & year<="2021-06-01"], na.rm=TRUE),
    ex_2021_lw = sum(excess.lwr[year>="2020-07-01" & year<="2021-06-01"], na.rm=TRUE),
    ex_2021_up = sum(excess.upr[year>="2020-07-01" & year<="2021-06-01"], na.rm=TRUE),
    ex_2122 = sum(excess[year>="2021-07-01" & year<="2022-06-01"], na.rm=TRUE),
    ex_2122_lw = sum(excess.lwr[year>="2021-07-01" & year<="2022-06-01"], na.rm=TRUE),
    ex_2122_up = sum(excess.upr[year>="2021-07-01" & year<="2022-06-01"], na.rm=TRUE)
  ) %>%
  ungroup() 

#####-------------- Table excess death in periods --------------#####

# Remember to change [t] in [c]

caption= "All-cause excess death in Europe during three pandemic periods (the first COVID-19 wave, and the first and second pandemic year after the shock)."
excess.EU[,c(1:4)] %>%
  mutate_all(linebreak, align = 'l') %>%
  kbl(longtable=T, booktabs=T, caption=caption, escape=F,
      col.names = c("Country","March to June 2020",
                    "July 2020 to June 2021","July 2021 to June 2022"),
      align="l", format="latex", label="excess_counts_EU") %>%
  kable_styling(position="center", font_size=10) %>%
  add_header_above(c(" "=1,"Excess death (prediction  intervals)"=3))

#####-------------- Plot excess death in periods --------------#####

plot1 <- excess.EU %>%
  ggplot(aes(x=country, y=ex_2020, width=0.65, fill=factor(country))) +
  geom_point(aes(fct_reorder(country,ex_2020), colour=ex_2020), stat="identity", size=2) +
  geom_errorbar(aes(ymin=ex_2020_lw, ymax=ex_2020_up, colour=ex_2020), width=.5) +
  scale_colour_gradient(low = "blue", high = "red") +
  coord_flip() +
  labs(title = "March to June 2020", x="", y="") +
  theme_minimal() +
  theme(
    # Axis 
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 15),
    # Legend
    legend.position = "none",
    legend.title = element_blank(),
    legend.text=element_text(size=15))

plot2 <- excess.EU %>%
  ggplot(aes(x=country, y=ex_2021, width=0.65, fill=factor(country))) +
  geom_point(aes(fct_reorder(country,ex_2021), colour=ex_2021), stat="identity", size=2) +
  geom_errorbar(aes(ymin=ex_2021_lw, ymax=ex_2021_up, colour=ex_2021), width=.5) +
  scale_colour_gradient(low = "blue", high = "red") +
  coord_flip() +
  labs(title="July 2020 to June 2021", x="", y="Excess death (counts)") +
  theme_minimal() +
  theme(
    # Axis 
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 15),
    # Legend
    legend.position = "none",
    legend.title = element_blank(),
    legend.text=element_text(size=15))

plot3 <- excess.EU %>%
  ggplot(aes(x=country, y=ex_2122, width=0.65, fill=factor(country))) +
  geom_point(aes(fct_reorder(country,ex_2122), colour=ex_2122), stat="identity", size=2) +
  geom_errorbar(aes(ymin=ex_2122_lw, ymax=ex_2122_up, colour=ex_2122), width=.5) +
  scale_colour_gradient(low = "blue", high = "red") +
  coord_flip() +
  labs(title = "July 2021 to June 2022", x="", y="") +
  theme_minimal() +
  theme(
    # Axis 
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 15),
    # Legend
    legend.position = "none",
    legend.title = element_blank(),
    legend.text=element_text(size=15))


plot.final <- grid.arrange(plot1, plot2, plot3, nrow = 1)
ggsave("fig/excess/FigEU 78.png", plot=plot.final, width=12, height=9, bg="white")



######################################################
##########          EXCESS ON RATES         ##########
######################################################
#####-------------- Directory and packages --------------#####

rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load data
load("fig/forecast/forecasts.EU.R.RData")

# devtools::install_github("teunbrand/ggh4x"force = TRUE)
library(tidyverse)
require(ggplot2)
library(splines)
library(zoo)
library(magic)
library(geofacet)
library(ggh4x)
library(colorspace)
library(knitr)
library(kableExtra)
library(gt)
library(magrittr)
library(gridExtra)
library(thematic)


#####-------------- Plot excess rates --------------#####

figEU81 <- db %>%
  filter(year>="2020-02-01", year<="2022-06-01") %>%
  filter(model=="Smooth trend model" & 
           country %in% 
           c("Austria","Bulgaria","Croatia","Czechia","Denmark","Estonia",
             "Finland","France","Germany","Greece","Hungary","Iceland",
             "Italy","Luxembourg","Netherlands","Portugal","Romania",
             "Slovenia","Spain","Sweden","Switzerland") |
           model=="Modulation model" & 
           country %in% c("Ireland","Lithuania","Norway","Poland")) %>%
  rename("observed"="rates", "observed.fit"="rates.fit") %>%
  mutate(expected = exp(eta), lwr = exp(eta.lo), upr = exp(eta.up)) %>%
  # Excess death by month
  mutate(excess = observed-expected,
         excess.lwr = observed-upr,
         excess.upr = observed-lwr) %>%
  select(country,year,observed,observed.fit,expected,lwr,upr) %>%
  
  ggplot(aes(x = year)) +
  # Differences
  stat_difference(aes(ymin=expected, ymax=observed), size=0.5, alpha=0.3) +
  # facet_grid(vars(country), scales = "free") +
  facet_wrap(vars(country), ncol = 3, scales = "free") +
  # Prediction intervals
  geom_smooth(aes(x=year, y=expected, ymin=lwr, ymax=upr), 
              stat='identity', linetype="solid", color="white", 
              size = 0.5, fill="white", alpha = 1) +
  geom_smooth(aes(x=year, y=expected, ymin=lwr, ymax=upr), 
              stat='identity', linetype="dashed", fill="lightgray", 
              size = 0.5, color="black", alpha = 1) +
  # Observed
  geom_line(aes(y = observed), size=0.5) +
  scale_x_date(date_labels="%Y", expand=c(0.05,0.05)) +
  scale_fill_manual(values = c(colorspace::lighten("#C32E5A"),
                               colorspace::lighten("#3D85F7"),"grey60"),
                    labels = c("excess", "deficit", "same")) +
  labs(title = "", x = "", y = "Death rates") +
  guides(color = guide_legend(order = 1), fill = guide_legend(order = 2)) +
  theme_minimal() +
  theme(
    # Axis
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 13, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 13),
    # Facet
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    # Legend
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    legend.position="none", # "bottom"
    legend.key.width=unit(2.5, "cm"),
    # Title
    plot.title = element_text(hjust = 0.5, size = 15))


ggsave("fig/excess/FigEU 81.png", plot=figEU81, width=12, height=17, bg="white")


#####-------------- Compute excess rates in periods --------------#####

excess.EU <- db %>%
  filter(year>="2020-02-01") %>%
  filter(model=="Smooth trend model" & 
           country %in% c("Austria","Bulgaria","Croatia","Czechia","Denmark",
                          "Estonia","Finland","France","Germany","Greece",
                          "Hungary","Iceland","Italy","Luxembourg",
                          "Netherlands","Portugal","Romania","Slovenia",
                          "Spain","Sweden","Switzerland") |
           model=="Modulation model" & 
           country %in% c("Ireland","Lithuania","Norway","Poland")) %>%
  # Transform back to original scale
  rename("observed"="rates") %>%
  mutate(expected = exp(eta),
         lwr = exp(eta.lo),
         upr = exp(eta.up)) %>%
  # Excess death by month
  mutate(excess = observed-expected,
         excess.lwr = observed-upr,
         excess.upr = observed-lwr) %>%
  select(country,year,excess,excess.lwr,excess.upr) %>%
  # Excess death by pandemic periods
  group_by(country) %>%
  summarize(
    excess_2020 = paste0(
      round(sum(excess[year>="2020-02-01" & year<="2020-06-01"], na.rm=TRUE)*1000,2), " (",
      round(sum(excess.lwr[year>="2020-02-01" & year<="2020-06-01"], na.rm=TRUE)*1000,2), " ; ",
      round(sum(excess.upr[year>="2020-02-01" & year<="2020-06-01"], na.rm=TRUE)*1000,2), ")"),
    excess_2021 = paste0(
      round(sum(excess[year>="2020-07-01" & year<="2021-06-01"], na.rm=TRUE)*1000,2), " (",
      round(sum(excess.lwr[year>="2020-07-01" & year<="2021-06-01"], na.rm=TRUE)*1000,2), " ; ",
      round(sum(excess.upr[year>="2020-07-01" & year<="2021-06-01"], na.rm=TRUE)*1000,2), ")"),
    excess_2122 = paste0(
      round(sum(excess[year>="2021-07-01" & year<="2022-06-01"], na.rm=TRUE)*1000,2), " (",
      round(sum(excess.lwr[year>="2021-07-01" & year<="2022-06-01"], na.rm=TRUE)*1000,2), " ; ",
      round(sum(excess.upr[year>="2021-07-01" & year<="2022-06-01"], na.rm=TRUE)*1000,2), ")"),
    # For the plot
    ex_2020 = sum(excess[year>="2020-02-01" & year<="2020-06-01"], na.rm=TRUE)*1000,
    ex_2020_lw = sum(excess.lwr[year>="2020-02-01" & year<="2020-06-01"], na.rm=TRUE)*1000,
    ex_2020_up = sum(excess.upr[year>="2020-02-01" & year<="2020-06-01"], na.rm=TRUE)*1000,
    ex_2021 = sum(excess[year>="2020-07-01" & year<="2021-06-01"], na.rm=TRUE)*1000,
    ex_2021_lw = sum(excess.lwr[year>="2020-07-01" & year<="2021-06-01"], na.rm=TRUE)*1000,
    ex_2021_up = sum(excess.upr[year>="2020-07-01" & year<="2021-06-01"], na.rm=TRUE)*1000,
    ex_2122 = sum(excess[year>="2021-07-01" & year<="2022-06-01"], na.rm=TRUE)*1000,
    ex_2122_lw = sum(excess.lwr[year>="2021-07-01" & year<="2022-06-01"], na.rm=TRUE)*1000,
    ex_2122_up = sum(excess.upr[year>="2021-07-01" & year<="2022-06-01"], na.rm=TRUE)*1000
  ) %>%
  ungroup() 


#####-------------- Table excess rates in periods --------------#####

caption= "Excess CDRs (x 1000) in 25 European countries. The excess is divided in three pandemic periods: the first COVID-19 wave, and the first and second pandemic year after the shock."
excess.EU[,c(1:4)] %>%
  mutate_all(linebreak, align = 'l') %>%
  kbl(longtable=T, booktabs=T, caption=caption, escape=F,
      col.names = c("Country","March to June 2020",
                    "July 2020 to June 2021","July 2021 to June 2022"),
      align="l", format="latex", label="excess_rates_EU") %>%
  kable_styling(position="center", font_size=10) %>%
  add_header_above(c(" "=1,"Excess death (prediction  intervals)"=3))


#####-------------- Plot excess rates in periods --------------#####

plot1 <- excess.EU %>%
  ggplot(aes(x=country, y=ex_2020, width=0.65, fill=factor(country))) +
  geom_point(aes(fct_reorder(country,ex_2020), colour=ex_2020), stat="identity", size=2) +
  geom_errorbar(aes(ymin=ex_2020_lw, ymax=ex_2020_up, colour=ex_2020), width=.5) +
  scale_colour_gradient(low = "blue", high = "red", limits = c(-5,62)) +
  coord_flip() +
  labs(title = "March to June 2020", x="", y="") +
  theme_minimal() +
  theme(
    # Axis 
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 15),
    # Legend
    legend.position = "none",
    legend.title = element_blank(),
    legend.text=element_text(size=15))

plot2 <- excess.EU %>%
  ggplot(aes(x=country, y=ex_2021, width=0.65, fill=factor(country))) +
  geom_point(aes(fct_reorder(country,ex_2021), colour=ex_2021), stat="identity", size=2) +
  geom_errorbar(aes(ymin=ex_2021_lw, ymax=ex_2021_up, colour=ex_2021), width=.5) +
  scale_colour_gradient(low = "blue", high = "red", limits = c(-5,62)) +
  coord_flip() +
  labs(title="July 2020 to June 2021", x="", y="Excess death (rates)") +
  theme_minimal() +
  theme(
    # Axis 
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 15),
    # Legend
    legend.position = "none",
    legend.title = element_blank(),
    legend.text=element_text(size=15))

plot3 <- excess.EU %>%
  ggplot(aes(x=country, y=ex_2122, width=0.65, fill=factor(country))) +
  geom_point(aes(fct_reorder(country,ex_2122), colour=ex_2122), stat="identity", size=2) +
  geom_errorbar(aes(ymin=ex_2122_lw, ymax=ex_2122_up, colour=ex_2122), width=.5) +
  scale_colour_gradient(low = "blue", high = "red", limits = c(-5,62)) +
  coord_flip() +
  labs(title = "July 2021 to June 2022", x="", y="") +
  theme_minimal() +
  theme(
    # Axis 
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 15),
    # Legend
    legend.position = "none",
    legend.title = element_blank(),
    legend.text=element_text(size=15))


plot.final <- grid.arrange(plot1, plot2, plot3, nrow = 1)
ggsave("fig/excess/FigEU 79.png", plot=plot.final, width=12, height=9, bg="white")

