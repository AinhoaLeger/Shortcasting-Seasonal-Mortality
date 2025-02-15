
######################################################
##########                DATA              ##########
######################################################
#####-------------- Directory and packages --------------#####

rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
require(ggplot2)
library(splines)
library(zoo)
library(data.table)
library(readxl)
library(foreign)
library(magic)

# Load function for forecasting
source("fun/ForFun total.R")

#####-------------- Load death series --------------#####

# Eurostat death data
death <- read_xlsx("data/deaths_EU.xlsx")

# Eurostat population data
pop <- readxl::read_xlsx("data/pop_1jan_EU.xlsx")

#####-------------- Prepare the data on deaths --------------#####

death <- death %>%
  filter(country %in% c("Bulgaria","Czechia","Denmark","Germany",
                        "Estonia","Ireland","Greece","Spain","France",
                        "Croatia","Italy","Lithuania","Luxembourg",
                        "Hungary","Netherlands","Austria","Poland",
                        "Portugal","Romania","Slovenia","Finland",
                        "Sweden","Iceland","Norway","Switzerland")) %>%
  pivot_longer(c(3:30), names_to = "year", values_to = "deaths")

death <- death %>%
  filter(TIME != "Total", TIME != "Unknown") %>%
  mutate(month = factor(TIME, 
                       levels=c("January","February","March","April","May",
                                "June","July","August","September","October",
                                "November","December"), 
                       labels=c("January"="01","February"="02","March"="03",
                                "April"="04","May"=05,"June"="06","July"="07",
                                "August"=08,"September"=09,"October"="10",
                                "November"="11","December"="12"))) %>%
  mutate(date = as.Date(as.yearmon(paste(year, month, sep = "-")))) %>%
  select(country, date, deaths) %>%
  arrange(country,date)

#####-------------- Prepare the data on population --------------#####

avg.last.2 <- 
  function(x) if (length(x) < 2) {rep(NA, length(x))} else {
    shift( rollmeanr(x, 2, fill=NA, align="right"), -1)}

exposure <- pop %>%
  pivot_longer(c(2:30), names_to = "year", values_to = "pop") %>%
  # Compute the average death count of the 5 preceding years
  group_by(country) %>%
  arrange(country,year) %>% 
  mutate(exposure = avg.last.2(pop)) %>%
  ungroup() %>%
  na.omit() %>% 
  mutate(exposure = exposure/12)

rates <- (death %>% mutate(year = as.character(year(date)))) %>% 
  left_join(exposure, by=c("country","year")) %>%
  na.omit() %>%
  mutate(rates=deaths/exposure) 

#####-------------- Save death and rates file ------------#####

# save(death, rates, file = "fig/data/series.EU.RData")


######################################################
##########               PLOTS              ##########
######################################################
#####-------------- Plots manuscript --------------#####

db.fig1 <- rates %>% 
  select(-c(pop, exposure, year)) %>%
  filter(country %in% c("Sweden"), date >= "2000-01-01", date <= "2022-06-01") %>% 
  rename(`Death counts` = deaths, `Death rates` = rates) %>%
  pivot_longer(c(3,4), names_to = "type", values_to = "value") 

db.fig2 <- rates %>% 
  filter(country %in% c("Sweden"), date >= "2000-01-01") %>% 
  mutate(date.fit = c(date[1:242],rep(NA,length(date)-242)),
         deaths.fit = c(deaths[1:242],rep(NA,length(deaths)-242)),
         date.fit = c(date[1:242],rep(NA,length(date)-242)),
         rates.fit = c(rates[1:242],rep(NA,length(rates)-242))) %>%
  select(-c(pop, exposure, year, deaths, rates)) %>%
  rename(`Death counts` = deaths.fit, `Death rates` = rates.fit) %>%
  pivot_longer(c(4,5), names_to = "type", values_to = "value")

fig13 <- ggplot() +
  geom_line(data = db.fig1, aes(x=date, y=value), color="red", size=1) +  
  geom_line(data = db.fig2, aes(x=date.fit, y=value), color="grey", size=1) +
  facet_grid(vars(type), scales = "free") + 
  scale_x_date(date_labels="%Y", date_breaks="2 year", expand=c(0.05,0.05)) +  
  labs(title= "Sweden", x = "", y = "") +
  theme_minimal() +
  theme(
    # Axis
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 15),
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

ggsave("fig/data/FigEU 13.png", plot=fig13, width=10, height=8, bg="white")

