#### Date: July 25th 2026
#### Purpose: Calculate mean temperature and daily light levels for experiments
#### Project: Seagrass Wasting Disease in Palau

##### Load Packages
library(tidyverse)

##### Load data 
kochs_raw <- read.csv("Data/Hobologgers/exp koch postulate 2025-12-11 16_23_33.csv")
exp1_raw <- read.csv("Data/Hobologgers/Exp-1 2025-04-28 21_37_38 PDT (Data PDT).csv") |> 
  mutate(sensor = "1")
exp2_raw <- read.csv("Data/Hobologgers/Exp-2 2025-04-28 21_32_54 PDT (Data PDT).csv") |> 
  mutate(sensor = "2")

##### Clean data ###############################################################
# Clean Koch's
kochs <- kochs_raw |> 
  rename(Date_time = Date.Time..Palau.Standard.Time., 
         Temp_C = Temperature.....C., 
         Light_lux = Light....lux.) |> 
  select(Date_time, Temp_C, Light_lux) |> 
  mutate(Date_time = as.POSIXct(Date_time, format = "%m/%d/%Y %H:%M:%S"))
  

# Remove first five minutes and last five minutes from data 
last_row <- nrow(kochs)

kochs_df <- kochs[7:(last_row-6),]

# Check that this is correct
head(kochs)
head(kochs_df)

tail(kochs)
tail(kochs_df)

# Clean 2025 Herbivory Experiment
herb_exp <- rbind(exp1_raw, exp2_raw) |> 
  rename(Date_time = Date.Time..PDT., 
         Temp_C = Temperature.....C., 
         Light_lux = Light....lux.) |> 
  select(Date_time, Temp_C, Light_lux, sensor) |> 
  mutate(Date_time = as.POSIXct(Date_time, format = "%m/%d/%Y %H:%M")) |> 
  arrange(Date_time)


###### Plotting ################################################################
## Kochs
# Temperature
# About 25.75 hours worth of data (~1 Day)

ggplot(data = kochs_df, aes(x = Date_time, y = Temp_C)) +
  geom_point() +
  theme_bw()+
  labs(y = "Temperature", x = "Date and time")

# Light
ggplot(data = kochs_df, aes(x = Date_time, y = Light_lux)) +
  geom_point() +
  theme_bw()+
  labs(y = "Light (lux)", x = "Date and time")

## Herbivory Experiment
# Temperature
ggplot(data = herb_exp, aes(x = Date_time, y = Temp_C, color = sensor)) +
  geom_point() +
  theme_bw()+
  labs(y = "Temperature", x = "Date and time") 
# Will need to remove the first and last bits when they were in the air conditioning

# Light
ggplot(data = herb_exp, aes(x = Date_time, y = Light_lux, color = sensor)) +
  geom_point() +
  theme_bw()+
  labs(y = "Light (lux)", x = "Date and time")


# Further cleaning of the herbivory experiment data based on time spent in air conditoning
head(herb_exp, n = 200)
tail(herb_exp, n = 100)
herb_exp_df <- herb_exp[92:17923, ] #the 92nd row if the first time one of the loggers goes up and then back down
# same but opposite for 17,923

# Temp cleaned
ggplot(data = herb_exp_df, aes(x = Date_time, y = Temp_C, color = sensor)) +
  geom_point() +
  theme_bw()+
  labs(y = "Temperature", x = "Date and time") 

# Light cleaned
ggplot(data = herb_exp_df, aes(x = Date_time, y = Light_lux, color = sensor)) +
  geom_point() +
  theme_bw()+
  labs(y = "Light (lux)", x = "Date and time")

##### Average Temp and Daytime light ###########################################
## Koch's 
# Descriptive statistics for temp
mean(kochs_df$Temp_C)
min(kochs_df$Temp_C)
max(kochs_df$Temp_C)
sd(kochs_df$Temp_C)

# Light
sd(kochs_df$Light_lux)
min(kochs_df$Light_lux)
max(kochs_df$Light_lux)

# Average daytime light
kochs_df |> 
  filter(Light_lux > 1) |> 
  summarize(mean_light = mean(Light_lux), 
            sd_light = sd(Light_lux))


## Herbivory experiment
# Descriptive statistics for temp
mean(herb_exp_df$Temp_C)
min(herb_exp_df$Temp_C)
max(herb_exp_df$Temp_C)
sd(herb_exp_df$Temp_C)

# Light
sd(herb_exp_df$Light_lux)
min(herb_exp_df$Light_lux)
max(herb_exp_df$Light_lux)

# Average daytime light
herb_exp_df |> 
  filter(Light_lux > 1) |> 
  summarize(mean_light = mean(Light_lux), 
            sd_light = sd(Light_lux))

