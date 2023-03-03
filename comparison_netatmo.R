####
# Script to evaluate palm results with netatmo data
####

# Set the working directory
setwd("your working directory")

# Load required packages
library(ggplot2)
library(lubridate)
library(reshape2)
library(hydroGOF)
library(data.table)
library(ggpubr)
library(dplyr)

### Data preparation and filtering ----
# Load the data
palm <- read.csv("palm_netatmo_parent.csv")
palm_child <- read.csv("palm_netatmo_child.csv")
netatmo <- read.csv("netatmo_city.csv")

# Convert the time columns to a POSIXct format
palm$time <- as.POSIXct(palm$time, format = "%Y-%m-%d %H:%M:%OS")
palm_child$time <- as.POSIXct(palm_child$time, format = "%Y-%m-%d %H:%M:%OS")
netatmo$date <- as.POSIXct(netatmo$date, format = "%Y/%m/%d %H:%M:%S")

# Round the time columns of the PALM data to the nearest hour
palm$time <- round_date(palm$time, unit = "hour")
palm_child$time <- round_date(palm_child$time, unit = "hour")

# Rename the first column of netatmo data
colnames(netatmo)[1] <- c("station_id")

# Rename the third column of the palm datasets
colnames(palm)[3] <- c("date")
colnames(palm_child)[3] <- c("date")

# Merge the palm and netatmo data, keeping all rows from the palm data and only matching rows from the netatmo data
compare <- merge(palm, netatmo, by = c("date","station_id"), all.x = T)

# Select only certain columns from the merged data
compare <- compare[c(1,2,4,6,7,15)]

# Repeat the above steps for the child domain
compare_child <- merge(palm_child, netatmo, by = c("date","station_id"), all.x = T)
compare_child <- compare_child[c(1,2,4,6,7,15)]

# Remove the last timestep from the data
end <- which(compare$date == as.POSIXct("2020-08-09 06:00:00"))
compare <- compare[-end,]

end_child <- which(compare_child$date == as.POSIXct("2020-08-09 06:00:00"))
compare_child <- compare_child[-end_child,]

# Define the number of time steps (n_timestep), here 29
n_timestep <- 29

# Count the number of NA values in the Netatmo data of for the parent domain per station
count_na <- aggregate(cbind(sumNA = is.na(temp))~station_id, compare, sum)

# Calculate the percentage of non missing values for each station
count_na$perc <- 1-(count_na$sumNA/n_timestep)

# Repeat the above steps for the child domain
count_na_c <- aggregate(cbind(sumNA = is.na(temp))~station_id, compare_child, sum)
count_na_c$perc <- 1-(count_na_c$sumNA/n_timestep)

# Identify stations in the parent domain with less than 80% of valid data
missing_data <- which(count_na$perc < 0.8)

# Extract the station IDs of the stations identified in the previous step
missing_stations <- count_na$station_id[missing_data]

# Repeat the above steps for the child domain
missing_data_c <- which(count_na_c$perc < 0.8)
missing_stations_c <- count_na_c$station_id[missing_data_c]

# Identify the rows in the parent domain that correspond to the stations identified as having less than 80% of valid data
missing_index = which(compare$station_id %in% missing_stations)

# Remove the identified rows
compare <- compare[-missing_index,]

# Repeat the above two steps for the "compare_child" data frame
missing_index_c = which(compare_child$station_id %in% missing_stations_c)
compare_child <- compare_child[-missing_index_c,]

# Write the data to csv files
write.csv(compare, "city_palm_netatmo.csv", row.names = F)
write.csv(compare_child, "city_palm_netatmo_child.csv", row.names = F)

### Plots ----
## Boxplot time series
# Melt the "compare" data frame to prepare it for plotting
compare_melt <- reshape2::melt(compare, measure.vars = c("palm_ta","temp"), value.name = "temperature")

# Rename the variable columns in the melted data frame
levels(compare_melt$variable) <- c("PALM","Netatmo")

# Repeat the same process for the "compare_child" data frame
compare_child_melt <- reshape2::melt(compare_child, measure.vars = c("palm_ta","temp"), value.name = "temperature")
levels(compare_child_melt$variable) <- c("PALM","Netatmo")

# Convert the date column to a factor variable
date_factor = as.factor(compare_melt$date)

# Define the date breaks and labels for the x-axis
date_breaks = c("2020-08-08 04:00:00","2020-08-08 10:00:00","2020-08-08 16:00:00","2020-08-08 22:00:00",
                "2020-08-09 04:00:00")
date_labels = c("04:00","10:00","16:00","22:00","04:00")

# Create the first box plot (parent)
p1 <- ggplot(compare_melt, aes(x = as.factor(date), y = temperature, fill = variable)) +
  # Add a box plot layer
  geom_boxplot(outlier.shape = NA) +
  # Label the x-axis
  xlab("time") +
  # Label the y-axis
  ylab("air temperature [°C]") +
  # Customise the x-axis
  scale_x_discrete(breaks=factor(date_breaks), labels = date_labels) +
  # Customise the y-axis
  scale_y_continuous(breaks = seq(15,45,5), limits = c(17,40)) +
  # Apply a black and white theme
  theme_bw() +
  # Customise the theme
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size=24))

# Create the second box plot (child)
p2 <- ggplot(compare_child_melt, aes(x = as.factor(date), y = temperature, fill = variable)) +
  geom_boxplot(outlier.shape = NA) +
  xlab("time") +
  ylab("air temperature [°C]") +
  scale_x_discrete(breaks=factor(date_breaks), labels = date_labels) +
  scale_y_continuous(breaks = seq(15,45,5), limits = c(17,40))
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size=24))

# Save the plot as a PNG image file
png("boxplot_ts_city.png", width = 1100, height = 600)
# Arrange the two plots (p1 and p2) side by side
# Include a common legend for both plots
ggarrange(p1, p2, ncol = 2, widths = c(2.1,2), align="h", common.legend = T, legend = "bottom")
# Close the PNG device
dev.off()
  
## Timeseries of individual stations
# Save the plot as a PNG image file
png("ts_palm_netatmo_city_parent.png", width = 1200, height = 1800)
# Create the plot
ggplot(compare, aes(x=date))+
  geom_line(aes(y=palm_ta, colour = "PALM"))+
  geom_point(aes(y=temp, colour = "Netatmo"))+
  # Manually set the colour
  scale_colour_manual(values = c("PALM"="#f75f55","Netatmo"="#00adb3"))+
  # Create subplots based on station_id with 10 rows
  facet_wrap(~station_id, nrow = 10)+
  # Label the x axis
  xlab("date")+
  # Label the y axis
  ylab("air temperature [°C]")+
  # Customise the x-axis
  scale_x_datetime(breaks = "4 hour", date_labels = "%H:%M")+
  # Customise the y-axis
  scale_y_continuous(breaks = seq(15,45,5))+
  # Apply a black and white theme
  theme_bw()+
  # Customise the theme
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title = element_text(size = 28),
        axis.title.y = element_text(margin = margin(r=10)),
        axis.text = element_text(size = 24),
        axis.text.x = element_blank(),
        strip.text = element_text(size=18),
        legend.text = element_text(size=28))+
  # Set the size of the legend points
  guides(colour=guide_legend(override.aes = list(size=2)))
# Close the PNG device
dev.off()

# Repeat the above steps for the child domain
png("ts_palm_netatmo_city_child.png", width = 1200, height = 1800)
ggplot(compare_child, aes(x=date))+
  geom_line(aes(y=palm_ta, colour = "PALM"))+
  geom_point(aes(y=temp, colour = "Netatmo"))+
  scale_colour_manual(values = c("PALM"="#f75f55","Netatmo"="#00adb3"))+
  facet_wrap(~station_id, nrow = 3)+
  xlab("date")+
  ylab("air temperature [°C]")+
  scale_x_datetime(breaks = "4 hour", date_labels = "%H:%M")+
  scale_y_continuous(breaks = seq(15,45,5))+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title = element_text(size = 28),
        axis.title.y = element_text(margin = margin(r=10)),
        axis.text = element_text(size = 24),
        axis.text.x = element_blank(),
        strip.text = element_text(size=18),
        legend.text = element_text(size=28))+
  guides(colour=guide_legend(override.aes = list(size=2)))
dev.off()

### Statistics ----
## For all stations, divided by domain
# Convert the data frame to a data table for the parent domain
compare_dt <- data.table(compare)

# Remove missing values from the data table
compare_dt <- na.omit(compare_dt)

# Calculate summary statistics
parent_stat <- compare_dt[, list(palm_mean = mean(palm_ta),
                                 netatmo_mean = mean(temp),
                                 palm_sd = sd(palm_ta), 
                                 netatmo_sd = sd(temp),
                                 # Calculate the Pearson correlation coefficient
                                 pearson = cor.test(palm_ta, temp, method = "pearson")$estimate,
                                 # Calculate the R-squared value of the linear regression 
                                 rsq = summary(lm(temp ~ palm_ta))$r.squared,
                                 # Calculate the slope of the linear regression
                                 slope = summary(lm(temp ~ palm_ta))$coefficients[2,1],
                                 # Calculate the intercept of the linear regression
                                 intercept = summary(lm(temp ~ palm_ta))$coefficients[1,1],
                                 # Calculate the root mean squared error (RMSE) between
                                 rmse = rmse(palm_ta, temp),
                                 # Calculate the mean squared error (MSE) between
                                 mse = mse(palm_ta, temp),
                                 # Calculate the index of agreement (IOA) between
                                 ioa = d(palm_ta, temp),
                                 # Calculate the bias between
                                 bias = Metrics::bias(temp, palm_ta),
                                 # Count the number of stations
                                 nr_stations = length(unique(station_id)),
                                 # Add a column indicating the domain
                                 domain = "parent")]

# Repeat the above steps for the child domain
compare_child_dt <- data.table(compare_child)
compare_child_dt <- na.omit(compare_child_dt)
child_stat <- compare_child_dt[, list(palm_mean = mean(palm_ta),
                                      netatmo_mean = mean(temp),
                                      palm_sd = sd(palm_ta), 
                                      netatmo_sd = sd(temp),
                                      pearson = cor.test(palm_ta, temp, method = "pearson")$estimate,
                                      rsq = summary(lm(temp ~ palm_ta))$r.squared,
                                      slope = summary(lm(temp ~ palm_ta))$coefficients[2,1],
                                      intercept = summary(lm(temp ~ palm_ta))$coefficients[1,1],
                                      rmse = rmse(palm_ta, temp),
                                      mse = mse(palm_ta, temp),
                                      ioa = d(palm_ta, temp),
                                      bias = Metrics::bias(temp, palm_ta),
                                      nr_stations = length(unique(station_id)),
                                      domain = "child")]


# Statistics for the parent domain reduced to child domain stations
# Filter the parent domain data frame to only include stations that are also in child domain data frame
stations_child <- unique(compare_child$station_id)
parent_filter <- subset(compare, compare$station_id %in% stations_child)

# Repeat the above steps for the statistics
parent_filter_dt <- data.table(parent_filter)
parent_filter_dt <- na.omit(parent_filter_dt)
parent_filter_stat <- parent_filter_dt[, list(palm_mean = mean(palm_ta),
                                              netatmo_mean = mean(temp),
                                              palm_sd = sd(palm_ta), 
                                              netatmo_sd = sd(temp),
                                              pearson = cor.test(palm_ta, temp, method = "pearson")$estimate,
                                              rsq = summary(lm(temp ~ palm_ta))$r.squared,
                                              slope = summary(lm(temp ~ palm_ta))$coefficients[2,1],
                                              intercept = summary(lm(temp ~ palm_ta))$coefficients[1,1],
                                              rmse = rmse(palm_ta, temp),
                                              mse = mse(palm_ta, temp),
                                              ioa = d(palm_ta, temp),
                                              bias = Metrics::bias(temp, palm_ta),
                                              nr_stations = length(unique(station_id)),
                                              domain = "parent_reduced")]

# Combine all summary statistics into one data frame
all_stats <- rbind(parent_stat, child_stat, parent_filter_stat)

# Save the statistics to a file
write.csv(all_stats, "evaluation_statistic.csv", row.names = F)

## Statistics for temporal clusters
# Filter the rows from compare_dt where date is between the specified range
tc_1 <- compare_dt[date %between% c("2020-08-08 00:00:00", "2020-08-08 05:00:00")]
tc_2 <- compare_dt[date %between% c("2020-08-08 06:00:00", "2020-08-08 09:00:00")]
tc_3 <- compare_dt[date %between% c("2020-08-08 10:00:00", "2020-08-08 13:00:00")]
tc_4 <- compare_dt[date %between% c("2020-08-08 14:00:00", "2020-08-08 17:00:00")]
tc_5 <- compare_dt[date %between% c("2020-08-08 18:00:00", "2020-08-08 21:00:00")]
tc_6 <- compare_dt[date %between% c("2020-08-08 22:00:00", "2020-08-09 01:00:00")]
tc_7 <- compare_dt[date %between% c("2020-08-09 02:00:00", "2020-08-09 05:00:00")]

# Calculate statistics
tc_1_stat <- tc_1[, list(# Root Mean Squared Error (RMSE)
                         rmse = rmse(palm_ta, temp),
                         # Mean Squared Error (MSE)
                         mse = mse(palm_ta, temp),
                         # Index of Agreement (IOA)
                         ioa = d(palm_ta, temp),
                         # Bias between the two values
                         bias = Metrics::bias(temp, palm_ta),
                         # Pearson correlation coefficient
                         pearson = cor.test(palm_ta, temp, method = "pearson")$estimate,
                         palm_mean = mean(palm_ta),
                         netatmo_mean = mean(temp),
                         palm_sd = sd(palm_ta),
                         netatmo_sd = sd(temp),
                         # Time slot number
                         tc = 1)]

# Repeat the above steps for all time clusters
tc_2_stat <- tc_2[, list(rmse = rmse(palm_ta, temp),
                         mse = mse(palm_ta, temp),
                         ioa = d(palm_ta, temp),
                         bias = Metrics::bias(temp, palm_ta),
                         pearson = cor.test(palm_ta, temp, method = "pearson")$estimate,
                         palm_mean = mean(palm_ta),
                         netatmo_mean = mean(temp),
                         palm_sd = sd(palm_ta),
                         netatmo_sd = sd(temp),
                         tc = 2)]

tc_3_stat <- tc_3[, list(rmse = rmse(palm_ta, temp),
                         mse = mse(palm_ta, temp),
                         ioa = d(palm_ta, temp),
                         bias = Metrics::bias(temp, palm_ta),
                         pearson = cor.test(palm_ta, temp, method = "pearson")$estimate,
                         palm_mean = mean(palm_ta),
                         netatmo_mean = mean(temp),
                         palm_sd = sd(palm_ta),
                         netatmo_sd = sd(temp),
                         tc = 3)]

tc_4_stat <- tc_4[, list(rmse = rmse(palm_ta, temp),
                         mse = mse(palm_ta, temp),
                         ioa = d(palm_ta, temp),
                         bias = Metrics::bias(temp, palm_ta),
                         pearson = cor.test(palm_ta, temp, method = "pearson")$estimate,
                         palm_mean = mean(palm_ta),
                         netatmo_mean = mean(temp),
                         palm_sd = sd(palm_ta),
                         netatmo_sd = sd(temp),
                         tc = 4)]

tc_5_stat <- tc_5[, list(rmse = rmse(palm_ta, temp),
                         mse = mse(palm_ta, temp),
                         ioa = d(palm_ta, temp),
                         bias = Metrics::bias(temp, palm_ta),
                         pearson = cor.test(palm_ta, temp, method = "pearson")$estimate,
                         palm_mean = mean(palm_ta),
                         netatmo_mean = mean(temp),
                         palm_sd = sd(palm_ta),
                         netatmo_sd = sd(temp),
                         tc = 5)]

tc_6_stat <- tc_6[, list(rmse = rmse(palm_ta, temp),
                         mse = mse(palm_ta, temp),
                         ioa = d(palm_ta, temp),
                         bias = Metrics::bias(temp, palm_ta),
                         pearson = cor.test(palm_ta, temp, method = "pearson")$estimate,
                         palm_mean = mean(palm_ta),
                         netatmo_mean = mean(temp),
                         palm_sd = sd(palm_ta),
                         netatmo_sd = sd(temp),
                         tc = 6)]

tc_7_stat <- tc_7[, list(rmse = rmse(palm_ta, temp),
                         mse = mse(palm_ta, temp),
                         ioa = d(palm_ta, temp),
                         bias = Metrics::bias(temp, palm_ta),
                         pearson = cor.test(palm_ta, temp, method = "pearson")$estimate,
                         palm_mean = mean(palm_ta),
                         netatmo_mean = mean(temp),
                         palm_sd = sd(palm_ta),
                         netatmo_sd = sd(temp),
                         tc = 7)]

# Combine all summary statistics into one data frame
tc_stat = rbind(tc_1_stat, tc_2_stat, tc_3_stat, tc_4_stat, tc_5_stat, tc_6_stat, tc_7_stat)

# Repeat the above steps for the child domain
tc_1_c <- compare_child_dt[date %between% c("2020-08-08 00:00:00", "2020-08-08 05:00:00")]
tc_2_c <- compare_child_dt[date %between% c("2020-08-08 06:00:00", "2020-08-08 09:00:00")]
tc_3_c <- compare_child_dt[date %between% c("2020-08-08 10:00:00", "2020-08-08 13:00:00")]
tc_4_c <- compare_child_dt[date %between% c("2020-08-08 14:00:00", "2020-08-08 17:00:00")]
tc_5_c <- compare_child_dt[date %between% c("2020-08-08 18:00:00", "2020-08-08 21:00:00")]
tc_6_c <- compare_child_dt[date %between% c("2020-08-08 22:00:00", "2020-08-09 01:00:00")]
tc_7_c <- compare_child_dt[date %between% c("2020-08-09 02:00:00", "2020-08-09 05:00:00")]

tc_1_c_stat <- tc_1_c[, list(rmse = rmse(palm_ta, temp),
                             mse = mse(palm_ta, temp),
                             ioa = d(palm_ta, temp),
                             bias = Metrics::bias(temp, palm_ta),
                             pearson = cor.test(palm_ta, temp, method = "pearson")$estimate,
                             palm_mean = mean(palm_ta),
                             netatmo_mean = mean(temp),
                             palm_sd = sd(palm_ta),
                             netatmo_sd = sd(temp),
                             tc = 1)]

tc_2_c_stat <- tc_2_c[, list(rmse = rmse(palm_ta, temp),
                             mse = mse(palm_ta, temp),
                             ioa = d(palm_ta, temp),
                             bias = Metrics::bias(temp, palm_ta),
                             pearson = cor.test(palm_ta, temp, method = "pearson")$estimate,
                             palm_mean = mean(palm_ta),
                             netatmo_mean = mean(temp),
                             palm_sd = sd(palm_ta),
                             netatmo_sd = sd(temp),
                             tc = 2)]

tc_3_c_stat <- tc_3_c[, list(rmse = rmse(palm_ta, temp),
                             mse = mse(palm_ta, temp),
                             ioa = d(palm_ta, temp),
                             bias = Metrics::bias(temp, palm_ta),
                             pearson = cor.test(palm_ta, temp, method = "pearson")$estimate,
                             palm_mean = mean(palm_ta),
                             netatmo_mean = mean(temp),
                             palm_sd = sd(palm_ta),
                             netatmo_sd = sd(temp),
                             tc = 3)]

tc_4_c_stat <- tc_4_c[, list(rmse = rmse(palm_ta, temp),
                             mse = mse(palm_ta, temp),
                             ioa = d(palm_ta, temp),
                             bias = Metrics::bias(temp, palm_ta),
                             pearson = cor.test(palm_ta, temp, method = "pearson")$estimate,
                             palm_mean = mean(palm_ta),
                             netatmo_mean = mean(temp),
                             palm_sd = sd(palm_ta),
                             netatmo_sd = sd(temp),
                             tc = 4)]

tc_5_c_stat <- tc_5_c[, list(rmse = rmse(palm_ta, temp),
                             mse = mse(palm_ta, temp),
                             ioa = d(palm_ta, temp),
                             bias = Metrics::bias(temp, palm_ta),
                             pearson = cor.test(palm_ta, temp, method = "pearson")$estimate,
                             palm_mean = mean(palm_ta),
                             netatmo_mean = mean(temp),
                             palm_sd = sd(palm_ta),
                             netatmo_sd = sd(temp),
                             tc = 5)]

tc_6_c_stat <- tc_6_c[, list(rmse = rmse(palm_ta, temp),
                             mse = mse(palm_ta, temp),
                             ioa = d(palm_ta, temp),
                             bias = Metrics::bias(temp, palm_ta),
                             pearson = cor.test(palm_ta, temp, method = "pearson")$estimate,
                             palm_mean = mean(palm_ta),
                             netatmo_mean = mean(temp),
                             palm_sd = sd(palm_ta),
                             netatmo_sd = sd(temp),
                             tc = 6)]

tc_7_c_stat <- tc_7_c[, list(rmse = rmse(palm_ta, temp),
                             mse = mse(palm_ta, temp),
                             ioa = d(palm_ta, temp),
                             bias = Metrics::bias(temp, palm_ta),
                             pearson = cor.test(palm_ta, temp, method = "pearson")$estimate,
                             palm_mean = mean(palm_ta),
                             netatmo_mean = mean(temp),
                             palm_sd = sd(palm_ta),
                             netatmo_sd = sd(temp),
                             tc = 7)]

tc_c_stat = rbind(tc_1_c_stat, tc_2_c_stat, tc_3_c_stat, tc_4_c_stat, tc_5_c_stat, tc_6_c_stat, tc_7_c_stat)

# Save the temporal statistics to files
write.csv(tc_stat, "statistics_city_tc.csv")
write.csv(tc_c_stat, "statistics_city_tc_child.csv")