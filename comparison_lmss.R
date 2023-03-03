####
# Script to plot PALM and LMSS data
####

# Set the working directory
setwd("your working directory")

# Load required libraries
library(ggplot2)
library(lubridate)
library(ggpubr)

# Read the data
palm <- read.csv("palm_lmss.csv")
palm_ws <- read.csv("palm_lmss_ws.csv")
lmss <- read.csv("LMSS_080820_090820.csv")

# Convert the date column of the lmss data to POSIXct
lmss$date <- as.POSIXct(lmss$date)

# Convert the time column of the palm data to POSIXct 
palm$time <- as.POSIXct(palm$time, format = "%Y-%m-%d %H:%M:%OS")
palm_ws$date <- as.POSIXct(palm_ws$date, format = "%Y-%m-%d %H:%M:%OS")

# Round the time of the palm data to the nearest hour
palm$time <- round_date(palm$time, unit = "hour")
palm_ws$date <- round_date(palm_ws$date, unit = "hour")

# Change the name of the second column in the palm data
colnames(palm)[2] <- c("date")

# Merge the three data frames into one, joining on the date column
lmss_palm <- merge(lmss, palm, by = "date")
lmss_palm <- merge(lmss_palm, palm_ws, by = "date")

# Remove columns 5 and 8 from the data frame
lmss_palm <- lmss_palm[-c(5,8)]

# Round the value of the lmss wind speed
lmss_palm$WSMEAN <- round(lmss_palm$WSMEAN, 2)

# Calculate the difference between the palm and lmss air temperature 
lmss_palm$ta_diff <- lmss_palm$palm_ta - lmss_palm$TMEAN

# Create the first plot (temperature)
p1 <- ggplot(lmss_palm, aes(x=date))+
  # Add the data as lines
  geom_line(aes(y=TMEAN, colour = "LMSS"))+
  geom_line(aes(y=palm_ta, colour = "PALM"))+
  # Label the x-axis
  xlab("date")+
  # Label the y-axis
  ylab("ta [°C]")+
  # Set the colours for the lines and the labels in the legend
  scale_colour_manual(values = c("LMSS"="#00adb3", "PALM"="#f75f55"),
                      labels = c("LMSS","PALM"))+
  # Customise the x-axis
  scale_x_datetime(breaks = "7 hour", date_labels = "%d.%m %H:%M")+
  # Customise the y-axis
  scale_y_continuous(breaks = seq(15,45,5))+
  # Apply a black and white theme
  theme_bw()+
  # Customise the theme
  theme(axis.title = element_text(size = 28),
        axis.text = element_text(size = 24),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(r=10)),
        legend.title = element_blank(),
        legend.text = element_text(size=28))

# Create the second plot (humidity)
p2 <- ggplot(lmss_palm, aes(x=date))+
  geom_line(aes(y=HMEAN, colour = "LMSS"))+
  geom_line(aes(y=palm_rh, colour = "PALM"))+
  xlab("date")+
  ylab("rh [%]")+
  scale_colour_manual(values = c("LMSS"="#00adb3", "PALM"="#f75f55"),
                      labels = c("LMSS","PALM"))+
  scale_x_datetime(breaks = "7 hour", date_labels = "%d.%m %H:%M")+
  scale_y_continuous(breaks = seq(20,80,10))+
  theme_bw()+
  theme(axis.title = element_text(size = 28),
        axis.text = element_text(size = 24),
        axis.text.x = element_blank(),
        axis.title.y = element_text(margin = margin(r=10)),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=28))

# Create the third plot (wind speed)
p3 <- ggplot(lmss_palm, aes(x=date))+
  geom_line(aes(y=WSMEAN, colour = "LMSS"))+
  geom_line(aes(y=palm_ws, colour = "PALM"))+
  xlab("date")+
  ylab("ws [m/s]")+
  scale_colour_manual(values = c("LMSS"="#00adb3", "PALM"="#f75f55"),
                      labels = c("LMSS","PALM"))+
  scale_x_datetime(breaks = "7 hour", date_labels = "%d.%m %H:%M")+
  scale_y_continuous(breaks = seq(0,4,1))+
  theme_bw()+
  theme(axis.title = element_text(size = 28),
        axis.text = element_text(size = 24),
        axis.title.y = element_text(margin = margin(r=10)),
        legend.title = element_blank(),
        legend.text = element_text(size=28))

# Save the plot as a PNG image file
png("palm_lmss.png", width = 1200, height = 1200)
# Arrange the three plots (p1, p2, p3) on top of each other
# Include a common legend for the plots
ggarrange(p1, p2, p3, nrow = 3, heights = c(3,3,2), align = "v", common.legend = T, legend = "right")
# Close the PNG device
dev.off()