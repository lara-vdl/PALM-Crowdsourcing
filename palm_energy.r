####
# Script to plot PALM energy fluxes at the surface
####

# Set the working directory 
setwd("your working directory")

# Load the required libraries
library(ggplot2)
library(lubridate)

# Read the data
palm <- read.csv("palm_lmss_energy.csv")

# Convert the date column to a POSIXct object
palm$date <- as.POSIXct(palm$date, format = "%Y-%m-%d %H:%M:%OS")

# Round the date column to the nearest hour
palm$date <- round_date(palm$date, unit = "hour")

# Remove the last timestep from the data
palm <- palm[1:28,]

# Calculate total heat flux at the surface
palm$total_hf <- palm$sens_hf + palm$lat_hf + palm$ground_hf

# Create the plot
# Save the plot as a PNG image file
png("palm_energy.png", width = 1200, height = 900)
# Create the plot
ggplot(palm, aes(x=date)) + 
  # Add the heat fluxes as points and lines 
  geom_point(aes(y=total_hf, colour="total_hf"), size =2) +
  geom_line(aes(y=total_hf, colour="total_hf")) +
  geom_point(aes(y=sens_hf, colour="sens_hf"), size =2) +
  geom_line(aes(y=sens_hf, colour="sens_hf")) +
  geom_point(aes(y=lat_hf, colour="lat_hf"), size =2) +
  geom_line(aes(y=lat_hf, colour="lat_hf")) +
  geom_point(aes(y=ground_hf, colour="ground_hf"), size =2) +
  geom_line(aes(y=ground_hf, colour="ground_hf")) + 
  # Add a horizontal line at y = 0 
  geom_hline(yintercept = 0, colour ="#666666") +
  # Set the colours for the lines and the labels in the legend
  scale_colour_manual(values = c("total_hf"="black","sens_hf"="red","lat_hf"="blue","ground_hf"="orange"),
                      labels = c(expression(paste("Q","*")),expression(Q[H]),expression(Q[E]),expression(Q[G]))) +
  # Label the x axis
  xlab("date") +
  # Label the y axis
  ylab(expression(paste("energy fluxes [W ", m^-2,"]"))) +
  # customise the x-axis
  scale_x_datetime(breaks = "5 hour", date_labels = "%d.%m %H:%M") +
  # customise the x-axis
  scale_y_continuous(breaks = seq(-100,600,50)) +
  # Apply a black and white theme
  theme_bw() +
  # customise the theme
  theme(axis.title = element_text(size = 28),
        axis.text = element_text(size = 24),
        axis.title.y = element_text(margin = margin(r=10)),
        legend.title = element_blank(),
        legend.text = element_text(size=28),
        legend.text.align = 0)
# Close the PNG device
dev.off()