#this script produces a GIF file of viewing distances across time
library(LightLogR)
library(tidyverse)
library(gganimate)
library(magick)
library(legendry)
library(patchwork)
library(ggridges)
load("data/cleaned/data.RData")
load("scripts/color.Rdata")


# Set up spatial plot -----------------------------------------------------

#set visualization parameters
extras <- list(
  geom_tile(),
  scale_fill_viridis_c(breaks = c(0, 50, 100, 150, 200),
    direction = -1, limits = c(0, 200),
                       oob = scales::oob_squish_any,
    guide = "colbar"),
  theme_minimal(),
  guides(colour = "none"),
  coord_fixed(),
  theme(legend.position = "bottom"),
  labs(x = "X position (°)", y = "Y position (°)", 
       fill = "Distance (cm)", alpha = "Confidence (0-255)"))

# p <-
#   dataVEET3 |> 
#   # filter_Datetime(start = "2024-06-10 13:10:00", length = "20 mins") |> 
#   cut_Datetime(unit = "5 secs", New.colname = Datetime) |> 
#   filter_Time(start = "13:10:00", length = "00:19:59") |> 
#   add_Time_col() |> 
#   mutate(dist1 = ifelse(dist1 == 0, 5000, dist1)) |> #replace 0 distances with 5m
#   # filter(conf1 >= 0.1 | dist1 == Inf) |> #remove data that has less than 10% confidence
#   ggplot(aes(x=x.pos, y=y.pos, fill = dist1/10, group = seq_along(Time)))+ #plot the data
#   labs(title = "Spatially resolved distance from eye-level. Time: {frame_time}") +
#   extras + transition_time(Time) +
#   facet_wrap(~Date)



# # prepare distance data -------------------------------------------------
tz   <- "US/Central"   # Time zone in which device was recording (e.g., US Central Time)

path <- "data/02_VEET_L.csv.zip"

dataVEET3 <- LightLogR::import$VEET(path, tz = tz, modality = "TOF", manual.id = "VEET")

dataVEET3 <- dataVEET3 |>
  aggregate_Datetime(unit = "5 secs") |>     # aggregate to 5-second bins
  gap_handler(full.days = TRUE) |>           # explicit NA for any gaps in between
  add_Date_col(group.by = TRUE) |> 
  remove_partial_data(dist1_0, threshold.missing = "1 hour")

dataVEET3 <- 
  dataVEET3 |> 
  pivot_longer(
    cols = -c(Datetime, file.name, Id, is.implicit, time_stamp, modality, Date),
    names_to = c(".value", "position"),
    names_pattern = "(conf1|conf2|dist1|dist2)_(\\d+)"
  )

dataVEET3 <- 
  dataVEET3 |> 
  # ungroup() |> 
  mutate(position = as.numeric(position),
         y.pos = (position %/% 8)+1,
         y.pos = scale(y.pos, scale = FALSE)*52/8,
         x.pos = (position %% 8)+1,
         x.pos = scale(x.pos, scale = FALSE)*41/8,
         observation = cumsum(position == 0),
         across(starts_with("conf"), \(x) x/255)
  )

extras <- list(
  geom_tile(),
  scale_fill_viridis_c(breaks = c(0, 50, 100, 150, 200),
                       direction = -1, limits = c(0, 200),
                       oob = scales::oob_squish_any,
                       guide = "colbar"),
  theme_void(),
  theme(plot.title.position = "plot"),
  guides(colour = "none", fill = "none"),
  coord_fixed(expand = FALSE),
  theme(axis.title.y = element_text(angle = 90, size = 12),
        axis.text.y = element_text(angle = 90, size = 10),
        axis.ticks = element_line(),
        axis.ticks.length = unit(2, "pt")),
  labs(x = NULL, 
       # y = NULL,
       y = "Spatial distance",
       fill = "Distance (cm)", alpha = "Confidence (0-255)"))

dist_bar <-
dataVEET3 |> 
  filter_Datetime(start = "2025-06-18 16:15:00", length = "301 secs") |>
  cut_Datetime(unit = "5 secs", New.colname = Datetime) |>
  # filter_Time(start = "13:10:00", length = "00:19:59") |>
  add_Time_col() |>
    mutate(dist1 = ifelse(dist1 == 0, 5000, dist1)) |> #replace 0 distances with 5m
    # filter(conf1 >= 0.1 | dist1 == Inf) |> #remove data that has less than 10% confidence
    ggplot(aes(x=x.pos, y=y.pos, fill = dist1/10, group = seq_along(Time)))+ #plot the data
    transition_time(Time) +
  scale_y_continuous(breaks = c(-20, -10, 0, 10, 20), labels = c("-20°", "", "0°", "", "20°")) +
  # scale_x_continuous(breaks = c(-20, 0, 20), labels = \(x) character(3)) +
  # scale_y_continuous(breaks = c(-20, 0, 20), labels = \(x) character(3)) +
    extras 
  
# dist_bar <- animate(dist_bar, height = 160, width = 160, res = 110, fps = 5, duration = 1, renderer = magick_renderer())
dist_bar <- animate(dist_bar, height = 160, width = 180, res = 110, fps = 5, duration = 60, renderer = magick_renderer())
dist_bar

# prepare spectral data ---------------------------------------------------

#if a version of LightLogR ≤ 0.9.2 is used, this script needs to be imported, as the device data was collected with a newer firmware version that changes the output format.
source("scripts/VEET_import.R")
dataVEET2 <- import$VEET(path, tz = tz, modality = "PHO", manual.id = "VEET")

# Aggregate spectral data to 5-minute intervals and mark gaps
dataVEET2 <- dataVEET2 |>
  aggregate_Datetime(unit = "5 secs") |>     # aggregate to 5-second bins
  gap_handler(full.days = TRUE) |>           # explicit NA for any gaps in between
  add_Date_col(group.by = TRUE) |> 
  remove_partial_data(Clear, threshold.missing = "1 hour")

count.columns <- c("s415", "s445", "s480", "s515", "s555", "s590", "s630", 
                   "s680", "s910", "Dark", "Clear") #column names

gain.ratios <- #gain ratios as specified by the manufacturers reconstruction guide
  tibble(
    gain = c(0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512),
    gain.ratio = c(0.008, 0.016, 0.032, 0.065, 0.125, 0.25, 0.5, 1, 2, 3.95, 7.75)
  ) 

#normalize data
dataVEET2 <-
  dataVEET2 |> 
  mutate(across(c(s415:Clear), \(x) (x - Dark)/integration_time)) |>  #remove dark counts & scale by integration time
  normalize_counts( #function to normalize counts
    gain.columns = rep("Gain", 11), #all sensor channels share the gain value
    count.columns = count.columns, #sensor channels to normalize
    gain.ratios #gain ratios
  ) |> 
  select(-c(s415:Clear)) |> # drop original raw count columns
  rename_with(~ str_remove(.x, ".normalized"))


#import calibration matrix
calib_mtx <- 
  read_csv("data/VEET_calibration_matrix.csv", show_col_types = FALSE) |> 
  column_to_rownames("wavelength")

# Reconstruct spectral power distribution (SPD) for each observation
dataVEET2 <- dataVEET2 |> mutate(
  Spectrum = spectral_reconstruction(
    pick(s415:s910),   # pick the normalized sensor columns
    calibration_matrix = calib_mtx, 
    format = "long" # return a long-form list column (wavelength, intensity)
  )
)

dataVEET2 <- 
dataVEET2 |> 
  filter_Date(length = "1 day") |> #keep only observations for one day (from start)
  mutate( 
    Illuminance = Spectrum |> #Use the spectrum,...
      map_dbl(spectral_integration, #... call the function spectral_integration for each,...
              action.spectrum = "photopic", #... use the brightness sensitivity function,...
              general.weight = "auto") #... and apply the appropriate efficacy weight.
  )
  
ill_bar <-
dataVEET2 |>
  filter_Datetime(start = "2025-06-18 16:15:00", length = "300 secs") |>
  cut_Datetime(unit = "5 secs", New.colname = Datetime) |>
  # filter_Time(start = "13:10:00", length = "00:19:59") |>
  add_Time_col() |>
  drop_na(Illuminance) |> 
  ggplot(aes(x=1, y=Illuminance)) + #plot the data
  geom_col(fill = "gold1", position = "identity") +
  geom_hline(aes(yintercept = 0), linewidth = 0.25) +
  labs(title = "{frame_time}") +
  # scale_y_continuous(trans = "symlog") + 
  transition_time(Time) +
  labs(y = "Illuminance") +
  theme_void() +
  coord_cartesian(expand = FALSE)+
  theme(plot.title = element_text(hjust = 0),
        plot.title.position = "plot",
        axis.title.y = element_text(angle = 90, size = 12))

# ill_bar  <- animate(ill_bar, height = 160, width = 80, res = 110, fps = 5, duration = 1, renderer = magick_renderer())
ill_bar  <- animate(ill_bar, height = 160, width = 80, res = 110, fps = 5, duration = 60, renderer = magick_renderer())
ill_bar

spec_bar <-
dataVEET2 |> 
  filter_Datetime(start = "2025-06-18 16:15:00", length = "5 mins") |>
  cut_Datetime(unit = "5 secs", New.colname = Datetime) |>
  # filter_Time(start = "13:10:00", length = "00:19:59") |>
  add_Time_col() |>
  unnest(Spectrum) |> #create a long format of the data where the spectrum is unnested
  drop_na(irradiance) |> 
  ggplot(aes(x = wavelength)) +
  # scale_y_continuous(trans = "symlog") +
  # geom_line(aes(y = irradiance * 1000)) +
    ggridges::geom_ridgeline_gradient(ggplot2::aes(
      y = 0,
      height = irradiance * 1000,
      fill = wavelength
    ), linewidth = 0.25) +
  geom_hline(aes(yintercept = 0), linewidth = 0.25) +
  ggplot2::scale_fill_gradientn(colors = colorSchema, guide = "none") +
  theme_void() +
  labs(x = "Spectral irradiance") +
  coord_cartesian(
    # ylim = c(-0.2, NA), 
    expand = FALSE) +
  theme(
        axis.title.x = element_text(size = 12)) +
  transition_time(Time)

# spec_bar  <- animate(spec_bar , height = 100, width = 260, res = 110, fps = 5, duration = 1, renderer = magick_renderer())
spec_bar  <- animate(spec_bar , height = 100, width = 260, res = 110, fps = 5, duration = 60, renderer = magick_renderer())
spec_bar 

# Animation ---------------------------------------------------------------

new_gif <- image_append(c(spec_bar[1], image_append(c(ill_bar[1], dist_bar[1]))), stack = TRUE)
for(i in 1:300){
  combined <- image_append(c(spec_bar[i], image_append(c(ill_bar[i], dist_bar[i]))), stack = TRUE)
  new_gif <- c(new_gif, combined)
}

new_gif

image_write_gif(new_gif,"output/abstract_med.gif", delay = 1/5)


# p
# 
# Plot <- animate(p, height = 500, width = 500, res = 80, fps = 4, duration = 60, renderer = magick_renderer())

image_write_video(Plot, "output/abstract.mp4", framerate = 5)  
