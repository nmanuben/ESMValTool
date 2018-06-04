## Insurance products
args <- commandArgs(trailingOnly = TRUE)
params <- yaml::read_yaml(args[1])

plot_dir <- params$plot_dir
run_dir <- params$run_dir
work_dir <- params$work_dir
## Create working dirs if they do not exist
dir.create(plot_dir, recursive = TRUE)
dir.create(run_dir, recursive = TRUE)
dir.create(work_dir, recursive = TRUE)

input_files_per_var <- yaml::read_yaml(params$input_files)
var_names <- names(input_files_per_var)
input_files <- lapply(var_names, function(x) names(input_files_per_var[[x]]))
names(input_files) <- var_names
model_names <- lapply(input_files_per_var, function(x) unname(sapply(x, '[[', 'model')))

## Do not print warnings
#options(warn=-1)


#Var considered
var0 <- var_names[1]

#Region considered
lat.max <- params$lat_max
lat.min <- params$lat_min
lon.max <- params$lon_max
lon.min <- params$lon_min


#Start and end periods for the historical and projection periods
start_historical <- as.POSIXct(params$start_historical)
end_historical <- as.POSIXct(params$end_historical)
start_projection <- as.POSIXct(params$start_projection)
end_projection <- as.POSIXct(params$end_projection)

#Regime parameters
metric <- params$metric
rcp8.5 <- params$rcp8.5
rcp2.6 <- params$rcp2.6
rcp_scenario <- c(rcp8.5, rcp2.6)

library(s2dverification)
library(startR)
library(multiApply)
library(devtools)
library(climdex.pcic)
source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Climdex/R/Climdex.R')
source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Climdex/R/Threshold.R')
library(parallel)
source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Climdex/R/DTR.R')


historical_data <- Start(model = fullpath_filenames_historical,
                           var = var0,
                           var_var = 'var_names',
                           #  sdate = paste0(seq(1970, 2000, by = 1), "0101"),
                           time = values(list(as.POSIXct(start_historical), 
                                              as.POSIXct(end_historical))),
                           time_tolerance = as.difftime(0, units = 'days'), 
                           lat = values(list(lat.min, lat.max)),
                           lon = values(list(lon.min, lon.max)),
                           lon_var = 'lon',
                           lon_reorder = CircularSort(0, 360),
                           return_vars = list(time = 'model', lon = 'model', lat = 'model'),
                           retrieve = TRUE)
if (var0 == "tasmax") {
  op = ">"
  qtile = 0.9
} else {
  op = "<"
  qtile = 0.1
}
base_range <- c(as.numeric(substr(start_historical, 1, 4)), as.numeric(substr(end_historical, 1, 4)))
threshold <- Threshold(historical_data, base_range = base_range, qtiles = qtile, ncores = NULL)

for (i in 1 : length(fullpath_filenames_projection)){
  rcp_data <- Start(model = fullpath_filenames_projection[[i]],
                      var = var0,
                      var_var = 'var_names',
                      time = values(list(as.POSIXct(start_projection), 
                                         as.POSIXct(end_projection))),
                      time_tolerance = as.difftime(0, units = 'days'), 
                      lat = values(list(lat.min, lat.max)),
                      lon = values(list(lon.min, lon.max)),
                      lon_var = 'lon',
                      lon_reorder = CircularSort(0, 360),
                      return_vars = list(time = 'model', lon = 'model', lat = 'model'),
                      retrieve = TRUE)
  
 heatwave <- Heatwave(rcp_data, threshold, op = op, spell_length = 6, by.seasons = TRUE, ncores = NULL)
  ## Select summer season
  heatwave_summer <- heatwave$result[seq(2, dim(heatwave$result)[1] - 2, by = 4), 1, 1, , ]
  years <-  as.numeric(substr(start_projection, 1, 4)) : as.numeric(substr(end_projection, 1, 4))
  lat <- attr(rcp_data, "Variables")$dat1$lat
  lon <- attr(rcp_data, "Variables")$dat1$lon
  lon[lon > 180] <- lon[lon > 180] - 360
  lon_order <- sort(lon, index.return = TRUE)
  data <- heatwave_summer
  data <- aperm(data, c(3,2,1))
  names(dim(data)) <- c("lon", "lat", "time")
  data <- Subset(data, "lon", lon_order$ix)
  lon <- lon_order$x
  attributes(lon) <- NULL
  attributes(lat) <- NULL
  dim(lon) <-  c(lon = length(lon))
  dim(lat) <- c(lat = length(lat))
  #data <- Mean1Dim(dtr_indicator$indicator, 2)
  metadata <- list(index = list(dim = list(list(name='time', unlim = FALSE, prec='double'))))
  attr(data, 'variables') <- metadata
  time <- years
  attributes(time) <- NULL
  dim(time) <- c(time = length(time))
  metadata <- list(time = list(standard_name = 'time', long_name = 'time', units = 'years since 0-0-0 00:00:00', prec = 'double', dim = list(list(name='time', unlim = FALSE))))
  attr(time, "variables") <- metadata
  ArrayToNetCDF(list(metric= data, lat = lat, lon = lon, time = time), 
                paste0(var0, "_extreme_spell_duration", "_",model_names,"_", rcp_scenario[i], "_", start_projection, "_", end_projection, ".nc"))
  brks <- c(0 : 50, 5)
  title <- paste0("Number of days JJA ", var0, " ", substr(start_projection, 1, 4), "-", 
                  substr(end_projection, 1, 4), " ",op, " the ", substr(as.character(qtile), 3, 3),
                  "0th quantile for ",substr(start_historical, 1, 4), "-", substr(end_historical, 1, 4), 
                  " (",rcp_scenario[i], ")") 
  PlotEquiMap(Mean1Dim(data, 3), lon = lon, lat = lat, fill = FALSE, 
              brks = brks, color_fun = clim.palette("yellowred"), 
              units = paste0("Days" ), toptitle = title,
              fileout = paste0(var0, "_extreme_spell_duration", "_",model_names,"_", rcp_scenario[i], 
                               "_", start_projection, "_", end_projection, ".pdf"),
              title_scale = 0.5)
}


## Plots
source("https://earth.bsc.es/gitlab/rserrano/colorbar/raw/master/int_breaks.R")
new_breaks <- int_breaks(data[1, , ], method = "equal", zero_centered = FALSE, n=10)







