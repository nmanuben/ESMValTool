## Regimes namelist
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
start <- as.POSIXct(params$start)
end <- as.POSIXct(params$end)

#Regime parameters
ncenters <- params$ncenters
cluster_method <- params$cluster_method
EOFS <- params$EOFS


library(s2dverification)
library(startR)
library(multiApply)
library(devtools)
library(climdex.pcic)
source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Climdex/R/Climdex.R')
library(parallel)

fullpath_filenames <- input_files[[var0]]
historical_data <- Start(model = fullpath_filenames_historical,
              var = var0,
              var_var = 'var_names',
              #  sdate = paste0(seq(1970, 2000, by = 1), "0101"),
              time = values(list(as.POSIXct(start_historical), 
                                 as.POSIXct(end_historical))),
              time_tolerance = as.difftime(1, units = 'days'), 
              lat = values(list(lat.min, lat.max)),
              lon = values(list(lon.min, lon.max)),
              lon_var = 'lon',
              lon_reorder = CircularSort(0, 360),
              return_vars = list(time = 'model', lon = 'model', lat = 'model'),
              retrieve = TRUE)


lat <- attr(data, "Variables")$dat1$lat
lon <- attr(data, "Variables")$dat1$lon
time_dimension <- which(names(dim(data)) == "time")


###Compute the quantiles and standard deviation for the historical period.
if (var0 != "pr") {
  if (var0 == "tasmin") {
    metric <- "tx10p"
    quantile <- 0.1
  } else if (var0 == "tasmax") {
    metric <- "tx90p" 
    quantile <- 0.9
  } else if (var0 == "sfcWind") {
    metric <- "Wx"
    quantile <- 0.9
  }
  
  base_range <- as.numeric(c(substr(start_historical, 1, 4), substr(end_historical, 1, 4)))
  
  #Compute the 90th percentile for the historical period 
  thresholds <- Threshold(data, base_range=as.numeric(base_range),
                          qtiles = quantile, ncores = detectCores() -1)
  #Compute the relevant indice
  base_index <- Climdex(data = data, metric = metric, 
                               quantiles = thresholds, ncores = detectCores() - 1)
  base_sd <- Apply(list(base_index$result), target_dims = list(c(1)), AtomicFun = "sd")$output1
  base_sd_historical <- InsertDim(base_sd, 1, dim(base_index$result)[1])
  historical_index_standardized <- (base_index$result - 10) / base_sd_historical
  
} else {
  metric <- "cdd"
  base_index_cdd <- Climdex(data = data, metric = metric, ncores = detectCores() - 1)
  sd_cdd <- Apply(list(base_index_cdd$result), target_dims = list(c(1)), AtomicFun = "sd")$output1
  mean_cdd <- Apply(list(base_index_cdd$result), target_dims = list(c(1)), AtomicFun = "mean")$output1
  metric2 <- "rx5day"
  base_index_rx5day <- Climdex(data = data, metric = metric2, ncores = detectCores() - 1)
}

rm(historical_data)


#Compute the time series of the relevant index, using the quantiles and standard deviation from the index
for (i in 1 : length(rcp_scenario)) {
 projection_data <- Start(model = fullpath_filenames_projection[[i]],
                           var = var0,
                           var_var = 'var_names',
                           time = values(list(as.POSIXct(start_projection), 
                                              as.POSIXct(end_projection))),
                           time_tolerance = as.difftime(1, units = 'days'), 
                           lat = values(list(lat.min, lat.max)),
                           lon = values(list(lon.min, lon.max)),
                           lon_var = 'lon',
                           lon_reorder = CircularSort(0, 360),
                           return_vars = list(time = 'model', lon = 'model', lat = 'model'),
                           retrieve = TRUE)
  if (metric == "tx90p" || metrix = "tx10p" || metric = "Wx") {
    base_sd_proj <- InsertDim(base_sd, 1, dim(projection_data$result)[1])
    projection_index_standardized <- (projection_index$result - 10) / base_sd_proj
    rm(base_sd_proj)
  }
 attributes(lon) <- NULL
 attributes(lat) <- NULL
 attributes(years) <- NULL
 dim(years) <- c(length(years))
 dim(lon) <-  c(lon = length(lon))
 dim(lat) <- c(lat = length(lat))
 ArrayToNetCDF(list(tx90p = historical_index_standardized,  lat = lat, lon = lon), paste0(metric, "_","IPSL-CM5A-LR","_", start_historical, "_", end_historical, ".nc"))
 
}

  








