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

fullpath_filenames <- input_files[[var0]]

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

lat <- attr(historical_data, "Variables")$dat1$lat
lon <- attr(historical_data, "Variables")$dat1$lon
time_dimension <- which(names(dim(historical_data)) == "time")

attributes(lon) <- NULL
attributes(lat) <- NULL
# attributes(years) <- NULL
# dim(years) <- c(length(years))
dim(lon) <-  c(lon = length(lon))
dim(lat) <- c(lat = length(lat))

###Compute the quantiles and standard deviation for the historical period.

if (var0 == "tasmin") {
  metric <- "tx10p"
  quantile <- 0.1
} else if (var0 == "tasmax") {
  metric <- "tx90p" 
  quantile <- 0.9
} else if (var0 == "sfcWind") {
  metric <- "Wx"
  quantile <- 0.9
} else if (var0 == "pr") {
  metric <- c("cdd", "rx5day")
}
base_sd <- base_sd_historical <- base_mean <- list()
for (m in 1 : length(metric)) {
  base_range <- as.numeric(c(substr(start_historical, 1, 4), substr(end_historical, 1, 4)))
  
  #Compute the 90th percentile for the historical period
  if (var0 != "pr") {
    thresholds <- Threshold(historical_data, base_range=as.numeric(base_range),
                            qtiles = quantile, ncores = detectCores() -1)
    base_index <- Climdex(data = historical_data, metric = metric[m], 
                          quantiles = thresholds, ncores = detectCores() - 1)
  } else {
    base_index <- Climdex(data = historical_data, metric = metric[m], ncores = detectCores() - 1)
  }
  
  base_index$result <- Trend(base_index$result, posTR = 1, interval = 1)$detrended
  base_sd[[m]] <- Apply(list(base_index$result), target_dims = list(c(1)), AtomicFun = "sd")$output1
  base_sd_historical[[m]] <- InsertDim(base_sd[[m]], 1, dim(base_index$result)[1])
  
  if (var0 != "pr") {
    base_mean[[m]] <- 10
  } else {
    base_mean[[m]] <-  Apply(list(base_index$result), target_dims = list(c(1)), AtomicFun = "mean")$output1
    base_mean_historical <- InsertDim(base_mean[[m]], 1, dim(base_index$result)[1])
  }
  historical_index_standardized <- (base_index$result - base_mean_historical) / base_sd_historical[[m]]
  for (mod in 1 : dim(projection_data)[model_dim]) {
    ArrayToNetCDF(list(metric = historical_index_standardized[,mod,1, ,],  lat = lat, lon = lon), 
                  paste0(metric[m], "_",model_names[mod],"_", "historical", "_", start_historical, "_", end_historical, ".nc"))
  }
}  
  

#Compute the time series of the relevant index, using the quantiles and standard deviation from the index
for (i in 1 : length(fullpath_filenames_projection)) {
    projection_data <- Start(model = fullpath_filenames_projection[[i]],
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
  for (m in 1 : length(metric)) {
    
    if (var0 != "pr") {
      projection_index <- Climdex(data = projection_data, metric = metric[m], 
                                  quantiles = thresholds, ncores = detectCores() - 1)
      projection_mean <- 10
    } else {
      projection_index <- Climdex(data = projection_data, metric = metric[m], 
                                  ncores = detectCores() - 1)
      projection_mean <- InsertDim(base_mean[[m]], 1, dim(projection_index$result)[1])
    }
    
    base_sd_proj <- InsertDim(base_sd[[m]], 1, dim(projection_index$result)[1])
    projection_index_standardized <- (projection_index$result - projection_mean) / base_sd_proj
    model_dim <- which(names(dim(projection_index_standardized)) == "model")
    for (mod in 1 : dim(projection_data)[model_dim]) {
      ArrayToNetCDF(list(metric = projection_index_standardized[,mod,1, ,],  lat = lat, lon = lon), 
                    paste0(metric[m], "_",model_names[mod],"_", rcp_scenario[i], "_", start_projection, "_", end_projection, ".nc"))
      }
    }
}
    
    
      
     
      
  








