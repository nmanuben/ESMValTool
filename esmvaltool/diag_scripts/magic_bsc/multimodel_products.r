####REQUIRED SYSTEM LIBS
####ŀibssl-dev
####libnecdf-dev
####cdo

# conda install -c conda-forge r-ncdf4

Sys.setenv(TAR = '/bin/tar')
library(s2dverification)
library(startR, lib.loc='/home/Earth/ahunter/R/x86_64-unknown-linux-gnu-library/3.2/')
library(multiApply, lib.loc='/home/Earth/ahunter/R/x86_64-unknown-linux-gnu-library/3.2/')
library(ggplot2)
library(yaml)

##Until integrated into current version of s2dverification
source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-MagicWP5/R/AnoAgree.R')
source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Magic_WP6/R/WeightedMean.R')




#Parsing input file paths and creating output dirs
args <- commandArgs(trailingOnly = TRUE)
params <- read_yaml(args[1])
plot_dir <- params$plot_dir
run_dir <- params$run_dir
work_dir <- params$work_dir
## Create working dirs if they do not exist
dir.create(plot_dir, recursive = TRUE)
dir.create(run_dir, recursive = TRUE)
dir.create(work_dir, recursive = TRUE)

input_files_per_var <- yaml::read_yaml(params$input_files)
var_names <- names(input_files_per_var)
model_names <- lapply(input_files_per_var, function(x) x$model)
model_names <- unname(model_names)
var0 <- lapply(input_files_per_var, function(x) x$short_name)
fullpath_filenames <- names(var0)
var0 <- unname(var0)[1]
experiment <- lapply(input_files_per_var, function(x) x$exp)
experiment <- unlist(unname(experiment))

climatology_class <- params$climatology_class
anomaly_class <- params$anomaly_class
climatology_files <- which(unname(experiment) == as.character(climatology_class))
anomaly_files <- which(unname(experiment) == as.character(anomaly_class))

model_names <-  lapply(input_files_per_var, function(x) x$model)
model_names <- unlist(unname(model_names))[anomaly_files]

start_climatology <- lapply(input_files_per_var, function(x) x$start_year)
start_climatology <- c(unlist(unname(start_climatology))[climatology_files])[1]
end_climatology <- lapply(input_files_per_var, function(x) x$end_year)
end_climatology <- c(unlist(unname(end_climatology))[climatology_files])[1]

start_anomaly <- lapply(input_files_per_var, function(x) x$start_year)
start_anomaly <- c(unlist(unname(start_anomaly))[anomaly_files])[1]
end_anomaly <- lapply(input_files_per_var, function(x) x$end_year)
end_anomaly <- c(unlist(unname(end_anomaly))[anomaly_files])[1]

print(start_climatology)
print(end_anomaly)

## Do not print warnings
#options(warn=-1)

#Parameters for Season() function
monini <- params$monini
moninf <- params$moninf
monsup <- params$monsup

#Threshold for ensemble agreement
agreement_threshold <- params$agreement_threshold




### Load data and compute climatologies and anomalies

climatology_filenames <- fullpath_filenames[climatology_files]
reference_data <- Start(model = climatology_filenames,
                        var = var0,
                        var_var = 'var_names',
                        time = 'all',
                        lat = 'all',
                        lon = 'all',
                        lon_var = 'lon',
                        return_vars = list(time = 'model', lon = 'model', lat = 'model'),
                        retrieve = TRUE)

lat <- attr(reference_data, "Variables")$dat1$lat
lon <- attr(reference_data, "Variables")$dat1$lon

months <- paste0(month.abb[moninf],"-", month.abb[monsup])

attributes(lon) <- NULL
attributes(lat) <- NULL
dim(lon) <-  c(lon = length(lon))
dim(lat) <- c(lat = length(lat))

time_dim <- which(names(dim(reference_data)) == "time")
dims <- dim(reference_data)

dims <- append(dims, c(12, dims[time_dim] / 12), after = time_dim)
dims <- dims[-time_dim]
dim(reference_data) <- dims
names(dim(reference_data))[c(time_dim, time_dim + 1)] <- c("month", "year")

margins <- list(c(1 : length(dim(reference_data)))[-c(time_dim + 1)])
reference_seasonal_mean <- Season(reference_data, posdim = time_dim, monini = monini, moninf = moninf,
                                   monsup = monsup)

margins <- list(c(1 : length(dim(reference_seasonal_mean)))[-c(time_dim + 1)])
years_dim <- which(names(dim(reference_seasonal_mean)) == "year")
climatology <- Mean1Dim(reference_seasonal_mean, years_dim)

anomaly_filenames <- fullpath_filenames[anomaly_files]
rcp_data <- Start(model = anomaly_filenames,
                  var = var0,
                  var_var = 'var_names',
                  time = 'all',
                  lat = 'all',
                  lon = 'all',
                  lon_var = 'lon',
                  lon_reorder = CircularSort(0, 360),
                  return_vars = list(time = 'model', lon = 'model',
                                     lat = 'model'),
                  retrieve = TRUE)


time_dim <- which(names(dim(rcp_data)) == "time")
dims <- dim(rcp_data)
dims <- append(dims, c(12, dims[time_dim] / 12), after = time_dim)
dims <- dims[-time_dim]
dim(rcp_data) <- dims
names(dim(rcp_data))[c(time_dim, time_dim + 1)] <- c("month", "year")
#
attr <- ((attr(rcp_data,"Variables")$common))[2]
units <- attr[[1]]$units


proj_seasonal_mean <- Season(rcp_data, posdim = time_dim, monini = monini, moninf = moninf,
                             monsup = monsup)
#
years_dim <- which(names(dim(proj_seasonal_mean)) == "year")
climatology <- InsertDim(climatology, years_dim, lendim = dim(proj_seasonal_mean)[years_dim])
anomaly <- proj_seasonal_mean - climatology

time <- seq(start_anomaly, end_anomaly, by = 1)
month <- moninf + (moninf - monsup)
if (month <= 9) {
  month <- paste0(as.character(0), as.character(month))
}
month <- paste0("-", month, "-")
day <- "01"
time <- as.POSIXct(paste0(time, month, day), tz = "CET")
time <- julian(time, origin = as.POSIXct("1970-01-01"))


attributes(time) <- NULL
dim(time) <- c(time = length(time))
metadata <- list(time = list(standard_name = 'time', long_name = 'time', units = 'days since 1970-01-01 00:00:00', prec = 'double', dim = list(list(name='time', unlim = FALSE))))
attr(time, "variables") <- metadata
#Save the single model anomalies
for (mod in 1 : length(model_names)) {
  data <- anomaly[mod,1,1, , ,]
  data <- aperm(data, c(2,3,1))
  names(dim(data)) <- c("lat", "lon", "time")
  metadata <- list(variable = list(dim = list(list(name='time', unlim = FALSE))))
  attr(data, 'variables') <- metadata
  ArrayToNetCDF(list(variable = data, lat = lat, lon = lon, time = time),
                paste0(plot_dir,  "/", var0, "_", months, "_anomaly_",model_names[mod],"_", start_anomaly, "_", end_anomaly,"_", start_climatology, "_", end_climatology, ".nc"))
}


#Compute and save the multi-model anomalies
multi_model_anomaly <- Mean1Dim(anomaly, 1)
data <- multi_model_anomaly[1,1, , ,]
data <- aperm(data, c(2,3,1))
names(dim(data)) <- c("lat", "lon", "time")
metadata <- list(variable = list(dim = list(list(name='time', unlim = FALSE))))
attr(data, 'variables') <- metadata
model_names_filename <- paste(model_names, collapse = '_')
ArrayToNetCDF(list(variable = data, lon = lon, lat = lat, time = time),
              paste0(plot_dir, "/", var0, "_",months, "_multimodel-anomaly_", model_names_filename,"_", start_anomaly, "_", end_anomaly,"_", start_climatology, "_", end_climatology, ".nc"))


##Plots
data <- Mean1Dim(data,3)
colorbar_lim <- max(abs(max(data)),abs(min(data)))
brks <- seq(-colorbar_lim, colorbar_lim, length.out = 11)

title <- paste0(months, " ", var0, " anomaly (", end_anomaly, "-", start_anomaly,
       ") - (", end_climatology, "-", start_climatology, ")")
PlotEquiMap(data, lat = lat, lon = lon, brks = brks, units =units, toptitle = title, filled.continents = FALSE,
            fileout = paste0(plot_dir, '/anomaly_map.png'))
