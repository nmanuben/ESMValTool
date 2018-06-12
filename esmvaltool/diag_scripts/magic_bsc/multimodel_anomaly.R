####REQUIRED SYSTEM LIBS
####ŀibssl-dev
####libnecdf-dev
####cdo

# conda install -c conda-forge r-ncdf4

#R package dependencies installation script
#install.packages('yaml')
#install.packages('devtools')
#library(devtools)
#Sys.setenv(TAR = '/bin/tar')
#install_git('https://earth.bsc.es/gitlab/es/startR', branch = 'develop-hotfixes-0.0.2')
#install_git('https://earth.bsc.es/gitlab/es/easyNCDF', branch = 'master')


#Parsing input file paths and creating output dirs
#args <- c('/home/Earth/nmanuben/esmvaltool_output/namelist_anomaly_agreement_20180302_135018/run/anomaly_agreement/main/settings.yml')


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

#Var considered
var0 <- var_names[1]

#Region considered
lat.max <- params$lat_max
lat.min <- params$lat_min
lon.max <- params$lon_max
lon.min <- params$lon_min

#Start and end periods for the reference and projection periods
start_reference <- as.POSIXct(params$start_reference)
end_reference <- as.POSIXct(params$end_reference)
start_projection <- as.POSIXct(params$start_reference)
end_projection <- as.POSIXct(params$end_reference)

#Parameters for Season() function
monini <- params$monini
moninf <- params$moninf
monsup <- params$monsup

#Threshold for ensemble agreement
#agreement_threshold <- params$agreement_threshold


### Loading and computation begins from here

library(s2dverification)
library(startR)
library(multiApply)

fullpath_filenames <- input_files[[var0]]

#Load the reference period
reference_data <- Start(model = fullpath_filenames_reference,
                         var = var0,
                         var_var = 'var_names',
                         #  sdate = paste0(seq(1970, 2000, by = 1), "0101"),
                         time = values(list(as.POSIXct(start_reference), 
                                            as.POSIXct(end_reference))),
                         time_tolerance = as.difftime(15, units = 'days'), 
                         lat = values(list(lat.min, lat.max)),
                         lon = values(list(lon.min, lon.max)),
                         lon_var = 'lon',
                         lon_reorder = CircularSort(0, 360),
                         return_vars = list(time = 'model', lon = 'model', lat = 'model'),
                         retrieve = TRUE)
units <- attr(reference_data, 'Variables')$common[[var0]]$units
lat <- attr(reference_data, "Variables")$dat1$lat
lon <- attr(reference_data, "Variables")$dat1$lon

months <- paste0(month.abb[moninf],"-", month.abb[monsup])

time_dimension <- which(names(dim(reference_data)) == "time")
lon[lon > 180] <- lon[lon > 180] - 360
lon_order <- sort(lon, index.return = TRUE)
reference_data <- Subset(reference_data, "lon", lon_order$ix)
lon <- lon_order$x

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


rcp_data <- Start(model = fullpath_filenames_projection,
                  var = var0,
                  var_var = 'var_names',
                  time = values(list(as.POSIXct(start_projection),
                                     as.POSIXct(end_projection))),
                  time_tolerance = as.difftime(15, units = 'days'),
                  lat = values(list(lat.min, lat.max)),
                  lon = values(list(lon.min, lon.max)),
                  lon_var = 'lon',
                  lon_reorder = CircularSort(0, 360),
                  return_vars = list(time = 'model', lon = 'model',
                                     lat = 'model'),
                  retrieve = TRUE)

## Reorder the longitudes
rcp_data <- Subset(rcp_data, "lon", lon_order$ix)


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

time <- c(substr(start_projection, 1, 4), substr(end_projection, 1, 4))
time <- as.numeric(time[1]) : as.numeric(time[2])
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
                paste0(plot_dir, var0, "_", months, "_anomaly_",model_names[mod],"_", start_projection, "_", end_projection,"_", start_reference, "_", end_reference, ".nc"))
}

#Compute and save the multi-model anomalies
multi_model_anomaly <- Mean1Dim(anomaly, 1)
data <- multi_model_anomaly[1,1, , ,]
data <- aperm(data, c(2,3,1))
names(dim(data)) <- c("lat", "lon", "time")
metadata <- list(variable = list(dim = list(list(name='time', unlim = FALSE))))
attr(data, 'variables') <- metadata
model_names_filename <- paste(model_names, collapse = '_')
ArrayToNetCDF(list(variable = data, lat = lat, lon = lon, time = time), 
              paste0(plot_dir, var0, "_",months, "_multimodel-anomaly_", model_names_filename,"_", start_projection, "_", end_projection,"_", start_reference, "_", end_reference, ".nc"))


##Plots 
data <- Mean1Dim(data,3)
colorbar_lim <- max(abs(max(data)),abs(min(data)))
brks <- seq(-colorbar_lim, colorbar_lim, length.out = 11)


title <- paste0(months, " ", var0, " anomaly (", substr(end_projection, 1, 4), "-", substr(start_projection, 1, 4), 
       ") - (", substr(end_reference, 1, 4), "-", substr(start_reference, 1, 4), ")")
PlotEquiMap(data, lat = lat, lon = lon, brks = brks, units =units, toptitle = title, filled.continents = FALSE,
            fileout = paste0(plot_dir, '/anomaly_map.png'))





