####REQUIRED SYSTEM LIBS
####ŀibssl-dev
####libnecdf-dev
####cdo

#R package dependencies installation script
#install.packages('yaml')
#install.packages('devtools')
#library(devtools)
#Sys.setenv(TAR = '/bin/tar')
#install_git('https://earth.bsc.es/gitlab/es/s2dverification', branch = 'production')
#install_git('https://earth.bsc.es/gitlab/ces/multiApply', branch = 'develop-copyMargins')
#install_git('https://earth.bsc.es/gitlab/es/startR', branch = 'develop-chunking')
#install_git('https://earth.bsc.es/gitlab/es/easyNCDF', branch = 'master')





#Parsing input file paths and creating output dirs
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

#Parameters for Season() function
monini <- params$monini
moninf <- params$moninf
monsup <- params$monsup

#Threshold for ensemble agreement
agreement_threshold <- params$agreement_threshold





library(s2dverification)
library(startR)
library(ggplot2)
source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-MagicWP5/R/AnoAgree.R')
source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Magic_WP6/R/WeightedMean.R')
source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-debug-plot-ts/R/PlotTimeSeries.R')

fullpath_filenames <- input_files[[var0]]
historical_data <- Start(model = fullpath_filenames,
                         var = var0,
                         var_var = 'var_names',
                         #  sdate = paste0(seq(1970, 2000, by = 1), "0101"),
                         time = "all", 
                         lat = values(list(lat.min, lat.max)),
                         lon = values(list(lon.min, lon.max)),
                         lon_var = 'lon',
                         lon_reorder = CircularSort(0, 360),
                         return_vars = list(time = 'model', lon = 'model', lat = 'model'),
                         retrieve = FALSE)

forecast_time <- attr(historical_data, 'Variables')$dat1$time
time_1 <- which(substr(forecast_time, 1, 7) == substr(start_historical, 1, 7))
time_2 <- which(substr(forecast_time, 1, 7) == substr(end_historical, 1, 7))

historical_data <- Start(model = fullpath_filenames,
                         var = var0,
                         var_var = 'var_names',
                         #  sdate = paste0(seq(1970, 2000, by = 1), "0101"),
                         time = indices(time_1 : time_2),
                         lat = values(list(lat.min, lat.max)),
                         lon = values(list(lon.min, lon.max)),
                         lon_var = 'lon',
                         lon_reorder = CircularSort(0, 360),
                         return_vars = list(time = 'model', lon = 'model', lat = 'model'),
                         retrieve = TRUE)
#units <- attr(historical_data, 'Variables')[[var0]]$units
#print("AAAAAAAAAA")
#print(var0)
#print(str(attributes(historical_data)))

#londim <- which(names(dim(historical_data)) == 'lon') 
#latdim <- which(names(dim(historical_data)) == 'lat')
lat <- attr(historical_data, "Variables")$dat1$lat
lon <- attr(historical_data, "Variables")$dat1$lon
historical_data <- historical_data * 24 * 60 * 60
area_mean <- WeightedMean(historical_data, lon = lon, lat = lat)
names(dim(area_mean)) <- c('model', 'var', 'time')

a <- PlotTimeSeries(area_mean, 
               coords = list(model = model_names[[var0]], 
                             time = attr(historical_data, 'Variables')$dat1$time), 
               colour_along = 'model', 
               title = "RCP8.5: Monthly mean anomalies",
               title_sub = paste0("Projection period: (", substr(end_projection, 1, 4), "-", substr(start_projection, 1, 4), 
                                  ")\nReference period: (", substr(end_historical, 1, 4), "-", substr(start_historical, 1, 4), ")"), 
               title_y = paste0("Anomaly (", as.character(units), ")"), 
               title_x = "Year") + scale_x_datetime()
               #file_width = 12, file_height = 8, file_name = paste0(plot_dir, '/test.png'))
ggsave(width = 12, height = 8, filename = paste0(plot_dir, '/test.png'), plot = a)




  time_dim <- which(names(dim(historical_data)) == "time")
 dims <- dim(historical_data)
 dims <- append(dims, c(12, dims[time_dim] / 12), after = time_dim)
 dims <- dims[-time_dim]
 dim(historical_data) <- dims
 names(dim(historical_data))[c(time_dim, time_dim + 1)] <- c("month", "year")
 
# #Calculate the historical mean for each month
# 
 margins <- list(c(1 : length(dim(historical_data)))[-c(time_dim + 1)])
 historical_seasonal_mean <- Season(historical_data, posdim = time_dim, monini = monini, moninf = moninf,
                                    monsup = monsup)
# 
 margins <- list(c(1 : length(dim(historical_seasonal_mean)))[-c(time_dim + 1)])
 years_dim <- which(names(dim(historical_seasonal_mean)) == "year")
 climatology <- Mean1Dim(historical_seasonal_mean, years_dim)
# 
 rcp_data <- Start(dat = fullpath_filenames,
                   var = var0,
                   var_var = 'var_names',
                   time = "all", #values(list(first, second)),
                   lat = values(list(lat.min, lat.max)),
                   lon = values(list(lon.min, lon.max)),
                   lon_var = 'lon',
                   lon_reorder = CircularSort(0, 360),
                   return_vars = list(time = 'dat', lon = 'dat',
                                      lat = 'dat'),
                   retrieve = FALSE)
# 
 forecast_time <- attr(rcp_data, 'Variables')$dat1$time
 time_1 <- which(substr(forecast_time, 1, 7) == substr(start_projection, 1, 7))
 time_2 <- which(substr(forecast_time, 1, 7) == substr(end_projection, 1, 7))
# 
# 
 rcp_data <- Start(dat = fullpath_filenames,
                    var = var0,
                    var_var = 'var_names',
                   time = indices(time_1 : time_2), #values(list(first, second)),
                   lat = values(list(lat.min, lat.max)),
                   lon = values(list(lon.min, lon.max)),
                   lon_var = 'lon',
#                   #   time_var = 'time',
                   lon_reorder = CircularSort(0, 360),
                   return_vars = list(time = 'dat', lon = 'dat',
                                      lat = 'dat'),
                   retrieve = TRUE)
# 
# 
 time_dim <- which(names(dim(rcp_data)) == "time")
 dims <- dim(rcp_data)
 dims <- append(dims, c(12, dims[time_dim] / 12), after = time_dim)
 dims <- dims[-time_dim]
 dim(rcp_data) <- dims
 names(dim(rcp_data))[c(time_dim, time_dim + 1)] <- c("month", "year")
# 
 lat <- attr(rcp_data,"Variables")$dat1$lat
 lon <- attr(rcp_data,"Variables")$dat1$lon
# 
# 
 proj_seasonal_mean <- Season(rcp_data, posdim = time_dim, monini = monini, moninf = moninf,
                              monsup = monsup)
# 
 years_dim <- which(names(dim(proj_seasonal_mean)) == "year")
 climatology <- InsertDim(climatology, years_dim, lendim = dim(proj_seasonal_mean)[years_dim])
 anomaly <- proj_seasonal_mean - climatology
# 
# 
 multiyearmeananomaly <- Mean1Dim(anomaly, years_dim)
 dat_dim <- which(names(dim(multiyearmeananomaly)) == "dat")
 agreement <- AnoAgree(multiyearmeananomaly, members_dim = dat_dim)
 datmeananomaly <- Mean1Dim(multiyearmeananomaly,  dat_dim)
##

PlotEquiMap(datmeananomaly[1, 1, , ], lon = lon, lat = lat, color_fun = clim.palette("yellowred"), filled.continents = FALSE,
              dots = drop(agreement) > agreement_threshold, fileout = paste0(plot_dir, '/test2.png'))



####
####    figure_filename <- interface_get_figure_filename(diag_script_base,
####                                                 '',
####                                                 field_type0,
####                                                 '',
####                                                 model_idx)
####    figure_path <- file.path(plot_dir, diag_base, paste0(figure_filename, "map", '.', output_file_type))
####    info_output(figure_path, verbosity, 1)
####     units <- paste0("attr(historical_data, 'variables')[[1]]$",var0,"$units")
####  units <- eval(parse(text=units))
####  toptitle <- paste0("Anomaly RCP8.5:", "(", substr(end_projection, 1, 4), "-",
####                    substr(start_projection, 1, 4),")","-" , "(", substr(end_historical, 1, 4), "-", substr(start_historical, 1, 4) , ")")
####
####
####  PlotEquiMap(ensemblemeananomaly[1, 1, 1, , ], lon = lon, lat = lat, color_fun = clim.palette("yellowred"), filled.continents = FALSE,
####              dots = drop(agreement) > agreement_threshold, units = units, toptitle = toptitle,  fileout = figure_path)
####
####    londim <- which(names(dim(anomaly)) == "longitude")
####  latdim <- which(names(dim(anomaly)) == "latitude")
####  areamean_anomaly <- WeightedMean(data = anomaly, lon = lon, lat = lat, londim = londim, latdim = latdim)
####  #names(dim(areamean_anomaly)) <- names(dim(anomaly))[-c(londim,latdim)]
####  areamean_anomaly <- as.data.frame.table(areamean_anomaly)
####  names(areamean_anomaly) <- c(names(dim(anomaly))[-c(londim, latdim)], "values")
####  areamean_anomaly$year <- c(as.numeric(substr(start_projection, 1, 4)) : as.numeric(substr(end_projection, 1, 4)))
####
####  font_size <- 20
####  title_1 <- paste0("RCP 8.5:", "Monthly mean anomalies for ", "(", end_projection, "-", start_projection, ")", "-", "(", historical_end, "-", historical_start, ")")
####  time_series <- ggplot(data = areamean_anomaly, aes(x = year, y = values, group = ensemble, linetype = ensemble)) +
####    geom_line()   + ylab(paste("Anomaly RCP8.5:", "(", end_projection, "-", start_projection, ")", "(", end_historical, "-",start_historical , ")")) +  xlab("Year") + theme(text=element_text(size = font_size),legend.text=element_text(size = font_size),
####                                                                                                                                                                                              axis.title=element_text(size = font_size))
####  time_series <- ggplot() + stat_summary(data =  areamean_anomaly, fun.y= "mean", mapping = aes_string(x = areamean_anomaly$year, y = areamean_anomaly$values, group = areamean_anomaly$dat, color = areamean_anomaly$dat, fill = areamean_anomaly$dat),
####                                         geom = "line", size = 1, show.legend = FALSE) + ylab(paste("Anomaly", "(", as.character(units), ")")) +  xlab("Year") +
####    stat_summary(data = areamean_anomaly, geom = "ribbon", fun.ymin = "min", fun.ymax = "max", aes(x = areamean_anomaly$year, y = areamean_anomaly$values, group = areamean_anomaly$dat, color = areamean_anomaly$dat, fill = areamean_anomaly$dat),
####                 alpha = 0.3, show.legend = FALSE) + ggtitle(toptitle) +  theme(text=element_text(size = font_size),legend.text=element_text(size = font_size))
####
####
####  ### Time series plot
####    figure_path <- file.path(plot_dir, diag_base, paste0(figure_filename, "time_series", '.', output_file_type))
####  ggsave(time_series, device = NULL, path = figure_path_time_series)
