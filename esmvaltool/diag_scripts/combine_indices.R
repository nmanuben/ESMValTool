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


#Region considered
lat.max <- params$lat_max
lat.min <- params$lat_min
lon.max <- params$lon_max
lon.min <- params$lon_min


#Start and end periods for the historical and projection periods
startdate <- as.POSIXct(params$startdate)
enddate <- as.POSIXct(params$enddate)
first_year <- substr(startdate, 1, 4)
last_year <- substr(enddate, 1, 4)

#Regime parameters
metrics <- params$metric
rcp_scenario <- "rcp8.5"
weights <- params$weights

weight_name <- c() 
for (i in 1 : length(weights)) {
  if (!is.null(weights)) {
    weight_name[i] <- paste0("_", as.character(weights[i]), metrics[i])  
  } else {
  weight_name[i] <- ""
  }
}
#Plot parameters
font_size <- 12

library(startR)
library(s2dverification)
source("https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Magic_WP6/R/CombineIndices.R")
source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Magic_WP6/R/WeightedMean.R')

fullpath_filenames <- input_files[[var0]]

var0 <- "index"
data <- Start(model = fullpath_filenames,
              var = "index",
              var_var = 'var_names',
              #  sdate = paste0(seq(1970, 2000, by = 1), "0101"),
              time = values(list(as.POSIXct(startdate), 
                                 as.POSIXct(enddate))),
              time_tolerance = as.difftime(180, units = 'days'), 
              lat = values(list(lat.min, lat.max)),
              lon = values(list(lon.min, lon.max)),
              lon_var = 'lon',
              lon_reorder = CircularSort(0, 360),
              return_vars = list(time = 'model', lon = 'model', lat = 'model'),
              retrieve = TRUE)

indices_dim <- which(names(dim(data)) == "model")
indices <- list()
for (i in 1 : dim(data)[indices_dim]) {
  indices[[i]] <- Subset(data, along = indices_dim, indices = i)
}
if (!is.null(combine)) {
  aci <- CombineIndices(indices, weights = weights)
  
  time <- attributes(data)$Variables$dat1$time
  time <- as.POSIXct(paste0(time, month, day), tz = "CET")
  time <- julian(time, origin = as.POSIXct("1970-01-01"))
  
  attributes(time) <- NULL
  dim(time) <- c(time = length(time))
  metadata <- list(time = list(standard_name = 'time', long_name = 'time', units = 'days since 1970-01-01 00:00:00', prec = 'double', dim = list(list(name='time', unlim = FALSE))))
  attr(time, 'variables') <- metadata
  aci <- aci[1,1, , ,]
  aci <- aperm(aci, c(3,2,1))
  names(dim(aci)) <- c("lon", "lat", "time")
  metadata <- list(index = list(dim = list(list(name='time', unlim = FALSE, prec = 'double'))))
  attr(aci, 'variables') <- metadata
  lat <- attr(data, "Variables")$dat1$lat
  lon <- attr(data, "Variables")$dat1$lon
  attributes(lon) <- NULL
  attributes(lat) <- NULL
  dim(lon) <-  c(lon = length(lon))
  dim(lat) <- c(lat = length(lat))
  ArrayToNetCDF(list(metric= aci, time = time,  lat = lat, lon = lon), 
                paste0("aci", "_",model_name,"_", rcp_scenario, "_", startdate, "_", enddate, "_",paste(weight_name, collapse = '') , ".nc"))
}

if (!is.null(area_weights)) {
  area_mean <- list()
  for (i in 1 : length(indices)) {
    area_mean[[i]] <- WeightedMean(indices[[i]], lon = lon, lat = lat)
    attributes(area_mean[[i]]) <- NULL
    dim(area_mean[[i]]) <- c(area_mean = length(area_mean[[i]]))
    ArrayToNetCDF(list(metric= area_mean[[i]], time = time), 
                  paste0(metric[i], "_area_mean_",model_name,"_", rcp_scenario, "_", startdate, "_", enddate, "_", weight_name[i], ".nc"))
    
  }
}

if (!is.null(area_weights) && !is.null(combine)) {
  combined_area_average <-  CombineIndices(area_mean, weights = weights)
  attributes(combined_area_average) <- NULL
  dim(combined_area_average) <- c(combined_area_average = length(combined_area_average))
  ArrayToNetCDF(list(aci= combined_area_average, time = time), 
                paste0("aci", "_area_mean_",model_name,"_", rcp_scenario, "_", startdate, "_", enddate, "_", paste(weight_name, collapse = ''), ".nc"))
  years = as.numeric(first_year) : as.numeric(last_year)
  data_frame <- data.frame(aci = combined_area_average, cdd = area_mean[[1]], 
                           rx5day = area_mean[[2]], tx90p = area_mean[[3]], tx10p = area_mean[[4]],
                           wx = area_mean[[5]], years = years)
 
  
  xymelt <- melt(data_frame, id.vars = "years") 
  

  g <- ggplot(xymelt, aes(x = years, y = value, color = variable, linetype = variable, size = variable)) +
    theme_bw() +
    geom_line(aes(size = variable)) + ylab("Component") +  xlab("Year") + theme(text=element_text(size = font_size),legend.text=element_text(size = font_size),
                                                      axis.title=element_text(size = font_size)) + #scale_colour_manual(values=c("red", "black", "black", "black", "black", "black")) + 
   scale_size_manual(values = c(1, 0.5, 0.5, 0.5, 0.5, 0.5)) +  ggtitle(paste0("ACI components (",  model_name, " " ,rcp_scenario, " ", paste(weight_name, collapse = ''), ")"))
  ggsave(filename = paste0("aci", "_area_mean_",model_name,"_", rcp_scenario, "_", startdate, "_", enddate, "_", paste(weight_name, collapse = ''), ".pdf"), g, device = NULL)
}

###Plots




