##############################################################################################################################
# script to make maps of North American glacial front after the LGM
# also plots estimated shoreline to show sea level rise
# author: Dan MacGuigan
# dmacguig@buffalo.edu
##############################################################################################################################

library(RColorBrewer)
library(maps)
library(rgdal)
library(gridExtra)
library(maptools)
library(scales)
library(raster)
library(rgeos)
library(ggplot2)
library(grid)
library(mapplots)
library(dplyr)    
library(mapproj)
library(ggmap)
library(ggspatial)
library(sf)
library(ggrepel)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(cowplot)
library(shades)
library(gtools)
library(stringr)
library(sdmpredictors)
library(stars)
library(spatialEco)
library(smoothr)

##############################################################################################################################
# specify the following variables
##############################################################################################################################

wd <- "H:/NAGlacialFront/NorthAmericaLGMIceRetreat" # top level working directory
glacialShapefileName <- list.files(path=paste(wd, "/data/Dalton_etAl_2020_QuatSci_NAGlaciers/", sep=""), pattern="*.shp")
glacialShapefileName <- mixedsort(glacialShapefileName) # natural sort files
glacialTime <- str_remove(sapply(strsplit(glacialShapefileName, "_"), `[`, 2), "cal") # get dates
glacialTime_title <- formatC(as.numeric(glacialTime) * 1000, format="d", big.mark=",") # get dates in years
glacialTime_file <- as.numeric(glacialTime) * 1000
projection <- CRS("+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs ") # map projection, North America Lambert Conformal Conic
                                                              # see https://proj.org/operations/projections/lcc.html#parameters for more projections

setwd(wd)

##############################################################################################################################
# READ DATA
##############################################################################################################################

# RIVERS
# load shapefile for all US streams at 1:10 million scale
# data from https://www.naturalearthdata.com/downloads/10m-physical-vectors/
# I modified the below shapefile slightly to fix a problem with the Mississippi River delta
rivers <- readOGR(paste(wd, "/data/naturalEarthShapefiles/ne_10m_rivers_lake_centerlines_scale_rank.shp", sep=""))

# LAKES
lakes <- ne_download(scale = 10, type = 'lakes', category = 'physical')
# below code fixes some bad polygons in this shapefile
# simplify the polgons a tad (tweak 0.00001 to your liking)
lakes <- gSimplify(lakes, tol = 0.00001)
# this is a well known R / GEOS hack (usually combined with the above) to 
# deal with "bad" polygons
lakes <- gBuffer(lakes, byid=TRUE, width=0)
# any bad polys?
sum(gIsValid(lakes, byid=TRUE)==FALSE)

# BORDERS
borders <- ne_countries(type = 'countries', scale = 'medium')

# BATHEMETIC DATA
# from https://www.bio-oracle.org/code.php
# with help from https://medium.com/themarinedatascientist/producing-a-landmass-shapefile-of-the-last-glacial-maximum-in-r-1fb38b12b614
# I manually reprojected this raster layer with ArcGIS
bathymetry <- raster(paste(wd, "./data/BathymetryDepthMean_arcgis_reproj.tif", sep=""))

# SEA LEVEL DATA
# from Lambeck et al. 2014
# https://doi.org/10.1073/pnas.1411762111
# Table S3
sea_level <- read.csv(paste(wd, "./data/Lambeck-etAl_2014_PNAS_TableS3.csv", sep=""), header=TRUE)

###############################################################################################################################
# PLOT MAP
##############################################################################################################################

# convert to sf objects
rivers_sf <- st_as_sf(rivers)
lakes_sf <- st_as_sf(lakes)
borders_sf <- st_as_sf(borders)

# dummy df for adding annotation to map 
dum_df <- data.frame(Longitude=c(-120), Latitude=c(15)) # fiddle with these to get the annotaiton in the right spot
dum_df_sf <- st_as_sf(dum_df, coords=c("Longitude", "Latitude"), crs=4326, agr="constant") # convert df to sf object
dum_df_sf_reproj <- st_transform(dum_df_sf, crs = projection)
long <- st_coordinates(dum_df_sf_reproj)[1]
lat <- st_coordinates(dum_df_sf_reproj)[2]

# set extents for plot, play around with these values until the plot looks nice
# the bouding box polygon in long/lat projection, i.e. axis-aligned for plot
bb <- st_sfc(
  st_polygon(list(cbind(
    c(-130, -130, -60, -60, -130), # x-coordinates (longitudes) of points A,B,C,D
    c(10, 85, 85, 10, 10)     # y-coordinates (latitudes) of points A,B,C,D
  ))),
  crs = "+proj=longlat +datum=WGS84 +no_defs")

# now in in LAEA projection
laeabb <- st_transform(bb, crs = projection)

# the extent of the bounding box in the new projection
b <- st_bbox(laeabb)

for(i in 1:length(glacialShapefileName)){
  # read in and process glacial data
  # data from https://www.sciencedirect.com/science/article/abs/pii/S0277379119307619, appexdix C
  glacier <- readOGR(paste(wd, "/data/Dalton_etAl_2020_QuatSci_NAGlaciers/", glacialShapefileName[i], sep=""))
  glacier_sf <- st_as_sf(glacier)
  
  # get shoreline for the historical time period
  time <- as.numeric(glacialTime[i])
  sl <- sea_level$esl_m[which.min(abs(sea_level$time_ka-time))] # get the sea level for the closest time point in dataset
  # prep bathymetry data
  bathymetry_df <- as.data.frame(bathymetry, xy=TRUE) # convert to data frame
  bathymetry_df$BathymetryDepthMean_arcgis_reproj <- bathymetry_df$BathymetryDepthMean_arcgis_reproj - sl # adjust sea level for given time point
  bathymetry_df$BathymetryDepthMean_arcgis_reproj[bathymetry_df$BathymetryDepthMean_arcgis_reproj >= 0] <- NA # NA for cells with sea level > 0
  bathymetry_df$BathymetryDepthMean_arcgis_reproj <- -(bathymetry_df$BathymetryDepthMean_arcgis_reproj) # convert depth to positive values for plotting

  # make plot
  plot <- ggplot() +
    geom_raster(bathymetry_df, mapping = aes(x = x, y = y, fill=BathymetryDepthMean_arcgis_reproj)) +
    scale_fill_gradient(na.value="gray20", low="#000d42", high="#328aa8", name="ocean\ndepth\n(km)",
                        labels = c("0", "2", "4", "6", "8", "10"),
                        breaks = c(0, 2000, 4000, 6000, 8000, 10000),
                        trans = 'reverse') +
    geom_sf(data=borders_sf, fill="black", color=NA, size=0.2) +
    geom_sf(data=rivers_sf, color="#36c3ff", lineend = "round", aes(size=factor(strokeweig)), show.legend = FALSE,) +
    geom_sf(data=lakes_sf, fill="#36c3ff", color=NA) +
    geom_sf(data=glacier_sf, fill=alpha("white", 0.7), color="white", size=1) +
    scale_size_manual(values=seq(0.1,20,by=0.05)) +
    labs(title = paste(glacialTime_title[i], " years before present", sep="")) +
    annotate("text", x=long, y=lat, label=paste("sea level ", formatC(signif(-sl, 4), format='f', digits=1), " meters\nlower than present-day", sep=""), 
             color="white", hjust=0.5, vjust=0.5, size=5) +
    coord_sf(crs = projection, xlim = c(b["xmin"], b["xmax"]), ylim = c(b["ymin"], b["ymax"])) +
    ylab("") +
    xlab("") +
    theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.7), 
          panel.background = element_rect(fill = "#1c5163"),
          plot.title = element_text(color = "black", size = 24, hjust = 0.5),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank())
  
  # save each individual plot
  ggsave(paste(wd, "/plots/calibratedDates_seaLevel/", glacialTime_file[i], ".png", sep=""), plot=plot, units = "in", dpi=300, height=7, width=8)

}

for(i in 1:length(glacialShapefileName)){
  # read in and process glacial data
  # data from https://www.sciencedirect.com/science/article/abs/pii/S0277379119307619, appexdix C
  glacier <- readOGR(paste(wd, "/data/Dalton_etAl_2020_QuatSci_NAGlaciers/", glacialShapefileName[i], sep=""))
  glacier_sf <- st_as_sf(glacier)
  
  # get shoreline for the historical time period
  time <- as.numeric(glacialTime[i])
  sl <- sea_level$esl_m[which.min(abs(sea_level$time_ka-time))] # get the sea level for the closest time point in dataset
  sl <- sl * 3.28084 # convert to feet
  # prep bathymetry data
  bathymetry_df <- as.data.frame(bathymetry, xy=TRUE) # convert to data frame
  bathymetry_df$BathymetryDepthMean_arcgis_reproj <- bathymetry_df$BathymetryDepthMean_arcgis_reproj * 3.28084 # convert to feet
  bathymetry_df$BathymetryDepthMean_arcgis_reproj <- bathymetry_df$BathymetryDepthMean_arcgis_reproj - sl # adjust sea level for given time point
  bathymetry_df$BathymetryDepthMean_arcgis_reproj[bathymetry_df$BathymetryDepthMean_arcgis_reproj >= 0] <- NA # NA for cells with sea level > 0
  bathymetry_df$BathymetryDepthMean_arcgis_reproj <- -(bathymetry_df$BathymetryDepthMean_arcgis_reproj) # convert depth to positive values for plotting
  bathymetry_df$BathymetryDepthMean_arcgis_reproj <- bathymetry_df$BathymetryDepthMean_arcgis_reproj / 5280 # convert depth to miles for plotting
  
  # make plot
  plot <- ggplot() +
    geom_raster(bathymetry_df, mapping = aes(x = x, y = y, fill=BathymetryDepthMean_arcgis_reproj)) +
    scale_fill_gradient(na.value="gray20", low="#000d42", high="#328aa8", name="ocean\ndepth\n(miles)",
                        breaks = c(0, 1, 2, 3, 4, 5, 6, 7),
                        trans = 'reverse') +
    geom_sf(data=borders_sf, fill="black", color=NA, size=0.2) +
    geom_sf(data=rivers_sf, color="#36c3ff", lineend = "round", aes(size=factor(strokeweig)), show.legend = FALSE,) +
    geom_sf(data=lakes_sf, fill="#36c3ff", color=NA) +
    geom_sf(data=glacier_sf, fill=alpha("white", 0.7), color="white", size=1) +
    scale_size_manual(values=seq(0.1,20,by=0.05)) +
    labs(title = paste(glacialTime_title[i], " years before present", sep="")) +
    annotate("text", x=long, y=lat, label=paste("sea level ", formatC(signif(-sl, 4), format='f', digits=1), " feet\nlower than present-day", sep=""), 
             color="white", hjust=0.5, vjust=0.5, size=5) +
    coord_sf(crs = projection, xlim = c(b["xmin"], b["xmax"]), ylim = c(b["ymin"], b["ymax"])) +
    ylab("") +
    xlab("") +
    theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.7), 
          panel.background = element_rect(fill = "#1c5163"),
          plot.title = element_text(color = "black", size = 24, hjust = 0.5),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank())
  
  # save each individual plot
  ggsave(paste(wd, "/plots/calibratedDates_seaLevel_imperial/", glacialTime_file[i], ".png", sep=""), plot=plot, units = "in", dpi=300, height=7, width=8)
  
}
