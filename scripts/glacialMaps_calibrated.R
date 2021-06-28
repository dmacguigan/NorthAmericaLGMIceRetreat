##############################################################################################################################
# script to make maps of North American glacial front after the LGM
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

##############################################################################################################################
# specify the following variables
##############################################################################################################################

wd <- "H:/NAGlacialFront/NorthAmericaLGMIceRetreat" # top level working directory
glacialShapefileName <- list.files(path=paste(wd, "/data/Dalton_etAl_2020_QuatSci_NAGlaciers/", sep=""), pattern="*.shp")
glacialShapefileName <- mixedsort(glacialShapefileName) # natural sort files
glacialTime <- str_remove(sapply(strsplit(glacialShapefileName, "_"), `[`, 2), "cal") # get dates
glacialTime_title <- formatC(as.numeric(glacialTime) * 1000, format="d", big.mark=",") # get dates in years
glacialTime_file <- as.numeric(glacialTime) * 1000
projection <- CRS("+proj=lcc +lat_1=30 +lat_2=60 +lon_0=-100") # map projection, North America Lambert Conformal Conic
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

###############################################################################################################################
# PLOT MAP
##############################################################################################################################

# convert to sf objects
rivers_sf <- st_as_sf(rivers)
lakes_sf <- st_as_sf(lakes)
borders_sf <- st_as_sf(borders)

# set extents for plot
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

# main plot list
plot_list <- list()
#for(i in 1:length(glacialTime)){
for(i in 1:length(glacialShapefileName)){
  # read in and process glacial data
  # data from https://www.sciencedirect.com/science/article/abs/pii/S0277379119307619, appexdix C
  glacier <- readOGR(paste(wd, "./data/Dalton_etAl_2020_QuatSci_NAGlaciers/", glacialShapefileName[i], sep=""))
  glacier_sf <- st_as_sf(glacier)
  
  # make plot
  plot_list[[i]] <- ggplot() +
    geom_sf(data=borders_sf, fill="black", color=NA, size=0.2) +
    geom_sf(data=rivers_sf, color="#36c3ff", lineend = "round", aes(size=factor(strokeweig)), show.legend = FALSE,) +
    geom_sf(data=lakes_sf, fill="#36c3ff", color=NA) +
    geom_sf(data=glacier_sf, fill=alpha("white", 0.7), color="white", size=1) +
    scale_size_manual(values=seq(0.1,20,by=0.05)) +
    labs(title = paste(glacialTime_title[i], " years before present", sep="")) +
    coord_sf(crs = projection, xlim = c(b["xmin"], b["xmax"]), ylim = c(b["ymin"], b["ymax"])) +
    ylab("") +
    xlab("") +
    theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.7), 
          panel.background = element_rect(fill = "#1c5163"),
          plot.title = element_text(color = "black", size = 24, hjust = 0.5),
          axis.text.y = element_text(color = "black", size = 10),
          axis.text.x = element_text(color = "black", size = 10))
  
  # save each individual plot
  ggsave(paste(wd, "/plots/calibratedDates/", glacialTime_file[i], ".png", sep=""), plot=plot_list[[i]], units = "in", dpi=300, height=8, width=8)
}
