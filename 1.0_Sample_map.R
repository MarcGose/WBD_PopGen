# load packages

library(maps)
library(mapdata)
library(ggplot2)
library(dplyr)
library(rgdal)
library(raster)
library(ggsn)
library(rworldmap)
library(mapproj)
library(nord)

# define colour palette

col_palette <- c(nord("mountain_forms")[1],
                 nord("frost")[2],
                 nord("aurora")[4],
                 nord("aurora")[1],
                 nord("moose_pond")[3])


# load metadata with pop assinment, lat, long, ID

locations <- read.table("file_lists/meta.txt", header = T)

# load world map

world_map <- map_data("world", returnclass = "sf") %>%
  as_tibble() 

# plot

(sample_map <- ggplot() +
    geom_map(data = world_map, aes(map_id = region), map = world_map, fill = "lightgrey", inherit.aes = FALSE) +
    geom_point(data = locations, aes(x = Long, y = Lat, color = POPFINE), size = 2, alpha = 0.7) +
    scale_color_manual(values = col_palette) +
    coord_sf(crs = "+proj=moll") +
    expand_limits(x = world_map$long, y = world_map$lat) +
    theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white")) +
    coord_cartesian(xlim = c(-80, 50), ylim = c(45, 80)))

# change map projection

sample_map_ortho <- sample_map +
  coord_map("ortho", orientation = c(55, -20, 0))

# save figure

ggsave("figures/Figure_1a.pdf", sample_map_ortho, height = 16.5, width = 21, units = "cm")
