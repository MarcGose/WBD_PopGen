# Load packages

library(circlize)
library(tidyverse)
library(RColorBrewer)
library(nord)

# Define colour palettes

# broad

col_palette <- c(nord("mountain_forms")[1], #dark blue
                 nord("frost")[2], #light blue
                 nord("moose_pond")[3], #brown
                 nord("aurora")[4], #green
                 nord("aurora")[1]) #red

# fine

col_palette2 <- c(nord("mountain_forms")[1],
                 nord("aurora")[2],
                 nord("frost")[2],
                 nord("aurora")[3],
                 nord("aurora")[5] ,
                 nord("moose_pond")[3],
                 nord("moose_pond")[6],
                 nord("silver_mine")[2],
                 nord("aurora")[4],
                 nord("aurora")[1])

# Prepare BayesAss data

# Broad pop assignment

baysass <- tibble(Region = c("BAS","ICE","WSI","NS","WNA"),
                  BAS = c(0.689438,0.22145,0.023238,0.043625,0.022225),
                  ICE = c(0.014563,0.930063,0.014525,0.0265,0.014338),
                  WSI = c(0.015225,0.006475,0.675838,0.295988,0.006525),
                  NS = c(0.0112,0.014938,0.057688,0.911988,0.0042),
                  WNA = c(0.041525,0.166888,0.041263,0.041913,0.708413))

# Convert to matrix
baysass.mat = as.matrix(baysass[, 2:6])
baysass.mat

# Gene flow matrix
dimnames(baysass.mat) = list(source = baysass$Region, sink = baysass$Region)
baysass.mat

# Reset the circular layout parameters
circos.clear()

# Parameters for the circular layout
circos.par(start.degree = 30, gap.degree = 3)

# Plot chord diagram
chordDiagram(x = baysass.mat, grid.col = col_palette, grid.border = "black", transparency = 0.3,
             order = baysass$Region, directional = T, direction.type = "arrows",
             self.link = 1, preAllocateTracks = list(track.height = 0.01),
             annotationTrack = "grid", annotationTrackHeight = c(0.05, 0.05),
             link.border = "NA", link.sort = T, link.decreasing = T,
             link.arr.length = 0.4, link.arr.lty = 3, link.arr.col = "#252525", 
             link.largest.ontop = F)

chordDiagram(x = baysass.mat, grid.col = col_palette, transparency = 0.4,
             order = baysass$Region, directional = T, direction.type = "arrows",
             self.link = 1, preAllocateTracks = list(track.height = 0.05),
             annotationTrack = "grid", annotationTrackHeight = c(0.05, 0.05),
             link.border = "NA", link.sort = T, link.decreasing = T,
             link.arr.length = 0.2, link.arr.lty = 2, link.arr.col = "#252525", 
             link.largest.ontop = F)

BA3_broad <- circos.trackPlotRegion(track.index = 1,
                       bg.border = NA,
                       panel.fun = function(x, y) {
                         xlim = get.cell.meta.data("xlim")
                         sector.index = get.cell.meta.data("sector.index")
                         # Text direction
                         theta = circlize(mean(xlim), 1)[1, 1] %% 360
                         dd = ifelse(theta < 180 || theta > 360, "bending.inside", "bending.outside")
                         circos.text(x = mean(xlim), y = 0.8, labels = sector.index, facing = dd,
                                     niceFacing = TRUE, cex = 1, font = 2)
                       }
)

ggsave("figures/Figure_2c.pdf", BA3_broad, height = 16.5, width = 21, units = "cm")

# Fine pop assignment

baysass <- tibble(Region = c("BAS","W_SCOT","ICE","E_SCOT","DEN","IRE","GER","ENG","NEL", "WNA"),
                  BAS = c(0.6982,0.0227,0.1441,0.0158,0.0151,0.0157,0.0156,0.0249,0.0167,0.0154),
                  W_SCOT = c(0.0067,0.6928,0.0069,0.0066,0.0062,0.0068,0.0065,0.2481,0.0065,0.0064),
                  ICE = c(0.0115,0.0106,0.8810,0.0120,0.0111,0.0112,0.0121,0.0158,0.0110,0.0120),
                  E_SCOT = c(0.0064,0.0450,0.0270,0.6733,0.0064,0.0070,0.0124,0.2031,0.0064,0.0069),
                  DEN = c(0.0198,0.0202,0.0193,0.0203,0.6863,0.0206,0.0193,0.1345,0.0196,0.0206),
                  IRE = c(0.0251,0.0894,0.0237,0.0234,0.0231,0.6896,0.0237,0.0292,0.0247,0.0234),
                  GER = c(0.0134,0.0310,0.0145,0.0140,0.0143,0.0148,0.6821,0.1746,0.0136,0.0138),
                  ENG = c(0.0097,0.0330,0.0100,0.0099,0.0095,0.0097,0.0096,0.8799,0.0098,0.0096),
                  NEL = c(0.0196,0.0300,0.0191,0.0188,0.0188,0.0177,0.0188,0.1317,0.6867,0.0190),
                  WNA = c(0.0243,0.0230,0.0920,0.0245,0.0246,0.0237,0.0233,0.0236,0.0244,0.6919))

# Convert to matrix
baysass.mat = as.matrix(baysass[, 2:11])
baysass.mat

# Gene flow matrix
dimnames(baysass.mat) = list(source = baysass$Region, sink = baysass$Region)
baysass.mat

# Reset the circular layout parameters
circos.clear()

# Parameters for the circular layout
circos.par(start.degree = 30, gap.degree = 3)

# Plot chord diagram
chordDiagram(x = baysass.mat, grid.col = col_palette2, grid.border = "black", transparency = 0.3,
             order = baysass$Region, directional = T, direction.type = "arrows",
             self.link = 1, preAllocateTracks = list(track.height = 0.01),
             annotationTrack = "grid", annotationTrackHeight = c(0.05, 0.05),
             link.border = "NA", link.sort = T, link.decreasing = T,
             link.arr.length = 0.4, link.arr.lty = 3, link.arr.col = "#252525", 
             link.largest.ontop = F)

chordDiagram(x = baysass.mat, grid.col = col_palette2, transparency = 0.4,
             order = baysass$Region, directional = T, direction.type = "arrows",
             self.link = 1, preAllocateTracks = list(track.height = 0.05),
             annotationTrack = "grid", annotationTrackHeight = c(0.05, 0.05),
             link.border = "NA", link.sort = T, link.decreasing = T,
             link.arr.length = 0.2, link.arr.lty = 2, link.arr.col = "#252525", 
             link.largest.ontop = F)

BA3_fine <- circos.trackPlotRegion(track.index = 1,
                       bg.border = NA,
                       panel.fun = function(x, y) {
                         xlim = get.cell.meta.data("xlim")
                         sector.index = get.cell.meta.data("sector.index")
                         # Text direction
                         theta = circlize(mean(xlim), 1)[1, 1] %% 360
                         dd = ifelse(theta < 180 || theta > 360, "bending.inside", "bending.outside")
                         circos.text(x = mean(xlim), y = 0.8, labels = sector.index, facing = dd,
                                     niceFacing = TRUE, cex = 1, font = 2)
                       }
)

ggsave("figures/Figure_S8.pdf", BA3_fine, height = 16.5, width = 21, units = "cm")
