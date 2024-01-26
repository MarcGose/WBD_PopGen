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

col_palette2 <- c(nord("aurora")[2],
                nord("baie_mouton")[6],
                nord("aurora")[3])

# examine chain convergence from trace file

broad_trace <- fread("files/WBD_GLs.trace.txt", header = T)

(trace_fig <- ggplot(broad_trace) +
    aes(State, LogProb) +
    geom_line(linewidth = 1) +
    theme_minimal() + 
    scale_x_continuous() +
    scale_y_continuous())

ggsave("figures/FigS5.pdf", trace_fig, height = 16.5, width = 21, units = "cm")


# Prepare BayesAss data

# Broad pop assignment

baysass <- tibble(Region = c("BAS","ICE","WSI","NS","WNA"),
                  BAS = c(0.68892,0.22218,0.02228,0.04438,0.0222),
                  ICE = c(0.01512,0.93362,0.0152,0.02086,0.0152),
                  WSI = c(0.00682,0.00684,0.67348,0.30608,0.00678),
                  NS = c(0.00832,0.0156,0.00444,0.9679,0.00376),
                  WNA = c(0.04188,0.16656,0.04168,0.04152,0.70838))

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

ggsave("figures/Fig3c.pdf", BA3_broad, height = 16.5, width = 21, units = "cm")

# WSINS and E_SCOT

baysass <- tibble(Region = c("WSI","E_SCOT","NS"),
                  WSI = c(0.98582,0.0071,0.0071),
                  E_SCOT = c(0.31512,0.67574,0.00914),
                  NS = c(0.31972,0.00684,0.67348))
# Convert to matrix
baysass.mat = as.matrix(baysass[, 2:4])
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

BA3_WSINS <- circos.trackPlotRegion(track.index = 1,
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

ggsave("figures/Fig2c.pdf", BA3_WSINS, height = 16.5, width = 21, units = "cm")
