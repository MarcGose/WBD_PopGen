# Calculate Fst and plot Fst against distance for IBD

library(vcfR)
library(dartR)
library(RColorBrewer)
library(gplots)
library(scales)
library(tidyverse)
library(dplyr)
library(adespatial)
library(vegan)

# define colours for heatmap

heat_col <- colorRampPalette(brewer.pal(9, "Blues"))(100)

# define colour palette for IBD plot

col_palette <- c(nord("frost")[2], #lightblue
          nord("aurora")[4], #lightgreen
          nord("frost")[4], #darkblue
          nord("moose_pond")[4], #lightbrown
          nord("silver_mine")[2], #darkgreen
          nord("aurora")[1], #darkred
          nord("moose_pond")[3], #darkbrown
          nord("aurora")[2]) #lightred

# read in VCF
        
vcf <- read.vcfR("files/WBD_GLs.vcf")

# read in metadata with IDs and pop assignments

gen_pops <- read.table("file_lists/pop_file_broad_fine.info")

# calculate Fsts

genotypes <- vcfR2genlight(vcf)

genotypes@pop <- as.factor(gen_pops$V2)

fst_broad <- gl.fst.pop(genotypes, nboots = 10000, percent = 95, nclusters = 1) # examine fst_broad for pairwise FSTs and p-values of genetic pops

# calculate per sampling site Fst

genotypes@pop <- as.factor(gen_pops$V3)

fst_fine <- gl.fst.pop(genotypes, nboots = 10000, percent = 95, nclusters = 1)

Fsts <- fst_fine$Fsts

Fsts <- Fsts[-7,-7]

fst_df <- as.data.frame(Fsts)

fst_df <- fst_df[-7,-7]

# Reshape the data for plotting

fst_melted <- reshape2::melt(Fsts)
fst_melted$Var1 <- factor(fst_melted$Var1, levels = c("WNA", "ICE", "BAS", "W_SCOT", "IRE", "E_SCOT", "ENG", "GER", "DEN", "NEL"))
fst_melted$Var2 <- factor(fst_melted$Var2, levels = c("WNA", "ICE", "BAS", "W_SCOT", "IRE", "E_SCOT", "ENG", "GER", "DEN", "NEL"))

# Create the non-mirrored heatmap using ggplot2

(Fst_heat <- ggplot(fst_melted, aes(Var2, Var1)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradientn(colors = heat_col, 
                       values = rescale(c(0, 0.5, 1)),
                       na.value = "white",
                       guide = "legend", 
                       name = "Fst Values") +
  labs(x = "", y = "", title = "") +
  theme_minimal() +
  theme(plot.background = element_rect(fill = "white")))

# safe figure

ggsave("figures/Figure_2a.pdf", Fst_heat, height = 16.5, width = 21, units = "cm")

#IBD

source("https://raw.githubusercontent.com/jorgeassis/marineDistances/master/Script.R")

global.polygon <- "files/GSHHS_h_L1.shp"

contour(global.polygon = global.polygon, file = "files/sampling_sites.txt", file.sep = ";", file.dec = ".", file.strucutre = 2, file.header = FALSE, resolution = 0.1, buffer = c(1,1,1,1), export.file = TRUE)

distances <- read.delim("file_lists/Pairwise_Marine_Distances.txt", sep = ";")

Dist_Fst <- read.table("file_lists/DistvsFst.txt", header = T)

Dist <- as.matrix(Dist_Fst$Dist)

Fst <- as.matrix(Dist_Fst$Fst)

Fst_stand <- Fst / (1 - Fst)

df_dist <- data.frame(geodistance=as.vector(Dist),
                      gendistance=as.vector(Fst_stand),
                      Dist_Fst$Group2)
# RDA and ANOVA

rda <- rda(Fst_stand, Dist)
RsquareAdj(rda)
anova(rda, perm=1000)

# Plot

(ibd_plot <- ggplot(df_dist) + 
    aes(x = geodistance, y = gendistance, col = Dist_Fst.Group2) +
    geom_point(shape = 19, alpha = 0.6, stroke = 0.1, size = 8) +
    scale_color_manual(values = cols) +
    xlab("Geographic distance (km)") +
    ylab("Fst / (1 - Fst)") +
    theme_classic())

(ibd_plot_curve <- ggplot(df_dist) + 
    aes(x = geodistance, y = gendistance) +
    geom_point(shape = 19, alpha = 0.6, stroke = 0.1, size = 8) +
    geom_smooth(aes(x = geodistance, y = gendistance), method = "lm",
                se = TRUE, level= 0.95, alpha = 0.05, size = 1, col = "darkblue", fill = "darkblue") +
    xlab("Geographic distance (km)") +
    ylab("Fst / (1 - Fst)") +
    theme_classic())

# save figures

ggsave("figures/Figure_2b_1.pdf", ibd_plot, height = 16.5, width = 21, units = "cm")
ggsave("figures/Figure_2B_2.pdf", ibd_plot_curve, height = 16.5, width = 21, units = "cm")
