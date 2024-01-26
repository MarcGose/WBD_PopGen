# MLH and inbreeding

library(ggplot2)
library(inbreedR)
library(reshape)
library(dplyr)
library(data.table)
library(vcfR)
library(nord)

# define colour palette

col_palette <- c(nord("lumina")[2], #pink6
         nord("frost")[2], #lightblue1
         nord("aurora")[2], #orange
         nord("aurora")[3]) #yellow2 

# read in VCF

vcf <- read.vcfR("files/WBD_GLs.vcf")

# read in metadata with IDs and pop assignments

gen_pops <- read.table("file_lists/pop_file_broad_fine.info")

gt <- extract.gt(vcf)
gt <- as.data.frame(t(gt), stringsAsFactors = FALSE)
gt[gt == "."] <- NA
snp_geno <- do.call(cbind, apply(gt, 2, function(x) colsplit(x, "/", c("a","b"))))
WBD_SNPs <- convert_raw(snp_geno)
check_data(WBD_SNPs)
het <- MLH(WBD_SNPs)
barplot(het, ylab = "MLH", xlab = "ID")

df_het <- cbind(het, gen_pops)

mutate(df_het, across(where(is.character), as.factor))

# Plot MLH proportion per region
(MLH_per_region <- ggplot(df_het) +
  aes(x = het, y = V2) +
  labs(y = "Pop", x = "MLH") +
  geom_boxplot(fill = col_palette[1:4]) +
  scale_fill_manual(values = c(col_palette[1], col_palette[2], col_palette[3], col_palette[4])) +
  geom_point() +
  coord_flip() +
  theme_few())

(MLH_dist <- ggplot(df_het, aes(het)) +
  geom_histogram(aes(y=..density..),  position="identity", alpha=0.9, col = col_palette[1],
                 fill = col_palette[1]) +
  geom_density(alpha=0.6, col = col_palette[1], fill = col_palette[1]) +
  labs(y = "Density", x = "MLH") +
  theme_few() +
  xlim(c(0.05, 0.2)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  scale_fill_manual(values = col_palette[1]) +
  ggtitle("A"))

ggsave("figures/Fig3d.pdf", MLH_per_region, height = 16.5, width = 21, units = "cm")
ggsave("figures/FigS9.pdf", MLH_dist, height = 16.5, width = 21, units = "cm")
