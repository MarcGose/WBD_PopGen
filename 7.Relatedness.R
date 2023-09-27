# load packages

library(data.table)
library(dplyr)
library(ggplot2)
library(paletteer)

# define colour palette

col_palette <- c(paletteer_d("ggthemes::Red_Blue_Brown")[12], #red
                 paletteer_d("ggthemes::Summer")[7], #orange
                 paletteer_d("ggthemes::Green_Orange_Teal")[2], #green
                 paletteer_d("ggthemes::Red_Blue_Brown")[2], #blue
                 paletteer_d("ggthemes::Winter")[9]) #brown

# read in genome file from PLINK

gen <- fread("files/WBD_GLs.genome", header = T)

summary(gen$PI_HAT)

ngsrel <- fread("files/WBD_NGSrelate.res")

gen$R1 <- ngsrel$R1 
gen$R0 <- ngsrel$R0
gen$KING <- ngsrel$KING

#remove pairs with values close to thesholds

gen <- gen %>%
  mutate(kinship = (Z1 / 4) + (Z2 / 2)) %>%
  mutate(criteria_2 = case_when(kinship >= 1/2^(5/2) & kinship < 1/2^(3/2) & Z0 < 0.1 ~ "Parent-offspring",
                                kinship >= 1/2^(5/2) & kinship < 1/2^(3/2) & Z0 > 0.1 & Z0 < 0.365 ~ "Full-sibling",
                                kinship >= 1/2^(7/2) + 0.01 & kinship < 1/2^(5/2) + 0.01 & Z0 > 0.365 + 0.01 & Z0 < 1-(1/(2^(3/2))) + 0.01 ~ "Second-degree",
                                kinship >= 1/2^(9/2) + 0.01 & kinship < 1/2^(7/2) + 0.01 & Z0 > 1-(1/2^(3/2)) + 0.01 & Z0 < 1 -(1/2^(5/2)) + 0.01 ~ "Third-degree",
                                kinship < 1/2^(9/2) + 0.01 & Z0 > 1-(1/2^(5/2)) + 0.01 ~ "Unrelated",
                                TRUE ~ "Unknown"))
# plot

(related <- ggplot(filter(gen, criteria_2 != "Unknown"), aes(R1, KING)) +
  geom_point(size = 4, alpha = 0.5,
             aes(colour = factor(criteria_2, 
                                 levels = c("Parent-offspring",
                                            "Full-sibling",
                                            "Second-degree",
                                            "Third-degree",
                                            "Unrelated",
                                            "Unknown")))) +
  scale_colour_manual(values = col_palette) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.2),
        legend.text = element_text(size=13),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(alpha=1),
                               keywidth = 0.2,
                               keyheight = 0.2,
                               default.unit = "inch")) +
  ggtitle("B"))

(PI_HAT_dist <- ggplot(gen, aes(x=PI_HAT)) + 
  geom_histogram(col = col_palette[4], alpha = 0.9, fill = col_palette[4], binwidth = 0.01) + 
  scale_y_continuous(trans='identity') +
  scale_x_continuous(limits = c(-0.1, 0.5)) +
  xlab("Pairwise relatedness") + ylab("Count") +
  theme_classic())

ggsave("figures/Figure_S1.pdf", related, height = 16.5, width = 21, units = "cm")
ggsave("figures/Figure_S1b.pdf", PI_HAT_dist, height = 16.5, width = 21, units = "cm")
