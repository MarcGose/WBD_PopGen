# load packes

library(tidyverse)
library(purrr)
library(forcats)
library(readxl)
library(data.table)
library(wesanderson)
library(patchwork)
library(nord)
library(tidytext)

# define colour palette

col_palette <- c( nord("aurora")[1],
                  nord("victory_bonds")[3],
                  nord("aurora")[4],
                  nord("moose_pond")[3],
                  nord("frost")[3],
                  nord("baie_mouton")[7])

# load meta data

meta <- read.table("file_lists/pop_file_all_together.info", header = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~#
#     Ln and Delta K      #
#~~~~~~~~~~~~~~~~~~~~~~~~~#

data_path <- "C:/Users/MarcG/OneDrive/Desktop/WBD_PopGen/files"
files <- dir(data_path, pattern = "*.log")

df <- tibble(filename = files) %>%
  mutate(file_contents = map(filename,
                             ~ readLines(file.path(data_path, .))))

get_best_like <- function(df) {
  like <- df[grep(pattern = "best like=",df)]
  like <- substr(like,regexpr(pattern="=",like)[1]+1,regexpr(pattern="=",like)[1]+14)
}

lnl <- df %>%
  mutate(lnl = map(df$file_contents, get_best_like)) %>%
  select(-file_contents) %>%
  unnest(cols = c(lnl)) %>%
  #mutate(filename = gsub("ORYX_GL_HS_NS_NGSadmix_", "", filename)) %>%
  mutate(filename = gsub("WBD_GLs_NGSadmix_", "", filename)) %>%
  mutate(filename = gsub("_out.log", "", filename)) %>%
  separate(filename, c("K", "run"), sep = "_") %>%
  mutate(run = gsub("run", "", run)) %>%
  mutate(run = as.numeric(run)) %>%
  mutate(K = as.factor(K)) %>%
  mutate(lnl = as.numeric(lnl)) %>%
  group_by(K) %>%
  summarise(mean = mean(lnl),
            minlnl = min(lnl),
            maxlnl = max(lnl),
            sd = sd(lnl)) %>%
  mutate(K = as.numeric(gsub("K", "", K)))

lnl_plot <- ggplot(lnl, aes(K, mean)) +
  geom_point() +
  geom_pointrange(aes(ymin = minlnl, ymax = maxlnl),
                  size = 0.5) +
  theme_minimal() +
  ylab("Likelihood") +
  scale_y_continuous(labels = scales::scientific) +
  ggtitle("B")

lnl_plot

#~~ Delta K

deltaK <- lnl %>%
  mutate(LprimeK = c(NA,mean[-1]-mean[-length(mean)]),
         LdblprimeK = c(NA,LprimeK[-c(1,length(mean))]-(LprimeK)[-(1:2)],NA),
         delta = LdblprimeK/sd) 

deltak_plot <- ggplot(deltaK, aes(K, delta)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  scale_x_continuous(limits=c(2, 5)) +
  ylab("Delta K") +
  ggtitle("A") +
  scale_y_continuous(labels = scales::scientific)

deltak_plot

opt_k <- deltak_plot + lnl_plot

# safe figure

ggsave("figures/Figure_S2.pdf", opt_k, height = 16.5, width = 21, units = "cm")

#~~~~~~~~~~~~~~~~~~~~~~~~~~#
#         Plotting         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~#

files <- dir(data_path, pattern = "*.qopt")

df <- tibble(filename = files) %>%
  mutate(file_contents = map(filename,
                             ~ fread(file.path(data_path, .))))

# Unnest data

df <- unnest(df, cols = c(file_contents)) %>%
  mutate(ID = rep(meta$V1, 80)) %>%
  pivot_longer(cols = c(V1, V2, V3, V4, V5, V6, V7, V8)) %>%
  drop_na()

df <- left_join(df, meta, by = c("ID" = "V1"))

df <- df %>% arrange(ID, name == "V1")

# Split by K

K1 <- filter(df, grepl("K1_run1_out.qopt$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K2 <- filter(df, grepl("K2_run1_out.qopt$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K3 <- filter(df, grepl("K3_run1_out.qopt$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K4 <- filter(df, grepl("K4_run1_out.qopt$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K5 <- filter(df, grepl("K5_run1_out.qopt$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K6 <- filter(df, grepl("K6_run1_out.qopt$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

# Plot

k1_plot <- ggplot(K1, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", linewidth = 0.1) +
  facet_grid(~fct_inorder(V2), switch = "x", scales = "free", space = "free") +
  theme_minimal() + 
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1.25)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  ylab("K=1")

k1_plot

k2_plot <- ggplot(K2_test, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", linewidth = 0.1) +
  facet_grid(~fct_inorder(V2), switch = "x", scales = "free_x", space = "free") +
  theme_minimal() + 
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_reordered(expand = expansion(add = 1.25)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  ylab("K=2")

k2_plot

k3_plot <- ggplot(K3, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", linewidth = 0.1) +
  facet_grid(~fct_inorder(V2), switch = "x", scales = "free", space = "free") +
  theme_minimal() + 
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1.25)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(angle = 75, hjust = 1, vjust = -60),
        legend.position = "none") +
  ylab("K=3")

k3_plot 

k4_plot <- ggplot(K4, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", linewidth = 0.1) +
  facet_grid(~fct_inorder(V2), switch = "x", scales = "free", space = "free") +
  theme_minimal() + 
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1.25)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  ylab("K=4")

k4_plot

k5_plot <- ggplot(K5, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", linewidth = 0.1) +
  facet_grid(~fct_inorder(V2), switch = "x", scales = "free", space = "free") +
  theme_minimal() + 
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1.25)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  ylab("K=5")

k5_plot

k6_plot <- ggplot(K6, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", linewidth = 0.1) +
  facet_grid(~fct_inorder(V2), switch = "x", scales = "free", space = "free") +
  theme_minimal() +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1.25)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        #axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(angle = 75, hjust = 1, vjust = -60),
        legend.position = "none") +
  ylab("K=6")

k6_plot

ngsAdmix_combined <- k1_plot + k2_plot + k3_plot + k4_plot + k5_plot + k6_plot + plot_layout(ncol = 1)

ggsave("figures/Figure_S3.pdf", ngsAdmix_combined, height = 16.5, width = 21, units = "cm")

### Sorted plot ####

K3_meta <- read.table("files/K=3_sorted.txt")

df_long <- pivot_longer(K3_meta, cols = c(V1, V2, V3), names_to = "Cluster", values_to = "Admixture_Proportion")

plot3 <- ggplot(df_long, aes(factor(V5), Admixture_Proportion, fill = factor(Cluster))) +
  geom_col(color = "gray", linewidth = 0.1) +
  theme_minimal() +
  facet_grid(~fct_inorder(V4), switch = "x", scales = "free", space = "free") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1.25)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        #axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(angle = 75, hjust = 1, vjust = -60),
        legend.position = "none") +
  ylab("K=3")

ggsave("figures/Figure_1c.pdf", plot3, height = 16.5, width = 21, units = "cm")
