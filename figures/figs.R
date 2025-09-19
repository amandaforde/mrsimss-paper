## FIGURES:

## This script provides the code used to produce the figures in the manuscript
## concerning both the simulation study and empirical analyses.

## Outputs:
## Figs 1-2 (simulation study)
## Supplementary Figs 1-10 (simulation study)
## Figs 3-4 (empirical analyses)


## Load required packages:
library(ggplot2)
library(dplyr)
library(patchwork)


################################################################################

## Fig 1. Estimated causal effect for each method and simulation setting with
## sample sizes of 200,000 and zero overlap, averaged over 100 simulated pairs
## of exposure and outcome GWAS summary statistics for each setting.

res.all <- read.table("simulations/results/sims-res-200.csv",header=TRUE)
res.all <- res.all %>%
  mutate(method = recode(method,
                         "SimSS2-mr_ivw"  = "SimSS-2-IVW",
                         "SimSS2-mr_raps" = "SimSS-2-RAPS",
                         "SimSS3-mr_ivw"  = "SimSS-3-IVW",
                         "SimSS3-mr_raps" = "SimSS-3-RAPS",
                         "mr_ivw" = "IVW",
                         "mr_raps" = "RAPS",
                         "mr_divw" = "dIVW"
  ))
res.all$method <- factor(res.all$method, levels=c("SimSS-2-IVW","SimSS-3-IVW","SimSS-2-RAPS","SimSS-3-RAPS","IVW","RAPS","dIVW"))

my_palette <- c(
  "SimSS-2-IVW" = "darkred",
  "SimSS-3-IVW" = "tomato3",
  "SimSS-2-RAPS" = "#053061",
  "SimSS-3-RAPS" = "#4393c3",
  "IVW" = "#D95F02FF",
  "RAPS" = "#7570B3FF",
  "dIVW" = "#E6AB02FF"
)

plot_sim <- function(data, cor_label){
  ggplot(data, aes(x = method, y = b, color = method, fill = method)) +
    geom_boxplot(size=0.6,aes(fill=method, color=method, alpha=0.3),outlier.size = 1) +
    stat_summary(fun = mean, geom = "point", shape = 21, size = 1, fill = "white") +
    geom_hline(yintercept = 0.3, linetype = 1, color = "black", linewidth = 0.5) +
    facet_grid(h2~prop_effect, labeller="label_both") +
    coord_cartesian(ylim = c(0.2, 0.4)) +
    scale_color_manual(values = my_palette) +
    scale_fill_manual(values = my_palette) +
    labs(
      x = "MR Method",
      y = "Estimated Causal Effect",
      title = bquote(cor[xy]~"="~.(cor_label))
    ) +
    theme_minimal(base_size = 10) + theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 9),
      plot.title = element_text(size = 10, face = "bold"),
      strip.text = element_text(size = 8, face = "italic"),
      panel.grid = element_line(size = 0.2, color = "grey90")
    )
}


## Zero Overlap
res.0 <- res.all[res.all$frac_overlap == 1,]
plotA <- plot_sim(res.0[res.0$cor_xy == -0.1,], cor_label = -0.1)
plotB <- plot_sim(res.0[res.0$cor_xy == 0.1,],  cor_label = 0.1)
plotC <- plot_sim(res.0[res.0$cor_xy == 0.3,],  cor_label = 0.3)
plotD <- plot_sim(res.0[res.0$cor_xy == 0.5,],  cor_label = 0.5)

# Combine plots:
figure <- plotA + plotB + plotC + plotD +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A",title = "Zero Overlap"
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 10),
    plot.title = element_text(size = 12, face = "bold", hjust = 0)
  )

# Save figure
ggsave("Fig-1-SIMS-200-0.png", figure, width = 8, height = 10, dpi = 600)

################################################################################

## Fig 2. Estimated causal effect for each method and simulation setting with
## sample sizes of 200,000 and full overlap, averaged over 100 simulated pairs
## of exposure and outcome GWAS summary statistics for each setting.

## Full Overlap
res.1 <- res.all[res.all$frac_overlap == 1,]
plotA <- plot_sim(res.1[res.1$cor_xy == -0.1,], cor_label = -0.1)
plotB <- plot_sim(res.1[res.1$cor_xy == 0.1,],  cor_label = 0.1)
plotC <- plot_sim(res.1[res.1$cor_xy == 0.3,],  cor_label = 0.3)
plotD <- plot_sim(res.1[res.1$cor_xy == 0.5,],  cor_label = 0.5)

# Combine plots:
figure <- plotA + plotB + plotC + plotD +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A",title = "Full Overlap"
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 10),
    plot.title = element_text(size = 12, face = "bold", hjust = 0)
  )

# Save figure
ggsave("Fig-2-SIMS-200-1.png", figure, width = 8, height = 10, dpi = 600)

################################################################################

## Supplementary Fig 1. Estimated causal effect for each method and simulation
## setting with sample size = 200,000 and 25% overlap, averaged over 100
## simulated pairs of exposure and outcome GWAS summary statistics for each
## setting.

## 25% Overlap
res.25 <- res.all[res.all$frac_overlap == 0.25,]
plotA <- plot_sim(res.25[res.25$cor_xy == -0.1,], cor_label = -0.1)
plotB <- plot_sim(res.25[res.25$cor_xy == 0.1,],  cor_label = 0.1)
plotC <- plot_sim(res.25[res.25$cor_xy == 0.3,],  cor_label = 0.3)
plotD <- plot_sim(res.25[res.25$cor_xy == 0.5,],  cor_label = 0.5)

# Combine plots:
figure <- plotA + plotB + plotC + plotD +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A",title = "25% Overlap"
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 10),
    plot.title = element_text(size = 12, face = "bold", hjust = 0)
  )

# Save figure
ggsave("Supp-Fig-1-SIMS-200-25.png", figure, width = 8, height = 10, dpi = 600)

################################################################################

## Supplementary Fig 2. Estimated causal effect for each method and simulation
## setting with sample size = 200,000 and 50% overlap, averaged over 100
## simulated pairs of exposure and outcome GWAS summary statistics for each
## setting.

## 50% Overlap
res.5 <- res.all[res.all$frac_overlap == 0.5,]
plotA <- plot_sim(res.5[res.5$cor_xy == -0.1,], cor_label = -0.1)
plotB <- plot_sim(res.5[res.5$cor_xy == 0.1,],  cor_label = 0.1)
plotC <- plot_sim(res.5[res.5$cor_xy == 0.3,],  cor_label = 0.3)
plotD <- plot_sim(res.5[res.5$cor_xy == 0.5,],  cor_label = 0.5)

# Combine plots:
figure <- plotA + plotB + plotC + plotD +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A",title = "50% Overlap"
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 10),
    plot.title = element_text(size = 12, face = "bold", hjust = 0)
  )

# Save figure
ggsave("Supp-Fig-2-SIMS-200-50.png", figure, width = 8, height = 10, dpi = 600)

################################################################################

## Supplementary Fig 3. Estimated causal effect for each method and simulation
## setting with sample size = 200,000 and 75% overlap, averaged over 100
## simulated pairs of exposure and outcome GWAS summary statistics for each
## setting.

## Full Overlap
res.75 <- res.all[res.all$frac_overlap == 0.75,]
plotA <- plot_sim(res.75[res.75$cor_xy == -0.1,], cor_label = -0.1)
plotB <- plot_sim(res.75[res.75$cor_xy == 0.1,],  cor_label = 0.1)
plotC <- plot_sim(res.75[res.75$cor_xy == 0.3,],  cor_label = 0.3)
plotD <- plot_sim(res.75[res.75$cor_xy == 0.5,],  cor_label = 0.5)

# Combine plots:
figure <- plotA + plotB + plotC + plotD +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A",title = "75% Overlap"
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 10),
    plot.title = element_text(size = 12, face = "bold", hjust = 0)
  )

# Save figure
ggsave("Supp-Fig-3-SIMS-200-75.png", figure, width = 8, height = 10, dpi = 600)

################################################################################

## Supplementary Fig 4. Estimated causal effect for each method and simulation
## setting with sample size = 200,000, collapsed over heritability and
## proportion of true effect variants, averaged over 100 simulated pairs of
## exposure and outcome GWAS summary statistics for each setting.

plot.full <- ggplot(res.all, aes(x = method, y = b, color = method, fill = method)) +
  geom_boxplot(size=0.6,aes(fill=method, color=method, alpha=0.3),outlier.size = 1) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 1, fill = "white") +
  geom_hline(yintercept = 0.3, linetype = 1, color = "black", size = 0.5) +
  facet_grid(cor_xy~frac_overlap, labeller="label_both") +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  labs(
    x = "MR Method",
    y = "Estimated Causal Effect",
    title = "Sample size = 200,000"
  ) +
  theme_minimal(base_size = 10) + theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 9),
    plot.title = element_text(size = 10, face = "bold"),
    strip.text = element_text(size = 8, face = "italic"),
    panel.grid = element_line(size = 0.2, color = "grey90")
  )

# Save figure
ggsave("Supp-Fig-4-SIMS-200-ALL.png", plot.full, width = 8, height = 10, dpi = 600)

################################################################################

## Supplementary Fig 5. Estimated causal effect for each method and simulation
## setting with sample size = 500,000, collapsed over heritability and
## proportion of true effect variants, averaged over 100 simulated pairs of
## exposure and outcome GWAS summary statistics for each setting.

res.all <- read.table("simulations/results/sims-res-500.csv",header=TRUE)
res.all <- res.all %>%
  mutate(method = recode(method,
                         "SimSS2-mr_ivw"  = "SimSS-2-IVW",
                         "SimSS2-mr_raps" = "SimSS-2-RAPS",
                         "SimSS3-mr_ivw"  = "SimSS-3-IVW",
                         "SimSS3-mr_raps" = "SimSS-3-RAPS",
                         "mr_ivw" = "IVW",
                         "mr_raps" = "RAPS",
                         "mr_divw" = "dIVW"
  ))
res.all$method <- factor(res.all$method, levels=c("SimSS-2-IVW","SimSS-3-IVW","SimSS-2-RAPS","SimSS-3-RAPS","IVW","RAPS","dIVW"))

plot.full <- ggplot(res.all, aes(x = method, y = b, color = method, fill = method)) +
  geom_boxplot(size=0.6,aes(fill=method, color=method, alpha=0.3),outlier.size = 1) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 1, fill = "white") +
  geom_hline(yintercept = 0.3, linetype = 1, color = "black", size = 0.5) +
  facet_grid(cor_xy~frac_overlap, labeller="label_both") +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  labs(
    x = "MR Method",
    y = "Estimated Causal Effect",
    title = "Sample size = 500,000"
  ) +
  theme_minimal(base_size = 10) + theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 9),
    plot.title = element_text(size = 10, face = "bold"),
    strip.text = element_text(size = 8, face = "italic"),
    panel.grid = element_line(size = 0.2, color = "grey90")
  )

# Save figure
ggsave("Supp-Fig-5-SIMS-500-ALL.png", plot.full, width = 8, height = 10, dpi = 600)


################################################################################

## Supplementary Fig 6. Estimated causal effect for each method and simulation
## setting with sample size = 50,000, heritability = 0.3 and proportion of true
## effect variants = 0.001, across 100 simulated pairs of exposure and outcome
## GWAS summary statistics for each setting.

res.all <- read.table("simulations/results/sims-res-50.csv",header=TRUE)
res.all <- res.all %>%
  mutate(method = recode(method,
                         "SimSS2-mr_ivw"  = "SimSS-2-IVW",
                         "SimSS2-mr_raps" = "SimSS-2-RAPS",
                         "SimSS3-mr_ivw"  = "SimSS-3-IVW",
                         "SimSS3-mr_raps" = "SimSS-3-RAPS",
                         "mr_ivw" = "IVW",
                         "mr_raps" = "RAPS",
                         "mr_divw" = "dIVW"
  ))
res.all$method <- factor(res.all$method, levels=c("SimSS-2-IVW","SimSS-3-IVW","SimSS-2-RAPS","SimSS-3-RAPS","IVW","RAPS","dIVW"))

res.3.001 <- res.all[res.all$h2 == 0.3 & res.all$prop_effect ==  0.001,]
plot.3.001 <- ggplot(res.3.001, aes(x = method, y = b, color = method, fill = method)) +
  geom_boxplot(size=0.6,aes(fill=method, color=method, alpha=0.3),outlier.size = 1) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 1, fill = "white") +
  geom_hline(yintercept = 0.3, linetype = 1, color = "black", size = 0.5) +
  facet_grid(cor_xy~frac_overlap, labeller="label_both") +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  labs(
    x = "MR Method",
    y = "Estimated Causal Effect",
    title = "Sample size = 50,000",
    subtitle = "Heritability = 0.3, Proportion of effect SNPs = 0.001"
  ) +
  theme_minimal(base_size = 10) + theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 9),
    plot.title = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 9),
    strip.text = element_text(size = 8, face = "italic"),
    panel.grid = element_line(size = 0.2, color = "grey90")
  )

# Save figure
ggsave("Supp-Fig-6-SIMS-50-3-001.png", plot.3.001, width = 8, height = 10, dpi = 600)

################################################################################

## Supplementary Fig 7. Estimated causal effect for each method and simulation
## setting with sample size = 50,000, heritability = 0.3 and proportion of true
## effect variants = 0.01, across 100 simulated pairs of exposure and outcome
## GWAS summary statistics for each setting.

res.3.01 <- res.all[res.all$h2 == 0.3 & res.all$prop_effect ==  0.01,]
plot.3.01 <- ggplot(res.3.01, aes(x = method, y = b, color = method, fill = method)) +
  geom_boxplot(size=0.6,aes(fill=method, color=method, alpha=0.3),outlier.size = 1) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 1, fill = "white") +
  geom_hline(yintercept = 0.3, linetype = 1, color = "black", size = 0.5) +
  facet_grid(cor_xy~frac_overlap, labeller="label_both") +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  labs(
    x = "MR Method",
    y = "Estimated Causal Effect",
    title = "Sample size = 50,000",
    subtitle = "Heritability = 0.3, Proportion of effect SNPs = 0.01"
  ) +
  theme_minimal(base_size = 10) + theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 9),
    plot.title = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 9),
    strip.text = element_text(size = 8, face = "italic"),
    panel.grid = element_line(size = 0.2, color = "grey90")
  )

# Save figure
ggsave("Supp-Fig-7-SIMS-50-3-01.png", plot.3.01, width = 8, height = 10, dpi = 600)

################################################################################

## Supplementary Fig 8. Estimated causal effect for each method and simulation
## setting with sample size = 50,000, heritability = 0.7 and proportion of true
## effect variants = 0.001, across 100 simulated pairs of exposure and outcome
## GWAS summary statistics for each setting.

res.7.001 <- res.all[res.all$h2 == 0.7 & res.all$prop_effect ==  0.001,]
plot.7.001 <- ggplot(res.7.001, aes(x = method, y = b, color = method, fill = method)) +
  geom_boxplot(size=0.6,aes(fill=method, color=method, alpha=0.3),outlier.size = 1) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 1, fill = "white") +
  geom_hline(yintercept = 0.3, linetype = 1, color = "black", size = 0.5) +
  facet_grid(cor_xy~frac_overlap, labeller="label_both") +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  labs(
    x = "MR Method",
    y = "Estimated Causal Effect",
    title = "Sample size = 50,000",
    subtitle = "Heritability = 0.7, Proportion of effect SNPs = 0.001"
  ) +
  theme_minimal(base_size = 10) + theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 9),
    plot.title = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 9),
    strip.text = element_text(size = 8, face = "italic"),
    panel.grid = element_line(size = 0.2, color = "grey90")
  )

# Save figure
ggsave("Supp-Fig-8-SIMS-50-7-001.png", plot.7.001, width = 8, height = 10, dpi = 600)

################################################################################

## Supplementary Fig 9. Estimated causal effect for each method and simulation
## setting with sample size = 50,000, heritability = 0.7 and proportion of true
## effect variants = 0.01, across 100 simulated pairs of exposure and outcome
## GWAS summary statistics for each setting.

res.7.01 <- res.all[res.all$h2 == 0.7 & res.all$prop_effect ==  0.01,]
plot.7.01 <- ggplot(res.7.01, aes(x = method, y = b, color = method, fill = method)) +
  geom_boxplot(size=0.6,aes(fill=method, color=method, alpha=0.3),outlier.size = 1) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 1, fill = "white") +
  geom_hline(yintercept = 0.3, linetype = 1, color = "black", size = 0.5) +
  facet_grid(cor_xy~frac_overlap, labeller="label_both") +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  labs(
    x = "MR Method",
    y = "Estimated Causal Effect",
    title = "Sample size = 50,000",
    subtitle = "Heritability = 0.7, Proportion of effect SNPs = 0.01"
  ) +
  theme_minimal(base_size = 10) + theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 9),
    plot.title = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 9),
    strip.text = element_text(size = 8, face = "italic"),
    panel.grid = element_line(size = 0.2, color = "grey90")
  )

# Save figure
ggsave("Supp-Fig-9-SIMS-50-7-01.png", plot.7.01, width = 8, height = 10, dpi = 600)

################################################################################

## Supplementary Fig 10. Estimated causal effect for SimSS-3-RAPS with four
## different thresholds implemented for varying overlap fractions with sample
## size = 50,000, heritability = 0.3, proportion of true effect variants = 0.01
## and exposure-outcome correlation = 0.5 across 100 simulated pairs of exposure
## and outcome GWAS summary statistics for each setting.


res.all <- read.table("simulations/results/sims-res-thresh.csv",header=TRUE)
res.all$threshold <- rep(c("5x10-8","5x10-6","5x10-5","5x10-4"), nrow(res.all)/4)

res.all$threshold <- factor(res.all$threshold, levels=c("5x10-8","5x10-6","5x10-5","5x10-4"))

my_cols <- c("#8FD19E","#49C16D","#2FA35C","#0F6D3B")

plot_full <- ggplot(res.all, aes(x = threshold, y = b, color = threshold, fill = threshold)) +
  geom_boxplot(size=0.6,aes(fill=threshold, color=threshold, alpha=0.3),outlier.size = 1) +
  # add mean point
  stat_summary(fun = mean, geom = "point", shape = 21, size = 1, fill = "white") +
  geom_hline(yintercept = 0.3, linetype = 1, color = "black", linewidth = 0.5) +
  facet_grid(cor_xy~frac_overlap, labeller="label_both") +
  scale_color_manual(values = my_cols) +
  scale_fill_manual(values = my_cols) +
  ylim(-0.2, 0.8) +
  labs(
    x = "Threshold",
    y = "Estimated Causal Effect",
    title = "Sample size = 50,000",
    subtitle = "Heritability = 0.3, Proportion of effect SNPs = 0.01"
  ) +
  theme_minimal(base_size = 10) + theme_bw() +
  scale_x_discrete(
    name = "Threshold",
    labels = c(
      expression("5×10"^-8),
      expression("5×10"^-6),expression("5×10"^-5),expression("5×10"^-4))) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 9),
    plot.title = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 9),
    strip.text = element_text(size = 8, face = "italic"),
    panel.grid = element_line(size = 0.2, color = "grey90")
  )

# Save figure
ggsave("Supp-Fig-10.png", plot_full, width = 9.5, height = 4, dpi = 600)


################################################################################
################################################################################

## Fig 3. Boxplots of the estimated causal effects for each method resulting
## from the 20 same-trait BMI-BMI analyses.

res.all <- read.table("real_data/results/res-BMI-BMI.csv",header=TRUE)

my_palette <- c(
  "SimSS-2-IVW" = "darkred",
  "SimSS-3-IVW" = "tomato3",
  "SimSS-2-RAPS" = "#053061",
  "SimSS-3-RAPS" = "#4393c3",
  "IVW" = "#D95F02FF",
  "RAPS" = "#7570B3FF",
  "Egger" = "#CC79A7",
  "Weighted median" = "#762a83",
  "dIVW" = "#E6AB02FF"
)

res.all$method <- factor(res.all$method, levels = c(
  "SimSS-2-IVW", "SimSS-3-IVW", "SimSS-2-RAPS", "SimSS-3-RAPS", "IVW", "Egger", "Weighted median", "RAPS", "dIVW"
))

plot.BMI.BMI <- ggplot(res.all, aes(y = method, x = b, color = method, fill = method)) +
  geom_boxplot(size=0.6,aes(fill=method, color=method, alpha=0.3),outlier.size = 1) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 1, fill = "white") +
  geom_vline(xintercept = 1, linetype = 1, color = "black", linewidth = 0.5) +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  labs(
    y = "MR Method",
    x = "Estimated Causal Effect",
    title = "BMI→BMI: Estimates Across 20 MR Analyses"
  ) +
  theme_minimal(base_size = 10) + theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 9),
    plot.title = element_text(size = 10, face = "bold"),
    strip.text = element_text(size = 8, face = "italic"),
    panel.grid = element_line(size = 0.2, color = "grey90")
  )

# Save figure
ggsave("Fig-3-BMI-BMI.png", plot.BMI.BMI, width = 8, height = 3, dpi = 600)


################################################################################

## Fig 4. Causal effect estimates for BMI-T2D analyses, with varying sample
## overlap.

res.all <- read.table("real_data/results/res-BMI-T2D.csv",header=TRUE)

res.all$overlap <-  factor(results_all$overlap, levels=c("0%","25%","50%","75%","100%"))
res.all$method <- factor(res.all$method, levels = c(
  "SimSS-2-IVW", "SimSS-3-IVW", "SimSS-2-RAPS", "SimSS-3-RAPS", "IVW", "Egger", "Weighted median", "RAPS", "dIVW"
))

# Add vertical gridlines
minor_lines <- c(1.25,1.75,2.25)

plot.T2D.BMI <- ggplot(results_all, aes(x = b, y = method, xmin = CI.lower, xmax = CI.upper, color = method)) +
  geom_vline(xintercept = minor_lines, linetype = "dotted", color = "grey80", size = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_point(size = 2.2) +
  geom_errorbarh(height = 0.2, size = 0.5) +
  facet_wrap(~overlap, nrow = 1, strip.position = "top") +
  scale_color_manual(values = my_palette) +
  scale_x_continuous(breaks = c(1, 1.5, 2)) +
  labs(
    x = "Estimated Causal Effect",
    y = "MR Method",
    title = "BMI→T2D: Estimates Across Overlap Scenarios"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )

# Save figure
ggsave("Fig-4-T2D-BMI.png", plot.T2D.BMI, width = 10, height = 3, dpi = 600)

