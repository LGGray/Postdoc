library(ggplot2)

gnomAD <- read.delim('~/Downloads/gnomAD_v4.1.0_13-30388771-30404967_2026_01_07_13_34_21.csv', header = TRUE, sep = ',')

gnomAD$label <- ifelse(gnomAD$Allele.Frequency < 1e-4 & gnomAD$phylop >= 2, 1, 0)
gnomAD <- gnomAD[order(gnomAD$Allele.Frequency, decreasing = TRUE), ]

gnomAD <- subset(gnomAD, VEP.Annotation != '')


library(ggplot2)
library(scales)
pdf("LINKS_study/figures/gnomAD_plot.pdf")
gnomAD$VEP.Annotation <- factor(gnomAD$VEP.Annotation)

ggplot(gnomAD, aes(x = Allele.Frequency, y = phylop)) +

  # Background: all variants
  geom_point(
    colour = "grey80",
    alpha = 0.4,
    size = 1
  ) +

  # Highlight: rare + conserved variants
  geom_point(
    data = subset(gnomAD, Allele.Frequency <= 1e-4 & phylop >= 2),
    aes(colour = VEP.Annotation),
    size = 1.6,
    alpha = 0.9
  ) +

  scale_x_log10(labels = label_scientific()) +
  geom_vline(xintercept = 1e-4, linetype = "dashed") +
  geom_hline(yintercept = 2, linetype = "dashed") +

  theme_classic(base_size = 12) +
  labs(
    x = "Allele frequency (log scale)",
    y = "PhyloP score",
    colour = "VEP annotation"
  )

dev.off()

