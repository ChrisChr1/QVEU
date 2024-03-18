library(ggrepel)
IRES.markers$delabel <- NA
IRES.markers$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
IRES.markers$diffexpressed[IRES.markers$avg_log2FC > 0.6 & IRES.markers$p_val_adj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
IRES.markers$diffexpressed[IRES.markers$avg_log2FC < -0.6 & IRES.markers$p_val_adj < 0.05] <- "DOWN"
IRES.markers$delabel[IRES.markers$diffexpressed != "NO"] <- IRES.markers$gene[IRES.markers$diffexpressed != "NO"]
#Create a column for the shape manual
IRES.markers$clusters <- factor(IRES.markers$cluster)
p1 <- ggplot(IRES.markers, aes(avg_log2FC, -log(p_val_adj,10), shape = factor(cluster), label = delabel)) + # -log10 conversion 
  geom_point(size=2) +
  xlab(expression("log"[2]*"(Fold Change)")) + 
  ylab(expression("-log"[10]*"Pvalue"))+
  geom_vline(xintercept=c(-0.6, 0.6), col="grey", linetype = 'dashed') +
  geom_hline(yintercept=-log10(0.05), col="grey", linetype = 'dashed')+
  coord_cartesian(ylim = c(0, 50)) +
  scale_color_brewer(palette="Paired") + theme(
    legend.position = c(1, 1),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
    )+
  labs(shape = 'Clusters') + geom_point(data=subset(IRES.markers, diffexpressed == "NO"), color="grey") +
  scale_shape_manual(values=1:nlevels(IRES.markers$clusters))
options(ggrepel.max.overlaps = 4)
p1 +
  facet_wrap(~factor(cluster))+  
  ylim(0,70) + labs(title = "Differential Expression across Clusters") + geom_text_repel()
