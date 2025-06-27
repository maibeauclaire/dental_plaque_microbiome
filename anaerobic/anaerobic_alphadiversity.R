####Alpha Diversity####
#Load these libraries
library(ggplot2)
library(ggpubr)
library(patchwork)
library(ggsci)

#Load the alpha table (summary file from mothur)
alpha = read.table(file=file.choose(), header=TRUE, sep="\t")
#Remove the label column from the alpha table
alpha.clean = alpha[,-which(names(alpha) %in% c("label"))]

#Load meta table
meta <- read.table(file=file.choose(),
                   comment="",
                   header=TRUE,
                   sep="\t")
#Set seed (some processes will be relying on the random number generator.Set random seed to make the analysis reproducible)
set.seed(8765)


#Merge alpha with metadata
meta_alpha = merge(meta, alpha.clean, by="group")


#Plot
my_comparisons <- list(c("DMSO", "C6-HSL"), c("DMSO", "C12-HSL"), 
                       c("C6-HSL", "C12-HSL"))
#sobs
p1 <- ggplot(meta_alpha, aes(x=Treatments, y=sobs, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=3) +
  labs(x="Treatments", y="sobs") + 
  scale_x_discrete(labels=c("DMSO" = "Control", "C6-HSL"= "C6-HSL", 
                            "C12-HSL" = "C12-HSL")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=1500, size = 5, method="anova") + 
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        title = element_text(size=14),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none")

#shannon
p2 <- ggplot(meta_alpha, aes(x=Treatments, y=shannon, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=3) +
  labs(x="Treatments", y="shannon") + 
  scale_x_discrete(labels=c("DMSO" = "Control", "C6-HSL"= "C6-HSL", 
                            "C12-HSL" = "C12-HSL")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=3.52, size = 5, method="anova") + 
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        title = element_text(size=14),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none")

#coverage
p3 <- ggplot(meta_alpha, aes(x=Treatments, y=coverage, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=3) +
  labs(x="Treatments", y="coverage") + 
  scale_x_discrete(labels=c("DMSO" = "Control", "C6-HSL"= "C6-HSL", 
                            "C12-HSL" = "C12-HSL")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=0.993, size = 5, method="anova") + 
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        title = element_text(size=14),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none")

#ace
p4 <- ggplot(meta_alpha, aes(x=Treatments, y=ace, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=3) +
  labs(x="Treatments", y="ace") + 
  scale_x_discrete(labels=c("DMSO" = "Control", "C6-HSL"= "C6-HSL", 
                            "C12-HSL" = "C12-HSL")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=14000, size = 5, method="anova") + 
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        title = element_text(size=14),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none")

#chao
p5 <- ggplot(meta_alpha, aes(x=Treatments, y=chao, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=3) +
  labs(x="Treatments", y="chao") + 
  scale_x_discrete(labels=c("DMSO" = "Control", "C6-HSL"= "C6-HSL", 
                            "C12-HSL" = "C12-HSL")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=6500, size = 5, method="kruskal.test") + 
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        title = element_text(size=14),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none")
#simpson
p6 <- ggplot(meta_alpha, aes(x=Treatments, y=simpson, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=3) +
  labs(x="Treatments", y="simpson") + 
  scale_x_discrete(labels=c("DMSO" = "Control", "C6-HSL"= "C6-HSL", 
                            "C12-HSL" = "C12-HSL")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=0.075, size = 5, method="kruskal.test") + 
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        title = element_text(size=14),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none")

alpha_all <- ggarrange(p1, p2, p3, p4, p5, p6, labels=c("A", "B", "C", "D", "E", "F"),
                       ncol =3, nrow=2)
ggsave("alpha_all.jpg", width=14, height=12)

###Biofilm Samples
biofilm = read.table(file=file.choose(), header=TRUE, sep="\t")
#Remove the label column from the alpha table
biofilm.clean = biofilm[,-which(names(biofilm) %in% c("label"))]
#Merge alpha with metadata
bio = merge(meta, biofilm.clean, by="group")
bio$title = "Biofilm"


#Plot
#sobs
p7 <- ggplot(bio, aes(x=Treatments, y=sobs, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=3) +
  labs(x="Treatments", y="sobs") + 
  scale_x_discrete(labels=c("DMSO" = "Control", "C6-HSL"= "C6-HSL", 
                            "C12-HSL" = "C12-HSL")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=1600, size = 5, method="anova") + 
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        title = element_text(size=14),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none")

###Figure 4B###
#shannon
p8 <- ggplot(bio, aes(x=Treatments, y=shannon, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=5) +
  labs(x="Treatments", y="Shannon") + 
  scale_x_discrete(limit = c("DMSO", "C6-HSL", "C12-HSL"))+
  theme(panel.spacing.x=unit(0.1, "line"),
        plot.title=element_text(size=16),
        axis.text= element_text(size = 20, face="bold"),
        axis.title.x=element_blank(),
        axis.title=element_text(size=22, face="bold"),
        strip.text.x=element_text(size=20, face="bold"),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=1), legend.position = "none") +
  facet_grid(. ~ title) + stat_compare_means(comparisons = my_comparisons, 
                                             label= "p.signif", 
                                             method="wilcox.test",size=8) +
  ylim(3.35, 3.480) + 
  scale_color_manual(values=c("#FF66FF", "#00CCFF", "#666666"))

#coverage
p9 <- ggplot(bio, aes(x=Treatments, y=coverage, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=3) +
  labs(x="Treatments", y="coverage") + 
  scale_x_discrete(labels=c("DMSO" = "Control", "C6-HSL"= "C6-HSL", 
                            "C12-HSL" = "C12-HSL")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=0.9935, size = 5, method="anova") + 
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        title = element_text(size=14),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none")

#chao
p10 <- ggplot(bio, aes(x=Treatments, y=chao, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=3) +
  labs(x="Treatments", y="chao") + 
  scale_x_discrete(labels=c("DMSO" = "Control", "C6-HSL"= "C6-HSL", 
                            "C12-HSL" = "C12-HSL")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=6200, size = 5, method="anova") + 
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        title = element_text(size=14),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none")

#ace
p11 <- ggplot(bio, aes(x=Treatments, y=ace, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=3) +
  labs(x="Treatments", y="ace") + 
  scale_x_discrete(labels=c("DMSO" = "Control", "C6-HSL"= "C6-HSL", 
                            "C12-HSL" = "C12-HSL")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=15000, size = 5, method="anova") + 
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        title = element_text(size=14),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none")

#simpson
p12 <- ggplot(bio, aes(x=Treatments, y=simpson, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=3) +
  labs(x="Treatments", y="simpson") + 
  scale_x_discrete(labels=c("DMSO" = "Control", "C6-HSL"= "C6-HSL", 
                            "C12-HSL" = "C12-HSL")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=0.062, size = 5, method="anova") + 
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        title = element_text(size=14),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none")
alpha_biofilm <- ggarrange(p7, p8, p9, p10, p11, p12, 
                           labels=c("A", "B", "C", "D", "E", "F"),
                           ncol =3, nrow=2)
ggsave("alpha_biofilm.jpg", width=14, height=12)

###Planktonic
###Get planktonic Samples
planktonic = read.table(file=file.choose(), header=TRUE, sep="\t")
#Remove the label column from the alpha table
planktonic.clean = planktonic[,-which(names(planktonic) %in% c("label"))]

#Merge planktonic alpha with metadata
pla = merge(meta, planktonic.clean, by="group")
pla$title = "Planktonic"

#Plot
#sobs
p13 <- ggplot(pla, aes(x=Treatments, y=sobs, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=3) +
  labs(x="Treatments", y="sobs") + 
  scale_x_discrete(labels=c("DMSO" = "Control", "C6-HSL"= "C6-HSL", 
                            "C12-HSL" = "C12-HSL")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=1400, size = 5, method="anova") + 
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        title = element_text(size=14),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none")

###Figure 4D###
#shannon
p14 <- ggplot(pla, aes(x=Treatments, y=shannon, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=5) +
  labs(x="Treatments", y="Shannon") + 
  scale_x_discrete(limit = c("DMSO", "C6-HSL", "C12-HSL"))+
  theme(panel.spacing.x=unit(0.1, "line"),
        plot.title=element_text(size=16),
        axis.text= element_text(size = 20, face="bold"),
        axis.title.x=element_blank(),
        axis.title=element_text(size=22, face="bold"),
        strip.text.x=element_text(size=20, face="bold"),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=1), legend.position = "none") +
  facet_grid(. ~ title) + stat_compare_means(comparisons = my_comparisons, 
                                             label= "p.signif", 
                                             method="wilcox.test",size=8) +
  ylim(3.32, 3.51) +
  scale_color_manual(values=c( "#FF66FF", "#00CCFF", "#666666"))

shannon <- ggarrange(p8, p14, align="h", ncol=2, nrow=2)
ggsave("shannon.jpg", shannon, height=10, width=10, dpi=320)

#coverage
p15 <- ggplot(pla, aes(x=Treatments, y=coverage, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=3) +
  labs(x="Treatments", y="coverage") + 
  scale_x_discrete(labels=c("DMSO" = "Control", "C6-HSL"= "C6-HSL", 
                            "C12-HSL" = "C12-HSL")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=0.993, size = 5, method="anova") + 
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        title = element_text(size=14),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none")

#chao
p16 <- ggplot(pla, aes(x=Treatments, y=chao, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=3) +
  labs(x="Treatments", y="chao") + 
  scale_x_discrete(labels=c("DMSO" = "Control", "C6-HSL"= "C6-HSL", 
                            "C12-HSL" = "C12-HSL")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=4750, size = 5, method="anova") + 
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        title = element_text(size=14),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none")

#ace
p17 <- ggplot(pla, aes(x=Treatments, y=ace, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=3) +
  labs(x="Treatments", y="ace") + 
  scale_x_discrete(labels=c("DMSO" = "Control", "C6-HSL"= "C6-HSL", 
                            "C12-HSL" = "C12-HSL")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=11000, size = 5, method="anova") + 
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        title = element_text(size=14),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none")

#simpson
p18 <- ggplot(pla, aes(x=Treatments, y=simpson, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=3) +
  labs(x="Treatments", y="simpson") + 
  scale_x_discrete(labels=c("DMSO" = "Control", "C6-HSL"= "C6-HSL", 
                            "C12-HSL" = "C12-HSL")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=0.075, size = 5, method="anova") + 
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        title = element_text(size=14),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none")
alpha_planktonic <-ggarrange(p13, p14, p15, p16, p17, p18, 
                             labels=c("A", "B", "C", "D", "E", "F"),
                             ncol =3, nrow=2)
ggsave("alpha_planktonic.jpg", width=14, height=12)
