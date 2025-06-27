####Alpha Diversity####
#Load these libraries
library(ggplot2)
library(ggpubr)
library(patchwork)
library(ggsci)
library(cowplot)

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

#Set comparisons
my_comparisons <- list(c("5A8", "SsoPox"), c("5A8", "GcL"), c("SsoPox", "GcL"))

#Plot 
#Sobs
p1 <- ggplot(meta_alpha, aes(x=Treatments, y=sobs, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=3) +
  labs(x="Treatments", y="sobs") + 
  scale_x_discrete(labels=c("5A8" = "Control", "GcL"= "GcL", 
                             "SsoPox" = "SsoPox")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=420, size = 5, method="anova") + 
    theme(panel.spacing.x=unit(0.1, "line"),
          axis.text= element_text(size = 12),
          axis.title=element_text(size=14),
          title = element_text(size=14),
          panel.background=element_rect(fill="white", colour="grey",
                                        size=, linetype="solid"),
          panel.border = element_rect(colour="black", fill=NA,
                                      size=0.5), legend.position = "none")

#coverage
p2 <- ggplot(meta_alpha, aes(x=Treatments, y=coverage, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=3) +
  labs(x="Treatments", y="coverage") + 
  scale_x_discrete(labels=c("5A8" = "Control", "GcL"= "GcL", 
                            "SsoPox" = "SsoPox")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=1.01, size = 5, method="anova") + 
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        title = element_text(size=14),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none")

#Shannon
ggplot(meta_alpha, aes(x=Treatments, y=shannon, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=3) +
  labs(x="Treatments", y="shannon") + 
  scale_x_discrete(labels=c("5A8" = "Control", "GcL"= "GcL", 
                            "SsoPox" = "SsoPox")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=3.5, size = 5, method="kruskal.test") + 
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        title = element_text(size=14),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none") 
 
#chao
p4 <- ggplot(meta_alpha, aes(x=Treatments, y=chao, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=3) +
  labs(x="Treatments", y="chao") + 
  scale_x_discrete(labels=c("5A8" = "Control", "GcL"= "GcL", 
                            "SsoPox" = "SsoPox")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=2100, size = 5, method="kruskal.test") + 
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        title = element_text(size=14),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none") 

#ace
p5 <- ggplot(meta_alpha, aes(x=Treatments, y=ace, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=3) +
  labs(x="Treatments", y="ace") + 
  scale_x_discrete(labels=c("5A8" = "Control", "GcL"= "GcL", 
                            "SsoPox" = "SsoPox")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=5200, size = 5, method="kruskal.test") + 
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        title = element_text(size=14),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none") 

#simpson
p6<- ggplot(meta_alpha, aes(x=Treatments, y=simpson, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=3) +
  labs(x="Treatments", y="simpson") + 
  scale_x_discrete(labels=c("5A8" = "Control", "GcL"= "GcL", 
                            "SsoPox" = "SsoPox")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=0.5, size = 5, method="kruskal.test") + 
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
#Load Biofilm Alpha Diversity Table from Mothur
biofilm = read.table(file=file.choose(), header=TRUE, sep="\t")

#Remove the label column from the alpha table
biofilm.clean = biofilm[,-which(names(biofilm) %in% c("label"))]

#Merge alpha with metadata
bio = merge(meta, biofilm.clean, by="group")
bio$title = "Biofilm"


#Plot
#shannon
#Figure 3B##
p7 <- ggplot(bio, aes(x=Treatments, y=shannon, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=5) +
  labs(x="Treatments", y="Shannon") + 
  scale_x_discrete(labels=c("5A8" = "5A8", "SsoPox" = "Sso","GcL" = "GcL"),
                   limit = c("5A8", "SsoPox", "GcL"))+
  theme(panel.spacing.x=unit(0.1, "line"),
        plot.title=element_text(size=16),
        axis.text= element_text(size = 20, face="bold"),
        axis.title=element_text(size=22, face="bold"),
        strip.text.x=element_text(size=20, face="bold"),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=1), legend.position = "none") +
  facet_grid(. ~ title) + stat_compare_means(comparisons = my_comparisons, 
              label= "p.signif", method="wilcox.test",size=6, tip.length = 0) +
ylim(2.6, 3.2) + scale_color_manual(values=c("#666666", "#FF66FF","#00CCFF"))
ggsave("shannon_biofilm.jpg", p7, width=4, height=4, dpi=320)

#ace
p8<- ggplot(bio, aes(x=Treatments, y=ace, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=3) +
  labs(x="Treatments", y="ace") + 
  scale_x_discrete(labels=c("5A8" = "Control", "GcL"= "GcL", 
                            "SsoPox" = "SsoPox")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=8500, size = 5, method="anova") + 
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        title = element_text(size=14),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none")

#simpson
p9 <- ggplot(bio, aes(x=Treatments, y=simpson, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=3) +
  labs(x="Treatments", y="simpson") + 
  scale_x_discrete(labels=c("5A8" = "Control", "GcL"= "GcL", 
                            "SsoPox" = "SsoPox")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=0.18, size = 5, method="anova") + 
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        title = element_text(size=14),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none")

#sobs
p10 <- ggplot(bio, aes(x=Treatments, y=sobs, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=3) +
  labs(x="Treatments", y="sobs") + 
  scale_x_discrete(labels=c("5A8" = "Control", "GcL"= "GcL", 
                            "SsoPox" = "SsoPox")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=850, size = 5, method="kruskal.test") + 
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        title = element_text(size=14),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none")

#coverage
p11 <- ggplot(bio, aes(x=Treatments, y=coverage, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=3) +
  labs(x="Treatments", y="coverage") + 
  scale_x_discrete(labels=c("5A8" = "Control", "GcL"= "GcL", 
                            "SsoPox" = "SsoPox")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=0.995, size = 5, method="kruskal.test") + 
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        title = element_text(size=14),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none")

#chao
p12 <- ggplot(bio, aes(x=Treatments, y=chao, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=3) +
  labs(x="Treatments", y="chao") + 
  scale_x_discrete(labels=c("5A8" = "Control", "GcL"= "GcL", 
                            "SsoPox" = "SsoPox")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=4000, size = 5, method="kruskal.test") + 
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        title = element_text(size=14),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none")

###Planktonic
###Load planktonic alpha diversity table 
planktonic = read.table(file=file.choose(), header=TRUE, sep="\t")
#Remove the label column from the alpha table
planktonic.clean = planktonic[,-which(names(planktonic) %in% c("label"))]

#Merge planktonic alpha with metadata
pla = merge(meta, planktonic.clean, by="group")
pla$title = "Planktonic"


#Plot
p13 <- ggplot(pla, aes(x=Treatments, y=sobs, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=3) +
  labs(x="Treatments", y="sobs") + 
  scale_x_discrete(labels=c("5A8" = "Control", "GcL"= "GcL", 
                            "SsoPox" = "SsoPox")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=220, size = 5, method="anova") + 
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        title = element_text(size=14),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none")

##Figure 3D###
#shannon
p14 <- ggplot(pla, aes(x=Treatments, y=shannon, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=5) +
  labs(x="Treatments", y="Shannon") + 
  scale_x_discrete(labels=c("5A8" = "5A8", "SsoPox" = "Sso","GcL" = "GcL"),
                   limit = c("5A8", "SsoPox", "GcL"))+
  theme(panel.spacing.x=unit(0.1, "line"),
        plot.title=element_text(size=16),
        axis.text= element_text(size = 20, face="bold"),
        axis.title=element_text(size=22, face="bold"),
        strip.text.x=element_text(size=20, face="bold"),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=1), legend.position = "none") +
  facet_grid(. ~ title) + stat_compare_means(comparisons = my_comparisons, 
              label= "p.signif", method="wilcox.test", size=6, tip.length=0) + 
  ylim(1.4, 2.7) +
scale_color_manual(values=c("#666666", "#FF66FF","#00CCFF"))
ggsave("shannon_planktonic.jpg", p14, width=4, height=4, dpi=320)

#coverage
p15 <- ggplot(pla, aes(x=Treatments, y=coverage, color = Treatments)) +
  geom_boxplot(outlier.shape =NA) +
  geom_jitter(size=3) +
  labs(x="Treatments", y="coverage") + 
  scale_x_discrete(labels=c("5A8" = "Control", "GcL"= "GcL", 
                            "SsoPox" = "SsoPox")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=1, size = 5, method="kruskal.test") + 
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
  scale_x_discrete(labels=c("5A8" = "Control", "GcL"= "GcL", 
                            "SsoPox" = "SsoPox")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=550, size = 5, method="kruskal.test") + 
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
  scale_x_discrete(labels=c("5A8" = "Control", "GcL"= "GcL", 
                            "SsoPox" = "SsoPox")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=880, size = 5, method="kruskal.test") + 
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
  scale_x_discrete(labels=c("5A8" = "Control", "GcL"= "GcL", 
                            "SsoPox" = "SsoPox")) +
  stat_compare_means(comparisons = my_comparisons, label="p.signif", 
                     method="wilcox.test") +
  stat_compare_means(label.y=0.45, size = 5, method="kruskal.test") + 
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        title = element_text(size=14),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none")

alpha_planktonic <-ggarrange(p13, p14, p15, p16, p17, p18, labels=c("A", "B", "C", "D", "E", "F"),
                                               ncol =3, nrow=2)
ggsave("alpha_planktonic.jpg", width=14, height=12)

shannon <- ggarrange(p7, p14, align="h", ncol=2, nrow=2)
ggsave("shannon.jpg", shannon, height=10, width=10, dpi=320)
