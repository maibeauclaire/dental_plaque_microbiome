####Taxa Summary####
#Load these libraries
#Load these libraries
library(ape)
library(dplyr)
library(ggplot2)
library(lme4)
library(phangorn)
library(phyloseq)
library(plotly)
library(tidyr)
library(vegan)
library(ggpubr)
library(ggsci)
library(patchwork)
library(forcats)
library(ggh4x)
library(RColorBrewer)
library(ggtext)

#Load OTU table (shared file from mothur)
OTU = read.table(file=file.choose(), header=TRUE, sep="\t")
#Use the "Group" column as row names in OTU
row.names(OTU) = OTU$Group
#Remove "label", "numOTUs", and "Group"
OTU.clean = OTU[,-which(names(OTU) %in% c("label", "numOtus", "Group"))]
#Transpose OTU table
OTU.t=t(OTU.clean)

#Load taxonomy table (taxonomy file from mothur)
tax = read.table(file=file.choose(), header=TRUE, sep="\t")
#Use the "OTU" column as row names in taxa table
row.names(tax) = tax$OTU
#remove all the OTUs that don't occur in our OTU.clean dataset
tax.clean = tax[row.names(tax) %in% colnames(OTU.clean),]
#seperate the "taxonomy" column so that each level is in its own column
tax.clean = separate(tax.clean, Taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
                     sep=";")
#Remove size, OTU and Species columns in taxa table
tax.clean = tax.clean[,-which(names(tax.clean) %in% c("Size", "OTU", "Species"))]

#Load metadata table
meta <- read.table(file=file.choose(),
                   comment="",
                   header=TRUE,
                   sep="\t")

#Use the "group" column as row names in meta
row.names(meta) = meta$group
#Remove group column in metadata table
meta = meta[,-which(names(meta) %in% c("group"))]

###Create physeq object
#Tell R which tables are each type
OTU.UF = otu_table(as.matrix(OTU.t), taxa_are_rows=TRUE)
tax.UF = tax_table(as.matrix(tax.clean))
meta.UF = sample_data(meta)

#Merge these tables into an object of class phyloseq
physeq=phyloseq(OTU.UF, tax.UF, meta.UF)

#Keep taxa that has the mean greater than 0.1
physeq2 = filter_taxa(physeq, function(x) mean(x) > 0.1, TRUE)

#Transform into relative abundances
physeq3 = transform_sample_counts(physeq2, function(x) x / sum(x) )

###All
##Genus Level
#Agglomerate taxa at genus level
all <- tax_glom(physeq3, taxrank='Genus')
#Create data frame from phyloseq object
all <- psmelt(all)
#Rename genus with < 5% abundance
all$Genus[all$Abundance < 0.05] <- "Others"
#Remove (100) in Genus names
all$Genus<- gsub("[(100)]","", all$Genus)

#Italicize genus species
all$Genus <- with(all, ifelse(Genus == "Others", "Others", 
                              paste0("*", Genus, "*")))
#New names for treatments
treats <- c("C12-HSL", "C6-HSL", "DMSO")
names(treats) <- c("C12-HSL", "C6-HSL", "DMSO")

#Plot
###Figure S4B###
#Samples
p1<- ggplot(all,
       aes(x = Sample, y = Abundance, 
           fill=fct_reorder(Genus, Abundance, .desc=TRUE))) +
  geom_bar(stat="identity", position='stack') +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), legend.position="bottom",
        panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 16),
        axis.title=element_text(size=20),
        strip.text.x=element_text(size=18),
        legend.title=element_text(size=16),
        legend.text=element_text(size=12),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=1))+
  labs(y="Relative Abundance", fill="Genus") +
  scale_fill_npg() +
  facet_nested(~ Type + Treatments, scale= "free",
               labeller=labeller(Treatments= treats),
               nest_line=element_line(colour="black"))

ggsave("taxasum_samples_anaerobic.jpg", p1, width=9.75, height=5, dpi=320)

###Figure 4A###
#Treatments
p2 <- ggplot(all,
       aes(x = Treatments, y = Abundance, 
           fill=fct_reorder(Genus, Abundance, .desc=TRUE))) +
  geom_bar(stat="identity", position='fill') +
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 18, face="bold"),
        axis.title=element_text(size=20, face="bold"),
        axis.title.x=element_blank(),
        strip.text.x=element_text(size=18, face="bold"),
        legend.title=element_text(size=16),
        legend.text=element_text(size=12),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=1))+
  labs(x="Treatments", y="Relative Abundance", fill="Genus") +
  scale_fill_npg() +
  facet_wrap("Type") +
  scale_x_discrete(limit=c("DMSO", "C6-HSL", "C12-HSL"))
ggsave("taxa_summary.jpg", p2, height=5, width=10, dpi=320)

##Phylum Level
#Agglomerate taxa at phylum level
phylum <- tax_glom(physeq3, taxrank='Phylum')
#Create data frame from phyloseq object
phylum_all <- psmelt(phylum)
#Rename Phylum with < 5% abundance
phylum_all$Phylum[phylum_all$Abundance < 0.05] <- "Others"
#Remove (100) in Phylum names
phylum_all$Phylum<- gsub("[(100)]","", phylum_all$Phylum)

#Plot
#By Samples
p3 <- ggplot(phylum_all,
       aes(x = Sample, y = Abundance, 
           fill=fct_reorder(Phylum, Abundance, .desc=TRUE))) +
  geom_bar(stat="identity", position='stack') +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), legend.position="bottom",
        panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        strip.text.x=element_text(size=12),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=1))+
  labs(y="Relative Abundance", fill="Phylum") +
  scale_fill_npg() +
  facet_nested(~ Type + Treatments, scale= "free",labeller=
                 labeller(Treatments= treats),
               nest_line=element_line(colour="black"))
ggsave("p3.jpg", width=7.46, height=4)

#By treatments
p4 <- ggplot(phylum_all,
       aes(x = Treatments, y = Abundance, 
           fill=fct_reorder(Phylum, Abundance, .desc=TRUE))) +
  geom_bar(stat="identity", position='fill') +
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        strip.text.x=element_text(size=12),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=1))+
  labs(x="Treatments", y="Relative Abundance", fill="Phylum") +
  scale_fill_npg() +
  facet_wrap("Type") +
  scale_x_discrete(labels=c("C12-HSL" = "C12-HSL", "C6-HSL" = "C6-HSL",
                            "DMSO" = "Control"))
ggsave("p4.jpg", width=6.46)


##Family Level
#Agglomerate taxa at Family level
fam <- tax_glom(physeq3, taxrank='Family')
#Create data frame from phyloseq object
fam_all <- psmelt(fam)
#Rename Family with < 5% abundance
fam_all$Family[fam_all$Abundance < 0.05] <- "Others"
#Remove (100) in Family names
fam_all$Family<- gsub("[(100)]","", fam_all$Family)

#Plot
#By Samples
p5 <- ggplot(fam_all,
       aes(x = Sample, y = Abundance, 
           fill=fct_reorder(Family, Abundance, .desc=TRUE))) +
  geom_bar(stat="identity", position='stack') +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), legend.position="bottom",
        panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        strip.text.x=element_text(size=12),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=1))+
  labs(y="Relative Abundance", fill="Family") +
  scale_fill_npg() +
  facet_nested(~ Type + Treatments, scale= "free",labeller=
                 labeller(Treatments= treats),
               nest_line=element_line(colour="black"))
ggsave("p5.jpg", width=9, height=4)

#By treatments
p6 <- ggplot(fam_all,
             aes(x = Treatments, y = Abundance, 
                 fill=fct_reorder(Family, Abundance, .desc=TRUE))) +
  geom_bar(stat="identity", position='fill') +
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        strip.text.x=element_text(size=12),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=1))+
  labs(x="Treatments", y="Relative Abundance", fill="Family") +
  scale_fill_npg() +
  facet_wrap("Type") +
  scale_x_discrete(labels=c("C12-HSL" = "C12-HSL", "C6-HSL" = "C6-HSL",
                            "DMSO" = "Control"))
ggsave("p6.jpg", width=6.46)

##Class Level
#Agglomerate taxa at Class level
class <- tax_glom(physeq3, taxrank='Class')
#Create data frame from phyloseq object
class_all <- psmelt(class)
#Rename Class with < 1% abundance
class_all$Class[class_all$Abundance < 0.05] <- "Others"
#Remove (100) in Class names
class_all$Class<- gsub("[(100)]","", class_all$Class)

#Plot
#By Samples
p7 <- ggplot(class_all,
             aes(x = Sample, y = Abundance, 
                 fill=fct_reorder(Class, Abundance, .desc=TRUE))) +
  geom_bar(stat="identity", position='stack') +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), legend.position="bottom",
        panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        strip.text.x=element_text(size=12),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=1))+
  labs(y="Relative Abundance", fill="Class") +
  scale_fill_npg() +
  facet_nested(~ Type + Treatments, scale= "free",labeller=
                 labeller(Treatments= treats),
               nest_line=element_line(colour="black"))
ggsave("p7.jpg", width=7.46, height=4)

#By treatments
p8 <- ggplot(class_all,
       aes(x = Treatments, y = Abundance, 
           fill=fct_reorder(Class, Abundance, .desc=TRUE))) +
  geom_bar(stat="identity", position='fill') +
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        strip.text.x=element_text(size=12),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=1))+
  labs(x="Treatments", y="Relative Abundance", fill="Class") +
  scale_fill_npg() +
  facet_wrap("Type") +
  scale_x_discrete(labels=c("C12-HSL" = "C12-HSL", "C6-HSL" = "C6-HSL",
                            "DMSO" = "Control"))
ggsave("p8.jpg", width=6.46)

##Order Level
#Agglomerate taxa at Order level
order <- tax_glom(physeq3, taxrank='Order')
#Create data frame from phyloseq object
order_all <- psmelt(order)
#Rename Class with < 1% abundance
order_all$Order[order_all$Abundance < 0.05] <- "Others"
#Remove (100) in Class names
order_all$Order<- gsub("[(100)]","", order_all$Order)

#Plot
#By samples
p9 <- ggplot(order_all,
             aes(x = Sample, y = Abundance, 
                 fill=fct_reorder(Order, Abundance, .desc=TRUE))) +
  geom_bar(stat="identity", position='stack') +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), legend.position="bottom",
        panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        strip.text.x=element_text(size=12),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=1))+
  labs(y="Relative Abundance", fill="Order") +
  scale_fill_npg() +
  facet_nested(~ Type + Treatments, scale= "free",labeller=
                 labeller(Treatments= treats),
               nest_line=element_line(colour="black"))
ggsave("p9.jpg", width=7.46, height=4)

#By treatments
p10 <- ggplot(order_all,
       aes(x = Treatments, y = Abundance, 
           fill=fct_reorder(Order, Abundance, .desc=TRUE))) +
  geom_bar(stat="identity", position='fill') +
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        strip.text.x=element_text(size=12),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=1))+
  labs(x="Treatments", y="Relative Abundance", fill="Order") +
  scale_fill_npg() +
  facet_wrap("Type") +
  scale_x_discrete(labels=c("C12-HSL" = "C12-HSL", "C6-HSL" = "C6-HSL",
                            "DMSO" = "Control"))
ggsave("p10.jpg", width=6.46)

#Set comparisons pairs
my_comparisons = list(c("C12-HSL", "C6-HSL"), c("C12-HSL", "DMSO"), 
                      c("C6-HSL", "DMSO"))
##Plot significant taxa
#Subset Biofilm
biofilm <- subset(all, Type == "Biofilm")

#Plot with relative abundance
p11 <- ggplot(biofilm, aes(x=Treatments, y= Abundance, color= Treatments)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~Genus) +
  labs(x="Treatments", y="Relatve Abundance") +
  scale_x_discrete(labels=c("C12-HSL" = "C12-HSL", "C6-HSL" = "C6-HSL",
                            "DMSO" = "DMSO")) +
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        strip.text.x=element_text(size=12),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  stat_compare_means(label.y= 0.25)
ggsave("p11.jpg", width = 9, height = 9)

#Subset Fusobacterium
Fus <- subset(biofilm, Genus == "Fusobacterium")
#Plot
p12 <- ggplot(Fus, aes(x=Treatments, y= Abundance, color= Treatments)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=5) +
  labs(x="Treatments", y="Relatve Abundance", title="Fusobacterium") +
  scale_x_discrete(labels=c("C12-HSL" = "C12-HSL", "C6-HSL" = "C6-HSL",
                            "DMSO" = "DMSO"), 
                   limit=c("DMSO", "C6-HSL", "C12-HSL")) +
  theme(panel.spacing.x=unit(0.1, "line"),
        plot.title=element_text(size=20),
        axis.text= element_text(size = 20, face="bold"),
        axis.title=element_text(size=22, face="bold"),
        axis.title.x=element_blank(),
        strip.text.x=element_text(size=20, face="bold"),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=1), legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", size=8, method="wilcox.test")+
  ylim(0.157, 0.173) +
  scale_color_manual(values=c("#FF66FF", "#00CCFF", "#666666"))

#Subset Porphyromonas
Por <- subset(biofilm, Genus == "Porphyromonas")
#Plot
p13 <- ggplot(Por, aes(x=Treatments, y= Abundance, color= Treatments)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=5) +
  labs(x="Treatments", y="Relatve Abundance", title="Porphyromonas") +
  scale_x_discrete(labels=c("C12-HSL" = "C12-HSL", "C6-HSL" = "C6-HSL",
                            "DMSO" = "DMSO"), 
                   limit=c("DMSO", "C6-HSL", "C12-HSL")) +
  theme(panel.spacing.x=unit(0.1, "line"),
        plot.title=element_text(size=20),
        axis.text= element_text(size = 20, face="bold"),
        axis.title=element_text(size=22, face="bold"),
        axis.title.x=element_blank(),
        strip.text.x=element_text(size=20, face="bold"),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=1), legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", size=8) +
  ylim(0.065, 0.085)+
  scale_color_manual(values=c("#FF66FF", "#00CCFF", "#666666"))

#Subset Veillonella
Vei <- subset(biofilm, Genus == "Veillonella")
#Plot
p14 <- ggplot(Vei, aes(x=Treatments, y= Abundance, color= Treatments)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=5) +
  labs(x="Treatments", y="Relatve Abundance", title="Veillonella") +
  scale_x_discrete(labels=c("C12-HSL" = "C12-HSL", "C6-HSL" = "C6-HSL",
                            "DMSO" = "DMSO"), 
                   limit=c("DMSO", "C6-HSL", "C12-HSL")) +
  theme(panel.spacing.x=unit(0.1, "line"),
        plot.title=element_text(size=20),
        axis.text= element_text(size = 20, face="bold"),
        axis.title=element_text(size=22, face="bold"),
        axis.title.x=element_blank(),
        strip.text.x=element_text(size=20, face="bold"),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=1), legend.position = "none") +
  stat_compare_means(comparison = my_comparisons, label = "p.signif", size=8) +
  ylim(0.082, 0.094) +
    scale_color_manual(values=c("#FF66FF", "#00CCFF", "#666666"))

###Figure S6C###
#Combine mutiple ggplots on one page 
sig_taxa_bio_an <- ggarrange(p12, p13, p14, align="h", ncol=3)

#Save picture
ggsave("sig_taxa_bio_anaerobic.jpg", sig_taxa_bio_an, 
       width = 16, height = 5, dpi=320)

#Planktonic
pla <- subset(all, Type == "Planktonic")

#Plot
ggplot(pla, aes(x=Treatments, y= Abundance, color= Treatments)) +
  geom_boxplot() +
  geom_jitter(size=5) +
  facet_wrap(~Genus) +
  labs(x="Treatments", y="Relatve Abundance") +
  scale_x_discrete(labels=c("C12-HSL" = "C12-HSL", "C6-HSL" = "C6-HSL",
                            "DMSO" = "Control")) +
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        strip.text.x=element_text(size=12),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")

ggsave("p15.jpg", width = 9, height = 9)

#Subset Streptococcus
Strep <- subset(pla, Genus == "*Streptococcus*")

####Figure S6D###
#Plot
p16 <- ggplot(Strep, aes(x=Treatments, y= Abundance, color= Treatments)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=5) +
  labs(x="Treatments", y="Relatve Abundance", title="Streptococcus") +
  scale_x_discrete(labels=c("C12-HSL" = "C12-HSL", "C6-HSL" = "C6-HSL",
                            "DMSO" = "DMSO"), 
                   limit=c("DMSO", "C6-HSL", "C12-HSL")) +
  theme(panel.spacing.x=unit(0.1, "line"),
        plot.title=element_text(size=20, face="italic"),
        axis.text= element_text(size = 20, face="bold"),
        axis.title=element_text(size=22, face="bold"),
        axis.title.x=element_blank(),
        strip.text.x=element_text(size=20, face="bold"),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=1), legend.position = "none") +
  stat_compare_means(comparison = my_comparisons, label = "p.signif", size=8, tip.length = 0) +
  ylim(0.080, 0.1) +
  scale_color_manual(values=c("#159772", "#d7600b", "#716cb0"))

sig_taxa_an <- ggarrange(p12,p13,p14,p16)
ggsave("sig_taxa_an.jpg", width= 10, height =10)
ggsave("sig_taxa_anaerobic_planktonic.jpg", p16, width=6, height = 6, dpi=320)
