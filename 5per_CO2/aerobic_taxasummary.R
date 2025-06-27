####Taxa Summary####
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
library(tidyr)

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

#Remove (100)
tax.clean$Genus = gsub("[(100)]", "", tax.clean$Genus)

#Load metadata table
meta <- read.table(file=file.choose(),
                   comment="",
                   header=TRUE,
                   sep="\t") %>%
  drop_na(group)

#Use the "group" column as row names in meta
row.names(meta) = meta$group

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
glom <- tax_glom(physeq3, taxrank='Genus')
#Create data frame from phyloseq object
all <- psmelt(glom)
#Rename genus with < 5% abundance
all$Genus[all$Abundance < 0.05] <- "Others"
#Remove (100) in Genus names
all$Genus<- gsub("[(100)]","", all$Genus)
#Transform
all <- transform(all, Genus=as.character(Genus))

#Rename Schaalia to Actinomyces
all$Genus <- gsub("Schaalia", "Actinomyces", all$Genus)
#Rename Lactobacillales
all$Genus <-gsub("Lactobacillales_unclassified", "Lactobacillales (Order)", all$Genus)


#New names for treatments
treats <- c("5A8", "Ssopox", "GcL")
names(treats) <- c("5A8", "SsoPox", "GcL")

#Italicize genus species
all$Genus <- with(all, ifelse(Genus == "Others", "Others", 
                              paste0("*", Genus, "*")))
all$Genus[all$Genus == "*Lactobacillales (Order)*"] <- "Lactobacillales*"


##Plot
###Figure S4A###
p1 <- ggplot(all,
       aes(x = Sample, y = Abundance, 
           fill=fct_reorder(Genus, Abundance, .desc=TRUE))) +
  geom_bar(stat="identity", position='stack') +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), legend.position="bottom",
        panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 16, face="bold"),
        axis.title=element_text(size=18, face="bold"),
        strip.text.x=element_text(size=16, face="bold"),
        legend.title=element_text(size=14),
        legend.text=element_markdown(size = 12),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=1))+
  labs(y="Relative Abundance", fill="Species") +
  scale_fill_brewer(palette = "Dark2") +
  facet_nested(~ Type + Treatments, scale= "free",
               labeller=labeller(Treatments= treats),
               nest_line=element_line(colour="black"))
ggsave("taxasum_samples.jpg", p1, width=7.46, height=4, dpi=320)

###Figure 3A###
p2 <- ggplot(all,
            aes(x = Treatments, y = Abundance, 
                fill=fct_reorder(Genus, Abundance, .desc=TRUE))) +
  geom_bar(stat="identity", position='fill') +
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 16, face="bold"),
        axis.title.x = element_blank(),
        axis.title=element_text(size=18, face="bold"),
        strip.text.x=element_text(size=16, face="bold"),
        legend.title=element_text(size=14),
        legend.text=element_markdown(size = 12),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=1))+
  labs(y="Relative Abundance", fill="Taxa") +
  scale_fill_brewer(palette = "Dark2") +
  facet_wrap("Type") +
  scale_x_discrete(labels=c("5A8" = "5A8", "SsoPox" = "Sso","GcL" = "GcL"),
                   limit = c("5A8", "SsoPox", "GcL"))
ggsave("taxasum_treatments.jpg", p2, height=3.25, width=7.46, dpi=320)

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
#Samples
p3 <-ggplot(phylum_all,
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
  facet_nested(~ Type + Treatments, scale= "free", labeller=
                 labeller(Treatments= treats),
               nest_line=element_line(colour="black"))
ggsave("p3.jpg", width=7.46, height=4)

#Treatments
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
  scale_x_discrete(labels=c("5A8" = "Control", "GcL" = "GcL",
                            "SsoPox" = "SsoPox"))
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
#Samples
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
  facet_nested(~ Type + Treatments, scale= "free", labeller=
                 labeller(Treatments= treats),
               nest_line=element_line(colour="black"))
ggsave("p5.jpg", width=7.46, height=4)

#Treatments
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
  scale_x_discrete(labels=c("5A8" = "Control", "GcL" = "GcL",
                            "SsoPox" = "SsoPox"))
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
#Samples
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
  facet_nested(~ Type + Treatments, scale= "free", labeller=
                 labeller(Treatments= treats),
               nest_line=element_line(colour="black"))
ggsave("p7.jpg", width=7.46, height=4)

#Treatments
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
  scale_x_discrete(labels=c("5A8" = "Control", "GcL" = "GcL",
                            "SsoPox" = "SsoPox"))
ggsave("p8.jpg", width=6.46)

##Order Level
#Agglomerate taxa at Order level
order <- tax_glom(physeq3, taxrank='Order')
#Create data frame from phyloseq object
order_all <- psmelt(order)
#Rename Class with < 5% abundance
order_all$Order[order_all$Abundance < 0.05] <- "Others"
#Remove (100) in Class names
order_all$Order<- gsub("[(100)]","", order_all$Order)

#Plot
#Samples
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
  facet_nested(~ Type + Treatments, scale= "free", labeller=
                 labeller(Treatments= treats),
               nest_line=element_line(colour="black"))
ggsave("p9.jpg", width=7.46, height=4)

#Treatments
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
  scale_x_discrete(labels=c("5A8" = "Control", "GcL" = "GcL",
                            "SsoPox" = "SsoPox"))
ggsave("p10.jpg", width=6.46)

#Set comparisons pairs
my_comparisons = list(c("5A8", "SsoPox"), c("5A8", "GcL"), c("SsoPox", "GcL"))

##Plot significant taxa
#Subset Biofilm
biofilm <- subset(all, Growth == "Biofilm")

#Plot with relative abundance
ggplot(biofilm, aes(x=Treatment, y= Abundance, color= Treatment)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~Genus) +
  labs(x="Treatments", y="Relatve Abundance") +
  scale_x_discrete(labels=c("5A8" = "Control", "GcL" = "GcL",
                            "SsoPox" = "SsoPox")) +
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        strip.text.x=element_text(size=12),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(label.y= 1)

#Subset Lactobacillales
Lacto <- subset(biofilm, Genus == "Lactobacillales*")

#Plot Lactobicallales
p12 <- ggplot(Lacto, aes(x=Treatments, y= Abundance, color= Treatments)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=5) +
  labs(x="Treatments", y="Relatve Abundance", title="Lactobacillales") +
  scale_x_discrete(labels=c("5A8" = "5A8", "SsoPox" = "Sso","GcL" = "GcL"),
                   limit = c("5A8", "SsoPox", "GcL")) +
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
  stat_compare_means(comparisons = my_comparisons, label= "p.signif", size=8, method="wilcox.test") +
  ylim(0.16, 0.255) +
  scale_color_manual(values=c("#666666", "#FF66FF","#00CCFF"))

#Subset Schaalia
Scha <- subset(biofilm, Genus == "*Actinomyces*")

#Plot
p13 <- ggplot(Scha, aes(x=Treatments, y= Abundance, color= Treatments)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=5) +
  labs(x="Treatments", y="Relatve Abundance", title="Actinomyces") +
  scale_x_discrete(labels=c("5A8" = "5A8", "SsoPox" = "Sso","GcL" = "GcL"),
                   limit = c("5A8", "SsoPox", "GcL"))+
  theme(panel.spacing.x=unit(0.1, "line"),
        plot.title=element_markdown(size=20),
        axis.text= element_text(size = 20, face="bold"),
        axis.title=element_text(size=22, face="bold"),
        axis.title.x=element_blank(),
        strip.text.x=element_text(size=20, face="bold"),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=1), legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, label= "p.signif", size=8, method="wilcox.test") +
  ylim(0.07, 0.15) +
  scale_color_manual(values=c("#666666", "#FF66FF","#00CCFF"))

#Subset Streptococcus
strep <- subset(biofilm, Genus == "*Streptococcus*")

#Plot
p14 <- ggplot(strep, aes(x=Treatment, y= Abundance, color= Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=5) +
  labs(x="Treatments", y="Relatve Abundance", title="Streptococcus") +
  scale_x_discrete(labels=c("5A8" = "5A8", "SsoPox" = "Sso","GcL" = "GcL"),
                    limit = c("5A8", "SsoPox", "GcL")) +
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
  stat_compare_means(comparisons = my_comparisons, label= "p.signif", size=8) +
  ylim(0.45, 0.60) +
  scale_color_manual(values=c("#666666", "#FF66FF","#00CCFF"))


###Figure S6A###
sig_taxa_bio <- ggarrange(p12, p13, p14, align="h", ncol=3)
                                
#Save picture
ggsave("sig_taxa_bio.jpg", sig_taxa_bio, width = 16, height = 5, dpi=320)

#Planktonic
pla <- subset(all, Type == "Planktonic")

#Plot 
p15 <- ggplot(pla, aes(x=Treatments, y= Abundance, color= Treatments)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~Genus) +
  labs(x="Treatments", y="Relatve Abundance") +
  scale_x_discrete(labels=c("5A8" = "Control", "GcL" = "GcL",
                            "SsoPox" = "SsoPox")) +
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        axis.title.x=element_blank(),
        strip.text.x=element_text(size=12),
        panel.background=element_rect(fill="white", colour="grey",
                                      size=, linetype="solid"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5), legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  stat_compare_means(label.y= 1.1)
ggsave("p15.jpg", width = 8, height = 6)
#Subset Abiotrophia
abi <- subset(pla, Genus == "*Abiotrophia*")

#Subset Lactobacillales
Lacto_pla <- subset(pla, Genus == "Lactobacillales (Order)")

#Plot Lactobicallales
p16 <- ggplot(Lacto_pla, aes(x=Treatments, y= Abundance, color= Treatments)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=5) +
  labs(x="Treatments", y="Relatve Abundance", title="Lactobacillales") +
  scale_x_discrete(labels=c("5A8" = "5A8", "SsoPox" = "Sso","GcL" = "GcL"),
                   limit = c("5A8", "SsoPox", "GcL")) +
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
  stat_compare_means(comparisons = my_comparisons, label= "p.signif", size=8) +
  ylim(0.1, 0.42) +
  scale_color_manual(values=c("#666666", "#FF66FF","#00CCFF"))

#Subset Schaalia
Scha_pla <- subset(pla, Genus == "Actinomyces")

#Plot
p17 <- ggplot(Scha_pla, aes(x=Treatments, y= Abundance, color= Treatments)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=5) +
  labs(x="Treatments", y="Relatve Abundance", title="Actinomyces") +
  scale_x_discrete(labels=c("5A8" = "5A8", "SsoPox" = "Sso","GcL" = "GcL"),
                   limit = c("5A8", "SsoPox", "GcL"))+
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
  stat_compare_means(comparisons = my_comparisons, label= "p.signif", size=8) +
  ylim(0.05, 0.13) +
  scale_color_manual(values=c("#666666", "#FF66FF","#00CCFF"))


#Subset Streptococcus
strep_pla <- subset(pla, Genus == "Streptococcus")

#Plot
p18 <- ggplot(strep_pla, aes(x=Treatments, y= Abundance, color= Treatments)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=5) +
  labs(x="Treatments", y="Relatve Abundance", title="Streptococcus") +
  scale_x_discrete(labels=c("5A8" = "5A8", "SsoPox" = "Sso","GcL" = "GcL"),
                   limit = c("5A8", "SsoPox", "GcL")) +
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
  stat_compare_means(comparisons = my_comparisons, label= "p.signif", size=8) +
  ylim(0.35, 0.85) +
  scale_color_manual(values=c("#666666", "#FF66FF","#00CCFF"))

###Figure S6B###
sig_taxa_pla <- ggarrange(p16, p17, p18, align="h", ncol=3)

#Save picture
ggsave("sig_taxa_pla.jpg", sig_taxa_pla, width = 16, height = 5, dpi=320)


