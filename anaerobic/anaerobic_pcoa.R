####Beta Diversity
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
library(VennDiagram)
library(venneuler)
library(ggpubr)
library(ggsci)
library(patchwork)

###All
#Load pcoa table (pcoa axes file from mothur)
pcoa <- read.table(file=file.choose(), header=T,
                   comment="",
                   sep="\t", row.names = 1)

#Load metadata table 
meta <- read.table(file=file.choose(),
                   comment="",
                   header=TRUE,
                   sep="\t")

#Make a vector with place holders
new_names <- rep("", ncol(pcoa))
#Fill in first with PC followed by the number
for(i in 1:ncol(pcoa)){
  new_names[i] <- paste("PC", i, sep="")
}

#Replace the column names of pcoa
names(pcoa) <- new_names

#Add group column in pcoa table
pcoa$group = row.names(pcoa)

#merge the metadata and pcoa table
PCOA <- merge(pcoa, meta, by="group")

#Plot
ggplot(PCOA) +
  geom_point(size=4, aes(x=PC1, y =PC2, color=Treatments)) +
  labs(x="PC1 (80.3%)", y="PC2 (3.5%)") +
  stat_ellipse(alpha = 0.3, geom="polygon", linetype="blank",
               aes(x = PC1, y = PC2, fill = Type)) +
  scale_color_discrete(name = "Treatments",
                       labels = c("C12-HSL","C6-HSL", "Control")) +
  scale_fill_hue(h.start = 20,
                 name = "Type",
                 labels = c("Biofilm", "Planktonic")) +
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 12),
        axis.title=element_text(size=14),
        title = element_text(size=14),
        panel.background=element_rect(fill="white", colour="black",
                                      linetype="solid"),
        panel.grid.minor = element_line(size = 0.1, linetype="solid",
                                        colour="gray"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key= element_blank()) +
  guides(color = guide_legend(override.aes = list(fill = "grey", size = 4)),
         fill = guide_legend(override.aes=list(shape = NA)))


###Biofilm
#Load pcoa table for biofilm (pcoa axes file from mothur) 
pcoa_bio <- read.table(file=file.choose(), header=T,
                       sep="\t",
                       comment="", row.names=1)

#Make a vector with place holders
new_names1 <- rep("", ncol(pcoa_bio))

#Fill in first with PC followed by the number
for(i in 1:ncol(pcoa_bio)){
  new_names1[i] <- paste("PC", i, sep="")
}

#Replace the column names of pcoa
names(pcoa_bio) <- new_names1

#Add group column in pcoa table
pcoa_bio$group = row.names(pcoa_bio)

#Subset metadata to get biofilm
meta_bio <- subset(meta, Type == "Biofilm")

#merge the metadata and pcoa table
PCOA_bio <- merge(pcoa_bio, meta_bio, by="group")
PCOA_bio$title <- "Biofilm"

#Plot
###Figure 4C###
p2 <- ggplot(PCOA_bio) +
  geom_point(size=5, aes(x=PC1, y =PC2, color=Treatments)) +
  labs(x="PC1 (19.1%)", y="PC2 (13.1%)") +
  theme(panel.spacing.x=unit(0.1, "line"),
        plot.title=element_text(size=16),
        axis.text= element_text(size = 20, face="bold"),
        axis.title=element_text(size=22, face="bold"),
        strip.text.x=element_text(size=20, face="bold"),
        panel.background=element_rect(fill="white", colour="black",linetype="solid"),
        panel.grid.minor = element_line(size = 0.1, linetype="solid",
                                        colour="gray"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=1),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key= element_blank()) +
  guides(color = guide_legend(override.aes = list(fill = "white", size = 4)),
         fill = guide_legend(override.aes=list(shape = NA))) +
  facet_grid(. ~ title) +
  scale_color_manual(values=c("#FF66FF", "#00CCFF", "#666666"))

###Planktonic
#Load pcoa table for planktonic (pcoa axes file from mothur) 
pcoa_pl <- read.table(file=file.choose(), header=T,
                      sep="\t",
                      comment="", row.names=1)

#Make a vector with place holders
new_names2 <- rep("", ncol(pcoa_pl))

#Fill in first with PC followed by the number
for(i in 1:ncol(pcoa_pl)){
  new_names2[i] <- paste("PC", i, sep="")
}

#Replace the column names of pcoa
names(pcoa_pl) <- new_names2

#Add group column in pcoa table
pcoa_pl$group = row.names(pcoa_pl)

#Subset metadata to get biofilm
meta_pl <- subset(meta, Type == "Planktonic")

#merge the metadata and pcoa table
PCOA_pl <- merge(pcoa_pl, meta_pl, by="group")
PCOA_pl$title <- "Planktonic"

#Plot
###Figure 4E###
p3 <- ggplot(PCOA_pl) +
  geom_point(size=5, aes(x=PC1, y =PC2, color=Treatments)) +
  labs(x="PC1 (50.4%)", y="PC2 (12.8%)") +
  theme(panel.spacing.x=unit(0.1, "line"),
        plot.title=element_text(size=16),
        axis.text= element_text(size = 20, face="bold"),
        axis.title=element_text(size=22, face="bold"),
        strip.text.x=element_text(size=20, face="bold"),
        panel.background=element_rect(fill="white", colour="black",linetype="solid"),
        panel.grid.minor = element_line(size = 0.1, linetype="solid",
                                        colour="gray"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=1),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key= element_blank()) +
  guides(color = guide_legend(override.aes = list(fill = "white", size = 4)),
         fill = guide_legend(override.aes=list(shape = NA))) +
  facet_grid(. ~ title) +
  scale_color_manual(values=c("#FF66FF", "#00CCFF", "#666666"))

#Save pictures
PCOA <- ggarrange(p2, p3, ncol=2, nrow=1, align="h", common.legend=TRUE, 
                  legend="right")
ggsave("PCOA.jpg", PCOA, height=5, width=11, dpi=320)

