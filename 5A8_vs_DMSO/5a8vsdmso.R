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
library(venneuler)
library(ggpubr)
library(ggsci)
library(patchwork)
library(glue)
library(tidyverse)
library(ggtext)

######PCOA##########################################################
###Biofilm
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
PCOA$title <- "Biofilm"


#Plot
p1 <- ggplot(PCOA, aes(x=PC1, y =PC2, color=Condition, linetype=Condition)) +
  geom_point(size=5) + 
  stat_ellipse(linetype=2) +
  labs(x="PC1 (99.4%)", y="PC2 (0.22%)") +
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 20, face="bold"),
        axis.title=element_text(size=22, face="bold"),
        strip.text.x=element_text(size=20, face="bold"),
        panel.background=element_rect(fill="white", colour="black",
                                      linetype="solid"),
        panel.grid.minor = element_line(size = 0.1, linetype="solid",
                                        colour="gray"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=1),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key= element_blank()) +
  guides(color = guide_legend(override.aes = list(fill = "grey", size = 4)),
         fill = guide_legend(override.aes=list(shape = NA))) +
  facet_grid(.~title) +
  scale_color_manual(values=c("#666666", "#7030A0"))

###Planktonic
#Load pcoa table (pcoa axes file from mothur)
pcoa1 <- read.table(file=file.choose(), header=T,
                   comment="",
                   sep="\t", row.names = 1)

#Make a vector with place holders
new_names <- rep("", ncol(pcoa1))

#Fill in first with PC followed by the number
for(i in 1:ncol(pcoa1)){
  new_names[i] <- paste("PC", i, sep="")
}

#Replace the column names of pcoa
names(pcoa1) <- new_names

#Add group column in pcoa table
pcoa1$group = row.names(pcoa1)

#merge the metadata and pcoa table
PCOA1 <- merge(pcoa1, meta, by="group")
PCOA1$title <- "Planktonic"  

#Plot
p2 <- ggplot(PCOA1, aes(x=PC1, y=PC2, color=Condition, linetype=Condition)) +
  geom_point(size=5) + 
  stat_ellipse(linetype=2) +
  labs(x="PC1 (92.9%)", y="PC2 (5.07%)") +
  theme(panel.spacing.x=unit(0.1, "line"),
        axis.text= element_text(size = 20, face="bold"),
        axis.title=element_text(size=22, face="bold"),
        strip.text.x=element_text(size=20, face="bold"),
        panel.background=element_rect(fill="white", colour="black",
                                      linetype="solid"),
        panel.grid.minor = element_line(size = 0.1, linetype="solid",
                                        colour="gray"),
        panel.border = element_rect(colour="black", fill=NA,
                                    size=1),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key= element_blank()) +
  guides(color = guide_legend(override.aes = list(fill = "grey", size = 4)),
         fill = guide_legend(override.aes=list(shape = NA))) +
  facet_grid(.~title) +
  scale_color_manual(values=c("#666666", "#7030A0"))


###Figure 2A###
#Align plots
PCOA_fin <- ggarrange(p1, p2, ncol=2, nrow=1, align="h", common.legend=TRUE, 
                  legend="right")
ggsave("PCOA_controls.jpg", PCOA_fin, height=4, width=10, dpi=320)

###lefse
##Planktonic
#Load taxonomy file
pl_taxonomy <- read_tsv("Planktonic/dental.final.an.0.02.cons.taxonomy") %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy,
           into=c("kingdom", "phylum", "class", "order", "family", "genus"),
           sep=";") %>%
  mutate(pretty_otu = str_replace(string=otu,
                                  pattern="tu0*",
                                  replacement = "TU "),
         genus = str_replace(string=genus,
                             pattern="(.*)",
                             replacement="*\\1*"),
         genus = str_replace(string=genus,
                             pattern="\\*(.*)_unclassified\\*",
                             replacement="Unclassified<br>*\\1*"))

#Load lefse file
pl_lefse <- read_tsv("Planktonic/dental.final.an.0.02.subsample.0.02.lefse_summary") 

#Remove NA
pl_lefse <- pl_lefse %>%
  drop_na(LDA) %>%
  filter(LDA > 4) %>%
  inner_join(., pl_taxonomy, by=c("OTU" = "otu")) %>%
  mutate(LDA = if_else(Class == "5A8", -1* LDA, LDA),
         genus=fct_reorder(genus, LDA))



#Plot lefse
###Figure S3###
p3 <- ggplot(pl_lefse, aes(x=LDA, y=genus, fill=Class)) +
  geom_col(color="black", linewidth=1) +
  labs(y=NULL, x="LDA Score (log 10)") +
  scale_x_continuous(limits=c(-10, 10), breaks = seq(-10, 10, by =2)) +
  scale_fill_manual(name=NULL, breaks=c("5A8", "DMSO"),
                    labels=c("Aerobic", "Anaerobic"),
                    values=c("#666666", "#7030A0")) +
  theme_classic() +
  theme(axis.text.y = element_markdown(),
        axis.text= element_text(size = 16, face="bold"),
        axis.title=element_text(size=18, face="bold"),
        legend.text=element_markdown(size = 16))

ggsave("planktonic_lefse.jpg", p3, width=12, height=9, dpi=320)

###lefse
##Biofilm
#Load taxonomy file
bf_taxonomy <- read_tsv("Biofilm/dental.final.an.0.02.cons.taxonomy") %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy,
           into=c("kingdom", "phylum", "class", "order", "family", "genus"),
           sep=";") %>%
  mutate(pretty_otu = str_replace(string=otu,
                                  pattern="tu0*",
                                  replacement = "TU "),
         genus = str_replace(string=genus,
                             pattern="(.*)",
                             replacement="*\\1*"),
         genus = str_replace(string=genus,
                             pattern="\\*(.*)_unclassified\\*",
                             replacement="Unclassified<br>*\\1*"))
#Load lefse file
bf_lefse <- read_tsv("Biofilm/dental.final.an.0.02.subsample.0.02.lefse_summary") 

#Remove NA
bf_lefse <- bf_lefse %>%
  drop_na(LDA) %>%
  filter(LDA > 4) %>%
  inner_join(., bf_taxonomy, by=c("OTU" = "otu")) %>%
  mutate(LDA = if_else(Class == "5A8", -1* LDA, LDA),
         genus=fct_reorder(genus, LDA))

#Plot lefse
####Figure 2B###
p4 <- ggplot(bf_lefse, aes(x=LDA, y=genus, fill=Class)) +
  geom_col(color="black", linewidth=1) +
  labs(y=NULL, x="LDA Score (log 10)") +
  scale_x_continuous(limits=c(-10, 5), breaks = seq(-10, 5, by =2)) +
  scale_fill_manual(name=NULL, breaks=c("5A8", "DMSO"),
                    labels=c("Aerobic", "Anaerobic"),
                    values=c("#666666", "#7030A0")) +
  theme_classic() +
  theme(axis.text.y = element_markdown(),
        axis.text= element_text(size = 16, face="bold"),
        axis.title=element_text(size=18, face="bold"),
        legend.text=element_markdown(size = 16))

ggsave("biofilm_lefse.jpg", p4, width=12, height=9, dpi=320)

