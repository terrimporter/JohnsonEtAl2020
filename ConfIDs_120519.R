# Teresita M. Porter, Dec. 5, 2019

library(vegan)
library(reshape2)
library(ggplot2)
library(ggpubr) # for common legend

###################################################################

# Read infiles
X16S<-read.csv(file="16S_taxonomic_assignments.csv", head=TRUE)
# For COI save each amplicon separately into their own df
COI<-read.csv(file="GJ_COI_taxonomic_assignments_raw.csv", head=TRUE)
F230R <- COI[grepl("FolF",COI$Amplicon_GlobalESV),]
mljg <- COI[grepl("mlCOI",COI$Amplicon_GlobalESV),]
# work with the unpaired reads instead to account for diversity in the longer v4 regions (400-600bp)
X18S<-read.csv(file="taxonomic_assignments_18S_unpaired.csv", head=TRUE)
# Separate out the R1 and R2 reads
Uni18S <- X18S[grepl("Uni18SF", X18S$GlobalASV),]
Uni18SR <- X18S[grepl("Uni18SR", X18S$GlobalASV),]

# Select phylum Arthropoda for F230R
F230R.Arth<-F230R[F230R$Phylum=="Arthropoda",]
# Select Metazoa for Leray amplicon
mljg.Meta<-mljg[mljg$Kingdom=="Metazoa",]

# 16S all Bacteria
# 18S has Archaea, Bacteria, Eukaryota (primers designed for Crustacea, Mollusca, Tunicata)
# Select 18S for Eukaryota only
Uni18S.Euk <- Uni18S[Uni18S$SuperKingdom=="Eukaryota",]
Uni18SR.Euk <- Uni18SR[Uni18SR$SuperKingdom=="Eukaryota",]

# Calc total taxa and total taxa confidently id'd
# based on RDP classifier for seqs > 250 bp https://rdp.cme.msu.edu/classifier/classifier.jsp
X16S.family<-length(unique(X16S$Family))
X16S.family_good<-length(unique(X16S$Family[X16S$fBP>=0.80]))

# Calc total taxa and total taxa confidently id'd (99% correct, 200bp)
# based on CO1v3 cutoffs for seqs ~ 200 bp https://github.com/terrimporter/CO1Classifier
F230R.family<-length(unique(F230R.Arth$Family))
F230R.family_good<-length(unique(F230R.Arth$Family[F230R.Arth$fBP>=0.20]))
mljg.family<-length(unique(mljg.Meta$Family))
mljg.family_good<-length(unique(mljg.Meta$Family[mljg.Meta$fBP>=0.20]))

# Calc total taxa and total taxa confidently id'd (95% correct, 400 bp)
# based on 18Sv3 cutoffs for seqs ~ 200 bp https://github.com/terrimporter/18SClassifier
Uni18S.family<-length(unique(Uni18S.Euk$Family))
Uni18S.family_good<-length(unique(Uni18S.Euk$Family[Uni18S.Euk$fBP>=0.80]))

Uni18SR.family<-length(unique(Uni18SR.Euk$Family))
Uni18SR.family_good<-length(unique(Uni18SR.Euk$Family[Uni18SR.Euk$fBP>=0.80]))

# create df for ggplot
X16S.df<-data.frame("rank"=c(rep("family",2)),
                   "status"=c("all","good"),
                   "value"=c(X16S.family, X16S.family_good))
F230R.df<-data.frame("rank"=c(rep("family",2)),
                    "status"=c("all","good"),
                    "value"=c(F230R.family, F230R.family_good))
mljg.df<-data.frame("rank"=c(rep("family",2)),
                     "status"=c("all","good"),
                     "value"=c(mljg.family, mljg.family_good))
Uni18S.df<-data.frame("rank"=c(rep("family",2)),
                    "status"=c("all","good"),
                    "value"=c(Uni18S.family, Uni18S.family_good))
Uni18SR.df<-data.frame("rank"=c(rep("family",2)),
                    "status"=c("all","good"),
                    "value"=c(Uni18SR.family, Uni18SR.family_good))

# create factors
X16S.df$status <- factor(X16S.df$status, levels = c("all","good"),
                         labels=c("All taxa", "Confidently identified taxa"))
F230R.df$status <- factor(F230R.df$status, levels = c("all","good"),
                         labels=c("All taxa", "Confidently identified taxa"))
mljg.df$status <- factor(mljg.df$status, levels = c("all","good"),
                          labels=c("All taxa", "Confidently identified taxa"))
Uni18S.df$status <- factor(Uni18S.df$status, levels = c("all","good"),
                         labels=c("All taxa", "Confidently identified taxa"))
Uni18SR.df$status <- factor(Uni18SR.df$status, levels = c("all","good"),
                           labels=c("All taxa", "Confidently identified taxa"))

# create bar plot with two series
p1<-ggplot(X16S.df, aes(fill=status, y=value, x=rank)) +
  ggtitle("16Sv3v4 Bacteria") +
  geom_bar(position="dodge",stat="identity") +
  scale_x_discrete(limits = rev(levels(rank))) +
  labs(x="Rank", y="Unique families") +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title=element_blank(),
        legend.text=element_text(size=14),
        legend.position = "bottom",
        legend.spacing.x = unit(0.5, 'cm'),
        axis.text=element_text(size=14),
        axis.text.x = element_blank(),
        axis.title=element_text(size=14),
        axis.title.x=element_blank())

p2 <- ggplot(F230R.df, aes(fill=status, y=value, x=rank)) +
  ggtitle("COI (F230R) Arthropoda") +
  geom_bar(position="dodge",stat="identity") +
  scale_x_discrete(limits = rev(levels(rank))) +
  labs(x="Rank", y="Unique families") +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title=element_blank(),
        legend.text=element_text(size=14),
        legend.position = "bottom",
        legend.spacing.x = unit(0.5, 'cm'),
        axis.text=element_text(size=14),
        axis.text.x = element_blank(),
        axis.title=element_text(size=14),
        axis.title.x=element_blank())

p3 <- ggplot(mljg.df, aes(fill=status, y=value, x=rank)) +
  ggtitle("COI (ml-jg) Metazoa") +
  geom_bar(position="dodge",stat="identity") +
  scale_x_discrete(limits = rev(levels(rank))) +
  labs(x="Rank", y="Unique families") +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title=element_blank(),
        legend.text=element_text(size=14),
        legend.position = "bottom",
        legend.spacing.x = unit(0.5, 'cm'),
        axis.text=element_text(size=14),
        axis.text.x = element_blank(),
        axis.title=element_text(size=14),
        axis.title.x=element_blank())

p4<-ggplot(Uni18S.df, aes(fill=status, y=value, x=rank)) +
  ggtitle("18Sv4 (Uni18S) Eukaryota") +
  geom_bar(position="dodge",stat="identity") +
  scale_x_discrete(limits = rev(levels(rank))) +
  labs(x="Rank", y="Unique families") +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title=element_blank(),
        legend.text=element_text(size=14),
        legend.position = "bottom",
        legend.spacing.x = unit(0.5, 'cm'),
        axis.text=element_text(size=14),
        axis.text.x = element_blank(),
        axis.title=element_text(size=14),
        axis.title.x=element_blank())

p5<-ggplot(Uni18SR.df, aes(fill=status, y=value, x=rank)) +
  ggtitle("18Sv4 (Uni18SR) Eukaryota") +
  geom_bar(position="dodge",stat="identity") +
  scale_x_discrete(limits = rev(levels(rank))) +
  labs(x="Rank", y="Unique families") +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title=element_blank(),
        legend.text=element_text(size=14),
        legend.position = "bottom",
        legend.spacing.x = unit(0.5, 'cm'),
        axis.text=element_text(size=14),
        axis.text.x = element_blank(),
        axis.title=element_text(size=14),
        axis.title.x=element_blank())


g <- ggarrange(p1, NULL, p2, p3, p4, p5, nrow=3, ncol=2, common.legend = TRUE, legend="bottom")

ggsave("Confidentids_120519.pdf", g, width = 8, height = 8)

