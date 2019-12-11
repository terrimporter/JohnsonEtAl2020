# Teresita M. Porter, Dec.5, 2019

library(stringr)
library(reshape2)
library(vegan)
library(purrr) # for map_dfr
library(ggplot2)
library(scales)
library(data.table)
library(gridExtra)

###################################################
# 16Sv3v4
###################################################

# Presence absence plots
X16S <- read.table("16S_taxonomic_assignments.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

# parse molecule and Filter ID from sample name
X16S.1 <- data.frame(X16S, do.call(rbind, str_split(X16S$SampleName, "_")), stringsAsFactors=FALSE)
names(X16S.1)[33:39] <- c("Name","Substrate","Molecule","FilterID","PCRStep","IlluminaSample","Blank")

# Add site column based on filterID
X16S.1$Site <- ifelse(X16S.1$FilterID %in%  c("1","9","2","10","3","11"), "A", "")
X16S.1$Site <- ifelse(X16S.1$FilterID %in%  c("4","12","5","13","6","14"), "B", X16S.1$Site)

#################################################################
# Create box plot figure
#################################################################

# get matrix for rarefaction
X16S.esv<-reshape2::dcast(X16S.1, Molecule+Site+FilterID+PCRStep ~ X16S.1$GlobalASV, value.var = "ASVsize", fun.aggregate = sum)

# merge first 2 colmns into one, then remove them
X16S.esv$sample<-paste(X16S.esv$Molecule, X16S.esv$Site, X16S.esv$FilterID, X16S.esv$PCRStep, sep="_")
X16S.esv <- X16S.esv[,-c(1:4)]

# move sample to rownames then delete
rownames(X16S.esv) <- X16S.esv$sample
X16S.esv$sample <- NULL

# remove columns with only zeros
X16S.esv_notnull<-X16S.esv[,colSums(X16S.esv) !=0]

#remove rows with only zeros & edit rownames
X16S.esv_notnull2<-X16S.esv_notnull[rowSums(X16S.esv_notnull) !=0,]

# remove PCR controls, field negative, and DNA extraction negative before plotting
X16S.esv_notnull2 <- X16S.esv_notnull2[-c(1:4,17:20,33:37,50:51),]

# Remove DNAs treated RNA samples before calculating 15th percentile
X16S.esv_notnull2.2 <- X16S.esv_notnull2[-c(25:48),]

#calculate 15th percentile for rrarefy function
set.seed(1234)
X16S.esv_15percentile<-quantile(rowSums(X16S.esv_notnull2.2), prob=0.15)
# 15% 
# 3809.7 

# Rarefy the dataset down to the 15th percentile before plotting
set.seed(1234)
X16S.rarefied <- rrarefy(X16S.esv_notnull2, sample=X16S.esv_15percentile)

# Convert to presence-absence matrix
X16S.rarefied[X16S.rarefied > 0] <-1

# Calc lineage richness
X16S.richness <- rowSums(X16S.rarefied)

# Move rownames to first column
X16S.richness.1 <- setDT(as.data.frame(X16S.richness), keep.rownames = TRUE)[]

# Split rn into Molecule and PCRstep fields
X16S.rn <- str_split(X16S.richness.1$rn, "_")
df.rn <- data.frame(matrix(unlist(X16S.rn), nrow=length(X16S.rn), byrow=T))
names(df.rn)[1:4] <- c("Molecule", "Site", "FilterID", "PCRStep")

# remove rn
X16S.richness.1$rn <- NULL

# combine Molcule and PCRstep fields with the rest of the df
X16S.richness.2 <- cbind(df.rn, X16S.richness.1)
names(X16S.richness.2)[5] <- "Richness"

# create factors
X16S.richness.2$Molecule <- factor(X16S.richness.2$Molecule, levels = c("DNA", "cDNA", "RNAt", "RNAf"))
#X16S.richness.2$Site <- factor(X16S.richness.2$Site, levels = c("A", "B"))
X16S.richness.2$PCRStep <- factor(X16S.richness.2$PCRStep, levels = c("1STEP", "2STEP"), labels=c("1 step", "2 step"))

# create boxplot
bp_16S <- ggplot(X16S.richness.2, aes(x = Molecule, y = Richness)) +
  geom_boxplot() +
  ggtitle("16S V3-V4") +
  xlab("") +
  ylab("ESV Richness") +
  facet_wrap(~PCRStep) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

###################################################
# COI (F230R and ml-jg) present separately
###################################################

# Presence absence plots
COI <- read.table("GJ_COI_taxonomic_assignments_raw.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

# parse molecule and Filter ID from sample name
COI.1 <- data.frame(COI, do.call(rbind, str_split(COI$SampleName, "_")), stringsAsFactors=FALSE)
names(COI.1)[34:36] <- c("Molecule","FilterID","PCRStep")

# parse amplicon from Amplicon_GlobalESV field
COI.2 <- data.frame(COI.1, do.call(rbind, str_split(COI.1$Amplicon_GlobalESV, "_")), stringsAsFactors=FALSE)
names(COI.2)[37:38] <- c("Amplicon","OTU_dontuse")

# Add site column based on filterID
COI.2$Site <- ifelse(COI.2$FilterID %in%  c("1","9","2","10","3","11"), "A", "")
COI.2$Site <- ifelse(COI.2$FilterID %in%  c("4","12","5","13","6","14"), "B", COI.2$Site)

#################################################################
# Create box plot figure
#################################################################

# get matrix for rarefaction
COI.esv<-reshape2::dcast(COI.2, Molecule+Site+FilterID+PCRStep+Amplicon ~ COI.2$Amplicon_GlobalESV, value.var = "ESVsize", fun.aggregate = sum)

# merge first 2 colmns into one, then remove them
COI.esv$sample<-paste(COI.esv$Molecule, COI.esv$Site, COI.esv$FilterID, COI.esv$PCRStep, COI.esv$Amplicon, sep="_")
COI.esv <- COI.esv[,-c(1:5)]

# move sample to rownames then delete
rownames(COI.esv) <- COI.esv$sample
COI.esv$sample <- NULL

# remove columns with only zeros
COI.esv_notnull<-COI.esv[,colSums(COI.esv) !=0]

#remove rows with only zeros & edit rownames
COI.esv_notnull2<-COI.esv_notnull[rowSums(COI.esv_notnull) !=0,]

# remove PCR controls, lab negative, and DNA extraction negative before plotting
COI.esv_notnull2 <- COI.esv_notnull2[-c(1:4,29:33,58:67,83:84),]

# Remove DNAse treated RNA samples before calculating 15th percentile
COI.esv_notnull2.2 <- COI.esv_notnull2[-c(49:72),]

#calculate 15th percentile for rrarefy function
COI.esv_15percentile<-quantile(rowSums(COI.esv_notnull2.2), prob=0.15)
# 15% 
# 21233.9 

# Curve looks good so rarefy to 15th percentile
# Rarefy the dataset down to the 15th percentile
set.seed(1234)
COI.rarefied <- rrarefy(COI.esv_notnull2, sample=COI.esv_15percentile)

# Convert to presence-absence matrix
COI.rarefied[COI.rarefied > 0] <-1

# Calc lineage richness
COI.richness <- rowSums(COI.rarefied)

# Move rownames to first column
COI.richness.1 <- setDT(as.data.frame(COI.richness), keep.rownames = TRUE)[]

# Split rn into Molecule and PCRstep fields
COI.rn <- str_split(COI.richness.1$rn, "_")
df.rn <- data.frame(matrix(unlist(COI.rn), nrow=length(COI.rn), byrow=T))
names(df.rn)[1:5] <- c("Molecule", "Site", "FilterID", "PCRStep", "Amplicon")

# remove rn
COI.richness.1$rn <- NULL

# combine Molcule and PCRstep fields with the rest of the df
COI.richness.2 <- cbind(df.rn, COI.richness.1)
names(COI.richness.2)[6] <- "Richness"

# create factors
COI.richness.2$Molecule <- factor(COI.richness.2$Molecule, levels = c("DNA", "cDNA", "RNAt", "RNAf"))
#COI.richness.2$Site <- factor(COI.richness.2$Site, levels = c("A", "B"))
COI.richness.2$Amplicon <- factor(COI.richness.2$Amplicon, 
                                  levels = c("FolF", "mlCOI"),
                                  labels = c("F230R", "ml-jg"))
COI.richness.2$PCRStep <- factor(COI.richness.2$PCRStep, levels = c("1STEP", "2STEP"), labels=c("1 step", "2 step"))

# create boxplot
bp_COI <- ggplot(COI.richness.2, aes(x = Molecule, y = Richness)) +
  geom_boxplot(aes(fill=Amplicon)) +
  ggtitle("COI") +
  xlab("") +
  ylab("ESV Richness") +
  facet_wrap(~PCRStep) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "bottom")

###################################################
# 18Sv4 (present F230R + ml-jg separately)
###################################################

# Presence absence plots
X18S <- read.table("taxonomic_assignments_18S_unpaired.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

# parse molecule and Filter ID from sample name
X18S.1 <- data.frame(X18S, do.call(rbind, str_split(X18S$SampleName, "_")), stringsAsFactors=FALSE)
names(X18S.1)[33:41] <- c("Name","Substrate","Molecule","FilterID","PCRStep","IlluminaSample","Lane","Read","Lane2")

# parse amplicon from Amplicon_GlobalESV field
X18S.2 <- data.frame(X18S.1, do.call(rbind, str_split(X18S.1$GlobalASV, "_")), stringsAsFactors=FALSE)
names(X18S.2)[42:43] <- c("Primer","OTU_dontuse")

# Add site column based on filterID
X18S.2$Site <- ifelse(X18S.2$FilterID %in%  c("1","9","2","10","3","11"), "A", "")
X18S.2$Site <- ifelse(X18S.2$FilterID %in%  c("4","12","5","13","6","14"), "B", X18S.2$Site)

# Create box plot figure
#################################################################

# get matrix for rarefaction
X18S.esv<-reshape2::dcast(X18S.2, Molecule+Site+FilterID+PCRStep+Primer ~ X18S.2$GlobalASV, value.var = "ASVsize", fun.aggregate = sum)

# merge first 2 colmns into one, then remove them
X18S.esv$sample<-paste(X18S.esv$Molecule, X18S.esv$Site, X18S.esv$FilterID, X18S.esv$PCRStep, X18S.esv$Primer, sep="_")
X18S.esv <- X18S.esv[,-c(1:5)]

# move sample to rownames then delete
rownames(X18S.esv) <- X18S.esv$sample
X18S.esv$sample <- NULL

# remove columns with only zeros
X18S.esv_notnull<-X18S.esv[,colSums(X18S.esv) !=0]

#remove rows with only zeros & edit rownames
X18S.esv_notnull2<-X18S.esv_notnull[rowSums(X18S.esv_notnull) !=0,]

# remove PCR controls, lab negative and DNA extraction negative before plotting
X18S.esv_notnull2 <- X18S.esv_notnull2[-c(1:8,33:39,64:72,93:96),]

# Remove DNAse treated RNA samples before calculating 15th percentile
X18S.esv_notnull2.2 <- X18S.esv_notnull2[-c(49:87),]

#calculate 15th percentile for rrarefy function
set.seed(1234)
X18S.esv_15percentile<-quantile(rowSums(X18S.esv_notnull2.2), prob=0.15)
# 15% 
# 18982.35 

# Curve looks good so rarefy to 15th percentile
# Rarefy the dataset down to the 15th percentile
set.seed(1234)
X18S.rarefied <- rrarefy(X18S.esv_notnull2, sample=X18S.esv_15percentile)

# Convert to presence-absence matrix
X18S.rarefied[X18S.rarefied > 0] <-1

# Calc lineage richness
X18S.richness <- rowSums(X18S.rarefied)

# Move rownames to first column
X18S.richness.1 <- setDT(as.data.frame(X18S.richness), keep.rownames = TRUE)[]

# Split rn into Molecule and PCRstep fields
X18S.rn <- str_split(X18S.richness.1$rn, "_")
df.rn <- data.frame(matrix(unlist(X18S.rn), nrow=length(X18S.rn), byrow=T))
names(df.rn)[1:5] <- c("Molecule", "Site", "FilterID", "PCRStep", "Primer")

# remove rn
X18S.richness.1$rn <- NULL

# combine Molcule and PCRstep fields with the rest of the df
X18S.richness.2 <- cbind(df.rn, X18S.richness.1)
names(X18S.richness.2)[6] <- "Richness"

# create factors
X18S.richness.2$Molecule <- factor(X18S.richness.2$Molecule, levels = c("DNA", "cDNA", "RNAt", "RNAf"))
#X18S.richness.2$Site <- factor(X18S.richness.2$Site, levels = c("A", "B"))
X18S.richness.2$Primer <- factor(X18S.richness.2$Primer, 
                                  levels = c("Uni18SF", "Uni18SR"),
                                  labels = c("Uni18S", "Uni18SR"))
X18S.richness.2$PCRStep <- factor(X18S.richness.2$PCRStep, levels = c("1STEP", "2STEP"), labels=c("1 step", "2 step"))

# create boxplot
bp_18S <- ggplot(X18S.richness.2, aes(x = Molecule, y = Richness)) +
  geom_boxplot(aes(fill=Primer)) +
  ggtitle("18S V4") +
  xlab("") +
  ylab("ESV Richness") +
  facet_wrap(~PCRStep) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "bottom")

g <- grid.arrange(bp_16S, bp_COI, bp_18S)

ggsave("boxplot_120519.pdf", g, width = 8, height = 8)
