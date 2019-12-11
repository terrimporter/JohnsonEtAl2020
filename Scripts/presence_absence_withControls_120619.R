# Teresita M. Porter, Dec. 6, 2019

library(stringr)
library(reshape2)
library(vegan)
library(purrr) # for map_dfr
library(ggplot2)
library(scales)
library(data.table)
library(gridExtra)
library(ggpubr) # for common legend

##################################################
# 16Sv3v4
##################################################

# Presence absence plots
X16S <- read.table("16S_taxonomic_assignments.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

# parse molecule and Filter ID from sample name
X16S.1 <- data.frame(X16S, do.call(rbind, str_split(X16S$SampleName, "_")), stringsAsFactors=FALSE)
names(X16S.1)[33:39] <- c("Name","Substrate","Molecule","FilterID","PCRStep","IlluminaSample","Blank")

# Add site column based on filterID
X16S.1$Site <- ifelse(X16S.1$FilterID %in%  c("1","9","2","10","3","11"), "A", "")
X16S.1$Site <- ifelse(X16S.1$FilterID %in%  c("4","12","5","13","6","14"), "B", X16S.1$Site)
X16S.1$Site <- ifelse(X16S.1$FilterID %in%  c("7"), "Field nc", X16S.1$Site)
X16S.1$Site <- ifelse(X16S.1$FilterID %in%  c("8"), "Lab nc", X16S.1$Site)
X16S.1$Site <- gsub("^$", "Extraction nc", X16S.1$Site)

# filter for good family rank assignments
# Based on RDP classifier website default cutoff for seqs > 250 bp
# https://rdp.cme.msu.edu/classifier/classifier.jsp
X16S.f <- X16S.1[X16S.1$fBP >= 0.80,]

# get matrix for rarefaction
X16S.family<-reshape2::dcast(X16S.f, Molecule+Site+PCRStep ~ X16S.f$Family, value.var = "ASVsize", fun.aggregate = sum)

# merge first 2 colmns into one, then remove them
X16S.family$sample<-paste(X16S.family$Molecule, X16S.family$Site, X16S.family$PCRStep, sep="_")
X16S.family <- X16S.family[,-c(1:3)]

# move sample to rownames then delete
rownames(X16S.family) <- X16S.family$sample
X16S.family$sample <- NULL

# remove columns with only zeros
X16S.family_notnull<-X16S.family[,colSums(X16S.family) !=0]

#remove rows with only zeros & edit rownames
X16S.family_notnull2<-X16S.family_notnull[rowSums(X16S.family_notnull) !=0,]

# KEEP PCR controls, lab, and DNA Extraction ncs for plotting 
# but remove before calculating 15th percentile
X16S.family_notnull2.2 <- X16S.family_notnull2[-c(1:2,7:8,13:15,20),]

# Keep Remove DNAse treated RNA samples for plotting
# but before calculating 15th percentile
X16S.family_notnull2.2 <- X16S.family_notnull2.2[-c(9:16),]

#calculate 15th percentile for rrarefy function
set.seed(1234)
X16S.family_15percentile<-quantile(rowSums(X16S.family_notnull2.2), prob=0.15)
# 15% 
# 14613.3 

# Create presence absence figure
# Rarefy the dataset down to the 15th percentile before plotting
set.seed(1234)
X16S.rarefied <- rrarefy(X16S.family_notnull2, sample=X16S.family_15percentile)

# Convert to presence-absence matrix
X16S.rarefied[X16S.rarefied > 0] <-1

# Move rownames to first column
X16S.rarefied.1 <- setDT(as.data.frame(X16S.rarefied), keep.rownames = TRUE)[]

# Split rn into Molecule and PCRstep fields
X16S.rn <- str_split(X16S.rarefied.1$rn, "_")
df.rn <- data.frame(matrix(unlist(X16S.rn), nrow=length(X16S.rn), byrow=T))
names(df.rn)[1:3] <- c("Molecule", "Site", "PCRStep")

# remove rn
X16S.rarefied.1$rn <- NULL

# combine Molcule and PCRstep fields with the rest of the df
X16S.rarefied.2 <- cbind(df.rn, X16S.rarefied.1)

# make long format for ggplot
X16S.long <- reshape2::melt(X16S.rarefied.2, id.vars=c("Molecule", "Site", "PCRStep"))

# add an empty entry so that the nubmer of panels for 16S is same as for COI and 18S
X16S.empty <- data.frame(Molecule="PCRnc", PCRStep="2STEP", Site="A", value=0, variable="Acetobacteraceae")
X16S.long <- rbind(X16S.long, X16S.empty)

##################################################
# COI - ml-jg and F230R (use these combined since just presenting presence absence data)
##################################################

# Presence absence plots
COI <- read.table("GJ_COI_taxonomic_assignments_raw.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

# parse molecule and Filter ID from sample name
COI.1 <- data.frame(COI, do.call(rbind, str_split(COI$SampleName, "_")), stringsAsFactors=FALSE)
names(COI.1)[34:36] <- c("Molecule","FilterID","PCRStep")

# Add site column based on filterID
COI.1$Site <- ifelse(COI.1$FilterID %in%  c("1","9","2","10","3","11"), "A", "")
COI.1$Site <- ifelse(COI.1$FilterID %in%  c("4","12","5","13","6","14"), "B", COI.1$Site)
COI.1$Site <- ifelse(COI.1$FilterID %in%  c("7"), "Field nc", COI.1$Site)
COI.1$Site <- ifelse(COI.1$FilterID %in%  c("8"), "Lab nc", COI.1$Site)
COI.1$Site <- gsub("^$", "Extraction nc", COI.1$Site)

# filter for good family rank assignments
# Based on CO1 Classifier v3 cutoff (99% correct) for seqs ~ 200 bp
# https://github.com/terrimporter/CO1Classifier
COI.f <- COI.1[COI.1$fBP >= 0.20,]

# get matrix for rarefaction
COI.family<-reshape2::dcast(COI.f, Molecule+Site+PCRStep ~ COI.f$Family, value.var = "ESVsize", fun.aggregate = sum)

# merge first 2 colmns into one, then remove them
COI.family$sample<-paste(COI.family$Molecule, COI.family$Site, COI.family$PCRStep, sep="_")
COI.family <- COI.family[,-c(1:3)]

# move sample to rownames then delete
rownames(COI.family) <- COI.family$sample
COI.family$sample <- NULL

# remove columns with only zeros
COI.family_notnull<-COI.family[,colSums(COI.family) !=0]

#remove rows with only zeros & edit rownames
COI.family_notnull2<-COI.family_notnull[rowSums(COI.family_notnull) !=0,]

# KEEP PCR controls, lab, and DNA Extraction ncs for plotting 
# but remove before calculating 15th percentile
COI.family_notnull2.2<-COI.family_notnull2[-c(5:7,12:18,22,26:27),]

# Keep Remove DNAse treated RNA samples for plotting
# but before calculating 15th percentile
COI.family_notnull2.2 <- COI.family_notnull2.2[-c(9:14),]

#calculate 15th percentile for rrarefy function
COI.family_15percentile<-quantile(rowSums(COI.family_notnull2.2), prob=0.15)
# 15% 
# 90279.35 

#################################################################
# Create presence absence figure
#################################################################

# Rarefy the dataset down to the 15th percentile before plotting
set.seed(1234)
COI.rarefied <- rrarefy(COI.family_notnull2, sample=COI.family_15percentile)

# Convert to presence-absence matrix
COI.rarefied[COI.rarefied > 0] <-1

# Move rownames to first column
COI.rarefied.1 <- setDT(as.data.frame(COI.rarefied), keep.rownames = TRUE)[]

# Split rn into Molecule and PCRstep fields
COI.rn <- str_split(COI.rarefied.1$rn, "_")
df.rn <- data.frame(matrix(unlist(COI.rn), nrow=length(COI.rn), byrow=T))
names(df.rn)[1:3] <- c("Molecule", "Site", "PCRStep")

# remove rn
COI.rarefied.1$rn <- NULL

# combine Molcule and PCRstep fields with the rest of the df
COI.rarefied.2 <- cbind(df.rn, COI.rarefied.1)

# make long format for ggplot
COI.long <- reshape2::melt(COI.rarefied.2, id.vars=c("Molecule", "Site", "PCRStep"))

##################################################
# 18Sv4
##################################################

# Presence absence plots
X18S <- read.table("taxonomic_assignments_18S_unpaired.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

# parse molecule and Filter ID from sample name
X18S.1 <- data.frame(X18S, do.call(rbind, str_split(X18S$SampleName, "_")), stringsAsFactors=FALSE)
names(X18S.1)[33:41] <- c("Name","Substrate","Molecule","FilterID","PCRStep","IlluminaSample","Lane","Read","Lane2")

# Add site column based on filterID
X18S.1$Site <- ifelse(X18S.1$FilterID %in%  c("1","9","2","10","3","11"), "A", "")
X18S.1$Site <- ifelse(X18S.1$FilterID %in%  c("4","12","5","13","6","14"), "B", X18S.1$Site)
X18S.1$Site <- ifelse(X18S.1$FilterID %in%  c("7"), "Field nc", X18S.1$Site)
X18S.1$Site <- ifelse(X18S.1$FilterID %in%  c("8"), "Lab nc", X18S.1$Site)
X18S.1$Site <- gsub("^$", "Extraction nc", X18S.1$Site)

# filter for good family rank assignments
# based on 18S classifier v3 (95% correct) for 200 bp seqs https://github.com/terrimporter/18SClassifier
X18S.f <- X18S.1[X18S.1$fBP >= 0.80,]

# get matrix for rarefaction
X18S.family<-reshape2::dcast(X18S.f, Molecule+Site+PCRStep ~ X18S.f$Family, value.var = "ASVsize", fun.aggregate = sum)

# merge first 2 colmns into one, then remove them
X18S.family$sample<-paste(X18S.family$Molecule, X18S.family$Site, X18S.family$PCRStep, sep="_")
X18S.family <- X18S.family[,-c(1:3)]

# move sample to rownames then delete
rownames(X18S.family) <- X18S.family$sample
X18S.family$sample <- NULL

# remove columns with only zeros
X18S.family_notnull<-X18S.family[,colSums(X18S.family) !=0]

#remove rows with only zeros & edit rownames
X18S.family_notnull2<-X18S.family_notnull[rowSums(X18S.family_notnull) !=0,]

# KEEP PCR controls, lab, and DNA Extraction ncs for plotting 
# but remove before calculating 15th percentile
X18S.family_notnull2.2 <- X18S.family_notnull2[-c(5:7,12:18,23),]

# Keep Remove DNAse treated RNA samples for plotting
# but before calculating 15th percentile
X18S.family_notnull2.2 <- X18S.family_notnull2.2[-c(9:16),]

#calculate 15th percentile for rrarefy function
X18S.family_15percentile<-quantile(rowSums(X18S.family_notnull2.2), prob=0.15)
# 15% 
# 51836.3   

#################################################################
# Create presence absence figure
#################################################################

# Rarefy the dataset down to the 15th percentile before plotting
set.seed(1234)
X18S.rarefied <- rrarefy(X18S.family_notnull2, sample=X18S.family_15percentile)

# Convert to presence-absence matrix
X18S.rarefied[X18S.rarefied > 0] <-1

# Move rownames to first column
X18S.rarefied.1 <- setDT(as.data.frame(X18S.rarefied), keep.rownames = TRUE)[]

# Split rn into Molecule and PCRstep fields
X18S.rn <- str_split(X18S.rarefied.1$rn, "_")
df.rn <- data.frame(matrix(unlist(X18S.rn), nrow=length(X18S.rn), byrow=T))
names(df.rn)[1:3] <- c("Molecule", "Site", "PCRStep")

# remove rn
X18S.rarefied.1$rn <- NULL

# combine Molcule and PCRstep fields with the rest of the df
X18S.rarefied.2 <- cbind(df.rn, X18S.rarefied.1)

# make long format for ggplot
X18S.long <- reshape2::melt(X18S.rarefied.2, id.vars=c("Molecule", "Site", "PCRStep"))

# put into a single df for a cleaner plot
X16S.long$amplicon <- "16S" 
COI.long$amplicon <- "COI" 
X18S.long$amplicon <- "18S"
merged <- rbind(X16S.long, COI.long, X18S.long)

# create factors
merged$amplicon <- factor(merged$amplicon, levels=c("16S", "COI", "18S"),
                          labels=c("16S V3-V4", "COI", "18S V4"))
merged$Site <- factor(merged$Site, levels=c("A", "B", "Field nc", "Lab nc", "Extraction nc"),
                      labels=c("Site A", "Site B", "Field nc", "Lab nc", "Extraction nc"))
merged$variable <- factor(merged$variable, levels=rev(sort(unique(merged$variable))))
merged$value <- factor(merged$value, levels=unique(merged$value),
                       labels=c("Detected", "Not detected"))
merged$Molecule <- factor(merged$Molecule, levels = c("DNA", "cDNA", "RNAt", "RNAf","PCRnc", "PCRpc"))
merged$PCRStep <- factor(merged$PCRStep, levels = c("1STEP", "2STEP"), labels=c("1 step", "2 step"))

# create heatmap
p1 <- ggplot(merged, aes(x = Site, y = variable)) +
  geom_tile(aes(fill = value, color=value)) +
  xlab("Sites") +
  ylab("Families") +
  facet_grid(vars(amplicon), vars(Molecule, PCRStep), scales="free_y") +
  scale_fill_manual(values = c("black", "white")) +
  scale_color_manual(values = c("black", "white")) +
  theme_bw() +
  theme(panel.border = element_rect(fill=NA, color=alpha("black", 1), size=1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_text(angle = 0)) +
  guides(fill=FALSE, color=FALSE)

ggsave("presenceAbsence_withControls_120619.pdf", p1, width = 8, height = 11)
