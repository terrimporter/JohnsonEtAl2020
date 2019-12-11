# Teresita M. Porter, Dec. 6, 2019

library(stringr)
library(vegan)
library(ggplot2)
library(ggpubr) # for common legend

# to calculate venn counts
#source("http://www.bioconductor.org/biocLite.R")
#biocLite("limma")
library(limma)
library(ggforce) # to draw circles with ggplot

##################################################
# 16Sv3v4
##################################################

# Presence absence plots
X16S <- read.table("16S_taxonomic_assignments.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

# parse molecule and Filter ID from sample name
X16S.1 <- data.frame(X16S, do.call(rbind, str_split(X16S$SampleName, "_")), stringsAsFactors=FALSE)
names(X16S.1)[33:39] <- c("Name","Substrate","Molecule","FilterID","PCRStep","IlluminaSample","Blank")

# exclude RNAf (came from separate extraction)
X16S.1 <-X16S.1[X16S.1$Molecule != "RNAf",]

# Add site column based on filterID
X16S.1$Site <- ifelse(X16S.1$FilterID %in%  c("1","9","2","10","3","11"), "A", "")
X16S.1$Site <- ifelse(X16S.1$FilterID %in%  c("4","12","5","13","6","14"), "B", X16S.1$Site)

# Add bottle column based on filterID (only interested in bottles 1-6)
X16S.1$Bottle <- ifelse(X16S.1$FilterID %in%  c("1","9"), "1", NA)
X16S.1$Bottle <- ifelse(X16S.1$FilterID %in%  c("2","10"), "2", X16S.1$Bottle)
X16S.1$Bottle <- ifelse(X16S.1$FilterID %in%  c("3","11"), "3", X16S.1$Bottle)
X16S.1$Bottle <- ifelse(X16S.1$FilterID %in%  c("4","12"), "4", X16S.1$Bottle)
X16S.1$Bottle <- ifelse(X16S.1$FilterID %in%  c("5","13"), "5", X16S.1$Bottle)
X16S.1$Bottle <- ifelse(X16S.1$FilterID %in%  c("6","14"), "6", X16S.1$Bottle)
X16S.1$Bottle <- ifelse(X16S.1$FilterID %in%  c("7"), NA, X16S.1$Bottle)

# Only work with data from bottles 1-6
X16S.2 <- X16S.1[!is.na(X16S.1$Bottle),]

# get matrix for rarefaction
X16S.esv <- reshape2::dcast(X16S.2, Molecule+Bottle ~ X16S.2$GlobalASV, value.var = "ASVsize", fun.aggregate = sum)

# merge first 2 colmns into one, then remove them
X16S.esv$sample<-paste(X16S.esv$Molecule, X16S.esv$Bottle, sep="_")
X16S.esv <- X16S.esv[,-c(1:2)]

# move sample to rownames then delete
rownames(X16S.esv) <- X16S.esv$sample
X16S.esv$sample <- NULL

# remove columns with only zeros
X16S.esv_notnull<-X16S.esv[,colSums(X16S.esv) !=0]

#remove rows with only zeros & edit rownames
X16S.esv_notnull2<-X16S.esv_notnull[rowSums(X16S.esv_notnull) !=0,]

# remove PCR controls, lab, and DNA extraction negatives before plotting
X16S.esv_notnull2 <- X16S.esv_notnull2[-c(13:15),]

# Remove DNAse treated RNA samples before calculating 15th percentile
X16S.esv_notnull2.2 <- X16S.esv_notnull2[-c(13:18),]

#calculate 15th percentile for rrarefy function
set.seed(1234)
X16S.esv_15percentile<-quantile(rowSums(X16S.esv_notnull2.2), prob=0.15)
# 15% 
# 23098.75 

# Rarefy the dataset down to the 15th percentile before plotting
set.seed(1234)
X16S.rarefied <- rrarefy(X16S.esv_notnull2, sample=X16S.esv_15percentile)

# Convert to presence-absence matrix
X16S.rarefied[X16S.rarefied > 0] <-1

# Create venn

# subset by bottle
X16S.1 <- X16S.rarefied[grepl("_1", rownames(X16S.rarefied)),]
X16S.2 <- X16S.rarefied[grepl("_2", rownames(X16S.rarefied)),]
X16S.3 <- X16S.rarefied[grepl("_3", rownames(X16S.rarefied)),]
X16S.4 <- X16S.rarefied[grepl("_4", rownames(X16S.rarefied)),]
X16S.5 <- X16S.rarefied[grepl("_5", rownames(X16S.rarefied)),]
X16S.6 <- X16S.rarefied[grepl("_6", rownames(X16S.rarefied)),]

# calc total richness
X16S.1.tot.richness <- specnumber(colSums(X16S.1))
X16S.2.tot.richness <- specnumber(colSums(X16S.2))
X16S.3.tot.richness <- specnumber(colSums(X16S.3))
X16S.4.tot.richness <- specnumber(colSums(X16S.4))
X16S.5.tot.richness <- specnumber(colSums(X16S.5))
X16S.6.tot.richness <- specnumber(colSums(X16S.6))

# transpose for venn counts
X16S.1.t <- t(X16S.1)
X16S.2.t <- t(X16S.2)
X16S.3.t <- t(X16S.3)
X16S.4.t <- t(X16S.4)
X16S.5.t <- t(X16S.5)
X16S.6.t <- t(X16S.6)

# create counts for Venn
X1.venn <- vennCounts(X16S.1.t)
X2.venn <- vennCounts(X16S.2.t)
X3.venn <- vennCounts(X16S.3.t)
X4.venn <- vennCounts(X16S.4.t)
X5.venn <- vennCounts(X16S.5.t)
X6.venn <- vennCounts(X16S.6.t)

# ignore genera with zero counts
X1.venn[1,4] <- 0
X2.venn[1,4] <- 0
X3.venn[1,4] <- 0
X4.venn[1,4] <- 0
X5.venn[1,4] <- 0
X6.venn[1,4] <- 0

# draw venn diagram with ggplot instead
# df with coordinates for circles
df.venn <- data.frame(x = c(-0.866, 0.866, 0),
                      y = c(1, 1, -0.5),
                      labels = c('DNA', 'cDNA', 'RNAt'))
df.venn$labels <- factor(df.venn$labels, levels=c('DNA', 'cDNA', 'RNAt'))

# transform venncounts into df, add coordinates for labels
class(X1.venn) <- 'matrix'
df.1.venn <- as.data.frame(X1.venn)[-1,] %>%
  mutate(x = c(0, 1.2, 0.8, -1.2, -0.8, 0, 0),
         y = c(-0.8, 1.2, 0, 1.2, 0, 1.4, 0.4))
class(X2.venn) <- 'matrix'
df.2.venn <- as.data.frame(X2.venn)[-1,] %>%
  mutate(x = c(0, 1.2, 0.8, -1.2, -0.8, 0, 0),
         y = c(-0.8, 1.2, 0, 1.2, 0, 1.4, 0.4))
class(X3.venn) <- 'matrix'
df.3.venn <- as.data.frame(X3.venn)[-1,] %>%
  mutate(x = c(0, 1.2, 0.8, -1.2, -0.8, 0, 0),
         y = c(-0.8, 1.2, 0, 1.2, 0, 1.4, 0.4))
class(X4.venn) <- 'matrix'
df.4.venn <- as.data.frame(X4.venn)[-1,] %>%
  mutate(x = c(0, 1.2, 0.8, -1.2, -0.8, 0, 0),
         y = c(-0.8, 1.2, 0, 1.2, 0, 1.4, 0.4))
class(X5.venn) <- 'matrix'
df.5.venn <- as.data.frame(X5.venn)[-1,] %>%
  mutate(x = c(0, 1.2, 0.8, -1.2, -0.8, 0, 0),
         y = c(-0.8, 1.2, 0, 1.2, 0, 1.4, 0.4))
class(X6.venn) <- 'matrix'
df.6.venn <- as.data.frame(X6.venn)[-1,] %>%
  mutate(x = c(0, 1.2, 0.8, -1.2, -0.8, 0, 0),
         y = c(-0.8, 1.2, 0, 1.2, 0, 1.4, 0.4))

# draw venn with ggplot
X16S.1.tot.richness.comma <- format(X16S.1.tot.richness,big.mark=",", trim=TRUE)
X16S.2.tot.richness.comma <- format(X16S.2.tot.richness,big.mark=",", trim=TRUE)
X16S.3.tot.richness.comma <- format(X16S.3.tot.richness,big.mark=",", trim=TRUE)
X16S.4.tot.richness.comma <- format(X16S.4.tot.richness,big.mark=",", trim=TRUE)
X16S.5.tot.richness.comma <- format(X16S.5.tot.richness,big.mark=",", trim=TRUE)
X16S.6.tot.richness.comma <- format(X16S.6.tot.richness,big.mark=",", trim=TRUE)

X1 <- paste("16S V3-V4\nBottle 1 =",
           X16S.1.tot.richness.comma)
X2 <- paste("Bottle 2 =",
           X16S.2.tot.richness.comma)
X3 <- paste("Bottle 3 =",
            X16S.3.tot.richness.comma)
X4 <- paste("Bottle 4 =",
            X16S.4.tot.richness.comma)
X5 <- paste("Bottle 5 =",
            X16S.5.tot.richness.comma)
X6 <- paste("Bottle 6 =",
            X16S.6.tot.richness.comma)

v1 <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  ggtitle(paste0(X1)) +
  geom_circle(aes(fill=labels), size = 0.5) +
  geom_circle(size = 0.5, show.legend = FALSE, color = "black", alpha=0.3) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#33638DFF","#20A387FF", "#95D840FF")) +
  scale_color_manual(values=c("#33638DFF","#20A387FF", "#95D840FF"), guide = FALSE) +
  labs(fill = NULL) +
  annotate("text", x = df.1.venn$x, y = df.1.venn$y, label = df.1.venn$Counts, size = 3)

v2 <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  ggtitle(paste0(X2)) +
  geom_circle(aes(fill=labels), size = 0.5) +
  geom_circle(size = 0.5, show.legend = FALSE, color = "black", alpha=0.3) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#33638DFF","#20A387FF", "#95D840FF")) +
  scale_color_manual(values=c("#33638DFF","#20A387FF", "#95D840FF"), guide = FALSE) +
  labs(fill = NULL) +
  annotate("text", x = df.2.venn$x, y = df.2.venn$y, label = df.2.venn$Counts, size = 3)

v3 <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  ggtitle(paste0(X3)) +
  geom_circle(aes(fill=labels), size = 0.5) +
  geom_circle(size = 0.5, show.legend = FALSE, color = "black", alpha=0.3) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#33638DFF","#20A387FF", "#95D840FF")) +
  scale_color_manual(values=c("#33638DFF","#20A387FF", "#95D840FF"), guide = FALSE) +
  labs(fill = NULL) +
  annotate("text", x = df.3.venn$x, y = df.3.venn$y, label = df.3.venn$Counts, size = 3)

v4 <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  ggtitle(paste0(X4)) +
  geom_circle(aes(fill=labels), size = 0.5) +
  geom_circle(size = 0.5, show.legend = FALSE, color = "black", alpha=0.3) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#33638DFF","#20A387FF", "#95D840FF")) +
  scale_color_manual(values=c("#33638DFF","#20A387FF", "#95D840FF"), guide = FALSE) +
  labs(fill = NULL) +
  annotate("text", x = df.4.venn$x, y = df.4.venn$y, label = df.4.venn$Counts, size = 3)

v5 <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  ggtitle(paste0(X5)) +
  geom_circle(aes(fill=labels), size = 0.5) +
  geom_circle(size = 0.5, show.legend = FALSE, color = "black", alpha=0.3) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#33638DFF","#20A387FF", "#95D840FF")) +
  scale_color_manual(values=c("#33638DFF","#20A387FF", "#95D840FF"), guide = FALSE) +
  labs(fill = NULL) +
  annotate("text", x = df.5.venn$x, y = df.5.venn$y, label = df.5.venn$Counts, size = 3)

v6 <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  ggtitle(paste0(X6)) +
  geom_circle(aes(fill=labels), size = 0.5) +
  geom_circle(size = 0.5, show.legend = FALSE, color = "black", alpha=0.3) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#33638DFF","#20A387FF", "#95D840FF")) +
  scale_color_manual(values=c("#33638DFF","#20A387FF", "#95D840FF"), guide = FALSE) +
  labs(fill = NULL) +
  annotate("text", x = df.6.venn$x, y = df.6.venn$y, label = df.6.venn$Counts, size = 3)

##################################################
# COI - ml-jg and F230R (use these combined since just presenting presence absence data)
##################################################

# Presence absence plots
COI <- read.table("GJ_COI_taxonomic_assignments_raw.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

# parse molecule and Filter ID from sample name
COI.1 <- data.frame(COI, do.call(rbind, str_split(COI$SampleName, "_")), stringsAsFactors=FALSE)
names(COI.1)[34:36] <- c("Molecule","FilterID","PCRStep")

# exclude RNAf (came from separate extraction)
COI.1 <-COI.1[COI.1$Molecule != "RNAf",]

# Add site column based on filterID
COI.1$Site <- ifelse(COI.1$FilterID %in%  c("1","9","2","10","3","11"), "A", "")
COI.1$Site <- ifelse(COI.1$FilterID %in%  c("4","12","5","13","6","14"), "B", COI.1$Site)

# Add bottle column based on filterID (only interested in bottles 1-6)
COI.1$Bottle <- ifelse(COI.1$FilterID %in%  c("1","9"), "1", NA)
COI.1$Bottle <- ifelse(COI.1$FilterID %in%  c("2","10"), "2", COI.1$Bottle)
COI.1$Bottle <- ifelse(COI.1$FilterID %in%  c("3","11"), "3", COI.1$Bottle)
COI.1$Bottle <- ifelse(COI.1$FilterID %in%  c("4","12"), "4", COI.1$Bottle)
COI.1$Bottle <- ifelse(COI.1$FilterID %in%  c("5","13"), "5", COI.1$Bottle)
COI.1$Bottle <- ifelse(COI.1$FilterID %in%  c("6","14"), "6", COI.1$Bottle)
COI.1$Bottle <- ifelse(COI.1$FilterID %in%  c("7"), NA, COI.1$Bottle)

# Only work with data from bottles 1-6
COI.2 <- COI.1[!is.na(COI.1$Bottle),]

# get matrix for rarefaction
COI.esv<-reshape2::dcast(COI.1, Molecule+Bottle ~ COI.1$Amplicon_GlobalESV, value.var = "ESVsize", fun.aggregate = sum)

# merge first 2 colmns into one, then remove them
COI.esv$sample<-paste(COI.esv$Molecule, COI.esv$Bottle, sep="_")
COI.esv <- COI.esv[,-c(1:2)]

# move sample to rownames then delete
rownames(COI.esv) <- COI.esv$sample
COI.esv$sample <- NULL

# remove columns with only zeros
COI.esv_notnull<-COI.esv[,colSums(COI.esv) !=0]

#remove rows with only zeros & edit rownames
COI.esv_notnull2<-COI.esv_notnull[rowSums(COI.esv_notnull) !=0,]

# Remove PCR controls, lab negative, and extraction negative before plotting
COI.esv_notnull2<-COI.esv_notnull2[-c(15:17),]

# Remove DNAse treated RNA samples before calculating 15th percentile
COI.esv_notnull2.2 <- COI.esv_notnull2[-c(15:21),]

#calculate 15th percentile for rrarefy function
COI.esv_15percentile<-quantile(rowSums(COI.esv_notnull2.2), prob=0.15)
# 15% 
# 105839.8 

#################################################################
# Create presence absence figure
#################################################################

# Rarefy the dataset down to the 15th percentile before plotting
set.seed(1234)
COI.rarefied <- rrarefy(COI.esv_notnull2, sample=COI.esv_15percentile)

# Convert to presence-absence matrix
COI.rarefied[COI.rarefied > 0] <-1

# Create venn

# subset by bottle
COI.1 <- COI.rarefied[grepl("_1", rownames(COI.rarefied)),]
COI.2 <- COI.rarefied[grepl("_2", rownames(COI.rarefied)),]
COI.3 <- COI.rarefied[grepl("_3", rownames(COI.rarefied)),]
COI.4 <- COI.rarefied[grepl("_4", rownames(COI.rarefied)),]
COI.5 <- COI.rarefied[grepl("_5", rownames(COI.rarefied)),]
COI.6 <- COI.rarefied[grepl("_6", rownames(COI.rarefied)),]

# calc total richness
COI.1.tot.richness <- specnumber(colSums(COI.1))
COI.2.tot.richness <- specnumber(colSums(COI.2))
COI.3.tot.richness <- specnumber(colSums(COI.3))
COI.4.tot.richness <- specnumber(colSums(COI.4))
COI.5.tot.richness <- specnumber(colSums(COI.5))
COI.6.tot.richness <- specnumber(colSums(COI.6))

# transpose for venn counts
COI.1.t <- t(COI.1)
COI.2.t <- t(COI.2)
COI.3.t <- t(COI.3)
COI.4.t <- t(COI.4)
COI.5.t <- t(COI.5)
COI.6.t <- t(COI.6)

# create counts for Venn
X1.venn <- vennCounts(COI.1.t)
X2.venn <- vennCounts(COI.2.t)
X3.venn <- vennCounts(COI.3.t)
X4.venn <- vennCounts(COI.4.t)
X5.venn <- vennCounts(COI.5.t)
X6.venn <- vennCounts(COI.6.t)

# ignore genera with zero counts
X1.venn[1,4] <- 0
X2.venn[1,4] <- 0
X3.venn[1,4] <- 0
X4.venn[1,4] <- 0
X5.venn[1,4] <- 0
X6.venn[1,4] <- 0

# use df with coordinates for circles (from above)
class(X1.venn) <- 'matrix'
df.1.venn <- as.data.frame(X1.venn)[-1,] %>%
  mutate(x = c(0, 1.2, 0.8, -1.2, -0.8, 0, 0),
         y = c(-0.8, 1.2, 0, 1.2, 0, 1.4, 0.4))
class(X2.venn) <- 'matrix'
df.2.venn <- as.data.frame(X2.venn)[-1,] %>%
  mutate(x = c(0, 1.2, 0.8, -1.2, -0.8, 0, 0),
         y = c(-0.8, 1.2, 0, 1.2, 0, 1.4, 0.4))
class(X3.venn) <- 'matrix'
df.3.venn <- as.data.frame(X3.venn)[-1,] %>%
  mutate(x = c(0, 1.2, 0.8, -1.2, -0.8, 0, 0),
         y = c(-0.8, 1.2, 0, 1.2, 0, 1.4, 0.4))
class(X4.venn) <- 'matrix'
df.4.venn <- as.data.frame(X4.venn)[-1,] %>%
  mutate(x = c(0, 1.2, 0.8, -1.2, -0.8, 0, 0),
         y = c(-0.8, 1.2, 0, 1.2, 0, 1.4, 0.4))
class(X5.venn) <- 'matrix'
df.5.venn <- as.data.frame(X5.venn)[-1,] %>%
  mutate(x = c(0, 1.2, 0.8, -1.2, -0.8, 0, 0),
         y = c(-0.8, 1.2, 0, 1.2, 0, 1.4, 0.4))
class(X6.venn) <- 'matrix'
df.6.venn <- as.data.frame(X6.venn)[-1,] %>%
  mutate(x = c(0, 1.2, 0.8, -1.2, -0.8, 0, 0),
         y = c(-0.8, 1.2, 0, 1.2, 0, 1.4, 0.4))

# draw venn with ggplot
COI.1.tot.richness.comma <- format(COI.1.tot.richness,big.mark=",", trim=TRUE)
COI.2.tot.richness.comma <- format(COI.2.tot.richness,big.mark=",", trim=TRUE)
COI.3.tot.richness.comma <- format(COI.3.tot.richness,big.mark=",", trim=TRUE)
COI.4.tot.richness.comma <- format(COI.4.tot.richness,big.mark=",", trim=TRUE)
COI.5.tot.richness.comma <- format(COI.5.tot.richness,big.mark=",", trim=TRUE)
COI.6.tot.richness.comma <- format(COI.6.tot.richness,big.mark=",", trim=TRUE)

X1 <- paste("COI (F230R + ml-jg)\nBottle 1 =",
            COI.1.tot.richness.comma)
X2 <- paste("Bottle 2 =",
            COI.2.tot.richness.comma)
X3 <- paste("Bottle 3 =",
            COI.3.tot.richness.comma)
X4 <- paste("Bottle 4 =",
            COI.4.tot.richness.comma)
X5 <- paste("Bottle 5 =",
            COI.5.tot.richness.comma)
X6 <- paste("Bottle 6 =",
            COI.6.tot.richness.comma)

v7 <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  ggtitle(paste0(X1)) +
  geom_circle(aes(fill=labels), size = 0.5) +
  geom_circle(size = 0.5, show.legend = FALSE, color = "black", alpha=0.3) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#33638DFF","#20A387FF", "#95D840FF")) +
  scale_color_manual(values=c("#33638DFF","#20A387FF", "#95D840FF"), guide = FALSE) +
  labs(fill = NULL) +
  annotate("text", x = df.1.venn$x, y = df.1.venn$y, label = df.1.venn$Counts, size = 3)

v8 <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  ggtitle(paste0(X2)) +
  geom_circle(aes(fill=labels), size = 0.5) +
  geom_circle(size = 0.5, show.legend = FALSE, color = "black", alpha=0.3) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#33638DFF","#20A387FF", "#95D840FF")) +
  scale_color_manual(values=c("#33638DFF","#20A387FF", "#95D840FF"), guide = FALSE) +
  labs(fill = NULL) +
  annotate("text", x = df.2.venn$x, y = df.2.venn$y, label = df.2.venn$Counts, size = 3)

v9 <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  ggtitle(paste0(X3)) +
  geom_circle(aes(fill=labels), size = 0.5) +
  geom_circle(size = 0.5, show.legend = FALSE, color = "black", alpha=0.3) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#33638DFF","#20A387FF", "#95D840FF")) +
  scale_color_manual(values=c("#33638DFF","#20A387FF", "#95D840FF"), guide = FALSE) +
  labs(fill = NULL) +
  annotate("text", x = df.3.venn$x, y = df.3.venn$y, label = df.3.venn$Counts, size = 3)

v10 <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  ggtitle(paste0(X4)) +
  geom_circle(aes(fill=labels), size = 0.5) +
  geom_circle(size = 0.5, show.legend = FALSE, color = "black", alpha=0.3) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#33638DFF","#20A387FF", "#95D840FF")) +
  scale_color_manual(values=c("#33638DFF","#20A387FF", "#95D840FF"), guide = FALSE) +
  labs(fill = NULL) +
  annotate("text", x = df.4.venn$x, y = df.4.venn$y, label = df.4.venn$Counts, size = 3)

v11 <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  ggtitle(paste0(X5)) +
  geom_circle(aes(fill=labels), size = 0.5) +
  geom_circle(size = 0.5, show.legend = FALSE, color = "black", alpha=0.3) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#33638DFF","#20A387FF", "#95D840FF")) +
  scale_color_manual(values=c("#33638DFF","#20A387FF", "#95D840FF"), guide = FALSE) +
  labs(fill = NULL) +
  annotate("text", x = df.5.venn$x, y = df.5.venn$y, label = df.5.venn$Counts, size = 3)

v12 <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  ggtitle(paste0(X6)) +
  geom_circle(aes(fill=labels), size = 0.5) +
  geom_circle(size = 0.5, show.legend = FALSE, color = "black", alpha=0.3) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#33638DFF","#20A387FF", "#95D840FF")) +
  scale_color_manual(values=c("#33638DFF","#20A387FF", "#95D840FF"), guide = FALSE) +
  labs(fill = NULL) +
  annotate("text", x = df.6.venn$x, y = df.6.venn$y, label = df.6.venn$Counts, size = 3)

##################################################
# 18Sv4
##################################################

# Presence absence plots
X18S <- read.table("taxonomic_assignments_18S_unpaired.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

# parse molecule and Filter ID from sample name
X18S.1 <- data.frame(X18S, do.call(rbind, str_split(X18S$SampleName, "_")), stringsAsFactors=FALSE)
names(X18S.1)[33:41] <- c("Name","Substrate","Molecule","FilterID","PCRStep","IlluminaSample","Lane","Read","Lane2")

# exclude RNAf (came from separate extraction)
X18S.1 <-X18S.1[X18S.1$Molecule != "RNAf",]

# Add site column based on filterID
X18S.1$Site <- ifelse(X18S.1$FilterID %in%  c("1","9","2","10","3","11"), "A", "")
X18S.1$Site <- ifelse(X18S.1$FilterID %in%  c("4","12","5","13","6","14"), "B", X18S.1$Site)

# Add bottle column based on filterID (only interested in bottles 1-6)
X18S.1$Bottle <- ifelse(X18S.1$FilterID %in%  c("1","9"), "1", NA)
X18S.1$Bottle <- ifelse(X18S.1$FilterID %in%  c("2","10"), "2", X18S.1$Bottle)
X18S.1$Bottle <- ifelse(X18S.1$FilterID %in%  c("3","11"), "3", X18S.1$Bottle)
X18S.1$Bottle <- ifelse(X18S.1$FilterID %in%  c("4","12"), "4", X18S.1$Bottle)
X18S.1$Bottle <- ifelse(X18S.1$FilterID %in%  c("5","13"), "5", X18S.1$Bottle)
X18S.1$Bottle <- ifelse(X18S.1$FilterID %in%  c("6","14"), "6", X18S.1$Bottle)
X18S.1$Bottle <- ifelse(X18S.1$FilterID %in%  c("7"), NA, X18S.1$Bottle)

# Only work with data from bottles 1-6
X18S.2 <- X18S.1[!is.na(X18S.1$Bottle),]

# get matrix for rarefaction
X18S.esv<-reshape2::dcast(X18S.2, Molecule+Bottle ~ X18S.2$GlobalASV, value.var = "ASVsize", fun.aggregate = sum)

# merge first 2 colmns into one, then remove them
X18S.esv$sample<-paste(X18S.esv$Molecule, X18S.esv$Bottle, sep="_")
X18S.esv <- X18S.esv[,-c(1:2)]

# move sample to rownames then delete
rownames(X18S.esv) <- X18S.esv$sample
X18S.esv$sample <- NULL

# remove columns with only zeros
X18S.esv_notnull<-X18S.esv[,colSums(X18S.esv) !=0]

#remove rows with only zeros & edit rownames
X18S.esv_notnull2<-X18S.esv_notnull[rowSums(X18S.esv_notnull) !=0,]

# remove PCR controls, field negative, and DNA extraction negative before plotting
X18S.esv_notnull2 <- X18S.esv_notnull2[-c(13:14),]

# Remove DNAse treated RNA samples before calculating 15th percentile
X18S.esv_notnull2.2 <- X18S.esv_notnull2[-c(13:18),]

#calculate 15th percentile for rrarefy function
X18S.esv_15percentile<-quantile(rowSums(X18S.esv_notnull2.2), prob=0.15)
# 15% 
# 119048.4  

#################################################################
# Create presence absence figure
#################################################################

# Rarefy the dataset down to the 15th percentile before plotting
set.seed(1234)
X18S.rarefied <- rrarefy(X18S.esv_notnull2, sample=X18S.esv_15percentile)

# Convert to presence-absence matrix
X18S.rarefied[X18S.rarefied > 0] <-1

# Create venn

# subset by bottle
X18S.1 <- X18S.rarefied[grepl("_1", rownames(X18S.rarefied)),]
X18S.2 <- X18S.rarefied[grepl("_2", rownames(X18S.rarefied)),]
X18S.3 <- X18S.rarefied[grepl("_3", rownames(X18S.rarefied)),]
X18S.4 <- X18S.rarefied[grepl("_4", rownames(X18S.rarefied)),]
X18S.5 <- X18S.rarefied[grepl("_5", rownames(X18S.rarefied)),]
X18S.6 <- X18S.rarefied[grepl("_6", rownames(X18S.rarefied)),]

# calc total richness
X18S.1.tot.richness <- specnumber(colSums(X18S.1))
X18S.2.tot.richness <- specnumber(colSums(X18S.2))
X18S.3.tot.richness <- specnumber(colSums(X18S.3))
X18S.4.tot.richness <- specnumber(colSums(X18S.4))
X18S.5.tot.richness <- specnumber(colSums(X18S.5))
X18S.6.tot.richness <- specnumber(colSums(X18S.6))

# transpose for venn counts
X18S.1.t <- t(X18S.1)
X18S.2.t <- t(X18S.2)
X18S.3.t <- t(X18S.3)
X18S.4.t <- t(X18S.4)
X18S.5.t <- t(X18S.5)
X18S.6.t <- t(X18S.6)

# create counts for Venn
X1.venn <- vennCounts(X18S.1.t)
X2.venn <- vennCounts(X18S.2.t)
X3.venn <- vennCounts(X18S.3.t)
X4.venn <- vennCounts(X18S.4.t)
X5.venn <- vennCounts(X18S.5.t)
X6.venn <- vennCounts(X18S.6.t)

# ignore genera with zero counts
X1.venn[1,4] <- 0
X2.venn[1,4] <- 0
X3.venn[1,4] <- 0
X4.venn[1,4] <- 0
X5.venn[1,4] <- 0
X6.venn[1,4] <- 0

# use df with coordinates for circles (from above)

# transform venncounts into df, add coordinates for labels
class(X1.venn) <- 'matrix'
df.1.venn <- as.data.frame(X1.venn)[-1,] %>%
  mutate(x = c(0, 1.2, 0.8, -1.2, -0.8, 0, 0),
         y = c(-0.8, 1.2, 0, 1.2, 0, 1.4, 0.4))
class(X2.venn) <- 'matrix'
df.2.venn <- as.data.frame(X2.venn)[-1,] %>%
  mutate(x = c(0, 1.2, 0.8, -1.2, -0.8, 0, 0),
         y = c(-0.8, 1.2, 0, 1.2, 0, 1.4, 0.4))
class(X3.venn) <- 'matrix'
df.3.venn <- as.data.frame(X3.venn)[-1,] %>%
  mutate(x = c(0, 1.2, 0.8, -1.2, -0.8, 0, 0),
         y = c(-0.8, 1.2, 0, 1.2, 0, 1.4, 0.4))
class(X4.venn) <- 'matrix'
df.4.venn <- as.data.frame(X4.venn)[-1,] %>%
  mutate(x = c(0, 1.2, 0.8, -1.2, -0.8, 0, 0),
         y = c(-0.8, 1.2, 0, 1.2, 0, 1.4, 0.4))
class(X5.venn) <- 'matrix'
df.5.venn <- as.data.frame(X5.venn)[-1,] %>%
  mutate(x = c(0, 1.2, 0.8, -1.2, -0.8, 0, 0),
         y = c(-0.8, 1.2, 0, 1.2, 0, 1.4, 0.4))
class(X6.venn) <- 'matrix'
df.6.venn <- as.data.frame(X6.venn)[-1,] %>%
  mutate(x = c(0, 1.2, 0.8, -1.2, -0.8, 0, 0),
         y = c(-0.8, 1.2, 0, 1.2, 0, 1.4, 0.4))

# draw venn with ggplot
X18S.1.tot.richness.comma <- format(X18S.1.tot.richness,big.mark=",", trim=TRUE)
X18S.2.tot.richness.comma <- format(X18S.2.tot.richness,big.mark=",", trim=TRUE)
X18S.3.tot.richness.comma <- format(X18S.3.tot.richness,big.mark=",", trim=TRUE)
X18S.4.tot.richness.comma <- format(X18S.4.tot.richness,big.mark=",", trim=TRUE)
X18S.5.tot.richness.comma <- format(X18S.5.tot.richness,big.mark=",", trim=TRUE)
X18S.6.tot.richness.comma <- format(X18S.6.tot.richness,big.mark=",", trim=TRUE)

X1 <- paste("18S V4\nBottle 1 =",
            X18S.1.tot.richness.comma)
X2 <- paste("Bottle 2 =",
            X18S.2.tot.richness.comma)
X3 <- paste("Bottle 3 =",
            X18S.3.tot.richness.comma)
X4 <- paste("Bottle 4 =",
            X18S.4.tot.richness.comma)
X5 <- paste("Bottle 5 =",
            X18S.5.tot.richness.comma)
X6 <- paste("Bottle 6 =",
            X18S.6.tot.richness.comma)

v13 <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  ggtitle(paste0(X1)) +
  geom_circle(aes(fill=labels), size = 0.5) +
  geom_circle(size = 0.5, show.legend = FALSE, color = "black", alpha=0.3) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#33638DFF","#20A387FF", "#95D840FF")) +
  scale_color_manual(values=c("#33638DFF","#20A387FF", "#95D840FF"), guide = FALSE) +
  labs(fill = NULL) +
  annotate("text", x = df.1.venn$x, y = df.1.venn$y, label = df.1.venn$Counts, size = 3)

v14 <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  ggtitle(paste0(X2)) +
  geom_circle(aes(fill=labels), size = 0.5) +
  geom_circle(size = 0.5, show.legend = FALSE, color = "black", alpha=0.3) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#33638DFF","#20A387FF", "#95D840FF")) +
  scale_color_manual(values=c("#33638DFF","#20A387FF", "#95D840FF"), guide = FALSE) +
  labs(fill = NULL) +
  annotate("text", x = df.2.venn$x, y = df.2.venn$y, label = df.2.venn$Counts, size = 3)

v15 <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  ggtitle(paste0(X3)) +
  geom_circle(aes(fill=labels), size = 0.5) +
  geom_circle(size = 0.5, show.legend = FALSE, color = "black", alpha=0.3) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#33638DFF","#20A387FF", "#95D840FF")) +
  scale_color_manual(values=c("#33638DFF","#20A387FF", "#95D840FF"), guide = FALSE) +
  labs(fill = NULL) +
  annotate("text", x = df.3.venn$x, y = df.3.venn$y, label = df.3.venn$Counts, size = 3)

v16 <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  ggtitle(paste0(X4)) +
  geom_circle(aes(fill=labels), size = 0.5) +
  geom_circle(size = 0.5, show.legend = FALSE, color = "black", alpha=0.3) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#33638DFF","#20A387FF", "#95D840FF")) +
  scale_color_manual(values=c("#33638DFF","#20A387FF", "#95D840FF"), guide = FALSE) +
  labs(fill = NULL) +
  annotate("text", x = df.4.venn$x, y = df.4.venn$y, label = df.4.venn$Counts, size = 3)

v17 <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  ggtitle(paste0(X5)) +
  geom_circle(aes(fill=labels), size = 0.5) +
  geom_circle(size = 0.5, show.legend = FALSE, color = "black", alpha=0.3) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#33638DFF","#20A387FF", "#95D840FF")) +
  scale_color_manual(values=c("#33638DFF","#20A387FF", "#95D840FF"), guide = FALSE) +
  labs(fill = NULL) +
  annotate("text", x = df.5.venn$x, y = df.5.venn$y, label = df.5.venn$Counts, size = 3)

v18 <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  ggtitle(paste0(X6)) +
  geom_circle(aes(fill=labels), size = 0.5) +
  geom_circle(size = 0.5, show.legend = FALSE, color = "black", alpha=0.3) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#33638DFF","#20A387FF", "#95D840FF")) +
  scale_color_manual(values=c("#33638DFF","#20A387FF", "#95D840FF"), guide = FALSE) +
  labs(fill = NULL) +
  annotate("text", x = df.6.venn$x, y = df.6.venn$y, label = df.6.venn$Counts, size = 3)

g <- ggarrange(v1, v7, v13, v2, v8, v14, v3, v9, v15, v4, v10, v16, v5, v11, v17, v6, v12, v18, nrow=6, ncol=3, common.legend = TRUE, legend="bottom")
ggsave("venn_esv_120619.pdf", g, width = 8, height = 10)

