# Teresita M. Porter, Dec. 5, 2019

library(stringr)
library(reshape2)
library(vegan) # rarecurve
library(purrr) # for map_dfr
library(ggplot2)
library(scales)
library(data.table)
library(gridExtra)
#library(ggpubr) # for common legend
library(cowplot) # get legend

###################################################################
# Edit rarecurve function to remove the horizontal lines
rarecurve2 <- function (x, step = 1, sample, xlab = "Sample Size", ylab = "Species", 
                        label = TRUE, col, lty, ...) 
{
  x <- as.matrix(x)
  if (!identical(all.equal(x, round(x)), TRUE)) 
    stop("function accepts only integers (counts)")
  if (missing(col)) 
    col <- par("col")
  if (missing(lty)) 
    lty <- par("lty")
  tot <- rowSums(x)
  S <- specnumber(x)
  if (any(S <= 0)) {
    message("empty rows removed")
    x <- x[S > 0, , drop = FALSE]
    tot <- tot[S > 0]
    S <- S[S > 0]
  }
  nr <- nrow(x)
  col <- rep(col, length.out = nr)
  lty <- rep(lty, length.out = nr)
  out <- lapply(seq_len(nr), function(i) {
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) 
      n <- c(n, tot[i])
    drop(rarefy(x[i, ], n))
  })
  Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
  Smax <- sapply(out, max)
  plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = xlab, ylab = ylab, 
       type = "n", ...)
  if (!missing(sample)) {
    abline(v = sample) 
    rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"), 
                                           y = z, xout = sample, rule = 1)$y)
    #    abline(h = rare, lwd = 0.5) #turn off horizontal lines
  }
  for (ln in seq_along(out)) {
    N <- attr(out[[ln]], "Subsample")
    lines(N, out[[ln]], col = col[ln], lty = lty[ln], ...)
  }
  if (label) { 
    ordilabel(cbind(tot, S), labels = rownames(x), ...)
  }
  invisible(out)
}

###################################################################
# 16Sv3v4
###################################################################

# Presence absence plots
X16S <- read.table("16S_taxonomic_assignments.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

# parse molecule and Filter ID from sample name
X16S.1 <- data.frame(X16S, do.call(rbind, str_split(X16S$SampleName, "_")), stringsAsFactors=FALSE)
names(X16S.1)[33:39] <- c("Name","Substrate","Molecule","FilterID","PCRStep","IlluminaSample","Blank")

# Add site column based on filterID
X16S.1$Site <- ifelse(X16S.1$FilterID %in%  c("1","9","2","10","3","11"), "A", "")
X16S.1$Site <- ifelse(X16S.1$FilterID %in%  c("4","12","5","13","6","14"), "B", X16S.1$Site)

# filter for good family rank assignments
# Based on RDP website for sequences > 250bp (default)
# https://rdp.cme.msu.edu/classifier/classifier.jsp
X16S.f <- X16S.1[X16S.1$fBP >= 0.80,]

# get matrix for rarefaction for boxplot
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

# remove PCR controls, lab negative, and DNA extraction negative before plotting
X16S.esv_notnull2<-X16S.esv_notnull2[-c(1:4,17:20,33:37,50:51),]

# Remove DNAse treated RNA samples before calculating 15th percentile
X16S.esv_notnull2.2 <- X16S.esv_notnull2[-c(25:28),]

#calculate 15th percentile for rrarefy function (including PCRpc and PCRnc)
X16S.esv_15percentile<-quantile(rowSums(X16S.esv_notnull2.2), prob=0.15)
# 15% 
# 3809.7  

# Plot rarefaction curve (reflects sampling for box plot)
set.seed(1234)
X16S.esv_rarecurveout<-rarecurve2(X16S.esv_notnull2, 
                                      sample=X16S.esv_15percentile, 
                                      step=500, 
                                      label=T)

# Reformat vegan list as df (cols OTU, raw.read)
X16S.esv.rare <- lapply(X16S.esv_rarecurveout, function(x){
  b <- as.data.frame(x)
  b <- data.frame(OTU = b[,1], raw.read = rownames(b))
  b$raw.read <- as.numeric(gsub("N", "",  b$raw.read))
  return(b)
})

# Add sample names to vegan output (df) (rownames)
X16S.esv.sample_names <- rownames(X16S.esv_notnull2)
names(X16S.esv.rare) <- X16S.esv.sample_names

# Map rownames to vegan output (df)
X16S.esv.rare <- map_dfr(X16S.esv.rare, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "sample")

# Parse out metadata from sample
X16S.esv.rare<-data.frame(X16S.esv.rare, do.call(rbind, str_split(X16S.esv.rare$sample,"_")))
names(X16S.esv.rare)[4:7]<-c("Molecule","Site","FilterID","PCRStep")

# Create factors
X16S.esv.rare$Molecule <- factor(X16S.esv.rare$Molecule, 
                                     levels=c("PCRpc","cDNA", "DNA", "RNAf", "RNAt", "PCRnc"))
X16S.esv.rare$PCRStep <- factor(X16S.esv.rare$PCRStep, 
                                 levels=c("1STEP","2STEP"),
                                labels=c("1 step", "2 step"))

# This rarefaction curve reflects samples analyzed for the box plot 
# vertical line represents 15th percentile after excluding PCRpc and PCRnc
p1.tmp <- ggplot(data = X16S.esv.rare) +
  ggtitle("16S V3-V4") +
  labs(x="Reads", y="ESVs") +
  geom_point(aes(x = raw.read, y = OTU, color = Molecule, shape = PCRStep), size=2) +
  geom_vline(xintercept = X16S.esv_15percentile, linetype = "dashed") +
  scale_shape_manual(values = c(1,17)) +
  scale_x_continuous(label = comma) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10),
        legend.position = "right") +
  guides(shape=guide_legend(nrow=2), 
         color = guide_legend(nrow=2),
         byrow=TRUE)

l <- get_legend(p1.tmp)

p1<-ggplot(data = X16S.esv.rare) +
  ggtitle("16S V3-V4") +
  labs(x="Reads", y="ESVs") +
  geom_point(aes(x = raw.read, y = OTU, color = Molecule, shape = PCRStep), size=2) +
  geom_vline(xintercept = X16S.esv_15percentile, linetype = "dashed") +
  scale_shape_manual(values = c(1,17)) +
  scale_x_continuous(label = comma) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10),
        legend.position = "none") +
  guides(shape=guide_legend(nrow=2), 
         color = guide_legend(nrow=2),
         byrow=TRUE)

###################################################################
# COI (F230R and mlCOI plot separately)
###################################################################

# Presence absence plots
COI <- read.table("GJ_COI_taxonomic_assignments_raw.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

# parse molecule and Filter ID from sample name
COI.1 <- data.frame(COI, do.call(rbind, str_split(COI$SampleName, "_")), stringsAsFactors=FALSE)
names(COI.1)[34:36] <- c("Molecule","FilterID","PCRStep")

# parse amplicon from Amplicon_GlobalESV
COI.2 <- data.frame(COI.1, do.call(rbind, str_split(COI.1$Amplicon_GlobalESV, "_")), stringsAsFactors=FALSE)
names(COI.2)[37:38] <- c("Amplicon","GlobalESV_dontuse")

# Add site column based on filterID
COI.2$Site <- ifelse(COI.2$FilterID %in%  c("1","9","2","10","3","11"), "A", "")
COI.2$Site <- ifelse(COI.2$FilterID %in%  c("4","12","5","13","6","14"), "B", COI.2$Site)

# filter for good family rank assignments (99% correct)
# COI Classifier v3 cutoffs at https://github.com/terrimporter/CO1Classifier for seqs ~200 bp
COI.f <- COI.2[COI.2$fBP >= 0.20,]

# filter by amplicon
COI.F230 <- COI.f[COI.f$Amplicon =="FolF",]
COI.mlCOI <- COI.f[COI.f$Amplicon == "mlCOI",]

# get matrix for rarefaction for boxplot analysis
COI.esv.F230<-reshape2::dcast(COI.F230, Molecule+Site+FilterID+PCRStep ~ COI.F230$Amplicon_GlobalESV, value.var = "ESVsize", fun.aggregate = sum)
COI.esv.mlCOI<-reshape2::dcast(COI.mlCOI, Molecule+Site+FilterID+PCRStep ~ COI.mlCOI$Amplicon_GlobalESV, value.var = "ESVsize", fun.aggregate = sum)

# merge first 2 colmns into one, then remove them
COI.esv.F230$sample<-paste(COI.esv.F230$Molecule, COI.esv.F230$Site, COI.esv.F230$FilterID, COI.esv.F230$PCRStep, sep="_")
COI.esv.F230 <- COI.esv.F230[,-c(1:4)]

COI.esv.mlCOI$sample<-paste(COI.esv.mlCOI$Molecule, COI.esv.mlCOI$Site, COI.esv.mlCOI$FilterID, COI.esv.mlCOI$PCRStep, sep="_")
COI.esv.mlCOI <- COI.esv.mlCOI[,-c(1:4)]

# move sample to rownames then delete
rownames(COI.esv.F230) <- COI.esv.F230$sample
COI.esv.F230$sample <- NULL

rownames(COI.esv.mlCOI) <- COI.esv.mlCOI$sample
COI.esv.mlCOI$sample <- NULL

# remove columns with only zeros
COI.esv_notnull.F230 <- COI.esv.F230[,colSums(COI.esv.F230) !=0]
COI.esv_notnull.mlCOI <- COI.esv.mlCOI[,colSums(COI.esv.mlCOI) !=0]

#remove rows with only zeros & edit rownames
COI.esv_notnull2.F230 <- COI.esv_notnull.F230[rowSums(COI.esv_notnull.F230) !=0,]
COI.esv_notnull2.mlCOI <- COI.esv_notnull.mlCOI[rowSums(COI.esv_notnull.mlCOI) !=0,]

# Remove PCR controls, lab negative, DNA extraction negative before plotting
COI.esv_notnull2.F230 <- COI.esv_notnull2.F230[-c(1:3,16:18,31:36,45:46),]
COI.esv_notnull2.mlCOI <- COI.esv_notnull2.mlCOI[-c(13,26:29),]

# Remove DNAse treated RNA samples before calculating 15th percentile
COI.esv_notnull2.2.F230 <- COI.esv_notnull2.F230[-c(25:39),]
COI.esv_notnull2.2.mlCOI <- COI.esv_notnull2.mlCOI[-c(25:31),]

#calculate 15th percentile for rrarefy function (including PCRpc and PCRnc)
COI.esv_15percentile.F230 <- quantile(rowSums(COI.esv_notnull2.2.F230), prob=0.15)
# 15% 
# 11607

COI.esv_15percentile.mlCOI <- quantile(rowSums(COI.esv_notnull2.2.mlCOI), prob=0.15)
# 15% 
# 7829.7 

# Plot rarefaction curve (reflects sampling for box plot)
set.seed(1234)
COI.esv_rarecurveout.F230 <- rarecurve2(COI.esv_notnull2.F230, 
                                  sample=COI.esv_15percentile.F230, 
                                  step=500, 
                                  label=T)
COI.esv_rarecurveout.mlCOI <- rarecurve2(COI.esv_notnull2.mlCOI, 
                                        sample=COI.esv_15percentile.mlCOI, 
                                        step=500, 
                                        label=T)

# Reformat vegan list as df (cols OTU, raw.read)
COI.esv.rare.F230 <- lapply(COI.esv_rarecurveout.F230, function(x){
  b <- as.data.frame(x)
  b <- data.frame(OTU = b[,1], raw.read = rownames(b))
  b$raw.read <- as.numeric(gsub("N", "",  b$raw.read))
  return(b)
})
COI.esv.rare.mlCOI <- lapply(COI.esv_rarecurveout.mlCOI, function(x){
  b <- as.data.frame(x)
  b <- data.frame(OTU = b[,1], raw.read = rownames(b))
  b$raw.read <- as.numeric(gsub("N", "",  b$raw.read))
  return(b)
})

# Add sample names to vegan output (df) (rownames)
COI.esv.sample_names.F230 <- rownames(COI.esv_notnull2.F230)
names(COI.esv.rare.F230) <- COI.esv.sample_names.F230

COI.esv.sample_names.mlCOI <- rownames(COI.esv_notnull2.mlCOI)
names(COI.esv.rare.mlCOI) <- COI.esv.sample_names.mlCOI

# Map rownames to vegan output (df)
COI.esv.rare.F230 <- map_dfr(COI.esv.rare.F230, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "sample")
COI.esv.rare.mlCOI <- map_dfr(COI.esv.rare.mlCOI, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "sample")

# Parse out metadata from sample
COI.esv.rare.F230 <- data.frame(COI.esv.rare.F230, do.call(rbind, str_split(COI.esv.rare.F230$sample,"_")))
names(COI.esv.rare.F230)[4:7]<-c("Molecule","Site","FilterID","PCRStep")

COI.esv.rare.mlCOI <- data.frame(COI.esv.rare.mlCOI, do.call(rbind, str_split(COI.esv.rare.mlCOI$sample,"_")))
names(COI.esv.rare.mlCOI)[4:7]<-c("Molecule","Site","FilterID","PCRStep")

# Create factors
COI.esv.rare.F230$Molecule <- factor(COI.esv.rare.F230$Molecule, 
                                 levels=c("PCRpc","cDNA", "DNA", "RNAf", "RNAt", "PCRnc"))
COI.esv.rare.mlCOI$Molecule <- factor(COI.esv.rare.mlCOI$Molecule, 
                                     levels=c("PCRpc","cDNA", "DNA", "RNAf", "RNAt", "PCRnc"))

# color by study, plot with ggplot
# This curve reflects samples analyzed for the box plot 
# vertical line represents 15th percentile after excluding PCRpc and PCRnc
p2a<-ggplot(data = COI.esv.rare.F230) +
  ggtitle("COI (F230R)") +
  labs(x="Reads", y="ESVs") +
  geom_point(aes(x = raw.read, y = OTU, color = Molecule, shape = PCRStep), size=2) +
  geom_vline(xintercept = COI.esv_15percentile.F230, linetype = "dashed") +
  scale_shape_manual(values = c(1,17)) +
  scale_x_continuous(label = comma) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10),
        legend.position = "none") 

p2b <- ggplot(data = COI.esv.rare.mlCOI) +
  ggtitle("COI (ml-jg)") +
  labs(x="Reads", y="ESVs") +
  geom_point(aes(x = raw.read, y = OTU, color = Molecule, shape = PCRStep), size=2) +
  geom_vline(xintercept = COI.esv_15percentile.mlCOI, linetype = "dashed") +
  scale_shape_manual(values = c(1,17)) +
  scale_x_continuous(label = comma) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10),
        legend.position = "none") 

###################################################################
# 18Sv4
###################################################################

# Presence absence plots
X18S <- read.table("taxonomic_assignments_18S_unpaired.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

# parse molecule and Filter ID from sample name
X18S.1 <- data.frame(X18S, do.call(rbind, str_split(X18S$SampleName, "_")), stringsAsFactors=FALSE)
names(X18S.1)[33:41] <- c("Name","Substrate","Molecule","FilterID","PCRStep","IlluminaSample","Lane","Read","Lane2")

# parse amplicon from Amplicon_GlobalESV field
X18S.2 <- data.frame(X18S.1, do.call(rbind, str_split(X18S.1$GlobalASV, "_")), stringsAsFactors=FALSE)
names(X18S.2)[42:43] <- c("Amplicon","OTU_dontuse")

# Add site column based on filterID
X18S.2$Site <- ifelse(X18S.2$FilterID %in%  c("1","9","2","10","3","11"), "A", "")
X18S.2$Site <- ifelse(X18S.2$FilterID %in%  c("4","12","5","13","6","14"), "B", X18S.2$Site)

# filter for good class rank assignments
# Based on 18S v3 classifier (95% correct) for 200 bp fragments https://github.com/terrimporter/18SClassifier
X18S.f <- X18S.2[X18S.2$fBP >= 0.80,]

# filter by amplicon
X18S.F <- X18S.f[X18S.f$Amplicon =="Uni18SF",]
X18S.R <- X18S.f[X18S.f$Amplicon == "Uni18SR",]

# get matrix for rarefaction (for boxplot)
X18S.esv.F<-reshape2::dcast(X18S.F, Molecule+Site+FilterID+PCRStep ~ X18S.F$GlobalASV, value.var = "ASVsize", fun.aggregate = sum)
X18S.esv.R<-reshape2::dcast(X18S.R, Molecule+Site+FilterID+PCRStep ~ X18S.R$GlobalASV, value.var = "ASVsize", fun.aggregate = sum)

# merge first 2 colmns into one, then remove them
X18S.esv.F$sample<-paste(X18S.esv.F$Molecule, X18S.esv.F$Site, X18S.esv.F$FilterID, X18S.esv.F$PCRStep, sep="_")
X18S.esv.F <- X18S.esv.F[,-c(1:4)]

X18S.esv.R$sample<-paste(X18S.esv.R$Molecule, X18S.esv.R$Site, X18S.esv.R$FilterID, X18S.esv.R$PCRStep, sep="_")
X18S.esv.R <- X18S.esv.R[,-c(1:4)]

# move sample to rownames then delete
rownames(X18S.esv.F) <- X18S.esv.F$sample
X18S.esv.F$sample <- NULL

rownames(X18S.esv.R) <- X18S.esv.R$sample
X18S.esv.R$sample <- NULL

# remove columns with only zeros
X18S.esv_notnull.F<-X18S.esv.F[,colSums(X18S.esv.F) !=0]
X18S.esv_notnull.R<-X18S.esv.R[,colSums(X18S.esv.R) !=0]

#remove rows with only zeros & edit rownames
X18S.esv_notnull2.F<-X18S.esv_notnull.F[rowSums(X18S.esv_notnull.F) !=0,]
X18S.esv_notnull2.R<-X18S.esv_notnull.R[rowSums(X18S.esv_notnull.R) !=0,]

# Remove PCR controls, lab negative, DNA extraction negative before plotting
X18S.esv_notnull2.F<-X18S.esv_notnull.F[-c(1:4,16:17,30:33),]
X18S.esv_notnull2.R<-X18S.esv_notnull.R[-c(1:3,16,29:31),]

# Remove DNAse treated RNA samples before calculating 15th percentile
X18S.esv_notnull2.2.F <- X18S.esv_notnull2.F[-c(24:35),]
X18S.esv_notnull2.2.R <- X18S.esv_notnull2.R[-c(25:35),]

#calculate 15th percentile for rrarefy function (including PCRpc and PCRnc)
X18S.esv_15percentile.F<-quantile(rowSums(X18S.esv_notnull2.2.F), prob=0.15)
# 15% 
# 9841.1 
X18S.esv_15percentile.R<-quantile(rowSums(X18S.esv_notnull2.2.R), prob=0.15)
# 15% 
# 2628.25 

# Plot rarefaction curve (reflects sampling for box plot)
set.seed(1234)
X18S.esv_rarecurveout.F<-rarecurve2(X18S.esv_notnull2.F, 
                                  sample=X18S.esv_15percentile.F, 
                                  step=500, 
                                  label=T)
X18S.esv_rarecurveout.R<-rarecurve2(X18S.esv_notnull2.R, 
                                    sample=X18S.esv_15percentile.R, 
                                    step=500, 
                                    label=T)

# Reformat vegan list as df (cols OTU, raw.read)
X18S.esv.rare.F <- lapply(X18S.esv_rarecurveout.F, function(x){
  b <- as.data.frame(x)
  b <- data.frame(OTU = b[,1], raw.read = rownames(b))
  b$raw.read <- as.numeric(gsub("N", "",  b$raw.read))
  return(b)
})
X18S.esv.rare.R <- lapply(X18S.esv_rarecurveout.R, function(x){
  b <- as.data.frame(x)
  b <- data.frame(OTU = b[,1], raw.read = rownames(b))
  b$raw.read <- as.numeric(gsub("N", "",  b$raw.read))
  return(b)
})

# Add sample names to vegan output (df) (rownames)
X18S.esv.sample_names.F <- rownames(X18S.esv_notnull2.F)
names(X18S.esv.rare.F) <- X18S.esv.sample_names.F

X18S.esv.sample_names.R <- rownames(X18S.esv_notnull2.R)
names(X18S.esv.rare.R) <- X18S.esv.sample_names.R

# Map rownames to vegan output (df)
X18S.esv.rare.F <- map_dfr(X18S.esv.rare.F, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "sample")

X18S.esv.rare.R <- map_dfr(X18S.esv.rare.R, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "sample")

# Parse out metadata from sample
X18S.esv.rare.F<-data.frame(X18S.esv.rare.F, do.call(rbind, str_split(X18S.esv.rare.F$sample,"_")))
names(X18S.esv.rare.F)[4:7]<-c("Molecule","Site","FilterID","PCRStep")

X18S.esv.rare.R<-data.frame(X18S.esv.rare.R, do.call(rbind, str_split(X18S.esv.rare.R$sample,"_")))
names(X18S.esv.rare.R)[4:7]<-c("Molecule","Site","FilterID","PCRStep")

# Create factors
X18S.esv.rare.F$Molecule <- factor(X18S.esv.rare.F$Molecule, 
                                 levels=c("PCRpc","cDNA", "DNA", "RNAf", "RNAt", "PCRnc"))
X18S.esv.rare.R$Molecule <- factor(X18S.esv.rare.R$Molecule, 
                                   levels=c("PCRpc","cDNA", "DNA", "RNAf", "RNAt", "PCRnc"))

# color by study, plot with ggplot
# This curve reflects samples analyzed for the box plot 
# vertical line represents 15th percentile after excluding PCRpc and PCRnc
p3a<-ggplot(data = X18S.esv.rare.F) +
  ggtitle("18S V4 (Uni18S)") +
  labs(x="Reads", y="ESVs") +
  geom_point(aes(x = raw.read, y = OTU, color = Molecule, shape = PCRStep), size=2) +
  geom_vline(xintercept = X18S.esv_15percentile.F, linetype = "dashed") +
  scale_shape_manual(values = c(1,17)) +
  scale_x_continuous(label = comma) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10),
        legend.position = "none")

p3b<-ggplot(data = X18S.esv.rare.R) +
  ggtitle("18S V4 (Uni18SR)") +
  labs(x="Reads", y="ESVs") +
  geom_point(aes(x = raw.read, y = OTU, color = Molecule, shape = PCRStep), size=2) +
  geom_vline(xintercept = X18S.esv_15percentile.R, linetype = "dashed") +
  scale_shape_manual(values = c(1,17)) +
  scale_x_continuous(label = comma) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10),
        legend.position = "none")

lay <-rbind(c(1,2),
            c(3,4),
            c(5,6))

g <- grid.arrange(p1, l, p2a, p2b, p3a, p3b, layout_matrix = lay)
ggsave("RarefactionCurves_120519.pdf", g, height = 8, width = 8)
