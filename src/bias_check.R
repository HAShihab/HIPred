#!/usr/bin/Rscript

# benchmarks
df <- read.csv("./data/Petrovski_DatasetS1.csv", sep="\t", head=TRUE)

omim  <- subset(df, OMIM.Haploinsufficiency != "")
omimd <- subset(df, OMIM.de.novo...Haploinsuficciency != "")
mgil  <- subset(df, MGI.Lethality.orthologs != "")
mgis  <- subset(df, MGI.Seizure.orthologs != "")



# parse PubMed records
df <- read.csv("./tmp/PubMed_Counts.txt", sep=" ", head=FALSE)

omim  <- df[omim[, 'OMIM.disease.genes'], ]$V2
omimd <- df[omimd[, 'OMIM.disease.genes'], ]$V2
mgil  <- df[mgil[, 'OMIM.disease.genes'], ]$V2
mgis  <- df[mgis[, 'OMIM.disease.genes'], ]$V2

#

wilcox.test(mgil, omim, alternative="less")
wilcox.test(mgil, omimd, alternative="less")

wilcox.test(mgis, omim, alternative="less")
wilcox.test(mgis, omimd, alternative="less")


essm <- read.csv("./data/ASD2.txt", head=FALSE)
essm <- df[essm$V1, ]$V2

wilcox.test(essm, omim, alternative="less")
wilcox.test(essm, omimd, alternative="less")


