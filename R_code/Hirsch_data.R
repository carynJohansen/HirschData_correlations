#Purpose: the purpose of this script is to clean and process the Hirsch expression and phenotype data for the analysis

setwd("~/Dropbox/UCD/rotations/RRI/Hirsch_correlations")

##########################
# Expression data

#the main thing here is to remove all rows that have a zero - get one data set with genes that express in all varieties

# load the expression data

maize_exp <- read.table("data/processed/maize_503genotypes_raw_expression_values_FPKM_cleaned.txt", sep="\t", header=TRUE)

#where are the NA's? Are there patterns to where they are located?
exp.na <- maize_exp[!complete.cases(maize_exp),]
#genes with NA's, and number of NA's

na_counting_function <- function(v) {
  sum(!complete.cases(v))
}

na.info <- data.frame(gene =exp.na$gene, num.Na = apply(exp.na, 1, na_counting_function))
hist(na.info$num.Na, breaks=50)
summary(na.info)

na.info[na.info$num.Na == 503,]
#GRMZM2G070092

maize_exp[maize_exp$gene == "GRMZM2G070092",]

# genes with NA's
maize_exp <- maize_exp[complete.cases(maize_exp),]
dim(maize_exp)

#just the expression values
exp <- maize_exp[,c(1,5:507)]

exp_limtedset <- exp[rowSums(exp[,2:504]==0) == 0,]
dim(exp_limtedset)
exp_limtedset <- exp_limtedset[complete.cases(exp_limtedset),]
dim(exp_limtedset)

#why is this different from their dataset?
# they have 14,968 annotatedB73 reference benes and 1425 RTAs expression in every line
# they used cufflinks confidence interval boundaries to determine expression

min(exp_limtedset)
max(exp_limtedset)

y <- unlist(exp_limtedset)
mean(y)

#take out rows with really really low expression, that is still above zero
exp_limtedset2 <- exp_limtedset[rowSums(exp_limtedset[,2:504] <= 0.5) == 0,]
write.table(exp_limtedset2, "data/processed/maize_503genotypes_FPKM_expressed.txt", quote = F, row.names = F, col.names = T)
#unlist:
x <- unlist(exp_limtedset2)
mean(x)
min(x)

############################
# Phenotypic data

transition <- read.table("data/processed/tpc119982SupplementalDS6_cleaned.csv", sep=",", header=T)
phenotypes <- read.table("data/processed/WIDIV_2010_paper_phenotypic_data_cleaned.csv", sep=",", header=T, na.strings = ".")

# first, take average of the phenotype data, because it has reps
str(phenotypes)
phenotypes$Rep <- as.factor(as.character(phenotypes$Rep))

measures <- names(phenotypes)[7:20]

p.rep1 <- phenotypes[phenotypes$Rep == 1,c(3,7:20)]
p.rep2 <- phenotypes[phenotypes$Rep == 2,c(3,7:20)]

#create phenotype list
library(data.table)

#create averages across the two reps, by the genotype (in their data frame, this is "Entry"), ignoring NA's
p.average <- rbindlist(list(p.rep1, p.rep2))[,lapply(.SD, mean, na.rm=TRUE), list(Entry)]

# there are more genotypes in the phenotype data than in the expression data
#WHY????
e.gtypes <- names(maize_exp)[5:507]
p.gtypes <- levels(p.average$Entry)

# Create one data.frame with all the phenotype data:

# In phenotype data frame, remove genotypes that do not have expression data

#vector of the genotypes not in the expressoin:
p.notin.e <- setdiff(levels(p.average$Entry), names(maize_exp)[5:507])
p.average <- p.average[!(p.average$Entry %in% p.notin.e),]

#alter the expression, transition, and p.average data.frames to all have the same genotypes name colum
#then rbind.fill

m.exp <- maize_exp[,c(1,5:507)]
names(m.exp) <- c("trait", names(m.exp)[2:504])



#reshape Phenotype data
library(reshape2)

#transition data
m.transition <- melt(transition, id.vars = "Genotype")
w.transition <- reshape(m.transition, idvar = "variable", timevar = "Genotype", direction = "wide")
names(w.transition) <- gsub("value.", "", names(w.transition))
names(w.transition) <- c("trait", names(w.transition)[2:503])

#Physiological data
m.phenotype <- melt(p.average, id.vars = "Entry")
w.phenotype <- reshape(m.phenotype, idvar="variable", timevar="Entry", direction="wide")
names(w.phenotype) <- gsub("value.", "", names(w.phenotype))
names(w.phenotype) <- c("trait", names(w.phenotype)[2:416])

library(plyr)
all.phe <- rbind.fill(list(w.phenotype, w.transition, m.exp))
write.table(all.phe, "data/processed/maize_alltraits_G_nT.csv", sep=",", quote=FALSE, row.names = FALSE, col.names = TRUE)


#make with cleaned up maize expression - zero-expression removed

names(exp_limtedset2) <- c("trait", names(exp_limtedset2)[2:504])

all.phe <- rbind.fill(list(w.phenotype, w.transition, exp_limtedset2))
write.table(all.phe, "data/processed/maize_g_nT.csv", sep=",", quote=FALSE, row.names = FALSE, col.names = TRUE)
