library("DEP")
library("dplyr")
library("EnrichmentBrowser")
library("limma")
library("Biobase")
library("dendextend")
library("RSkittleBrewer")
library("preprocessCore")
library("vsn")
library("broom")
library("snpStats")
library("sva")
library("genefilter")
library("SummarizedExperiment")


#Read proteinGroups file

proteinGroups <- read.table(file = "proteinGroups.txt",
                            header = TRUE,sep = "\t",quote = "\"'",dec = ".",numerals = c("warn.loss"),
                            row.names = NULL,na.strings = c("NA","NaN","Infinite"))   


#Filter dataset

#Filter for contaminant proteins, reverse hits, and only identifed my site.

proteinGroups <- filter(proteinGroups, Reverse != "+", Potential.contaminant != "+", Only.identified.by.site != "+")



#Filter for protines with 2 or more unique peptides

proteinGroups <- filter(proteinGroups, Unique.peptides > 1)



#Filter for protines with 2 or more unique peptides

proteinGroups <- filter(proteinGroups, Reporter.intensity.count.1 > 11)



## Make row ID unique name

# Are there any duplicated gene names?

proteinGroups$Gene.names %>% duplicated() %>% any() 



# Make unique names using the annotation in the "Gene.names" column as primary names and the annotation in "Protein.IDs" as name for those that do not have an gene name.

proteinGroups <- make_unique(proteinGroups, "Gene.names", "Protein.IDs", delim = ";")

rownames(proteinGroups) <- proteinGroups[,"name"]



#Build an ExpressionSet from proteinGroups data.frame



#1 Assay data



#create a minimal ExpressionSet object with just the assay data wanted

assaydata <- data.matrix(proteinGroups[,c(32:40,42)]) #channel 10 removed, outlier in PCA

DSM2020_O <- ExpressionSet(assayData=assaydata)



#2 Phenotypic data.  Data describing treatment, control, batch, and other covariates.



#Make a phenotype table describing your experiment.
#phenotable <- matrix(c("C","D","C","D","C","D","C","D","C","D"))
phenotable <- matrix(c("C","D","C","C","D","C","D","D","C","D"))


colnames(phenotable) <- c("Treatment")

rownames(phenotable) <- colnames(exprs(DSM2020_O))

phenotable <- as.data.frame(phenotable)

phenotable <- new("AnnotatedDataFrame", data = phenotable)



# Verify row names of phenotable are the same as the column names of the assay data of the expression set.

all(rownames(phenotable)==colnames(exprs(DSM2020_O)))

DSM2020_O <- ExpressionSet(assayData=assaydata, phenoData=phenotable)



#3 Feature data.

featuretable <- proteinGroups[c(1:31,43:68)]

featuretable <- as.data.frame(featuretable)

featuretable <- new("AnnotatedDataFrame", data = featuretable)



# Verify row names of featuretable are the same as the row names of the assay data of the expression set.

all(rownames(featuretable)==rownames(exprs(DSM2020_O)))

DSM2020_O <- ExpressionSet(assayData=assaydata, phenoData=phenotable, featureData=featuretable)

#ExpressionSet is ready



#Load dataset

pdata=pData(DSM2020_O)
edata=exprs(DSM2020_O)
fdata = fData(DSM2020_O)



#Log2 transform all edata

edata = log2(edata)



#Quantile normalization
#norm_edata = normalize.quantiles(as.matrix(edata))
#colnames(norm_edata) <- colnames(edata)
#row.names(norm_edata) <- row.names(edata)
#norm_edata <- normalizeBetweenArrays(edata, method = "cyclicloess", cyclic.method="affy")
norm_edata <- justvsn(edata)

# By default calculates the euclidean distance between rows.

dist1 = dist(t(norm_edata))



#Heatmap of euclidean distance between samples.

colramp = colorRampPalette(c(3,"white",2))(9)
heatmap(as.matrix(dist1),col=colramp,Colv=NA,Rowv=NA)



#hierarchical clustering of samples.

hclust1 = hclust(dist1)
plot(hclust1,hang=-1)



#Remove unknown batch variables with surrogate variable analysis.

mod = model.matrix(~Treatment,data=pdata)
mod0 = model.matrix(~1, data=pdata)
sva1 = sva(norm_edata,mod,mod0) #n.sv=2



#Add the surrogate variables to the model matrix and perform the model fit.

modsv = cbind(mod,sva1$sv)
fitsv = lm.fit(modsv,t(norm_edata))



#statistics with limma

fit_limma = lmFit (norm_edata, modsv)
ebayes_limma = eBayes(fit_limma)
head(ebayes_limma)

names(ebayes_limma)
hist(ebayes_limma$p.value[,"TreatmentD"],col=2)


 
#adjust P-values for multiple hypothesis testing.  Benjamini Hochberg

qstats_obj <- p.adjust(ebayes_limma$p.value[,"TreatmentD"], method = "BH", length(ebayes_limma$p.value[,"TreatmentD"]))
sum(qstats_obj < 0.05)
hist(qstats_obj,col=2)


#Make a proteingroups object with limma p values, q values, and expression diffrences

expdiff <- as.matrix(rowMeans(norm_edata[,c(2,5,7,8,10)])-rowMeans(norm_edata[,c(1,3,4,6,9)]))
colnames(expdiff) <- "log2 fold change"

p.value <- ebayes_limma$p.value[,"TreatmentD"]
p.value <- as.data.frame(p.value)
names(p.value)[1] <- "p value"

qstats_obj <- as.matrix(qstats_obj)
qstats_obj <- as.data.frame(qstats_obj)
names(qstats_obj)[1] <- "q value"

negLogq = -log2(qstats_obj)
negLogq <- as.data.frame(negLogq)
names(negLogq)[1] <- "#-log2(q)"


stats_table <- cbind( expdiff, p.value, qstats_obj, negLogq)

proteinGroups[,c(32:40,42)] <- norm_edata[,1:10]
proteinGroups_limma <- merge(proteinGroups, stats_table, by=0)

write.csv(proteinGroups_limma, file ="proteinGroups_limma_vehicle_SVA.csv")







