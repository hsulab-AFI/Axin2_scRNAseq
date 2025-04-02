# Microarray analysis
BiocManager::install("GEOquery")
BiocManager::install("limma")
BiocManager::install("affy")
BiocManager::install("gcrma")

#Load the necessary libraries
library(GEOquery)
library(Biobase)
library(limma)
library(affy)
library(gcrma)

# load series and platform data from GEO
gse <- getGEO("GSE50796", GSEMatrix =TRUE, getGPL=TRUE)
if (length(gse) > 1) idx <- grep("GPL6883", attr(gse, "names")) else idx <- 1
gse <- gse[[idx]]
class(gse)

## check how many platforms used
length(gse)

pData(gse) ## print the sample information
fData(gse) ## print the gene annotation
exprs(gse) ## print the expression data

# assume the data are on a log2 scale; typically in the range of 0 to 16.

# Check the normalisation and scales used
# if the values go beyond 16, so we will need to perform a log2 transformation

## exprs get the expression levels as a data frame and get the distribution
summary(exprs(gse))
exprs(gse) <- log2(exprs(gse))
boxplot(exprs(gse),outline=FALSE)

# Inspect the clinical variables
library(dplyr)
sampleInfo <- pData(gse)
sampleInfo

# Let's pick just the columns contain information we need

sampleInfo <- select(sampleInfo, title, characteristics_ch1)

## Optionally, rename to more convenient column names
sampleInfo <- rename(sampleInfo,group = characteristics_ch1, patient = title)
Patient_id = c("Fused_1","Fused_2","Fused_3","Fused_4","Fused_5","Fused_6","Fused_7","Patent_1",
               "Patent_2","Patent_3","Patent_4","Patent_5","Patent_6","Patent_7")
sampleInfo$Patient_id = Patient_id


library(pheatmap)
## argument use="c" stops an error if there are any missing data points

corMatrix <- cor(exprs(gse),use="c")
pheatmap(corMatrix)  

## Print the rownames of the sample information and check it matches the correlation matrix
rownames(sampleInfo)
colnames(corMatrix)
## If not, force the rownames to match the columns
rownames(sampleInfo) <- colnames(corMatrix)
pheatmap(corMatrix,
         labels_row = Patient_id,
         labels_col = Patient_id,
         annotation_col=sampleInfo)


## compute PCA plot on the samples
library(ggplot2)
library(ggrepel)
## MAKE SURE TO TRANSPOSE THE EXPRESSION MATRIX

pca <- prcomp(t(exprs(gse)))

## Join the PCs to the sample information
cbind(sampleInfo, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=group,label=paste(Patient_id))) + geom_point() + geom_text_repel()


# Differential Expression
library(limma)
design <- model.matrix(~0+sampleInfo$group)
design
## the column names are a bit ugly, so we will rename
colnames(design) <- c("Patent","Synostosed")

summary(exprs(gse))

## calculate median expression level
cutoff <- median(exprs(gse))

## TRUE or FALSE for whether each gene is "expressed" in each sample
is_expressed <- exprs(gse) > cutoff

## Identify genes expressed in more than 2 samples

keep <- rowSums(is_expressed) > 2

## check how many genes are removed / retained.
table(keep)

## subset to just those expressed genes
gse <- gse[keep,]

fit <- lmFit(exprs(gse), design)
head(fit$coefficients)

contrasts <- makeContrasts(Synostosed - Patent, levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
topTable(fit2)
topTable(fit2, coef=1)

# If we want to know how many genes are differentially-expressed overall we can use the decideTests function.
decideTests(fit2)

# Further processing
anno <- fData(gse)
anno

anno <- select(anno,Symbol,Entrez_Gene_ID,Chromosome,Cytoband)
fit2$genes <- anno
topTable(fit2)

full_results <- topTable(fit2, number=Inf)
full_results
full_results <- tibble::rownames_to_column(full_results,"ID")

## Make sure you have ggplot2 loaded
library(ggplot2)
library(ggrepel)
p_cutoff <- 0.05
fc_cutoff <- 1
topN <- 100

full_results %>% 
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  mutate(Rank = 1:n(), Label = ifelse(Rank < topN, Symbol,"")) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant,label=Label)) + geom_point() + geom_text_repel(col="black")


# We can also filter according to p-value (adjusted) and fold-change cut-offs
p_cutoff <- 0.05
fc_cutoff <- 1
filter(full_results, P.Value < 0.05)

# These results can be exported with the write_csv function
library(readr)
write_csv(full_results, "full_results.csv")

filtered_results <- filter(full_results, P.Value < 0.05)
write_csv(filtered_results, "filtered_results_p_0.05_FC_1.csv")



# most differentially-expressed genes
## Use to top 10 genes for illustration
topN <- 10
ids_of_interest <- mutate(filtered_results, Rank = 1:n()) %>% 
  filter(Rank < topN) %>% 
  pull(ID)
gene_names <- mutate(filtered_results, Rank = 1:n()) %>% 
  filter(Rank < topN) %>% 
  pull(Symbol) 


## Get the rows corresponding to ids_of_interest and all columns
gene_matrix <- exprs(gse)[ids_of_interest,]
pheatmap(gene_matrix,
         labels_row = gene_names,
         scale = "row")

# sort filtered results in descending order
up = filtered_results[ which(filtered_results$P.Value < 0.05 & filtered_results$logFC > 0),]
down = filtered_results[ which(filtered_results$P.Value < 0.05 & filtered_results$logFC < 0),]
uptop10 = top_n(up, 10, logFC)
downtop10 = top_n(down, 10, -logFC)
sig = rbind(uptop10, downtop10)

topN <- 100
##
ids_of_interest <- mutate(sig, Rank = 1:n()) %>% 
  filter(Rank < topN) %>% 
  pull(ID)
gene_names <- mutate(sig, Rank = 1:n()) %>%
  filter(Rank < topN) %>% 
  pull(Symbol)
## Get the rows corresponding to ids_of_interest and all columns
gene_matrix <- exprs(gse)[ids_of_interest,]
Patients = c("Fused_7","Fused_4","Fused_2","Fused_1","Fused_3","Fused_5","Fused_6","Patent_4",
             "Patent_5","Patent_2","Patent_3","Patent_7","Patent_1","Patent_6")
pheatmap(gene_matrix,
         labels_row = gene_names,
         labels_col = Patient_id,
         scale = "row")


# Draw heatmap of enriched genes of signaling pathways
BMPgenes <- c("GDF10", "TMEM100", "TWSG1", "SFRP2", "TGFB3", "DLX5", "BMP8B", "SOST")
interested_pathways <- data.frame(BMPgenes)

## BMP
ids_of_BMPs <- mutate(filtered_results, Rank = 1:n()) %>% 
  filter(filtered_results$Symbol %in% interested_pathways$BMPgenes2) %>% 
  pull(ID)
gene_names_BMP <- mutate(filtered_results, Rank = 1:n()) %>% 
  filter(filtered_results$Symbol %in% interested_pathways$BMPgenes2) %>% 
  pull(Symbol) 
 
### Get the rows corresponding to ids_of_interest and all columns
BMP_matrix <- exprs(gse)[ids_of_BMPs,]
callback_rev <- function(hc, mat){
  hc$order <- rev(hc$order)
  hc
}
pheatmap(BMP_matrix,
         labels_row = gene_names_BMP,
         clustering_callback = callback_rev,
         scale = "row")


# Pathway analysis with DAVID tools results
dat1 = read.table("downregulated_GO_BP.txt", sep='\t', header=T)
#find p<0.05 GO term
dat = dat1[which(dat1$PValue < 0.05),]

#show top 20
show = dat[c(1:20),]

#order in counts
show = show[order(show$Count,decreasing=T),]
nrows = nrow(show)

#draw plot
ggplot(show, aes(x=Term,y=Count,fill=PValue)) + geom_bar(width=0.8,stat="identity") + coord_flip() + scale_fill_gradient(low = "red", high = "green") + scale_x_discrete(limits=show$Term[30:1])

