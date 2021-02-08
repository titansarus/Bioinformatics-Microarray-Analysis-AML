library(GEOquery)
library(limma)
library(Biobase)
library(pheatmap)
library(reshape2)
library(plyr)
library(ggplot2)
library(stringr)

seriesName <- "GSE48558"
platformName <- "GPL6244"

gset <-
  getGEO(
    seriesName,
    GSEMatrix = TRUE,
    AnnotGPL = TRUE ,
    destdir = "Data"
  )
if (length(gset) > 1) {
  idx <- grep(platformName, attr(gset, "names"))
} else {
  idx <- 1
}
gset <- gset[[idx]]

gset <-
  gset[, which(gset$source_name_ch1 == "AML Patient" |
                 gset$`phenotype:ch1` == "Normal")]


func <- function(x) {
  if (gset$source_name_ch1[x] == "AML Patient") {
    return("Test")
  } else {
    spll <- strsplit2(gset$source_name_ch1[x] , "\\+")[1, 1]
    return(paste0("Normal_" , spll))
  }
}


gr <- sapply(1:length(gset$`phenotype:ch1`) , func)


expr <- exprs(gset)
print(paste0("Max Expr: " , max(expr)))
print(paste0("Min Expr: " , min(expr)))
# expr <- log2(1 + expr)
# exprs(gset)<- expr


pdf("Results/boxplot.pdf" , width = 32)
boxplot(expr)
dev.off()


#expr <- normalizeQuantiles(expr)
#exprs(gset) <- expr

#pdf("Results/boxplot_afterNormalize.pdf" , width = 32)
#boxplot(expr)
#dev.off()



pca <- prcomp(expr)
pdf("Results/pca.pdf" , width = 10, height = 10)
plot(pca)
plot(pca$x[, 1:2])
dev.off()


expr.scaled <- t(scale(t(expr) , scale = FALSE))



pca <- prcomp(expr.scaled)
pdf("Results/pca_scaled.pdf" ,
    width = 10,
    height = 10)
plot(pca)
plot(pca$x[, 1:2])
dev.off()


pcar <- data.frame(pca$rotation[, 1:3] , group = gr)
pdf("Results/PCA_Samples.pdf" ,
    width = 15 ,
    height = 15)
ggplot(pcar , aes(PC1 , PC2 , color = group , size = 4)) + geom_point() + theme_bw()
dev.off()



pdf("Results/corheat.pdf" , width = 30, height = 30)
pheatmap(cor(expr) , labels_row =  gr , labels_col =  gr )
dev.off()



##Analyze

gr2 <- gr
gr2[which((pcar$PC2 > 0.11 & pcar$group == "Test"))] <-
  "AML_Near_CD34"


gr2 <- factor(gr2)
gset$description <- gr2

onehot <- model.matrix(~ description + 0, gset)
colnames(onehot) <- levels(gr2)

colnames(onehot) <-
  c(sapply(colnames(onehot) , function(x)
    str_replace(x, " ", "_")))

fit <- lmFit(gset , onehot)

contrast <-
  makeContrasts(AML_Near_CD34 - Normal_CD34 , levels = onehot)

fit2 <- contrasts.fit(fit , contrast)
fit2 <- eBayes(fit2 , 0.01)
tT <- topTable(fit2,
               adjust = "fdr",
               sort.by = "B" , number = Inf)

tT2 <- subset(tT , select = c("Gene.symbol" , "Gene.title", "adj.P.Val"  , "logFC"))

tT2 <- subset(tT , select = c("Gene.symbol" , "Gene.title", "adj.P.Val"  , "logFC"))
write.table(
  tT2 ,
  file = "Results/aml_cd34.txt" ,
  sep = "\t" ,
  row.names = FALSE ,
  quote = FALSE
)

aml.up <- subset(tT2, logFC > 1 & adj.P.Val < 0.05)
aml.up.genes <-unique( as.character(
  strsplit2( (aml.up$Gene.symbol),"///")))
write.table(
  aml.up.genes,
  "Results/amd_cd34_up_genes.txt" ,
  quote = FALSE,
  row.names = FALSE ,
  col.names = FALSE
)



aml.down <- subset(tT2, logFC < -1 & adj.P.Val < 0.05)
aml.down.genes <- unique(as.character(
  strsplit2( (aml.down$Gene.symbol),"///")))
write.table(
  aml.down.genes,
  "Results/amd_cd34_down_genes.txt" ,
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
