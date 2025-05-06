library(limma)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Methylation matrix (rows = CpG, columns = samples)
meth_matrix <- read.csv("bulk_samples.csv", row.names=1)
colnames(meth_matrix) <- sapply(strsplit(colnames(meth_matrix), "_"), `[`, 1)
meth_matrix

clusters = read.csv("sample_clusters.csv", row.names=1)
clusters

# Let's delete the outlier
dim(meth_matrix)
meth_matrix$GSM1052037 = NULL
dim(meth_matrix)

dim(clusters)
clusters = clusters[clusters$Sample!="GSM1052037",]
dim(clusters)

# Reorder the columns to fit the clusters order
meth_matrix = meth_matrix[clusters$Sample]

# Sample metadata with groups
groups <- factor(clusters$Cluster)
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)

# Fit model
fit <- lmFit(meth_matrix, design)
contrast.matrix <- makeContrasts(c("cluster_1-cluster_2", "cluster_1-cluster_3", "cluster_2-cluster_3"), levels=design)
contrast.matrix <- makeContrasts(
  cluster_1_vs_2 = cluster_1 - cluster_2,
  cluster_1_vs_3 = cluster_1 - cluster_3,
  cluster_2_vs_3 = cluster_2 - cluster_3,
  levels = design
)
fit2 <- contrasts.fit(fit, contrast.matrix)

# Get the differentially methylated CpGs
fit4 <- treat(fit2, lfc = 0, trend = TRUE, robust = TRUE)
treat1 <- topTreat(fit4, coef = 1, number = Inf, p.value = 0.05, adjust.method = "BH")
dim(treat1)
treat2 <- topTreat(fit4, coef = 2, number = Inf, p.value = 0.05, adjust.method = "BH")
dim(treat2)
treat3 <- topTreat(fit4, coef = 3, number = Inf, p.value = 0.05, adjust.method = "BH")
dim(treat3)

# Annotation
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(anno[, c("Name", "UCSC_RefGene_Name", "UCSC_RefGene_Group")])


