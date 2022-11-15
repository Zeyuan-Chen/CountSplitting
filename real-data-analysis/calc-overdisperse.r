library(Seurat)
library(countsplit)
library(caret)
library(pracma)
set.seed(1234)

args = commandArgs(trailingOnly=TRUE)

# choose from "hpc", "amy"
region = as.character(args[1])
ann.col = as.character(args[2])
L.version = as.character(args[3]) 
gene.version = as.character(args[4]) 

print(paste0("region: ", region))
print(paste0("ann.col: ", ann.col))
print(paste0("L.version: ", L.version))
print(paste0("gene.version: ", gene.version))

# region  = "hpc" #hpc, amy, pbmc
# ann.col = "seurat_clusters"
# L.version = "pcs" # bool pcs
# gene.version = "random" # random or hvg

base.dir = file.path("/u/project/halperin/johnsonc/Integrate_CCA/Notebook/PoissonThin")
res.dir  = file.path(base.dir, "res")
fig.dir  = file.path(base.dir, "fig")

if (region != "pbmc"){
    rna.obj = readRDS(paste0("../../Data/RNA/blk6-dba_", region, ".seu.rds"))
}else{
    print("load and preprocess pbmc 3k")
    data(pbmc.counts, package="countsplit")
    rownames(pbmc.counts) <- sapply(rownames(pbmc.counts), function(u) stringr::str_replace_all(u, "_","-"))
    rna.obj = CreateSeuratObject(counts = pbmc.counts, min.cells = 3, min.features = 200)

    rna.obj[["percent.mt"]] <- PercentageFeatureSet(rna.obj, pattern = "^MT-")

    rna.obj <- subset(rna.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
    rna.obj <- NormalizeData(rna.obj, normalization.method = "LogNormalize", scale.factor = 10000)
    rna.obj <- FindVariableFeatures(rna.obj, selection.method = "vst", nfeatures = 2000)
    rna.obj <- ScaleData(rna.obj, features = VariableFeatures(object = rna.obj))
    rna.obj <- RunPCA(rna.obj, features = VariableFeatures(object = rna.obj))
    rna.obj <- FindNeighbors(rna.obj, dims = 1:10)
    rna.obj <- FindClusters(rna.obj, resolution = 0.5)
    rna.obj <- RunUMAP(rna.obj, dims = 1:10)
}


if (L.version == "bool"){
    print(paste0("converting celltype annotation to one hot encoded version as L.hat"))
    L.df = rna.obj@meta.data[,"seurat_clusters", drop =F]
    dummy <- dummyVars(" ~ .", data=L.df)
    L.df <- data.frame(predict(dummy, newdata = L.df)) 
    #avoid singular
    print(paste0("Drop the first cluster column that has: ",sum(L.df[,1])/nrow(L.df) * 100, " percent of the cells"))
    L.df = L.df[,2:ncol(L.df),drop = F]
}else{
    print("directly using top 10 PCs as L.hat")
    L.df = rna.obj@reductions$ pca@ cell.embeddings[, 1:10]
}


p = 2000 
if(gene.version == "hvg"){
    rna.obj <- FindVariableFeatures(rna.obj, selection.method = "vst", nfeatures = p)
    fit.genes = VariableFeatures(rna.obj)
    print(paste0("working on top ", p, " most variable genes"))
}else{
    fit.genes = sample(rownames(rna.obj), p)
    print(paste0("working on random ", p, " genes"))
}

if(region != "pbmc"){
    print("randomely down sampling to 200 genes for mouse dataset (save computational time)")
    p = 200
    fit.genes = sample(fit.genes, p)
}

# cells by highly variable genes 
X <- t(as.matrix(rna.obj@ assays$ RNA @ counts [fit.genes, ]))
#library size
sfs = rna.obj@meta.data$nCount_RNA
sfs = sfs/mean(sfs)

f = as.formula(paste0("X ~ ",paste(colnames(L.df), collapse = ' + ')," + offset(log(sfs))"))
print(f)

overdisps.file  = file.path(res.dir, paste0("overdisps.",  region, ".", L.version, ".", gene.version, ".rds"))
pred_means.file = file.path(res.dir, paste0("pred_means.", region, ".", L.version, ".", gene.version, ".rds"))

if(T){
#if (!file.exists(overdisps.file) | !file.exists(pred_means.file)){
    overdisps <- data.frame(matrix(NA, nrow=p, ncol=2))
    names(overdisps) <- c("means",  "nb1")
    pred_means <- X[,1:p]

    df = as.data.frame(L.df)
    for (i in 1:p) {
        df$X = X[,i]
        print(i)
        overdisps[i,1] <- mean(X[,i])
        ### Based on this stack overflow thread, I am almost positive that theta is b, not
        ### 1/b. https://stats.stackexchange.com/questions/10419/what-is-theta-in-a-negative-binomial-regression-fitted-with-r
        
        try1 <- try(mod <- MASS::glm.nb(f, data = df))
        if (class(try1) != "try-error") {
            overdisps[i,2] <- mod$theta
            pred_means[,i] <- predict(mod, type="response")
        }
    }
    saveRDS(overdisps,  overdisps.file)
    saveRDS(pred_means, pred_means.file)
}


overdisps  = readRDS(overdisps.file)
pred_means = readRDS(pred_means.file)


over_dis = as.matrix(t(pred_means/repmat(t(as.matrix(overdisps$nb1)), nrow(pred_means),1)))

frac1  = sum(over_dis < 1, na.rm=TRUE)/sum(!is.na(over_dis))
frac10 = sum(over_dis < 10, na.rm=TRUE)/sum(!is.na(over_dis))

### visualization
### Use random indices to avoid crashing the computer. 
randomindices <- sample(1:length(as.numeric(over_dis)), 
                        size=min(5e6, length(as.numeric(over_dis))))

over_dis = as.numeric(t(over_dis))[randomindices]
over_dis = over_dis[!is.na(over_dis)]
over_dis = over_dis[!is.infinite(over_dis)]
over_dis = log10(over_dis)

g = ggplot(data=NULL, aes(x=over_dis, y=..density..))+
    geom_histogram(binwidth = 0.5) + 
    #scale_x_log10() + 
    #original value set to over_dis = 1, since we are in log10, this is set to 0. also plot 1 which is 10 in original space
    geom_vline(xintercept=0, col="red")+
    geom_vline(xintercept=1, col="cyan")+

    theme_classic()+
    xlab("Overdispersion (log10 scale)")+
    ylab("Density")+
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    ggtitle(paste0("Histogram of Estimated Overdispersion Values ", region, "\n", 
                   "frac < 1: ", round(frac1, 3), "    frac < 10: ", round(frac10, 3)))+ 
    theme(plot.title = element_text(hjust = 0.5)) 

ggsave(file.path(fig.dir, paste0("overdisps.", region, ".", L.version, ".", gene.version, ".png")), g, 
       width = 6, height = 6, dpi = 500)
