library(DESeq2)
countdata <- read.table("/Users/rajsuba/Desktop/temp/latent.lytic.repeat.copy.txt", row.names=1, header = TRUE, sep = '\t')
countdata <- countdata[ ,6:ncol(countdata)]
countdata <- as.matrix(countdata)
coldata <- read.csv('/Users/rajsuba/Desktop/temp/labels.csv', row.names=1)
coldata <- coldata[,c("condition","sample")]
colnames(countdata)
colnames(countdata) <- gsub("X.home.suba.editing.aligned_with_rawreads.", "", colnames(countdata))
colnames(countdata) <- gsub("Aligned.sortedByCoord.out.bam", "", colnames(countdata))
colnames(countdata)
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
keep <- rowSums(counts(dds)) >= 1
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds)
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
write.csv(resdata, file="/Users/rajsuba/Desktop/temp/latent.lytic.repeat.csv")

plot_pca <- function(
  pca_data, 
  ylim=20, 
  xlim=100, 
  point_size=3,
  expand_circ=0.08, 
  s_shape=0.8, 
  col_alpha=0.25, 
  edge_width=1,
  color_by=NULL)
{
  
  # calculate variance
  percentVar <- round(100 * attr(pca_data, "percentVar"), digits = 2)
  
  # pca plot object
  pcaplot <- ggplot2::ggplot(pca_data, 
                             ggplot2::aes(x = PC1, y = PC2, color = color_by)) +
    ggplot2::geom_point(size=point_size) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ylim(-ylim, ylim) +
    xlim(-xlim, xlim) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.x = element_text(face = "bold", size = 16),
                   axis.text.x = element_text(size = 16)) +
    ggplot2::theme(axis.title.y = element_text(face = "bold", size = 16),
                   axis.text.y = element_text(size = 16)) +
    ggplot2::theme(legend.position="right", legend.text = element_text(size = 14),
                   legend.title = element_text(face = "bold", size = 16)) +
    ggalt::geom_encircle(ggplot2::aes(group=color_by, fill=color_by), 
                         s_shape=s_shape, 
                         expand=expand_circ, 
                         alpha=col_alpha,
                         size=edge_width)
  
  return(pcaplot)
}
