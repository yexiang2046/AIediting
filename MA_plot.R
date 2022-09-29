working_dir <- "/Users/rajsuba/Desktop/temp"
setwd(working_dir)
input <- read.csv("BC3.latent.lytic.csv", header=TRUE, sep = ",")
col_names <- colnames(input)
repSC <- 9 # Column ID where replicates begin.

wt_rep_no <- 2 # Number of wt replicates.
mut_rep_no <- 2 # Number of mut replicates.

padj_cutoff <- 0.05 # padj cutoff.
l2fc_cutoff <- 0 # log2FoldChange cutoff.
meanC_cutoff_up <- 1.5 # Mean counts ratio upper cutoff.
meanC_cutoff_down <- 0.66 # Mean counts ratio lower cutoff.

dotploter <- function(in_data, x, y, z, y_int, width, y_lim, x_lim, x_lab, y_lab){
  x_var <- enquo(x)
  y_var <- enquo(y)
  z_var <- enquo(z)
  map <- ggplot(data=in_data, aes(x=!!x_var, y=!!y_var)) + 
    geom_point(aes(color=!!z_var), alpha=0.4, size=1.5) +
    scale_color_manual(values = c("gray55", "#F8766D", "#00BFC4")) + 
    geom_hline(aes(yintercept = y_int), colour = "blue", size = width, linetype="dashed") +
    ylim(y_lim) + xlim(x_lim) +
    xlab(x_lab) + 
    ylab(y_lab) + 
    theme(axis.title.x = element_text(face = "bold", size = 14),
          axis.text.x = element_text(size = 10)) +
    theme(axis.title.y = element_text(face = "bold", size = 14),
          axis.text.y = element_text(size = 10)) +
    theme(legend.title = element_blank())
}

#### MA Plot -------------

input$Threshold <- "NS"
input$Threshold[(input$log2FoldChange >= 0 & input$padj < padj_cutoff)] <- "Up"
input$Threshold[(input$log2FoldChange < 0 & input$padj < padj_cutoff)] <- "Down"
#head(input)

input$Threshold <- factor(input$Threshold, levels=c("NS","Up","Down"))

ma_data <- input[,c("baseMean", "log2FoldChange", "Threshold")]
#head(ma_data)

map <- dotploter(ma_data, x=baseMean, y=log2FoldChange, Threshold, 0, 0.5,
                 c(-7.5, 7.5), c(0, 8e+04), "Mean expression", "Log2 Fold Change")
map
pvalue_cutoff <- 0.05
input$Threshold <- "NS"
input$Threshold[(input$log2FoldChange >= 0.5 & input$pvalue < pvalue_cutoff)] <- "Up"
input$Threshold[(input$log2FoldChange < -0.5 & input$pvalue < pvalue_cutoff)] <- "Down"
