library(ggplot2)
csv_data<-read.csv("/Users/rajsuba/Desktop/temp/BC3_total.csv")
x <- csv_data$log2fold
y <- csv_data$Lytic_Latent
hist(y,ylim=c(0,10000))

group <- ifelse(x > 0.5 & y > 0.2, "A",
                               ifelse(x > 0.5 & y < -0.2, "B",
                                     ifelse(x < -0.5 & y > 0.2, "C",
                                            ifelse(x < -0.5 & y < -0.2, "D", "E"))))

df <- data.frame(x = x, y = y, group = group)
custom_colors <- c("indianred2", "orchid2", "peachpuff2", "steelblue2", "gray90")
ggplot(data=df, aes(x = x, y = y)) +
  geom_point(aes(colour = group), size = 0.05, show.legend = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = custom_colors)+ xlim(-3,3)+ylim(-1,1)+ ggtitle("BC3_Editing vs expression") +
  xlab("log2fold") + ylab("Editing change")



csv_data<-read.csv("/Users/rajsuba/Desktop/temp/BC3_pvalue.csv")
library(GeneOverlap)
data(GeneOverlap)
go.obj <- newGeneOverlap(csv_data$downreg, 
                         csv_data$increased, 
                         3034)
print(go.obj)  # not tested yet.
go.obj <- testGeneOverlap(go.obj)
print(go.obj)

emp.data <- data.frame(
  a = c (38091,884), 
  b = c (979,46))

# Print the data frame.			
print(emp.data)
fisher.test(emp.data)

library(ggplot2)
csv_data<-read.csv("/Users/rajsuba/Desktop/temp/BC3_decreased_editing.csv")
x <- csv_data$log2fold
y <- csv_data$p_value

group <- ifelse(x > 0.5 & y < 0.05, "A",
                       ifelse(x < -0.5 & y < 0.05, "B","C"))

df <- data.frame(x = x, y = -log2(y), group = group)
custom_colors <- c("indianred2", "steelblue2", "gray28")
ggplot(data=df, aes(x = x, y = y)) +
  geom_point(aes(colour = group), size = 0.05, show.legend = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = custom_colors)+ ggtitle("BC3-Transcripts with decreased editing in lytic") +
  xlab("log2fold") + ylab("log10(pvalue)")+ xlim(-2,2)+ylim(0,20)
