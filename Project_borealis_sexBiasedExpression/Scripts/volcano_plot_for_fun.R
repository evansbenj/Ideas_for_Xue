library(tidyverse)


result_table <- read_csv("input_file")

result_volcano <- (result_table
                   %>% mutate(logP = (-log(pvalue)))
                    %>% mutate (case_cols = case_when(
                         logP >= 2 & abs(log2FoldChange) >= 5  ~ "Large FC",
                         logP >= 2 & abs(log2FoldChange) >= 3 ~ "Medium Significant",
                         logP >= 2 & abs(log2FoldChange) > 0 ~ "Significant",
                         logP < 2 ~ "NS"
                       ))
)

volcano_cols <- c("Large FC" = "red1", "Significant" = "steelblue3", "Medium Significant" ="black", "NS" = "grey")


volcano_plot <- (ggplot(result_volcano) + 
  geom_point(
    aes(x=log2FoldChange, y =logP, col=case_cols), 
    size=3, alpha=0.6
  ) +
  xlab("log2FC") + 
  ylab("-log (P-value)") +
  guides(alpha=FALSE, size = FALSE) + #don't include alpha and size infor on legend
  scale_colour_manual( #manipulate the color of the points and legends
    values = volcano_cols,
    name = "Transcript", #name of the legend
  )+
  theme_classic() +
  theme(
    axis.title.y = element_text(size=20, margin=margin(0,20,0,0)),
    axis.title.x = element_text(size=20, margin=margin(20,0,0,0)),
    axis.text = element_text(size = 18),
    legend.title = element_text(colour="black", size=13, face="bold"),
    legend.text = element_text(colour="black", size = 13),
    legend.justification = "top",
    legend.position = c(0.85, 1),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black")
  )
)

volcano_plot
