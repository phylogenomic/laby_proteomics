library(pacman)


pacman::p_load(tidyverse, data.table,plotly,here,grid,patchwork,ggrepel,
               RColorBrewer,ggVennDiagram,factoextra,FactoMineR,furrr,ggplotify)# CRAN

pacman::p_load(Biostrings,DEP,sva,ComplexHeatmap,SummarizedExperiment)

both <- 1:2551

mmetsp_det_M <- paste0("M",1:743)
mmetsp_det_M <- c(mmetsp_det_M,both)
jgi_det_J <- paste0("J",1:489)
jgi_det_J <- c(jgi_det_J,both)

y <- list("MMETSP" = mmetsp_det_M,
          "JGI" = jgi_det_J)
# 4D Venn diagram
venn_detected <- ggVennDiagram(y, lwd = 0.8, lty = 1) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  ggtitle("Detected proteins")


both <- 1:302

mmetsp_det_M <- paste0("M",1:228)
mmetsp_det_M <- c(mmetsp_det_M,both)
jgi_det_J <- paste0("J",1:93)
jgi_det_J <- c(jgi_det_J,both)

y <- list("MMETSP" = mmetsp_det_M,
          "JGI" = jgi_det_J)
# 4D Venn diagram
venn_significant <- ggVennDiagram(y, lwd = 0.8, lty = 1) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  ggtitle("Significant proteins")

venn_plot <- wrap_plots(venn_detected, venn_significant) +
  plot_annotation(tag_levels = "A", tag_suffix = ".") &
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 15, hjust = 0, vjust = 0))


scale_x <- 4
ggsave(
  filename = paste0("proteomics/img_IndependentSets/indset.fig2_t0_vennplot.pdf"),
  plot = venn_plot,
  width = 3.5*scale_x,
  height = 2.2*scale_x)
