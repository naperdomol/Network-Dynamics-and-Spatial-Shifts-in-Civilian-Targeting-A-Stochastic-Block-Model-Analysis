output_path <- ""   

time_periods <- list("78_81" = "1978-1981", 
                     "82_84" = "1982-1984", 
                     "85_87" = "1985-1987", 
                     "88_90" = "1988-1990", 
                     "91_93" = "1991-1993",
                     "94_96" = "1994-1996",
                     "97_99" = "1997-1999",
                     "00_01" = "2000-2001",
                     "02_04" = "2002-2004",
                     "05_07" = "2005-2007")

# Iterate over each time period and create three plots for each period
for (suffix in names(time_periods)) { 
  period <- time_periods[[suffix]]
  
  # File paths for saving the plots
  bipartite_plot_path <- paste0(output_path,"/bipartite_graph_", suffix, ".png")
  armed_units_plot_path <- paste0(output_path,"/armed_units_graph_", suffix, ".png")
  municipalities_plot_path <- paste0(output_path,"/municipalities_graph_", suffix, ".png")
  
  # Plot bipartite graph
  graph <- get(paste0("graph_", suffix))
  png(bipartite_plot_path, width = 450, height = 450)
  par(mar = c(0.1, 0.1, 0.1, 0.1), mfrow = c(1, 1))
  set.seed(326)
  print(
    plot(graph, 
         layout = layout_with_fr(graph, niter = 1000),
         vertex.label = NA, 
         vertex.color = V(graph)$color,
         vertex.shape = V(graph)$shape, 
         vertex.size = log(igraph::degree(graph)+1)*2, 
         vertex.frame.color = ifelse(V(graph)$type == F, 
                                     ifelse(V(graph)$armed_group_type == "guerrilla group", "#AA35B3", 
                                            ifelse(V(graph)$armed_group_type == "paramilitary group", "royalblue", 
                                                   ifelse(V(graph)$armed_group_type == "organized crime", "#FF4C00", "#3AB795"))),
                                     adjustcolor("gray", 0.8)),
         edge.color = adjustcolor("gray", 0.8),
         edge.lty = 1)
  )
  dev.off()
  
  # Plot armed units' network
  armed_units_graph <- get(paste0("armed_units_graph_", suffix))
  png(armed_units_plot_path, width = 450, height = 450)
  par(mar = c(0.1, 0.1, 0.1, 0.1), mfrow = c(1, 1))
  set.seed(498)
  print(plot(armed_units_graph,
             layout = layout_with_fr(armed_units_graph), 
             vertex.label = V(armed_units_graph)$armed_unit_label, 
             vertex.label.cex = 0.65,
             #vertex.label = NA,
             vertex.color = V(armed_units_graph)$color, 
             vertex.size = 2*(log(strength(armed_units_graph)+1)) + 2, 
             vertex.frame.color = NA,
             edge.color = adjustcolor("lightgray", 0.9),
             edge.width = (E(armed_units_graph)$weight)*0.6,
             edge.lty = 1))
  dev.off()
  
  # Plot municipalities' network
  municipalities_graph <- get(paste0("municipalities_graph_", suffix))
  png(municipalities_plot_path, width = 450, height = 450)
  unique_reg <- sort(unique(V(municipalities_graph)$region))
  color_palette <- viridis(length(unique_reg), option = "D")
  color_palette <- adjustcolor(color_palette, 0.7)
  reg_colors <- setNames(color_palette, unique_reg)
  V(municipalities_graph)$color <- reg_colors[as.character(V(municipalities_graph)$region  )]
  
  par(mar = c(0.1, 0.1, 0.1, 0.1), mfrow = c(1, 1))
  set.seed(713)
  print(
    plot(municipalities_graph,
         layout = layout_with_fr(municipalities_graph),
         vertex.label = NA,  
         vertex.color = V(municipalities_graph)$color,
         vertex.size = 2*log(strength(municipalities_graph) + 1), 
         edge.color = adjustcolor("gray", 0.8),
         edge.lty = 1)
  )
  legend("bottomleft", legend = names(reg_colors), col = reg_colors, pch = 15, pt.cex = 2, cex = 1.2, bty = "n")
  dev.off()
}