#Descriptive analysis of bipartite graph characteristics

#-------------------------------------------------------------------------------
# Periods
t <- c("78_81", "82_84", "85_87", "88_90", "91_93", "94_96", "97_99", "00_01", "02_04", "05_07")
periods <- c("1978-1981", "1982-1984", "1985-1987", "1988-1990", "1991-1993", "1994-1996", "1997-1999", "2000-2001", "2002-2004", "2005-2007")

municipalities_vertices <- rep(0, length(t))
armed_units_vertices <- rep(0, length(t))
armed_units_per_municipality <- rep(0, length(t))
edges <- rep(0, length(t))
density <- rep(0, length(t))
connected_components_count <- rep(0, length(t))
giant_component_size <- rep(0, length(t))
municipalities_count <- rep(0, length(t))
units_count <- rep(0, length(t))

for (i in seq_along(t)) {
  graph_name <- paste0("graph_", t[i])
  graph_i <- get(graph_name)
  
  municipalities_vertices[i] <- sum(V(graph_i)$type == TRUE)  # Municipalities vertices
  armed_units_vertices[i] <- sum(V(graph_i)$type == FALSE)  # Armed units vertices
  armed_units_per_municipality[i] <- armed_units_vertices[i] / municipalities_vertices[i]  # Armed units per municipality
  municipalities_count[i] <- length(V(graph_i)[type == TRUE & degree(graph_i) > 2])   # Count municipalities with presence of three or more armed units
  units_count[i] <- length(V(graph_i)[type == FALSE & degree(graph_i) > 2])    # Count armed units operating in three or more municipalities
  edges[i] <- ecount(graph_i)  # Edges
  density[i] <- edges[i] / (municipalities_vertices[i] * armed_units_vertices[i])  # Density, bipartite graphs
  connected_components_count[i] <- components(graph_i)$no  # Number of connected components
  giant_component_size[i] <- max(sapply(X = decompose(graph_i), FUN = vcount))  # Size of the largest component (giant component)
  
}

df_charact <- data.frame(
  period = periods,
  edges = edges,
  armed_units_vertices = armed_units_vertices,
  municipalities_vertices = municipalities_vertices,
  municipalities_count = municipalities_count,
  units_count = units_count,
  armed_units_per_municipality = armed_units_per_municipality,
  density = density,
  connected_components_count = connected_components_count,
  giant_component_size = giant_component_size
)

df_charact$period <- factor(df_charact$period, levels = periods)

# Degree distribution
armed_units_ds <- data.frame(period = character(), vertex = character(), degree = numeric())

for (i in seq_along(t)) {
  graph <- paste0("graph_", t[i])
  graph_i <- get(graph)
  degree_dist <- degree(graph_i,v=V(graph_i)[type==F])
  temp <- data.frame(period = periods[i], vertex = V(graph_i)$name[V(graph_i)[type==F]] , degree = degree_dist)
  armed_units_ds <- rbind(armed_units_ds, temp)
  rm(temp)
}
armed_units_ds$period <- factor(armed_units_ds$period, levels = periods)

municipalities_data <- data.frame(period = character(), degree = numeric())

for (i in seq_along(t)) {
  graph <- paste0("graph_", t[i])
  graph_i <- get(graph)
  degree_dist <- degree(graph_i,v=V(graph_i)[type==T])
  temp <- data.frame(period = periods[i], vertex = V(graph_i)$name[V(graph_i)[type==T]] , degree = degree_dist)
  municipalities_data <- rbind(municipalities_data, temp)
  rm(temp)
}
municipalities_data$period <- factor(municipalities_data$period, levels = periods)

#-------------------------------------------------------------------------------
#Plots
output_path <- "C:/Users/nperd/Documents/UN/Estadística/Redes/Proyecto/plots/"     

# Armed substructures (units) vertices
png(filename = paste0(output_path, "armed_units_vertices.png"), width = 650, height = 450)
ggplot(df_charact, aes(x = period)) +
  geom_line(aes(y = armed_units_vertices, color = "Armed units vertices"), size = 1.2, group = 1) +
  geom_line(aes(y = units_count, color = "Armed units in 3 or more municipalities"), size = 1.2, linetype = "dashed", group = 1) +
  scale_y_continuous(name = "Number of armed units",
    breaks = seq(0, max(df_charact$armed_units_vertices), by = 20)) +
  scale_color_manual(name = "", values = c("Armed units vertices" = "black", "Armed units in 3 or more municipalities" = "grey40")) +
  labs(x = "", color = "") +
  theme_minimal(base_family = "serif") +  # Serif
  theme(axis.text.x = element_text(size = 16, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.position = "bottom",
    legend.text = element_text(size = 16),  
    legend.key.width = unit(2, "cm"))
dev.off()

armed_units_summary <- armed_units_ds %>%
  group_by(period) %>%
  summarise(
    mean_degree = mean(degree),
    sd_degree = sd(degree),
    skewness = skewness(degree),
    q3_degree = quantile(degree, 0.75),
    max_degree = max(degree)
  )

# Municipalities vertices
png(filename = paste0(output_path, "municipalities_vertices.png"), width = 650, height = 450)
ggplot(df_charact, aes(x = period)) +
  geom_line(aes(y = municipalities_vertices, color = "Municipalities vertices"), size = 1.2, group = 1) +
  geom_line(aes(y = municipalities_count, color = "Municipalities with 3 or more armed substructures"), size = 1.2, linetype = "dashed", group = 1) +
  scale_y_continuous(
    name = "Number of municipalities",
    breaks = seq(0, max(df_charact$municipalities_vertices), by = 100)) +
  scale_color_manual(name = "", values = c("Municipalities vertices" = "black", "Municipalities with 3 or more armed substructures" = "grey40")) +
  labs( x = "", color = "") +
  theme_minimal(base_family = "serif") +  # Serif
  theme(axis.text.x = element_text(size = 16, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.position = "bottom",
    legend.text = element_text(size = 16),  
    legend.key.width = unit(2, "cm"))
dev.off()

municipalities_summary <- municipalities_data %>%
  group_by(period) %>%
  summarise(
    mean_degree = mean(degree),
    sd_degree = sd(degree),
    skewness = skewness(degree),
    q3_degree = quantile(degree, 0.75),
    max_degree = max(degree)
  )

# Edges and density
png(filename = paste0(output_path, "edges_density_bipartite.png"), width = 650, height = 450)
ggplot(df_charact, aes(x = periods)) +
  geom_line(aes(y = edges, color = "Edges"), size = 1.2, group = 1) +
  scale_y_continuous(
    name = "Edges",
    limits = c(0, max(df_charact$edges, na.rm = TRUE) * 1.1),
    breaks = seq(0, max(df_charact$edges, na.rm = TRUE), by = 300),
    sec.axis = sec_axis(~ . * max(df_charact$density, na.rm = TRUE) / max(df_charact$edges, na.rm = TRUE),
                        name = "Density",
                        labels = scales::number_format(accuracy = 0.01))) +
  geom_line(aes(y = density * max(edges) / max(density), color = "Density"), size = 1.2, group = 1, linetype = "dashed") +
  scale_color_manual(name = "", values = c("Edges" = "black", "Density" = "grey40")) +
  theme_minimal(base_family = "serif") +  # Serif
  labs(x = "") +
  theme(plot.title = element_text(hjust = 0.5, family = "serif", size = 14),
    axis.text.x = element_text(size = 16, angle = 90, hjust = 1, family = "serif"),
    axis.text.y = element_text(size = 16, family = "serif"),
    axis.title.x = element_text(size = 16, family = "serif"),
    axis.title.y = element_text(size = 16, family = "serif"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 16),  
    legend.key.width = unit(2, "cm"))
dev.off()

# Connected components and giant component size
png(filename = paste0(output_path, "components_bipartite.png"), width = 650, height = 450)
ggplot(df_charact, aes(x = period)) +
  geom_line(aes(y = connected_components_count, color = "Connected components"), size = 1.2, group = 1) +
  geom_line(aes(y = giant_component_size * max(connected_components_count) / max(giant_component_size), color = "Giant component size"), size = 1.2, linetype = "dashed", group = 1) +
  scale_y_continuous(name = "Connected components",
    breaks = seq(0, max(df_charact$connected_components_count), by = 5),
    sec.axis = sec_axis(~ . * max(df_charact$giant_component_size) / max(df_charact$connected_components_count), 
                        name = "Number of vertices",
                        breaks = seq(0, max(df_charact$giant_component_size), by = 100))) +
  scale_color_manual(name = "", values = c("Connected components" = "black", "Giant component size" = "grey40")) +
  labs(color = "", x = "") +
  theme_minimal(base_family = "serif") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),  
    plot.subtitle = element_text(size = 16, color = "gray50", hjust = 0.5),
    axis.text.x = element_text(size = 16, angle = 90, hjust = 1, family = "serif"),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.position = "bottom",
    legend.text = element_text(size = 16),  
    legend.key.width = unit(2, "cm"))
dev.off()

# Degree distribution
png(filename = paste0(output_path, "degree_dist_armed_units.png"), width = 650, height = 450)
ggplot(armed_units_ds, aes(x = period, y = degree)) +
  geom_boxplot(fill = adjustcolor("royalblue", 0.5), color = "black") +
  labs(x = "", y = "Vertex degree") +
  theme_minimal(base_family = "serif") +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray50"),
    axis.text.x = element_text(size = 16, angle = 90, hjust = 1, family = "serif"),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16))
dev.off()

png(filename = paste0(output_path, "degree_dist_municipalities.png"), width = 650, height = 450)
ggplot(municipalities_data, aes(x = period, y = degree)) +
  geom_boxplot(fill = "grey40", color = "black") +
  labs(x = "", y = "Vertex degree") +
  theme_minimal(base_family = "serif") +
  theme(axis.text.x = element_text(size = 16, angle = 90, hjust = 1, family = "serif"),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16))
dev.off()


#Top 5 vertices: degree
top_5_nodes <- list()
for (i in t) {
  graph_name <- paste0("graph_", i)
  graph <- get(graph_name)
  
  top5_degree_municipalities<- V(graph)[type == FALSE][order(degree(graph, v = V(graph)[type == FALSE]), decreasing = TRUE)[1:5]]$armed_unit_label
  top5_degree_armed_units <- V(graph)[type == TRUE][order(degree(graph, v = V(graph)[type == TRUE]), decreasing = TRUE)[1:5]]$name
  
  top_5_nodes[[i]] <- list(
    municipalities_deg = top5_degree_municipalities,
    armed_units_deg = top5_degree_armed_units
  )
}
rm(list=c("top5_degree_municipalities","top5_degree_armed_units"))

# Create a named vector to match t codes to period labels
names(periods) <- t  

top_q3_to_max_nodes <- list()

for (i in t) {
  graph <- get(paste0("graph_", i))
  period_label <- periods[i]
  
  deg_mun <- degree(graph, v = V(graph)[type == TRUE])
  deg_arm <- degree(graph, v = V(graph)[type == FALSE])
  
  mun_logical <- deg_mun > municipalities_summary$q3_degree[municipalities_summary$period == period_label] &
    deg_mun <= municipalities_summary$max_degree[municipalities_summary$period == period_label]
  arm_logical <- deg_arm > armed_units_summary$q3_degree[armed_units_summary$period == period_label] &
    deg_arm <= armed_units_summary$max_degree[armed_units_summary$period == period_label]
  
  top_mun_vertices <- V(graph)[type == TRUE][mun_logical]
  top_mun_degrees <- deg_mun[mun_logical]
  top_mun_vertices <- top_mun_vertices[order(top_mun_degrees, decreasing = TRUE)]
  
  top_arm_vertices <- V(graph)[type == FALSE][arm_logical]
  top_arm_degrees <- deg_arm[arm_logical]
  top_arm_vertices <- top_arm_vertices[order(top_arm_degrees, decreasing = TRUE)]

    top_q3_to_max_nodes[[i]] <- list(
    municipalities_deg = top_mun_vertices$name,
    municipalities_subregion = top_mun_vertices$subregion,
    armed_units_deg = top_arm_vertices$armed_unit_label
  )
}