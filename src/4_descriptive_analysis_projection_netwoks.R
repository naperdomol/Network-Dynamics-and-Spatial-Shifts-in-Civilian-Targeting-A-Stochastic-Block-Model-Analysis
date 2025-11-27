# Descriptive analysis of projections characteristics
#-------------------------------------------------------------------------------
t <- c("78_81", "82_84", "85_87", "88_90", "91_93", "94_96", "97_99", "00_01", "02_04", "05_07")
periods <- c("1978-1981", "1982-1984", "1985-1987", "1988-1990", "1991-1993", "1994-1996", "1997-1999", "2000-2001", "2002-2004", "2005-2007")


# Municipalities projection networks

# Global structural properties
for (i in seq_along(t)) {
  graph_name <- paste0("municipalities_graph_", t[i])
  M_name <- paste0("M_", t[i])
  graph <- get(graph_name)
  B <- get(M_name)
  giant_component <- get(paste0("giant_comp_municipalities_", t[i]))
  A = B%*%t(B)
  A = ifelse(A != 0, 1, 0)
  diag(A) = 0
  A_S = A%*%A
  lower_tri_A = A_S
  lower_tri_A[!lower.tri(lower_tri_A)]  =  0
  
  vertices_municipalities[i] = vcount(graph)
  edges[i] = ecount(graph)
  degree_mean[i] = mean(degree(graph))
  density[i] = edge_density(graph)
  global_transitivity[i] = transitivity(graph, type="global")
  triangles[i] = sum(count_triangles(graph))/3
  triads[i] = sum(lower_tri_A)
  components[i]<-components(graph)$no
  giant_comp_size[i]<-max(sapply(X = decompose(graph), FUN = vcount))
  avg_path[i] = mean_distance(graph, directed = FALSE, unconnected = TRUE)
  diameter[i] = diameter(graph, directed = FALSE, unconnected = TRUE)
  degree_centralization[i] = centr_degree(graph = graph, loops = FALSE, normalized = TRUE)$centralization
  closeness_centralization[i] = centr_clo(giant_component, normalized = TRUE)$centralization
  betweenness_centralization[i] = centr_betw(graph, normalized = TRUE)$centralization
  cut_vertices[i] = length(articulation_points(graph))
  clique_number[i] = clique_num(graph)
  maxcliques_count[i] = count_max_cliques(giant_component)

  region <- as.factor(V(giant_component)$region)
  assort[i] <- assortativity_nominal(graph = giant_component, types = region, directed = FALSE)
  rm(A, A_S, lower_tri_A, region)
}

glob_char_municipalities <- data.frame(period = periods, edges = edges, vertices_municipalities = vertices_municipalities, degree_mean = degree_mean , 
                             density = density, global_transitivity=global_transitivity, triangles=triangles, triads=triads, 
                             components = components, giant_comp_size=giant_comp_size, avg_path=avg_path, diameter=diameter, 
                             degree_centralization=degree_centralization, closeness_centralization= closeness_centralization, betweenness_centralization=betweenness_centralization, 
                             cut_vertices=cut_vertices, clique_number=clique_number, maxcliques_count=maxcliques_count, assort=assort)
glob_char_municipalities$period <- factor(glob_char_municipalities$period, levels = periods)

# Local structural properties
loc_char_municipalities <- data.frame(period = character(), vertex = character(), 
                                     degree = numeric(), strength = numeric(), betweenness_centrality = numeric())
for (i in seq_along(t)) {
  graph_name <- paste0("municipalities_graph_", t[i])
  graph <- get(graph_name)
  degree <- degree(graph)
  strength <- strength(graph)
  betweenness_centrality <- betweenness(graph, normalized = TRUE)
  temp_df <- data.frame(period = periods[i], vertex = V(graph)$municipality_code, degree = degree, strength = strength, betweenness_centrality = betweenness_centrality)
  loc_char_municipalities <- rbind(loc_char_municipalities, temp_df)
  rm(temp_df)
}

loc_char_municipalities$period <- factor(loc_char_municipalities$period, levels = periods)

# Maximal cliques size distribution
maxcliques_municipalities <- data.frame(period = character(), clique_size = numeric())
for (i in seq_along(t)) {
  graph_name <- paste0("giant_comp_municipalities_", t[i])
  graph <- get(graph_name)
  max_cliques <- max_cliques(graph)
  cliques_size <- sapply(max_cliques, length)
  temp_df <- data.frame(period = periods[i], clique_size = cliques_size)
  maxcliques_municipalities <- rbind(maxcliques_municipalities, temp_df)
  rm(temp_df)
}
maxcliques_municipalities$period <- factor(maxcliques_municipalities$period, levels = periods)


#-------------------------------------------------------------------------------
#Plots
output_path <- "C:/Users/nperd/Documents/UN/Estadística/Redes/Proyecto/plots/"        

# Degree distribution
png(filename = paste0(output_path, "degree_dist_municipalities_proj.png"), , width = 650, height = 450)
ggplot(loc_char_municipalities, aes(x = period, y = degree)) +
  geom_boxplot(fill = "grey", color = "black") +
  labs(x = "", y = "Vertex degree") +
  theme_minimal(base_family = "serif") +
  theme(axis.text.x = element_text(size = 16, angle = 90, hjust = 1, family = "serif"),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))
dev.off()

degree_summary_municipalities <- loc_char_municipalities %>%
  group_by(period) %>%
  summarize(
    Min = min(degree, na.rm = TRUE),
    Q1 = quantile(degree, 0.25, na.rm = TRUE),
    Median = median(degree, na.rm = TRUE),
    Mean = mean(degree, na.rm = TRUE),
    Q3 = quantile(degree, 0.75, na.rm = TRUE),
    Max = max(degree, na.rm = TRUE),
    Sd = sd(degree, na.rm = TRUE),
    Skewness = skewness(degree, na.rm = TRUE)
  )

# Strength distribution
png(filename = paste0(output_path, "str_dist_municipalities_proj.png"), , width = 650, height = 450)
ggplot(loc_char_municipalities, aes(x = period, y = strength)) +
  geom_boxplot(fill = "grey", color = "black") +
  labs(x = "", y = "Vertex strength") +
  scale_y_continuous(breaks = seq(0, max(loc_char_municipalities$strength, na.rm = TRUE), by = 100)) + 
  theme_minimal(base_family = "serif") +
  theme(axis.text.x = element_text(size = 16, angle = 90, hjust = 1, family = "serif"),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))
dev.off()

summary_strength_municipalities <- loc_char_municipalities %>%
  group_by(period) %>%
  summarize(
    Min = min(strength, na.rm = TRUE),
    Q1 = quantile(strength, 0.25, na.rm = TRUE),
    Median = median(strength, na.rm = TRUE),
    Mean = mean(strength, na.rm = TRUE),
    Q3 = quantile(strength, 0.75, na.rm = TRUE),
    Max = max(strength, na.rm = TRUE),
    SD = sd(strength, na.rm = TRUE),
    Skewness = skewness(strength, na.rm = TRUE)
  )

# Density y global transitivity
png(filename = paste0(output_path, "dens_trans_municipalities_proj.png"), , width = 650, height = 450)
ggplot(glob_char_municipalities, aes(x = periods)) +
  geom_line(aes(y = density, linetype = "density", group = 1), size = 1.2, color = adjustcolor("black", 0.8)) +  
  geom_line(aes(y = global_transitivity * max(density, na.rm = TRUE) / max(global_transitivity, na.rm = TRUE), linetype = "Global transitivity", group = 1), size = 1.2, color = adjustcolor("grey", 0.8)) + 
  scale_y_continuous(name = "Density", sec.axis = sec_axis(~ . * max(glob_char_municipalities$global_transitivity, na.rm = TRUE) / max(glob_char_municipalities$density, na.rm = TRUE), 
                        name = "Global transitivity")) +
  scale_linetype_manual(values = c("density" = "solid", "Global transitivity" = "dashed")) +
  labs(x = "", y = "Density", linetype = "") +
  theme_minimal(base_family = "serif") +
  theme(axis.text.x = element_text(size = 16, angle = 90, hjust = 1, family = "serif"),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 16),  
    legend.key.width = unit(2, "cm"))
dev.off()

# Connected components
png(filename = paste0(output_path, "components_municipalities_proj.png"), , width = 650, height = 450)
ggplot(glob_char_municipalities, aes(x = period)) +
  geom_line(aes(y = components, linetype = "Connected components", group = 1), size = 1.2, color = adjustcolor("black", 0.9)) +
  geom_line(aes(y = giant_comp_size * max(components) / max(giant_comp_size), linetype = "Giant component size", group = 1), size = 1.2, color = adjustcolor("grey", 0.7)) +
  scale_y_continuous(name = "Connected components", breaks = seq(0, max(glob_char_municipalities$components), by = 5),
    sec.axis = sec_axis(~ . * max(glob_char_municipalities$giant_comp_size) / max(glob_char_municipalities$components), 
                        name = "Giant component size",
                        breaks = seq(0, max(glob_char_municipalities$giant_comp_size), by = 100))) +
  scale_linetype_manual(values = c("Connected components" = "solid", "Giant component size" = "dashed")) +
  labs(x = "", linetype = "") +
  theme_minimal(base_family = "serif") +
  theme(axis.text.x = element_text(size = 16, angle = 90, hjust = 1, family = "serif"),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.text = element_text(size = 16),  
        legend.key.width = unit(2, "cm"))
dev.off()


# Network level indices of centrality
png(filename = paste0(output_path, "centralization_indices_municipalities_proj.png"), , width = 650, height = 450)
ggplot(glob_char_municipalities, aes(x = period)) +
  geom_line(aes(y = degree_centralization, linetype = "Degree centr.", color = "Degree centr.", group = 1), size = 1.2) +
  geom_line(aes(y = closeness_centralization, linetype = "Closeness centr.", color = "Closeness centr.", group = 1), size = 1.2) +
  geom_line(aes(y = betweenness_centralization, linetype = "Betweenness centr.", color = "Betweenness centr.", group = 1), size = 1.2) +
  scale_y_continuous(name = "", breaks = seq(0, max(c(glob_char_municipalities$degree_centralization, glob_char_municipalities$closeness_centralization, glob_char_municipalities$betweenness_centralization)), by = 0.1)) +
  scale_linetype_manual(values = c("Degree centr." = "solid", "Closeness centr." = "dashed", "Betweenness centr." = "dotted")) +
  scale_color_manual(values = c("Degree centr." = "black", "Closeness centr." = "black", "Betweenness centr." = "grey")) +
  labs(x = "", linetype = "", color = "") +
  theme_minimal(base_family = "serif") +
  theme(axis.text.x = element_text(size = 16, angle = 90, hjust = 1, family = "serif"),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.text = element_text(size = 16),
        legend.key.width = unit(2, "cm"))
dev.off()


# Betweenness
png(filename = paste0(output_path, "betweenness_cent_municipalities_proj.png"), , width = 650, height = 450)
ggplot(loc_char_municipalities, aes(x = period, y = betweenness_centrality)) + 
  geom_boxplot(fill = "grey", color = "black") + 
  labs(x = "", y = "Betweenness centrality") + 
  theme_minimal(base_family = "serif") + 
  theme(axis.text.x = element_text(size = 16, angle = 90, hjust = 1, family = "serif"), 
    axis.text.y = element_text(size = 16), 
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16))
dev.off()

# Count of maximal cliques and clique number
png(filename = paste0(output_path, "cliques_municipalities_proj.png"), , width = 650, height = 450)
ggplot(glob_char_municipalities, aes(x = period)) +
  geom_line(aes(y = clique_number, linetype = "Clique number", group = 1), size = 1.2, color = adjustcolor("black", 0.8)) +
  geom_line(aes(y = maxcliques_count * (250 / 1500), linetype = "Number of maximal cliques", group = 1), size = 1.2, color = adjustcolor("grey", 0.7)) +  
  scale_y_continuous(name = "Clique number", limits = c(0, 250), sec.axis = sec_axis(~ . * (1500 / 250), name = "Number of maximal cliques")) +
  scale_linetype_manual(values = c("Clique number" = "solid", "Number of maximal cliques" = "dashed")) +  
  labs(x = "", linetype = "") +
  theme_minimal(base_family = "serif") +
  theme(axis.text.x = element_text(size = 16, angle = 90, hjust = 1, family = "serif"),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 16),  
    legend.key.width = unit(2, "cm"))
dev.off()

# Maximal cliques size distribution
cliques_summary_municipalities <- maxcliques_municipalities %>%
  group_by(period) %>%
  summarize(
    Min = min(clique_size, na.rm = TRUE),
    Q1 = quantile(clique_size, 0.25, na.rm = TRUE),
    Median = median(clique_size, na.rm = TRUE),
    Mean = mean(clique_size, na.rm = TRUE),
    Q3 = quantile(clique_size, 0.75, na.rm = TRUE),
    Max = max(clique_size, na.rm = TRUE),
    Sd = sd(clique_size, na.rm = TRUE),
    Skewness = skewness(clique_size, na.rm = TRUE))

png(filename = paste0(output_path, "cliques_dist_municipalities_proj.png"), , width = 650, height = 450)
ggplot(maxcliques_municipalities, aes(x = period, y = clique_size)) +
  geom_boxplot(fill = "grey", color = "black") +
  labs(x = "", y = "Maximal clique size") +
  theme_minimal(base_family = "serif") +
  theme(axis.text.x = element_text(size = 16, angle = 90, hjust = 1, family = "serif"),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16))
dev.off()


# Proportion of vertices in more than one, ten and thirty cliques
df <- data.frame()

for (i in seq_along(t)) {
  graph_name <- paste0("giant_comp_municipalities_", t[i])
  graph <- get(graph_name)
  
  max_cliques_list <- max_cliques(graph)
  cliques_df <- data.frame(vertex = unlist(max_cliques_list),
                           clique_id = rep(1:length(max_cliques_list),
                                           sapply(max_cliques_list, length)))
  vertex_cliques <- split(cliques_df$clique_id, cliques_df$vertex)
  
  
  largest_cliques_ids <- order(sapply(max_cliques_list, length), decreasing = TRUE)[1:5]
  vertex_largest_cliques <- lapply(vertex_cliques, function(cliques) 
    sum(cliques %in% largest_cliques_ids))
  
  temp <- data.frame(period = periods[i],
                     vertex_id = names(vertex_largest_cliques),
                     vertex = V(graph)$name,
                     cliques_count = sapply(vertex_cliques, length),
                     top5_cliques_count = unlist(vertex_largest_cliques))
  df <- rbind(df, temp)
  rm(temp, graph_name, max_cliques_list, cliques_df, vertex_cliques,vertex_largest_cliques, largest_cliques_ids)
}

proportions_df <- data.frame(period = unique(df$period),
                             p_1plus = numeric(length(unique(df$period))),
                             p_10plus = numeric(length(unique(df$period))),
                             p_30plus = numeric(length(unique(df$period))))

for (p in unique(df$period)) {
  period_data <- df[df$period == p, ]
  vertices <- nrow(period_data)

  cliques_1plus <- sum(period_data$cliques_count > 1)
  p_1plus <- round(cliques_1plus / vertices,2)
  cliques_10plus <- sum(period_data$cliques_count > 10)
  p_10plus <- round(cliques_10plus / vertices,2)
  cliques_30plus <- sum(period_data$cliques_count > 30)
  p_30plus <- round(cliques_30plus / vertices,2)
  
  proportions_df[proportions_df$period == p, ] <- c(p, p_1plus, p_10plus, p_30plus)
}

print(proportions_df)