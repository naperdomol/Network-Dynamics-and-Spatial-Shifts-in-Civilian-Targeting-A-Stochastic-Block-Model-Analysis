# Community Detection: Fast-greddy clustering
output_path <- ""   

#-------------------------------------------------------------------------------
# Armed units projection networks
periods<- c( "1982-1984", "1985-1987", "1988-1990", "1991-1993", 
             "1994-1996", "1997-1999", "2000-2001", "2002-2004", "2005-2007")
t <- c(#"78_81", 
  "82_84", "85_87", "88_90", "91_93", "94_96", "97_99", "00_01", "02_04", "05_07")

clusters <-  list()
mod_values <-  list()

for (period in t) {
  graph <- get(paste0("giant_comp_armed_unit_", period))
  
  #Fast greedy clustering
  cfg <-  cluster_fast_greedy(graph, weights = E(graph)$weight)
  V(graph)$membership_cfg <- cfg$membership   # Membership attribute
  assign(paste0("giant_comp_armed_unit_", period), graph)
  
  # Number of clusters and modularity
  clusters[period] <- length(unique(cfg$membership))
  mod_values[period] <- modularity(cfg)
  
  # Average strength within and between communities
  graph <- get(paste0("giant_comp_armed_unit_", period))
  membership <- V(graph)$membership_cfg
  communities <- sort(unique(membership))
  
  avg_strength_intra <- list()
  avg_strength_inter <- list()
  community_count <- list()
  for (i in communities) {
    avg_strength <- (sum(E(graph)[which(membership == i) %--% which(membership == i)]$weight)/(sum(membership == i)*(sum(membership == i)-1)))*2
    avg_strength_intra[[paste0("Community_", i)]] <- avg_strength
  }
  
  community_pairs <- combn(communities, 2)
  for (i in 1:ncol(community_pairs)) {
    comm1 <- community_pairs[1, i]
    comm2 <- community_pairs[2, i]
    edges_between <- E(graph)[which(membership == comm1) %--% which(membership == comm2)]
    avg_strength <- sum(E(graph)[edges_between]$weight)/(sum(membership == comm1)*sum(membership == comm2))
    avg_strength <- ifelse(is.nan(avg_strength), 0, avg_strength)
    avg_strength_inter[[paste0("Community_", comm1, "_", comm2)]] <- avg_strength
  }
  
  cfg_strg_mat <- matrix(0, nrow = length(communities), ncol = length(communities), dimnames = list(communities, communities))
  for (i in communities) {
    cfg_strg_mat[i, i] <- avg_strength_intra[[paste0("Community_", i)]]
  }
  
  for (i in 1:ncol(community_pairs)) {
    comm1 <- community_pairs[1, i]
    comm2 <- community_pairs[2, i]
    cfg_strg_mat[comm2, comm1] <- avg_strength_inter[[paste0("Community_", comm1, "_", comm2)]]
    cfg_strg_mat[comm1, comm2] <- avg_strength_inter[[paste0("Community_", comm1, "_", comm2)]]
  }
  assign(paste0("cfg_avg_strg_mat", period), cfg_strg_mat)    # Average edge strength within and between communities
}

df_s<- data.frame(period = periods, cluster=unlist(clusters), modularity=unlist(mod_values))

# Municipalities projection networks

t <- c(#"78_81", 
  "82_84", "85_87", "88_90", "91_93", "94_96", "97_99", "00_01", "02_04", "05_07")

clusters <-  list()
mod_values <-  list()


for (i in seq_along(t)) {
  period<-t[i]
  graph <- get(paste0("giant_comp_municipalities_", period))
  
  # Fast greedy clustering
  cfg <-  cluster_fast_greedy(graph, weights = E(graph)$weight)
  V(graph)$membership_cfg <- cfg$membership # Membership attribute
  assign(paste0("giant_comp_armed_unit_", period), graph)
  
  # Number of clusters and modularity
  clusters[period] <- length(unique(cfg$membership))
  mod_values[period] <- modularity(cfg)
  
  # Average strength within and between communities
  graph <- get(paste0("giant_comp_municipalities_", period))
  membership <- V(graph)$membership_cfg
  communities <- sort(unique(membership))
  
  avg_strength_intra <- list()
  avg_strength_inter <- list()
  community_count <- list()
  for (i in communities) {
    avg_strength <- (sum(E(graph)[which(membership == i) %--% which(membership == i)]$weight)/(sum(membership == i)*(sum(membership == i)-1)))*2
    avg_strength_intra[[paste0("Community_", i)]] <- avg_strength
  }
  
  community_pairs <- combn(communities, 2)
  for (i in 1:ncol(community_pairs)) {
    comm1 <- community_pairs[1, i]
    comm2 <- community_pairs[2, i]
    edges_between <- E(graph)[which(membership == comm1) %--% which(membership == comm2)]
    avg_strength <- sum(E(graph)[edges_between]$weight)/(sum(membership == comm1)*sum(membership == comm2))
    avg_strength <- ifelse(is.nan(avg_strength), 0, avg_strength)
    avg_strength_inter[[paste0("Community_", comm1, "_", comm2)]] <- avg_strength
  }
  
  cfg_strg_mat <- matrix(0, nrow = length(communities), ncol = length(communities), dimnames = list(communities, communities))
  for (i in communities) {
    cfg_strg_mat[i, i] <- avg_strength_intra[[paste0("Community_", i)]]
  }
  
  for (i in 1:ncol(community_pairs)) {
    comm1 <- community_pairs[1, i]
    comm2 <- community_pairs[2, i]
    cfg_strg_mat[comm2, comm1] <- avg_strength_inter[[paste0("Community_", comm1, "_", comm2)]]
    cfg_strg_mat[comm1, comm2] <- avg_strength_inter[[paste0("Community_", comm1, "_", comm2)]]
  }
  assign(paste0("cfg_m_avg_strg_mat", period), cfg_strg_mat)   # Average edge strength within and between clusters
}
df_m<- data.frame(period = periods, cluster=unlist(clusters), modularity=unlist(mod_values))


#-------------------------------------------------------------------------------
# Communities characterization

# Community size (fast-greedy algorithm)
df <- left_join(df_m, df_s, by = "period", suffix = c("_municipalities", "_armed_units"))

png(filename = paste0(output_path, "fgc_communities.png"), width = 650, height = 450)
ggplot(df, aes(x = period)) +
  geom_line(aes(y = cluster_municipalities, linetype = "Municipalities", group = 1), size = 1.2, color = adjustcolor("gold", 0.7)) +
  geom_line(aes(y = cluster_armed_units, linetype = "Armed structures", group = 1), size = 1.2, color = adjustcolor("royalblue", 0.6)) +
  scale_y_continuous(name = "Number of clusters",
                     limits = c(0, max(c(df$cluster_municipalities, df$cluster_armed_units)) + 3),  
                     breaks = seq(0, max(c(df$cluster_municipalities, df$cluster_armed_units)) + 3, by = 3)) +
  scale_linetype_manual(values = c("Municipalities" = "solid", "Armed structures" = "dashed")) +  
  labs(x = "", linetype = "") +
  theme_minimal(base_family = "serif") +
  theme(
    axis.text.x = element_text(size = 16, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 16),  
    legend.key.width = unit(2, "cm"))
dev.off()


# Modularity
png(filename = paste0(output_path, "modularity.png"), width = 650, height = 450)
ggplot(df, aes(x = period)) +
  geom_line(aes(y = modularity_municipalities, linetype = "Municipalities", group = 1), size = 1.2, color = adjustcolor("gold", 0.7)) +
  geom_line(aes(y = modularity_armed_units, linetype = "Armed structures", group = 1), size = 1.2, color = adjustcolor("royalblue", 0.6)) +
  scale_y_continuous(name = "Modularity",
                     limits = c(0, max(c(df$modularity_municipalities, df$modularity_armed_units)) + 0.5),  
                     breaks = seq(0, max(c(df$modularity_municipalities, df$modularity_armed_units)) + 0.5, by = 0.2)) +
  scale_linetype_manual(values = c("Municipalities" = "solid", "Armed structures" = "dashed")) +  
  labs(x = "", linetype = "") +
  theme_minimal(base_family = "serif") +
  theme(
    axis.text.x = element_text(size = 16, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 16),  
    legend.key.width = unit(2, "cm"))
dev.off()

rm(df, df_m, df_s)

# Intra- and inter-cluster average edge strength 

# Armed strutuctures networks
periods <- c("82_84", "85_87", "88_90", "91_93", "94_96", "97_99", "00_01", "02_04", "05_07")
data_list <- list() 
for (period in periods) {
  matrix_name <- paste0("cfg_avg_strg_mat", period)
  mat <- get(matrix_name)
  
  diagonal <- data.frame(
    avg_edge_strg = diag(mat),
    intra_cluster = 1,
    period = period
  )
  
  triangular <- which(upper.tri(mat, diag = FALSE), arr.ind = TRUE)
  upper_tri <- data.frame(
    avg_edge_strg = mat[triangular],
    intra_cluster = 0,
    period = period
  )
  
  combined <- rbind(diagonal, upper_tri)
  data_list[[period]] <- combined
  rm(mat, matrix_name, diagonal, triangular, upper_tri, combined)
}
data_combined <- do.call(rbind, data_list)

o <- c("82_84", "85_87","88_90","91_93","94_96","97_99","00_01","02_04","05_07")
n <- c("1982-1984", "1985-1987","1988-1990","1991-1993","1994-1996","1997-1999","2000-2001","2002-2004","2005-2007")
period_map <- setNames(n, o)
data_combined$period <- period_map[data_combined$period]

png(filename = paste0(output_path, "cluster_au_avg_strg.png"), width = 650, height = 450)
ggplot(data_combined, aes(x = period, y = avg_edge_strg, shape = factor(intra_cluster), color = factor(intra_cluster))) +
  geom_point(size = 4, alpha = 0.6) + 
  scale_shape_manual(name = " ",
                     values = c("1" = 16, "0" = 8), 
                     labels = c("Inter-cluster", "Intra-cluster")) +
  scale_color_manual(
    name = " ",
    values = c("1" = "forestgreen", "0" = "#D883C6"), 
    labels = c("Inter-cluster", "Intra-cluster")) +
  labs(x = "", y = "Avg. Edge Strength") + 
  theme_minimal(base_family = "serif") +
  theme(axis.text.x = element_text(size = 16, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 14, color = "gray50", hjust = 0.5))
dev.off()

# Municipalities networks

periods <- c("82_84", "85_87", "88_90", "91_93", "94_96", "97_99", "00_01", "02_04", "05_07")
data_list <- list() 
for (period in periods) {
  matrix_name <- paste0("cfg_m_avg_strg_mat", period)
  mat <- get(matrix_name)
  
  diagonal <- data.frame(
    avg_edge_strg = diag(mat),
    intra_cluster = 1,
    period = period
  )
  
  triangular <- which(upper.tri(mat, diag = FALSE), arr.ind = TRUE)
  upper_tri <- data.frame(
    avg_edge_strg = mat[triangular],
    intra_cluster = 0,
    period = period
  )
  
  combined <- rbind(diagonal, upper_tri)
  data_list[[period]] <- combined
  rm(mat, matrix_name, diagonal, triangular, upper_tri, combined)
}
data_combined <- do.call(rbind, data_list)

o <- c("82_84", "85_87","88_90","91_93","94_96","97_99","00_01","02_04","05_07")
n <- c("1982-1984", "1985-1987","1988-1990","1991-1993","1994-1996","1997-1999","2000-2001","2002-2004","2005-2007")
period_map <- setNames(n, o)
data_combined$period <- period_map[data_combined$period]

png(filename = paste0(output_path, "cluster_m_avg_strg.png"), width = 650, height = 450)
ggplot(data_combined, aes(x = period, y = avg_edge_strg, shape = factor(intra_cluster), color = factor(intra_cluster))) +
  geom_point(size = 4, alpha = 0.7) + 
  scale_shape_manual(name = " ",
                     values = c("1" = 16, "0" = 8), 
                     labels = c("Inter-cluster", "Intra-cluster")) +
  scale_color_manual( name = " ",
                      values = c("1" = "orange", "0" = "steelblue"), 
                      labels = c("Inter-cluster", "Intra-cluster")) +
  labs(x = "", y = "Avg. Edge Strength") + 
  theme_minimal(base_family = "serif") +
  theme(axis.text.x = element_text(size = 16, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 14, color = "gray50", hjust = 0.5))
dev.off()

# Cliques and communities
max_cliques_list <- max_cliques(giant_comp_municipalities_94_96)
cliques_df <- data.frame(vertex = unlist(max_cliques_list), 
                         clique_id = rep(1:length(max_cliques_list), 
                                         sapply(max_cliques_list, length)))
vertex_cliques <- split(cliques_df$clique_id, cliques_df$vertex)

clique_sizes <- sapply(max_cliques_list, length)
largest_cliques_ids <- order(clique_sizes, decreasing = TRUE)[1:5]
vertex_largest_cliques <- lapply(vertex_cliques, function(cliques) 
  sum(cliques %in% largest_cliques_ids))

df <- data.frame(vertex = names(vertex_largest_cliques),
                 membership = V(giant_comp_municipalities_94_96)$membership,
                 num_cliques = sapply(vertex_cliques, length),
                 num_largest_cliques = unlist(vertex_largest_cliques))

vertices_in_singletons<-V(giant_comp_municipalities_94_96)[c3$membership %in% which(sizes(cluster_walktrap(giant_comp_municipalities_94_96)) == 1)]
vertex_cliques <- lapply(vertices_in_singletons, function(vertex) {
  which(sapply(max_cliques_list, function(clique) vertex %in% clique))
})


# Assortativity

# Armed structures
ggplot(dfsubestrc_caract, aes(x = periodo)) +
  geom_point(aes(y = rand, shape = "Índice Rand"), size = 4, color = adjustcolor("#D883C6", 0.8)) +
  geom_point(aes(y = (assort + 1) / 2, shape = "Asortatividad"), size = 4, color = adjustcolor("royalblue", 0.7)) +  
  scale_y_continuous(
    name = "Índice Rand",
    limits = c(0, 1),  # Ajuste del eje para el Índice Rand
    sec.axis = sec_axis(~ . * 2 - 1, name = "Asortatividad", breaks = seq(-1, 1, by = 0.5))  # Ajuste para Asortatividad
  ) +
  scale_shape_manual(values = c("Índice Rand" = 18, "Asortatividad" = 18)) +  
  labs(x = "",
       shape = "") +
  theme_minimal(base_family = "serif") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),  
    plot.subtitle = element_text(size = 10, color = "gray50", hjust = 0.5),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 9),
    legend.position = "bottom"
  )

# Municipalities
ggplot(dfmpios_caract, aes(x = periodo)) +
  geom_point(aes(y = rand, shape = "Índice Rand"), size = 4, color = adjustcolor("steelblue", 0.8)) +
  geom_point(aes(y = (assort + 1) / 2, shape = "Asortatividad"), size = 4, color = adjustcolor("orange", 0.7)) +  
  scale_y_continuous(
    name = "Índice Rand",
    limits = c(0, 1),  # Ajuste del eje para el Índice Rand
    sec.axis = sec_axis(~ . * 2 - 1, name = "Asortatividad", breaks = seq(-1, 1, by = 0.5))  # Ajuste para Asortatividad
  ) +
  scale_shape_manual(values = c("Índice Rand" = 18, "Asortatividad" = 18)) +  
  labs(x = "",
       shape = "") +
  theme_minimal(base_family = "serif") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),  
    plot.subtitle = element_text(size = 12, color = "gray50", hjust = 0.5),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "bottom"
  )

