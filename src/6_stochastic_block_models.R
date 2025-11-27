#-------------------------------------------------------------------------------
# Stochastic Block Models
output_path <- "C:/Users/nperd/Documents/UN/Estadística/Redes/Proyecto/plots/"   


#-------------------------------------------------------------------------------
# Municipalities projection networks

periods <- list("78_81" = municipalities_graph_78_81,
                "82_84" = municipalities_graph_82_84,
                "85_87" = municipalities_graph_85_87,
                "88_90" = municipalities_graph_88_90,
                "91_93" = municipalities_graph_91_93,
                "94_96" = municipalities_graph_94_96,
                "97_99" = municipalities_graph_97_99,
                "00_01" = municipalities_graph_00_01,
                "02_04" = municipalities_graph_02_04,
                "05_07" = municipalities_graph_05_07)

set.seed(42)

for (period in names(periods)) {
  # Convert adjacency matrix to numeric matrix
  adj_matrix <- as.matrix(as_adjacency_matrix(periods[[period]], attr = "weight"))   
  assign(paste0("adj_matrix_m_", period), adj_matrix)
  

    # Fit SBM to the period's adjacency matrix
  model_sbm <- blockmodels::BM_poisson(membership_type = "SBM_sym", adj = adj_matrix, verbosity = 0, plotting = "")
  model_sbm$estimate()
  assign(paste0("model_sbm_m_", period), model_sbm)
  # Extract ICLs and find optimal number of partitions
  #model_sbm<-get(paste0("model_sbm_m_", period))
  ICLs <- model_sbm$ICL
  assign(paste0("ICLs_m_", period), ICLs)
  Q <- which.max(ICLs)
  assign(paste0("Q_m_", period), Q)
  
  # Save ICL values and plot results
  png(filename = paste0(output_path, "ICLs_m_", period,".png"), width = 400, height = 300)
  par(mfrow = c(1,1), mar = c(2.75,2.75,1.5,0.5), mgp = c(1.7,0.7,0))
  plot(ICLs, xlab = "Q", ylab = "ICL", type = "b", pch = 16)
  lines(x = c(Q, Q), y = c(min(ICLs), max(ICLs)), col = "red", lty = 2)
  dev.off()
  
  # Assign vertices to communities and calculate alpha
  Z <- model_sbm$memberships[[Q]]$Z
  assign(paste0("Z_m_", period), Z)
  labs <- apply(X = Z, MARGIN = 1, FUN = which.max)
  assign(paste0("labs_m_", period), labs)
  alpha <- table(labs) / vcount(periods[[period]])
  assign(paste0("alpha_m_", period), round(alpha, 3))
  graph <- periods[[period]]
  V(graph)$membership_sbm <- labs
  assign(paste0("municipalities_graph_", period), graph)
  alpha <- table(labs) / vcount(graph)
  assign(paste0("alpha_m_", period), round(alpha, 3))
  
  # Communities connection probabilities π_ij (π_ij)
  Lambda <- model_sbm$model_parameters[[Q]]$lambda
  assign(paste0("Lambda_m_", period), Lambda)
  
  # Network representing communities connection probabilities (π_ij) 
  g_c <- graph_from_adjacency_matrix(round(Lambda, 3), mode = "undirected", weighted = TRUE, diag = FALSE)
  V(g_c)$Lambda_intra_q <- diag(Lambda)
  V(g_c)$vert_prop <- alpha # vertices proportion
  assign(paste0("g_comm_m_", period), g_c)
}
rm(adj_matrix, model_sbm, Q, ICLs, Z, labs, alpha, Lambda, g_c)

#Switching labels

#1978-1981
V(g_comm_m_78_81)$name<-c(1,2,3,4)

#1982-1984
V(municipalities_graph_82_84)$membership_sbm <- dplyr::recode(
  V(municipalities_graph_82_84)$membership_sbm,
  `1` = 1, `2` = 2, `3` = 6, `4` = 4, `5` = 5,
  `6` = 7,  `7` = 3
)

V(g_comm_m_82_84)$name<-c(1,2,6,4,5,7,3)


#1985-1987
V(municipalities_graph_85_87)$membership_sbm <- dplyr::recode(
  V(municipalities_graph_85_87)$membership_sbm,
  `1` = 1, `2` = 4, `3` = 5, `4` = 2, `5` = 9,
  `6` = 7,  `7` = 6,  `8` = 8,  `9` = 3
)

V(g_comm_m_85_87)$name<-c(1,4,5,2,9,7,6,8,3)

#1988-1990
V(municipalities_graph_88_90)$membership_sbm <- dplyr::recode(
  V(municipalities_graph_88_90)$membership_sbm,
  `1` = 1, `2` = 2, `3` = 3, `4` = 11, `5` = 4,
  `6` = 5,  `7` = 6,  `8` = 10,  `9` = 8,  `10` = 7,
  `11` = 9
)

V(g_comm_m_88_90)$name<-c(1, 2, 3, 11, 4, 5, 6, 10, 8, 7, 9)

#1991-1993
V(municipalities_graph_91_93)$membership_sbm <- dplyr::recode(
  V(municipalities_graph_91_93)$membership_sbm,
  `1` = 8, `2` = 15, `3` = 6, `4` = 1, `5` = 2,
  `6` = 14,  `7` = 9,  `8` = 5,  `9` = 12,  `10` = 10,
  `11` = 7, `12` = 3, `13` = 13, `14` = 4, `15` = 11
)

V(g_comm_m_91_93)$name<-c(8,15,6,1,2,14,9,5,12,10,7,3,13,4,11)


#1994-1996
V(municipalities_graph_94_96)$membership_sbm <- dplyr::recode(
  V(municipalities_graph_94_96)$membership_sbm,
  `1` = 2, `2` = 4, `3` = 5, `4` = 3, `5` = 11,
  `6` = 6,  `7` = 8,  `8` = 9,  `9` = 10, `10` = 12,
  `11` = 7, `12` = 1
)

V(g_comm_m_94_96)$name<-c(2,4,5,3,11,6,8,9,10,12,7,1)

#1997-1999
V(municipalities_graph_97_99)$membership_sbm <- dplyr::recode(
  V(municipalities_graph_97_99)$membership_sbm,
  `1` = 4, `2` = 13, `3` = 20, `4` = 11, `5` = 10,
  `6` = 7,  `7` = 5,  `8` = 14,  `9` = 6,  `10` = 22,
  `11` = 16, `12` = 18, `13` = 3, `14` = 15, `15` = 19, 
  `16` = 17, `17` = 9, `18` = 12, `19` = 8, `20` = 1, 
  `21` = 21, `22` = 2
)

V(g_comm_m_97_99)$name<-c(4,13,20,11,10,7,5,14,6,22,16,18,3,15,19,17,9,12,8,1,21,2)

#2000-2001
V(municipalities_graph_00_01)$membership_sbm <- dplyr::recode(
  V(municipalities_graph_00_01)$membership_sbm,
  `1` = 8, `2` = 27, `3` = 18, `4` = 4, `5` = 10,
  `6` = 16,  `7` = 17,  `8` = 15,  `9` = 6,  `10` = 5,
  `11` = 2, `12` = 14, `13` = 26, `14` = 3, `15` = 21, 
  `16` = 25, `17` = 11, `18` = 24, `19` = 13, `20` = 9, 
  `21` = 23, `22` = 7, `23` = 20, `24` = 1, `25` = 12, 
  `26` = 19, `27` = 22
)

V(g_comm_m_00_01)$name<-c(8,27,18,4,10,16,17,15,6,5,2,14,26,3,21,25,11,24,
                          13,9,23,7,20,1,12,19,22)

#2002-2004
V(municipalities_graph_02_04)$membership_sbm <- dplyr::recode(
  V(municipalities_graph_02_04)$membership_sbm,
  `1` = 9, `2` = 27, `3` = 5, `4` = 18, `5` = 21,
  `6` = 4,  `7` = 28,  `8` = 17,  `9` = 6,  `10` = 14,
  `11` = 26, `12` = 10, `13` = 13, `14` = 25, `15` = 12, 
  `16` = 7, `17` = 29, `18` = 19, `19` = 2, `20` = 24, 
  `21` = 22, `22` = 1, `23` = 23, `24` = 3, `25` = 8, 
  `26` = 11, `27` = 15, `28` = 16, `29` = 20
)

V(g_comm_m_02_04)$name<-c(9,27,5,18,21,4,28,17,6,14,26,10,13,25,12,7,29,19,
                          2,24,22,1,23,3,8,11,15,16,20)

#2005-2007
V(municipalities_graph_05_07)$membership_sbm <- dplyr::recode(
  V(municipalities_graph_05_07)$membership_sbm,
  `1` = 4, `2` = 1, `3` = 15, `4` = 14, `5` = 9,
  `6` = 10,  `7` = 13,  `8` = 17,  `9` = 2,  `10` = 16,
  `11` = 11, `12` = 6, `13` = 18, `14` = 7, `15` = 12, 
  `16` = 8, `17` = 5, `18` = 3
)

V(g_comm_m_05_07)$name<-c(4,1,15,14,9,10,13,17,2,16,11,6,18,7,12,8,5,3)

#-------------------------------------------------------------------------------
# Sankey diagram

periods <- list(
  "1978-1981" = municipalities_graph_78_81,
  "1982-1984" = municipalities_graph_82_84,
  "1985-1987" = municipalities_graph_85_87,
  "1988-1990" = municipalities_graph_88_90,
  "1991-1993" = municipalities_graph_91_93,
  "1994-1996" = municipalities_graph_94_96,
  "1997-1999" = municipalities_graph_97_99,
  "2000-2001" = municipalities_graph_00_01,
  "2002-2004" = municipalities_graph_02_04,
  "2005-2007" = municipalities_graph_05_07
)

# 1: Extract memberships from each period
membership_df <- data.frame()

for (i in seq_along(periods)) {
  print(i)
  period_name <- names(periods)[i]
  graph <- periods[[i]]
  membership_df <- bind_rows(membership_df, 
                             data.frame(
                               municipality_code = V(graph)$name,  # municipality_code identifier
                               community = as.factor(V(graph)$membership_sbm),  # Community ID
                               period = period_name))
}

# 2: Track transitions between periods
full_grid <- expand.grid(municipality_code = unique(membership_df$municipality_code), 
                         period = unique(membership_df$period))

membership_data <- full_grid %>%
  left_join(membership_df, by = c("municipality_code", "period")) %>%
  arrange(municipality_code, period) %>%
  group_by(municipality_code) %>%
  mutate(
    next_community = lead(community),
    next_period = lead(period),
    prev_community = lag(community),
    entry = if_else((is.na(prev_community) & !is.na(community) & period!="1982-1984" ), "Entry", NA_character_),
    Inactive = if_else((is.na(next_community) & !is.na(community) & period!="2005-2007"), "Inactive", NA_character_)
  ) %>%
  mutate(prev_community = if_else(is.na(prev_community), entry, prev_community),
         next_community = if_else(is.na(next_community), Inactive, next_community),
         entry = coalesce(entry, if_else(!is.na(community) & period == "1982-1984", "Entry", NA_character_)),
         Inactive = coalesce(Inactive, if_else(!is.na(community) & period == "2005-2007", "Inactive", NA_character_))) %>%
  mutate(entry = as.integer(!is.na(entry)),
         Inactive = as.integer(!is.na(Inactive))) %>%
  ungroup()
rm(full_grid)

membership_data <- membership_data %>%
  group_by(municipality_code) %>%
  mutate(
    first_entry_period = min(period[entry == 1], na.rm = TRUE),
    last_inactive_period = max(period[Inactive == 1], na.rm = TRUE),
    community = coalesce(community, if_else(period >= first_entry_period & period <= last_inactive_period, "Inactive", NA_character_)),
    next_community = coalesce(next_community, if_else(period >= first_entry_period & period < last_inactive_period, "Inactive", NA_character_))
  ) %>%
  mutate(
    community = case_when(
      community != "Inactive" & nchar(as.character(community)) == 1 ~ paste0("0", community),
      TRUE ~ as.character(community) # Keep "Inactive" unchanged
    ),
    next_community = case_when(
      next_community != "Inactive" & nchar(as.character(next_community)) == 1 ~ paste0("0", next_community),
      TRUE ~ as.character(next_community)
    )
  ) %>%
  mutate(
    community = factor(community, 
                       levels = c(sort(unique(community[community != "Inactive"])), "Inactive")),
    next_community = factor(next_community, 
                            levels = c(sort(unique(next_community[next_community != "Inactive"])), "Inactive"))
  ) %>% 
  filter(!is.na(community) & period!="2005-2007") %>%
  ungroup()

# 3: Create transition dataset
transitions <- membership_data %>%
  group_by(period, community, next_period, next_community) %>%
  summarise(weight = n(), .groups = "drop") %>%
  mutate(source = paste(community, period, sep = "_"),
         target = paste(next_community, next_period, sep = "_")) %>%
  select(source, target, weight)

# 4: Define a set with the nodes
nodes <- data.frame(name = unique(c(transitions$source, transitions$target)))
nodes <- nodes %>%
  separate(name, into = c("community", "period"), sep = "_", remove = FALSE) %>%
  arrange(period,community)
nodes$id <- seq_along(nodes$name)

transitions$source_id <- match(transitions$source, nodes$name) - 1
transitions$target_id <- match(transitions$target, nodes$name) - 1
transitions<-as.data.frame(transitions)

# 5: plot
library(htmlwidgets)
library(htmltools)

period_labels <- '
<div style="display: grid; grid-template-columns: repeat(10, 1fr); text-align: center;
            margin-bottom: -20px; font-size: 14px; font-weight: bold;">
  <div>1978-1981</div> <div>1982-1984</div> <div>1985-1987</div> <div>1988-1990</div> 
  <div>1991-1993</div> <div>1994-1996</div> <div>1997-1999</div> <div>2000-2001</div> 
  <div>2002-2004</div> <div>2005-2007</div>
</div>
'

p <- sankeyNetwork(
  Links  = transitions,
  Nodes  = nodes,
  Source = "source_id",
  Target = "target_id",
  Value  = "weight",
  NodeID = "community",
  units  = "TWh",
  fontSize  = 12,
  nodeWidth = 10,
  iterations = 0
)

p <- htmlwidgets::prependContent(p, htmltools::HTML(period_labels))

# Use TRUE explicitly, not T
htmlwidgets::saveWidget(p, file = "sankey_diagram.html", selfcontained = TRUE)
browseURL("sankey_diagram.html")
webshot("sankey_diagram.html", file = paste0(output_path,"sankey_diagram.png"), vwidth = 1200, vheight = 800, zoom = 2)

###############################################################################
# Characterizing communities

summary(Z_m_85_87[cbind(1:vcount(giant_comp_municipalities_85_87), labs_m_85_87)])
summary(Z_m_88_90[cbind(1:vcount(giant_comp_municipalities_88_90), labs_m_88_90)])
summary(Z91_93[cbind(1:vcount(giant_comp_municipalities_91_93), labs91_93)])
summary(Z94_96[cbind(1:vcount(giant_comp_municipalities_94_96), labs94_96)])
summary(Z97_99[cbind(1:vcount(giant_comp_municipalities_97_99), labs97_99)])
summary(Z00_01[cbind(1:vcount(giant_comp_municipalities_00_01), labs00_01)])
summary(Z02_04[cbind(1:vcount(giant_comp_municipalities_02_04), labs02_04)])
summary(Z05_07[cbind(1:vcount(giant_comp_municipalities_05_07), labs05_07)])

# Distribution of vertex community assignments
round(alpha_m_78_81[order(alpha_m_78_81, decreasing = T)], 3)
round(alpha_m_82_84[order(alpha_m_82_84, decreasing = T)], 3)
round(alpha_m_85_87[order(alpha_m_85_87, decreasing = T)], 3)
round(alpha_m_88_90[order(alpha_m_88_90, decreasing = T)], 3)
round(alpha_m_91_93[order(alpha_m_91_93, decreasing = T)], 3)
round(alpha_m_94_96[order(alpha_m_94_96, decreasing = T)], 3)
round(alpha_m_97_99[order(alpha_m_97_99, decreasing = T)], 3)
round(alpha_m_00_01[order(alpha_m_00_01, decreasing = T)], 3)
round(alpha_m_02_04[order(alpha_m_02_04, decreasing = T)], 3)
round(alpha_m_05_07[order(alpha_m_05_07, decreasing = T)], 3)

# Plots
t <- c("78_81", "82_84", "85_87", "88_90", "91_93", "94_96", "97_99", "00_01", "02_04", "05_07")
seeds <- c(235, 54, 145, 658, 645, 654, 68, 654, 972,456)

for (i in seq_along(t)) {
  period <-t[10]
  seed <- seeds[10]
  g_comm <- get(paste0("g_comm_m_", period))
  #V(g_comm)$label <- paste0("Q", seq_along(V(g_comm)))
  V(g_comm)$label <- paste0("Q", V(g_comm)$name)
  
  color_function <- colorRampPalette(c("grey93", "grey20"))
  min_value <- min(c(V(g_comm)$Lambda_intra_q, E(g_comm)$weight))
  max_value <- max(c(V(g_comm)$Lambda_intra_q, E(g_comm)$weight))
  vertex_colors <- color_function(100)[as.numeric(cut(V(g_comm)$Lambda_intra_q, breaks = seq(min_value, max_value, length.out = 101)))]
  edge_colors <- color_function(100)[as.numeric(cut(E(g_comm)$weight, breaks = seq(min_value, max_value, length.out = 101)))]
  V(g_comm)$color <- vertex_colors
  E(g_comm)$color <- edge_colors
  
  mean_lambda_i <- mean(V(g_comm)$Lambda_intra_q)
  cv_lambda_i <- sd(V(g_comm)$Lambda_intra_q) / mean_lambda_i
  mean_lambda_ij <- mean(E(g_comm)$weight)
  cv_lambda_ij <- sd(E(g_comm)$weight) / mean_lambda_ij
  
  #png(filename = paste0(output_path, "/sbm_communities_mpal_", period, ".png"),width = 465, height = 545)
      set.seed(89)
      par(mar = c(0.1, 0.1, 1.3 , 2.5))
      plot(g_comm, layout = layout_with_fr(g_comm), 
           vertex.size = ((V(g_comm)$vert_prop)* 100 + 1),
          # vertex.label = V(g_comm)$label,
           vertex.label.size = 0.8,
           vertex.label.degree = (3/4)*pi,
           vertex.label.dist = 1,
           vertex.label.color = "black",
           edge.width = log(E(g_comm)$weight)+1,
           vertex.frame.color = "black",
           vertex.color = V(g_comm)$color,
           edge.color = E(g_comm)$color,
           edge.lty = 1)
      legend("topleft",
             legend = c(
               substitute(bar(lambda[i]) == val, list(val = round(mean_lambda_i, 2))),
               substitute(CV(lambda[i]) == val, list(val = round(cv_lambda_i, 2))),
               substitute(bar(lambda[ij]) == val, list(val = round(mean_lambda_ij, 2))),
               substitute(CV(lambda[ij]) == val, list(val = round(cv_lambda_ij, 2)))
             ),
             bty = "n", cex = 0.8, text.width = 2)
      image.plot(legend.only = TRUE, zlim = c(min_value, max_value), col = color_function(100),
                 legend.args = list(text = "Avg. Intensity (λ)", side = 4, line = 2, cex = 0.8),
                 axis.args = list(cex.axis = 0.8))
  dev.off()
}


#-------------------------------------------------------------------------------
# Mapping communities

custom_colors <- c(
  "1" = "#3584BB", "2" = "#B6CCEA", "3" ="#ff7518", "4" ="#fee999",
  "5" = "#34a532", "6" = "#a0dc92", "7" = "#DA3C3D", "8" = "#7E5E9C",
  "9" = "#c2a8db", "10" = "#8d4f42", "11" = "#c9a39a", "12" = "#df72ba",
  "13" = "#f3cbe6", "14" = "#828282", "15" = "#c9c9c9", "16" = "#c4c42d",
  "17" = "#e5e09a", "18" = "#4ec4cf", "19" = "#9ed8e2", "20" = "#4b6e99",
  "21" = "#d84315", "22" = "#ffab40", "23" = "#5e35b1", "24" = "#ab47bc",
  "25" = "#009688", "26" = "#4cb362", "27" = "#9bbb11", "28" = "#fdd835",
  "29" = "#ff7043"
)

# List of periods
periods <- c("1978-1981", "1982-1984", "1985-1987", "1988-1990", "1991-1993",
             "1994-1996", "1997-1999", "2000-2001", "2002-2004")

setwd("C:/Users/nperd/Documents/UN/Estadística/Redes/Proyecto/data/shapefiles")
mpioshp <- st_read("MGN_MPIO_POLITICO.shp",quiet=TRUE)   
deptoshp <- st_read("MGN_DPTO_POLITICO.shp",quiet=TRUE)   
mundoshp <- st_read("admin00.shp",quiet=TRUE)
mundocol <- mundoshp %>% 
  filter(CNTRY_NAME %in% c("Peru","Brazil","Venezuela","Ecuador","Panama")) %>% 
  st_as_sf()

for (p in 1:length(periods)) {
  temp <- membership_data %>%
    filter(period == periods[9] & community != "Inactive") %>%
    select(municipality_code, community)
  
  mapmpios <- mpioshp %>% 
    left_join(temp, by = c("MPIO_CCNCT" = "municipality_code")) %>%
    select(MPIO_CCNCT, community, geometry)
  mapmpios$community <- factor(as.numeric(as.character(mapmpios$community)), 
                               levels = sort(unique(as.numeric(as.character(mapmpios$community)))))
  
  mapmpios <- st_as_sf(mapmpios)
  colors_subset <- custom_colors[names(custom_colors) %in% levels(mapmpios$community)]
  #filename <- paste0("C:/Users/nperd/Documents/UN/Estadística/Redes/Proyecto/plots/municipalitiesmap_communities_", periods[p], ".png")
  
  #png(filename, width = 568, height = 568, res = 100)
  print(
    ggplot() +
      geom_sf(data = mapmpios, aes(fill = community), color = NA, lwd = 0.2) +
      geom_sf(data = mundocol, fill = "#DBDBDB", col = "darkgrey") +
      geom_sf(data = deptoshp, fill = NA, col = "darkgrey") +
      scale_fill_manual(values = colors_subset, name = "Q", na.value = "grey95") +
      coord_sf(xlim = c(-79.5, -66.5), ylim = c(-4.5, 13), expand = FALSE) +
      scale_x_continuous(breaks = seq(-79.5, -66.5, by = 4)) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
  )
  dev.off()
}

#"2005-2007"
temp <- membership_data %>%
  filter(period == "2002-2004" & next_community != "Inactive") %>%
  select(municipality_code, next_community)

mapmpios <- mpioshp %>% 
  left_join(temp, by = c("MPIO_CCNCT" = "municipality_code")) %>% 
  select(MPIO_CCNCT, next_community, geometry)
mapmpios$community <- factor(as.numeric(as.character(mapmpios$next_community)), 
                             levels = sort(unique(as.numeric(as.character(mapmpios$next_community)))))

mapmpios <- st_as_sf(mapmpios)
colors_subset <- custom_colors[names(custom_colors) %in% levels(mapmpios$community)]
filename <- paste0("C:/Users/nperd/Documents/UN/Estadística/Redes/Proyecto/plots/municipalitiesmap_communities_2005-2007.png")

png(filename, width = 650, height = 545, res = 100)
print(ggplot() +
        geom_sf(data = mapmpios, aes(fill = community), color = NA, lwd = 0.2) +
        geom_sf(data = mundocol, fill = "#DBDBDB", col = "darkgrey") +
        geom_sf(data = deptoshp, fill = NA, col = "darkgrey") +
        scale_fill_manual(values = colors_subset, name = "Q", na.value = "grey95") +
        coord_sf(xlim = c(-79.5, -66.5), ylim = c(-4.5, 13), expand = FALSE) +
        scale_x_continuous(breaks = seq(-79.5, -66.5, by = 4)) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
)
dev.off()

#-------------------------------------------------------------------------------
# Network plot by community

for (suffix in t) {
  municipalities_graph <- get(paste0("municipalities_graph_", suffix))
  
  png(paste0(output_path,"/graph_municipalities_", suffix, ".png"), width = 615, height = 545)
  
  unique_reg <- sort(unique(V(municipalities_graph)$membership_sbm))
  color_palette <- custom_colors[names(custom_colors) %in% unique_reg]
  color_palette <- adjustcolor(color_palette, 0.7)
  reg_colors <- setNames(color_palette, unique_reg)
  
  V(municipalities_graph)$color <- reg_colors[as.character(V(municipalities_graph)$membership_sbm)]
  
  par(mar = c(0.1, 0.1, 0.1, 0.1), mfrow = c(1, 1))
  set.seed(73)
  plot(
    municipalities_graph,
    layout = layout_with_fr(municipalities_graph),
    vertex.label = NA,
    vertex.color = V(municipalities_graph)$color,
    vertex.size = (1.5 * log(strength(municipalities_graph) + 1) + 2),
    edge.color = adjustcolor("gray", 0.95),
    edge.lty = 1
  )
  
  dev.off()
}
