setwd("C:/Users/nperd/Documents/UN/Estadística/Redes/Proyecto/Datos")

# Packages
packages <- c("igraph", "blockmodels", "sbm", "readxl","writexl", 
                   "ggplot2", "stringr", "tidyr", "dplyr", "viridis", 
                   "RColorBrewer", "e1071", "fields", "networkD3", "htmltools",
                    "htmlwidgets", "sf", "poweRlaw")

for (i in packages) {
  if (!require(i, character.only = TRUE)) {
    suppressWarnings(suppressMessages(install.packages(i)))
    suppressWarnings(suppressMessages(library(i, character.only = TRUE)))
  }
}


# Data

#data <- read_excel("data.xlsx", col_types = "text")
rm(list=ls())
colombian_municipalities <- read_excel("colombian_municipalities.xlsx", col_types = "text")
armed_units_metadata <- read_excel("armed_units_metadata.xlsx", col_types = "text")

  
#-------------------------------------------------------------------------------
# Bipartite graph (affiliation network)
  time_periods <- list(
    "78_81" = "1978-1981", 
    "82_84" = "1982-1984",
    "85_87" = "1985-1987",
    "88_90" = "1988-1990",
    "91_93" = "1991-1993",
    "94_96" = "1994-1996",
    "97_99" = "1997-1999",
    "00_01" = "2000-2001",
    "02_04" = "2002-2004",
    "05_07" = "2005-2007"
  )
for (suffix in names(time_periods)) {
  period <- time_periods[[suffix]]    
  # Subset data for the period
  assign(paste0("vic_", suffix), data[data$period == period, ])   
  vic_period <- get(paste0("vic_", suffix))
  
  # Select the military units active for each period
  for (i in setdiff(names(vic_period), c("region", "subregion", "municipality_code", "period", "victims_farc", "victims_paramilitary", "victims_post_demobilization"))) {
    sum_value <- sum(vic_period[, i])  # Calculate the sum of each column
    if (sum_value == 0) {  # Remove column if the sum is zero
      vic_period <- vic_period[, !names(vic_period) %in% i]
    }
  }
  
  # Incidence matrix
  assign(paste0("vic_", suffix), vic_period)
  assign(paste0("M_", suffix), as.matrix(vic_period[!colnames(vic_period) %in% c("region", "subregion", "municipality_code", "period", "victims_farc", "victims_paramilitary", "victims_post_demobilization")]))
  
  # Bipartite graph from the incidence matrix (rows: municipalities; columns: military units)
  assign(paste0("graph_", suffix), graph_from_biadjacency_matrix(
    get(paste0("M_", suffix)), 
    mode = "all", 
    weighted = NULL
  ))
  graph_period <- get(paste0("graph_", suffix))
  
  # Bipartite graph: Vertex characteristics
  
  armed_units <- data.frame(armed_unit_code = colnames(vic_period)[!colnames(vic_period) %in% c("region", "subregion", "municipality_code", "period", "victims_farc", "victims_paramilitary", "victims_post_demobilization")])
  armed_units <- left_join(armed_units, armed_units_metadata, by = "armed_unit_code")
  armed_units <- armed_units %>% arrange(armed_unit_code)
  municipalities <- vic_period[colnames(vic_period) %in% c("region", "subregion", "municipality_code")]

  V(graph_period)$name <- c(as.character(municipalities$municipality_code), armed_units$armed_unit_code)    # Set vertex names
  V(graph_period)$type <- c(rep(TRUE, nrow(municipalities)), rep(FALSE, nrow(armed_units)))    # Define vertex type: TRUE for municipalities, FALSE for military units
  V(graph_period)$armed_unit_label <- c(rep(NA, nrow(municipalities)), armed_units$armed_unit)    # Label military units
  V(graph_period)$armed_group_type <- c(rep(NA, nrow(municipalities)), armed_units$armed_group_type)    # Assign armed group type to military units
  V(graph_period)$subregion <- c(municipalities$subregion, rep(NA, nrow(armed_units)))    # Assign subregion to municipalities
  V(graph_period)$region <- c(municipalities$region, rep(NA, nrow(armed_units)))    # Assign region to municipalities
  V(graph_period)$color <- ifelse(V(graph_period)$type == F,    # Set vertex colors based on armed group type and vertex type
                                  ifelse(V(graph_period)$armed_group_type == "guerrilla group", adjustcolor("#AA35B3", 0.8),
                                         ifelse(V(graph_period)$armed_group_type == "paramilitary group", adjustcolor("royalblue", 0.7),
                                                ifelse(V(graph_period)$armed_group_type == "organized crime", adjustcolor("#FF4C00", 0.7),
                                                       adjustcolor("#3AB795", 0.8)))),
                                  adjustcolor("gold", 0.7))
  V(graph_period)$shape <- ifelse(V(graph_period)$type == FALSE, "circle", "square")   # Set vertex shapes based on vertex type
  
  assign(paste0("graph_", suffix), graph_period)
  
#-------------------------------------------------------------------------------
# Projections 
  
  projections <- bipartite_projection(graph_period, multiplicity = TRUE)     # 'multiplicity = T' option adds edge weights to the resulting projections, representing the number of common neighbors between the vertices in the original bipartite graph
  assign(paste0("armed_units_graph_", suffix), projections[[1]])      # Military units' network projection
  assign(paste0("municipalities_graph_", suffix), projections[[2]])     # Municipalities' network projection
  
  armed_units_graph <- get(paste0("armed_units_graph_", suffix))
  municipalities_graph <- get(paste0("municipalities_graph_", suffix))
  
  # Armed units network vertices characteristics  
  V(armed_units_graph)$armed_unit_label <- armed_units$armed_unit
  V(armed_units_graph)$armed_group_type <- armed_units$armed_group_type
  V(armed_units_graph)$armed_group <- armed_units$armed_group_code
  V(armed_units_graph)$color <- ifelse(V(armed_units_graph)$armed_group_type == "guerrilla group", adjustcolor("#AA35B3", 0.8),
                                          ifelse(V(armed_units_graph)$armed_group_type == "paramilitary group", adjustcolor("royalblue", 0.7),
                                                 ifelse(V(armed_units_graph)$armed_group_type == "organized crime", adjustcolor("#FF4C00", 0.7),
                                                        adjustcolor("#3AB795", 0.8))))

  assign(paste0("armed_units_graph_", suffix), armed_units_graph)

  # Municipalities network vertices characteristics
  V(municipalities_graph)$shape <- "circle"
  V(municipalities_graph)$municipality_code <- municipalities$municipality_code
  V(municipalities_graph)$region <- municipalities$region
  V(municipalities_graph)$subregion <- municipalities$subregion
  V(municipalities_graph)$vic_farc <- vic_period$victims_farc
  V(municipalities_graph)$vic_paramilitary <- vic_period$victims_paramilitary
  V(municipalities_graph)$vic_postdem <- vic_period$victims_post_demobilization
  
  assign(paste0("municipalities_graph_", suffix), municipalities_graph)

#-------------------------------------------------------------------------------
# Subgrpahs
  
  # Giant component, bipartite graph
  assign(paste0("giant_comp_", suffix), induced_subgraph(graph_period, which(components(graph_period)$membership == which.max(components(graph_period)$csize))))  
  
  # Giant component, military units' network  
  assign(paste0("giant_comp_armed_unit_", suffix), induced_subgraph(armed_units_graph, which(components(armed_units_graph)$membership == which.max(components(armed_units_graph)$csize))))  
  
  # Giant component, municipalities' network  
  assign(paste0("giant_comp_municipalities_", suffix), induced_subgraph(municipalities_graph, which(components(municipalities_graph)$membership == which.max(components(municipalities_graph)$csize))))  
}

rm(list = c("graph_period", "armed_units", "municipalities", "vic_period", "projections", "municipalities_graph", "armed_units_graph"))    

