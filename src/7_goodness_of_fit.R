# Goodness-of-fit
output_path <- ""   

#-------------------------------------------------------------------------------
# Monte Carlo simulations

periods <- c("78_81", "82_84", "85_87", "88_90", "91_93", "94_96", "97_99", "00_01", "02_04", "05_07")
ntrials <- 10000
set.seed(42)

simulation_results <- list()
observed_results <- list()
ppp_values <- list()


for (period in periods) {
  lambda.mat <- get(paste0("Lambda_m_", period))
  graph <- get(paste0("municipalities_graph_", period))
  nv <- vcount(graph)
  alpha <- get(paste0("alpha_m_", period))
  
  density <- vector("list", ntrials)
  glo.trans <- vector("list", ntrials)
  mean_strength <- vector("list", ntrials)
  sd_strength <- vector("list", ntrials)
  betw_centr <- vector("list", ntrials)
  max_coreness <- vector("list", ntrials)
  n_components <- vector("list", ntrials)
  
  for (i in 1:ntrials) {
    sbm <- sampleSimpleSBM(nv, alpha, list(mean = round(lambda.mat,5)), model = "poisson")
    adj_mat <- sbm$rNetwork()$networkData
    g.sbm <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", weighted = TRUE, diag = FALSE)
    
    comps <- components(g.sbm)
    giant_component <- induced_subgraph(g.sbm, which(comps$membership == which.max(comps$csize)))
    
    # Network statistics
    density[[i]] <- edge_density(g.sbm)
    glo.trans[[i]] <- transitivity(g.sbm, type = "global")
    mean_strength[[i]] <- mean(strength(giant_component))
    sd_strength[[i]] <- sd(strength(giant_component))
    betw_centr[[i]] <- centr_betw(giant_component, normalized = TRUE)$centralization
    max_coreness[[i]] <- max(coreness(giant_component))
    n_components[[i]] <- comps$no
  }
  
  simulation_df <- data.frame(
    density = unlist(density),
    glo.trans = unlist(glo.trans),
    mean_strength = unlist(mean_strength),
    sd_strength = unlist(sd_strength),
    betw_centr = unlist(betw_centr),
    max_coreness = unlist(max_coreness),
    components = unlist(n_components),
    period = period
  )
  simulation_results[[period]] <- simulation_df
  
  # Observed values
  comps_obs <- components(graph)
  giant_component_obs <- induced_subgraph(graph, which(comps_obs$membership == which.max(comps_obs$csize)))
  
  density_obs <- edge_density(graph)
  glo.trans_obs <- transitivity(graph, type = "global")
  mean_strength_obs <- mean(strength(giant_component_obs))
  sd_strength_obs <- sd(strength(giant_component_obs))
  betw_centr_obs <- centr_betw(giant_component_obs, normalized = TRUE)$centralization
  max_coreness_obs <- max(coreness(giant_component_obs))
  components_obs <- comps_obs$no
  
  observed_values <- data.frame(
    statistic = c("density", "glo.trans", "mean_strength", "sd_strength", "betw_centr", "max_coreness", "components"),
    value = c(density_obs, glo.trans_obs, mean_strength_obs, sd_strength_obs, betw_centr_obs, max_coreness_obs, components_obs),
    period = period
  )
  observed_results[[period]] <- observed_values
  
  # Posterior predictive p-values (ppp)
  obs_vals <- observed_values$value
  names(obs_vals) <- observed_values$statistic
  
  ppp_period <- sapply(names(obs_vals), function(stat) {
    sims <- simulation_df[[stat]]
    mean(sims > obs_vals[stat], na.rm = TRUE)  # one-sided
  })
  
  ppp_values[[period]] <- ppp_period
}

#-------------------------------------------------------------------------------
# Merge simulation results and prepare for plotting
simulation_results <- do.call(rbind, simulation_results)
observed_results <- do.call(rbind, observed_results)

o <- c("78_81", "82_84", "85_87", "88_90", "91_93", "94_96", "97_99", "00_01", "02_04", "05_07")
n <- c("1978-1981", "1982-1984", "1985-1987", "1988-1990", "1991-1993", "1994-1996", "1997-1999", "2000-2001", "2002-2004", "2005-2007")
period_map <- setNames(n, o)
names(ppp_values) <- period_map[names(ppp_values)]


ppp_df <- do.call(rbind, ppp_values)
ppp_df <- as.data.frame(ppp_df)
ppp_df$period <- factor(rownames(ppp_df), levels = n)

library(reshape2)
ppp_long <- melt(ppp_df, id.vars = "period", variable.name = "statistic", value.name = "ppp")

#-------------------------------------------------------------------------------
# Boxplots for ppp values, one plot per statistic
output_path <- "C:/Users/nperd/Documents/UN/Estadística/Redes/Proyecto/plots/"

statistics <- unique(ppp_long$statistic)

for (stat in statistics) {
  stat_data <- subset(ppp_long, statistic == stat)
  png(filename = paste0(output_path, "ppp_boxplot_", stat, ".png"), width = 600, height = 450)
  par(mar = c(6.5, 5, 0.5, 1))
  boxplot(ppp ~ period, data = stat_data,
          col = "gray90",
          main = "",
          xlab = "",
          ylab = "ppp",
          outline = TRUE,
          las = 2,
          ylim = c(0, 1),
          cex.axis = 1.2,  
          cex.lab = 1.5)   
  dev.off()
}
