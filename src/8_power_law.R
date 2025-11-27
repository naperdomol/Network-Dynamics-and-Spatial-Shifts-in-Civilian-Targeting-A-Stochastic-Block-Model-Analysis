install.packages("poweRlaw")
library(poweRlaw)



periods <- c("85_87", "88_90", "91_93", "97_99", "00_01", "02_04", "05_07")

pl_results <- list()
set.seed(1234)
for (suffix in periods) {
  g <- get(paste0("giant_comp_municipalities_", suffix))
  x <- strength(g, mode = "all")
  
  # Fit discrete power law
  m_pl <- displ$new(x)
  est_pl <- estimate_xmin(m_pl)
  m_pl$setXmin(est_pl)
  m_pl$setPars(estimate_pars(m_pl))
  p_pl <- bootstrap_p(m_pl, no_of_sims = 5000, threads = 2)
  
  # Log-normal
  m_ln <- dislnorm$new(x)
  m_ln$setXmin(m_pl$getXmin())
  m_ln$setPars(estimate_pars(m_ln))
  res_ln <- compare_distributions(m_pl, m_ln)
  
  # Exponential
  m_exp <- disexp$new(x)
  m_exp$setXmin(m_pl$getXmin())
  m_exp$setPars(estimate_pars(m_exp))
  res_exp <- compare_distributions(m_pl, m_exp)
  
  # Poisson
  m_pois <- dispois$new(x)
  m_pois$setXmin(m_pl$getXmin())
  m_pois$setPars(estimate_pars(m_pois))
  res_pois <- compare_distributions(m_pl, m_pois)
  
  pl_results[[suffix]] <- list(
    n = length(x),
    xmin = m_pl$getXmin(),
    alpha = round(m_pl$pars, 3),
    ntail = sum(x >= m_pl$getXmin()),
    p = round(p_pl$p, 3),
    gof = round(p_pl$gof, 5),
    lognorm = res_ln,
    gof = p_pl$gof,
    LR_lognorm = res_ln$test_statistic,
    p_lognorm = res_ln$p_one_sided,
    LR_exp = res_exp$test_statistic,
    p_exp = res_exp$p_one_sided,
    LR_pois = res_pois$test_statistic,
    p_pois = res_pois$p_one_sided
  )
}


results_table <- data.frame(
  period = character(),
  power_law_p = numeric(),
  LR_lognorm = numeric(),
  p_lognorm = numeric(),
  LR_exp = numeric(),
  p_exp = numeric(),
  LR_pois = numeric(),
  p_pois = numeric()
)

for (suffix in names(pl_results)) {
  res <- pl_results[[suffix]]
  results_table <- rbind(results_table, data.frame(
    period = suffix,
    power_law_p = res$p,
    LR_lognorm = round(res$LR_lognorm, 3),
    p_lognorm = round(res$p_lognorm, 3),
    LR_exp = round(res$LR_exp, 3),
    p_exp = round(res$p_exp, 3),
    LR_pois = round(res$LR_pois, 3),
    p_pois = round(res$p_pois, 3)
  ))
}

# Print or export
print(results_table)

summary_table <- data.frame(
  period = character(),
  n = integer(),
  ntail = integer(),
  alpha = numeric(),
  xmin = integer(),
  ks_stat = numeric(),
  p_value = numeric()
)

for (suffix in names(pl_results)) {
  res <- pl_results[[suffix]]
  summary_table <- rbind(summary_table, data.frame(
    period = suffix,
    n = res$n,
    ntail = res$ntail,
    alpha = res$alpha,
    xmin = res$xmin,
    ks_stat = round(res$gof, 5),
    p_value = res$p
  ))
}

print(summary_table)

