library(tidytext)
library(tidyverse)
library(Matrix)
library(tidyla)

set.seed(8675309)

# prep parallel framework
library(furrr)

plan(multisession, workers = parallel::detectCores() - 1)

# load data and sample rows
dtm <- read_rds("data-derived/sbir-dtm.rds")

rows <-
  sample(1:nrow(dtm), 1000)

# function to buld over a range of topics
estimate_num_topics <- function(dtm, k_range) {
  
  est_k <- 
    future_map(
      k_range, 
      function(k) {
        m <- 
          tidylda(
            data = dtm,
            k = k,
            iterations = 200,
            burnin = 150,
            calc_r2 = TRUE
          )
        
        tibble(
          k = k,
          mean_coherence = try(mean(m$summary$coherence)),
          var_coherence = try(var(m$summary$coherence)),
          pct_high_coherence = try(sum(m$summary$coherence > 0.01) / nrow(m$summary)),
          r2 = try(m$r2),
          likelihood = try(mean(m$log_likelihood$log_likelihood[51:200]))
        )
      }, .options = furrr_options(
        packages = c("tidylda", "tidyverse", "Matrix"), 
        seed = TRUE),
      .progress = TRUE
    ) |> 
    bind_rows()
  
  # use loess to get optimals by variable
  
  optimal_k <- try({
    
    opt_c <- loess(mean_coherence ~ k, data = est_k)
    
    opt_r2 <- loess(r2 ~ k, data = est_k)
    
    opt_ll <- loess(likelihood ~ k, data = est_k)
    
    
    tibble(
      coherence = opt_c$x[which.max(opt_c$fitted)],
      r2 = opt_r2$x[which.max(opt_r2$fitted)],
      likelihood = opt_ll$x[which.max(opt_ll$fitted)]
    )
  })
  
  
  
  list(
    est_k = est_k,
    optimal_k = optimal_k
  )
  
}

dtm_samp <- dtm[rows, ]

dtm_samp <- dtm_samp[, colSums(dtm_samp) > 0]

est_k <- estimate_num_topics(dtm = dtm_samp, k_range = seq(50, 300, by = 10))

write_rds(
  est_k,
  "data-derived/sbir-model-selection.rds"
)


