

library(tidylda)

library(tidyverse)

library(Matrix)

library(mvrsquared)

library(furrr)

plan(multisession, workers = parallel::detectCores() - 1)


models <- read_rds("data-derived/sbir-topic-models.rds")

dtm <- read_rds("data-derived/sbir-dtm.rds")

calc_partial_r2 <- function(model, dtm, threads = parallel::detectCores() - 1) {
  
  SSR_full <- calc_rsquared(
    y = dtm,
    yhat = list(
      x = Matrix::rowSums(dtm) * model$theta,
      w = model$beta
    ),
    return_ss_only = TRUE,
    ybar = colMeans(dtm),
    threads = threads
  )
  
  SSR_partial <- 
    future_map(
      seq_len(nrow(model$beta)),
      function(j) {
        
        beta_new <- model$beta[- j, ]
        
        theta_new <- model$theta[, - j]
        
        theta_new <- theta_new / rowSums(theta_new)
        
        theta_new[is.na(theta_new)] <- 0
        
        result <- calc_rsquared(
          y = dtm,
          yhat = list(
            x = Matrix::rowSums(dtm) * theta_new,
            w = beta_new
          ),
          return_ss_only = TRUE,
          ybar = colMeans(dtm),
          threads = 1
        )
      }, .options = furrr_options(
        packages = c("tidylda", "Matrix", "mvrsquared"), 
        seed = TRUE)
    )
  
  result <- 
    1 - SSR_full[1] / unlist(SSR_partial |> map(function(x) x[1]))
  
  result
}

models <- 
  models |>
  map(
    function(m){
      m$summary$part_r2 <- 
        calc_partial_r2(model = m, dtm = dtm)
      m
    }
  )


model_summaries <- 
  map(models, function(x){
    x$summary
  }) |>
  bind_rows() |>
  mutate(
    model = c(rep(1, 100), rep(2, 150))
  )

write_rds(model_summaries, "data-derived/sbir-model-summaries.rds")
