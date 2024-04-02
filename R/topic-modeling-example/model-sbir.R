### Preamble ----
# This runs the SBIR topic models and gets their R-squared values

### Load libraries ----
library(tidyverse)
library(tidytext)
library(tidylda)
library(janitor)
library(Matrix)

### Load the SBIR data and curate it ----
sbir <- read_csv("data-raw/sbir/award_data.csv") |>
  clean_names()

sbir <- 
  sbir |>
  filter(award_year %in% 2020:2022) |>
  mutate(
    id = 1:n(),
    abstract = str_conv(abstract, "UTF8")
  )

# make a DTM
sbir_dtm <-
  sbir |>
  select(
    id,
    text = abstract
  ) |>
  unnest_tokens(
    output = word,
    input = text,
    stopwords = stop_words$word,
    token = "ngrams",
    n_min = 1, 
    n = 1
  ) |>
  count(id, word) |>
  # filter words that are just numbers or or appear 2 or fewer times in a document
  filter(n > 2 | ! stringr::str_detect(word, "^[0-9]+$")) |> 
  cast_sparse(id, word, n)

# count term frequencies
tf <- tibble(
  word = colnames(sbir_dtm), 
  freq = colSums(sbir_dtm), 
  docs = colSums(sbir_dtm > 0)
)

# remove words where they appear in fewer than 5 documents
vocab <- tf$word[tf$docs >= 5]

sbir_dtm <- sbir_dtm[, vocab]

write_rds(
  sbir_dtm,
  "data-derived/sbir-dtm.rds"
)

### fit models with 100 and 150 topics ----
models <- 
  c(100, 150) |>
  parallel::mclapply(
    function(k) {
      tidylda(
        data = sbir_dtm,
        k = k,
        iterations = 200,
        burnin = 150,
        calc_r2 = TRUE
      )
    },
    mc.cores = 2
  )


### save the results ----
write_rds(
  models,
  file = "data-derived/sbir-topic-models.rds"
)

write_rds(
  tibble(
    k = models |> map(function(x) nrow(x$beta)) |> unlist(),
    r2 = models |> map(function(x) x$r2) |> unlist()
  ),
  file = "data-derived/sbir-topic-models-r2.rds"
)
