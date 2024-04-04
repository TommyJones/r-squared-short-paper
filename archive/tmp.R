library(httr2)
library(rvest)

req <- request("https://fraser.stlouisfed.org/api") |>
  req_headers(`X-API-Key` = "b6cb6e135c542ebe963890c917c8005f")


fomc <- 
  req |>
  req_url_path("/theme") |> # req_url_path("/title/605825/items") |> # 
  req_perform() 

fomc_body <-
  fomc$body |>
  rawToChar()

path <- 
  "https://fraser.stlouisfed.org/title/federal-open-market-committee-meeting-minutes-transcripts-documents-677/meeting-september-21-22-2021-605825"

out <-
  path |>
  read_html() |>
  html_elements("a") |>
  html_attr("href") |>
  unique()

