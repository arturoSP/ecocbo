# Read the data ----
dataFish <- readr::read_csv(
  "data-raw/PRCRMP_Fish_Biomass.csv",
  locale = readr::locale(encoding = "ISO-8859-1")
) |>
  filter(YEAR == 2021) |>
  #filter(YEAR == 2019)|>
  filter(`DEPTH ZONE` == "very shallow") |>
  mutate(
    REGION = stringr::str_replace_all(
      REGION,
      c(
        "Mona/Desecheo" = "WEST",
        "Vieques/Culebra" = "EAST",
        "Southeast" = "EAST"
      )
    )
  ) |>
  mutate(
    REGION = ifelse(
      LOCATION == "Fajardo",
      "EAST",
      ifelse(
        LOCATION == "Vega Baja",
        "NORTH",
        ifelse(
          LOCATION == "La Parguera",
          "SOUTH",
          ifelse(LOCATION == "San Juan", "NORTH", REGION)
        )
      )
    )
  ) |>
  filter(REGION != "Southeast") |>
  mutate(REGION = stringr::str_to_upper(REGION)) |>
  select(-c(YEAR, `DEPTH ZONE`, `SITE NAME`, TRANSECT)) |>
  mutate(across(everything(), ~ replace_na(.x, 0))) |>
  filter(REGION != "SOUTHWEST") |>
  as.data.frame() |>
  select(-1) |>
  mutate(REGION = factor(REGION), LOCATION = factor(LOCATION))

summary(dataFish)

usethis::use_data(dataFish, overwrite = TRUE)
