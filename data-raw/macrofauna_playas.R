# Read the data ----
data0 <- readxl::read_xls("./data-raw/macrofauna_playas.xls")

macrofDat <- data0 |>
  filter(Locality == "Progreso" | Locality == "Sisal") |>
  filter(Zone == "Mesolittoral") |>
  relocate(Locality, Site, .before= core) |>
  select(-core, -Zone) |>
  as.data.frame()

macrofDat[c(8:10),2] <- 3
macrofDat <- macrofDat[-4,]
macrofDat <- macrofDat[-c(20, 26),]
macrofDat <- macrofDat[-30,]
macrofDat <- macrofDat[-c(35, 37,39),]

usethis::use_data(macrofDat, overwrite=TRUE)

## Working parameters ----
cases = 4
N = 100
M = 20
n = 15
m = 10
k = 50
transformation = "square root"
method = "bray"
dummy = T
useParallel = T
model = "nested.symmetric"
type = "counts"
Sest.method = "average"

# Start the analysis ----
simResultsNested <- prep_data(macrofDat, type, Sest.method, cases, N, M, n, m, k,
                        transformation, method, dummy, useParallel, model)

usethis::use_data(simResultsNested2, overwrite = TRUE)

betaNested2 <- sim_beta(simResultsNested2, 0.05)
usethis::use_data(betaNested2, overwrite=TRUE)

plot_power(betaNested2)
