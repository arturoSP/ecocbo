epibionts <- read_csv("./data-raw/epibionts.csv")[,-1]

# After checking for the sites that are less similar it was found that E-2, M-1 and I-3 will be used.
 epiDat <- epibionts[epibionts$sector=="E"&epibionts$site==2 |
             epibionts$sector=="M"&epibionts$site==1 |
             epibionts$sector=="I"&epibionts$site==3,]
 epiDat[,"site"] <- paste0(epiDat$sector, epiDat$site)
 epiDat <- epiDat[,-1]

usethis::use_data(epiDat, overwrite = TRUE)
