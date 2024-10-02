## code to prepare `pseudo_log` dataset goes here

t <- seq(8192)
# pseudo_log_data = ts( sign ( cos(t*0.13423167+0.00) + 1.3 * cos(t*0.049432+2.314) + 0.7))

pseudo_log_data = ts( sign( cos(t*0.13423167+0.00) + 1.3 * cos(t*0.119432+2.314) + 0.134994 + 0.4*cos(t*0.653167) + 0.11 * cos(t*0.78913498) + 0.8) , start=0, deltat=1) 
pseudo_log <- list()


pseudo_log_spectrum <- list(
  Amp = c(1, 1.3, 0.134994, 0.4, 0.11), 
  Freq = c(0.13423167, 0.119432, 0, 0.653167, 0.78913498),
  Phases = c(0, 2.314, 0, 0.653167, 0))


attr(pseudo_log_spectrum,"class") <- 'mfft_deco'

pseudo_log$data <- pseudo_log_data
pseudo_log$spectrum <- pseudo_log_spectrum


usethis::use_data(pseudo_log, overwrite = TRUE)
