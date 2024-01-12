## code to prepare `harmonic_sample` dataset goes here

t <- seq(1024)
harmonic_sample_data = ts( cos(t*0.13423167+0.00) + 1.3 * cos(t*0.119432+2.314) + 0.134994 + 0.4*cos(t*0.653167) + 0.11 * cos(t*0.78913498) + rnorm(1024)*0.12 , start=0, deltat=1) 

harmonic_sample <- list()

harmonic_sample_spectrum <- list(
  Amp = c(1, 1.3, 0.134994, 0.4, 0.11), 
  Freq = c(0.13423167, 0.119432, 0, 0.653167, 0.78913498),
  Phase = c(0, 2.314, 0, 0.653167, 0))

attr(harmonic_sample_spectrum,"class") <- 'mfft_deco'

harmonic_sample$data <- harmonic_sample_data
harmonic_sample$spectrum <- harmonic_sample_spectrum



usethis::use_data(harmonic_sample, overwrite = TRUE)
