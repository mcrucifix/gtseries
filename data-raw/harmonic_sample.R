## code to prepare `harmonic_sample` dataset goes here

t <- seq(1024)
harmonic_sample = ts( cos(t*0.13423167+0.00) + 1.3 * cos(t*0.119432+2.314) + 0.134994 + 0.4*cos(t*0.653167) + 0.11 * cos(t*0.78913498) + rnorm(1024)*0.12 , start=0, deltat=1) 


usethis::use_data(harmonic_sample, overwrite = TRUE)
