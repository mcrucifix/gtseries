require(palinsol) 

insol_scaled <- ts(scale(sapply(seq(4096)*1000, function(t) Insol(ber90(t)))))
insol_truncated <- sign(insol_scaled+1.2)


usethis::use_data(insol_scaled, overwrite = TRUE, compress="xz")
usethis::use_data(insol_truncated, overwrite = TRUE, compress="xz")
