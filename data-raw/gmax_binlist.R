# Initializing a list of 1-kbp bins into which reads will be binned
gmax_binlist <- list()

# Generating breakpoints for the bins of every chromosome
# Here we consider 1-kbp bins ; this could be made longer or shorter
for(i in 1:nrow(delgbs::gmax_chr_sizes)) {
  gmax_binlist[[i]] <-
    (0:ceiling(delgbs::gmax_chr_sizes[["length"]][i] / 1000)) * 1000
}

# Naming each element of the list according to its chromosome
names(gmax_binlist) <- delgbs::gmax_chr_sizes$chr

# Saving to file using the devtools utility
devtools::use_data(gmax_binlist, compress = "xz")
