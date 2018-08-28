# Generating a data.frame of chromosome sizes
chr_sizes <- read.table("gmax_chr_size.txt",
                        head = TRUE,
                        stringsAsFactors = FALSE)

# Initializing a list of 5-kbp bins into which reads will be binned
Gmax_binlist_5kb <- list()

# Generating breakpoints for the bins of every chromosome
# Here we consider 1-kbp bins ; this could be made longer or shorter
for(i in 1:nrow(chr_sizes)) {
  Gmax_binlist_5kb[[i]] <- (0:ceiling(chr_sizes[["length"]][i] / 5000)) * 5000
}

# Naming each element of the list according to its chromosome
names(Gmax_binlist_5kb) <- chr_sizes$chr

save(Gmax_binlist_5kb, file = "data/Gmax_binlist_5kb.RData")
