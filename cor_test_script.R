library(here) 

# this needs cnasimtools
if(!require(cnasimtools)){
  remotes::install_github("vpparkkinen/cnasimtools")
}
library(cnasimtools)
source(here("functions.R"))

# create noisy datasets from random 1-asf targets.
# amount of noise chosen arbitrarily here,
# could also try ideal data with one irrelevant factor added first.
# any case, 
# change n = 4 to some reasonable number (100? 1000? 10000?) 
dats <- replicate(n = 1000, 
                  noisyDat(set_N = 100, n.asf = 1, noisefraction = .1), 
                  simplify = FALSE)

# store the DGSs for each data set
targets <- lapply(dats, function(x) attributes(x)$target)
targets <- lapply(targets, function(x) gsub("^\\(|\\)$", "", x))
targets <- lapply(targets, cna::noblanks)
outcomes <- lapply(targets, rhs)

# add irrelevant factor to datasets
datsX <- lapply(dats, function(x) {x$X <- rbinom(nrow(x), 1, .5); return(x)})

# add irrelevant factor to the DGS models
# modelsX <- lapply(targets, add_factor)

cor_matrices <- lapply(datsX, cor)
cor_matrices_outcome <- vector("list", length(cor_matrices))

for(i in seq_along(cor_matrices)){
  cor_matrices_outcome[[i]] <- cor_matrices[[i]][
    ,which(dimnames(cor_matrices[[i]])[[2]] == outcomes[[i]])
    ]
  cor_matrices_outcome[[i]] <- cor_matrices_outcome[[i]][
    -which(names(cor_matrices_outcome[[i]]) == outcomes[[i]])
    ]
}

cor_matrices_outcome <- lapply(cor_matrices_outcome, abs)

cor_matrices_outcome_nox <- lapply(cor_matrices_outcome, 
                                   function(x)
                                     x[-which(names(x) == "X")])

nox_grand_mean <- mean(unlist(cor_matrices_outcome_nox))

x_cors <- lapply(cor_matrices_outcome,
                 function(x) x[which(names(x) == "X")])

x_mean <- mean(unlist(x_cors))

# cor_matrices_differences <- lapply(cor_matrices_outcome, 
#                                    function(x) 
#                                      x[
#                                        -which(names(x) == "X")
#                                        ] - x[which(names(x) == "X")])
# 
# mean_diffs <- lapply(cor_matrices_differences, mean)
# grand_mean_diff <- mean(unlist(mean_diffs))

x_mean / nox_grand_mean

