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
dats <- replicate(n = 4, 
                  noisyDat(n.asf = 1, noisefraction = .15), 
                  simplify = FALSE)

# store the DGSs for each data set
targets <- lapply(dats, function(x) attributes(x)$target)

# add irrelevant factor to datasets
datsX <- lapply(dats, function(x) {x$X <- rbinom(nrow(x), 1, .5); return(x)})

# add irrelevant factor to the DGS models
modelsX <- lapply(targets, add_factor)

# datasets with irrelevant factor X are now in a list datsX. Corresponding
# models with one irrelevant factor value X or x added to the DGS 
# are in a list modelsX.
# The test loop should go through these such
# that for index i we do 
# cna_cortest(model = modelsX[[i]], data = datsX[[i]], ...), and 
# then calculate whatever we want to calculate from the results.


