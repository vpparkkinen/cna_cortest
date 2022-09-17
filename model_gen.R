# this needs cnasimtools
if(!require(cnasimtools)){
  remotes::install_github("vpparkkinen/cnasimtools")
}

library(cnasimtools)

# create noisy datasets from random 1-asf targets.
# amount of noise chosen arbitrarily here,
# could also try ideal data with one irrelevant factor added first.
# any case, assuming that this works, 
# change n = 4 to some reasonable number (100? 1000? 10000?) 
dats <- replicate(n = 4, 
                  noisyDat(n.asf = 1, noisefraction = .15), 
                  simplify = FALSE)

# store the DGSs for each data set
targets <- lapply(dats, function(x) attributes(x)$target)

# add irrelevant factor to datasets
datsX <- lapply(dats, function(x) {x$X <- rbinom(nrow(x), 1, .5); return(x)})

# function that adds a literal to a model, either as a disjunct or
# a conjunct
# 'fname' is the factor, not factor value, to be added to the model
# that is, value of 'fname' should be a capital letter
add_factor <- function(model, fname = "X"){ 
  fname <- toupper(fname) 
  model <- cna::noblanks(model)
  model <- gsub("^\\(|\\)$", "", model)
  lhs <- lhs(model) #grab the left hand side
  outcome <- rhs(model) # and the right hand side
  disjuncts <- unlist(strsplit(lhs, "\\+")) # lhs disjuncts
  as_disjunct <- sample(c(TRUE, FALSE), 1) # as its own disjunct, or not
  neg_irfac <- sample(c(TRUE, FALSE), 1) # to negate, or not
  if(neg_irfac){
    fname <- tolower(fname)
  }
  if(as_disjunct){
    newmodel_lhs <- paste0(lhs, "+", fname)
    newmodel <- paste0(newmodel_lhs, "<->", outcome)
  } else {
    pick_disj <- sample(seq_along(disjuncts), 1)
    new_disjunct <- paste0(disjuncts[pick_disj], "*", fname)
    new_lhs <- paste0(c(disjuncts[-pick_disj], new_disjunct), collapse = "+")
    newmodel <- paste0(new_lhs, "<->", outcome)
  }
  if(is.inus(newmodel)){
    return(newmodel)
  } else {
    add_factor(model, fname = fname)
  }
}

# add irrelevant factor to the DGS models
modelsX <- lapply(targets, add_factor)

# datasets with irrelevant factor X are now in a list datsX. Corresponding
# models with one irrelevant factor value X or x added to the DGS 
# are in a list modelsX.
# The test loop should go through these such
# that for index i we do 
# cna_cortest(model = modelsX[[i]], data = datsX[[i]], ..., and 
# then calculate whatever we want to calculate from the results.
