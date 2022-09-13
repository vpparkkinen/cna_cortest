library(frscore)
model <- "A*c+A*D+B*C<->E"
model <- cna::noblanks(model)
data <- d.error 
# DGS is A + B*C <-> E but above model con/cov=1
# looking at how individual factors correlate with outcome
# will not help detect the model is incorrect
lhss <- lhs(model) #grab the left hand side
rhs <- rhs(model) # and the right hand side
disjuncts <- unlist(strsplit(lhss, "\\+")) # separate disjuncts
facs <- lapply(disjuncts, cna:::tryparse) # and then the conjuncts
facs <- lapply(facs, all.vars) # and then the conjuncts


res <- vector("list", length(disjuncts))
names(res) <- disjuncts

for(dis in seq_along(disjuncts)){
  #get the alternative disjuncts that you want to suppress
  alt_disjuncts <- disjuncts[-which(disjuncts == disjuncts[[dis]])]
  facres <- vector("numeric", length(facs[[dis]]))
  names(facres) <- facs[[dis]]
  for(factor_index in seq_along(facs[[dis]])){
    test_factor <- facs[[dis]][factor_index]
    test_factor <- toupper(test_factor) # just for the correlation test
    cofactors <- facs[[dis]][-factor_index]
    # "cofactors present" expression:
    cofactors_present <- paste0(cofactors, sep = "*", collapse = "")
    # remove possible trailing "*":
    cofactors_present <- gsub("\\*$", "", cofactors_present)
    # wrap in parentheses to be extra careful
    cofactors_present <- paste0("(", cofactors_present, ")")
    # now create the expression to negate alternative disjuncts
    alt_disjuncts_suppressed <- paste0(alt_disjuncts, sep = "+", collapse = "")
    # remove possible trailing "+":
    alt_disjuncts_suppressed <- gsub("\\+$", "", alt_disjuncts_suppressed)
    alt_disjuncts_suppressed <- paste0("!(", alt_disjuncts_suppressed, ")")
    # now we can create the expression that goes to selectCases()
    ctrl_expression <- paste0(cofactors_present, "*", alt_disjuncts_suppressed)
    # now subset the data
    test_data <- ct2df(selectCases(ctrl_expression, data))
    ## NOW CALCULATE THE CORRELATION 
    facres[[factor_index]] <- cor(test_data[test_factor], test_data[toupper(rhs)])
    ## AND STORE INTO YOUR RESULTS  OBJECT FOR THE INNER LOOP
  }
  res[[dis]] <- facres
}

