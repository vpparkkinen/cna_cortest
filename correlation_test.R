library(frscore)
 


cna_cortest <- function(model, data){
  model <- cna::noblanks(model)
  lhss <- lhs(model) #grab the left hand side
  rhs <- rhs(model) # and the right hand side
  disjuncts <- unlist(strsplit(lhss, "\\+")) # separate disjuncts
  facs <- lapply(disjuncts, cna:::tryparse) # and then the conjuncts
  facs <- lapply(facs, all.vars) # and then the conjuncts
  
  
  res <- vector("list", length(disjuncts))
  names(res) <- disjuncts
  
  #calculate unconditional correlation if lhs is single disjunct with single conjunct
  if(length(disjuncts) == 1 & nchar(disjuncts[1]) == 1){
    facres <- vector("numeric", length(disjuncts))
    names(facres) <- disjuncts[[1]]
    facres[[1]] <- cor(data[toupper(disjuncts)], data[toupper(rhs)])
    res[[1]] <- facres
  } else {
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
        # addition to make it work with single disjunct
        if (alt_disjuncts_suppressed == "!()"){
          ctrl_expression <- cofactors_present
        }
        # addition to make it work with single conjunct
        if (cofactors_present == "()"){
          ctrl_expression <- alt_disjuncts_suppressed
        }
        # now subset the data
        test_data <- ct2df(selectCases(ctrl_expression, data))
        ## NOW CALCULATE THE CORRELATION 
        facres[[factor_index]] <- cor(test_data[test_factor], test_data[toupper(rhs)])
        ## AND STORE INTO YOUR RESULTS  OBJECT FOR THE INNER LOOP
      }
      res[[dis]] <- facres
    }
  }
  attr(res, "model") <- model
  return(res)
}


model <- "A<->E"

data <- d.error
cna_cortest(model, data)


# DGS is A + B*C <-> E but above model con/cov=1
# looking at how individual factors correlate with outcome
# will not help detect the model is incorrect
