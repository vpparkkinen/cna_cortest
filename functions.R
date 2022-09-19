
# function that adds a factor to a model, either as a disjunct, or
# an extra conjunct to an existing disjunct.
# 'fname' is the factor, not factor value, to be added to the model,
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


cna_cortest <- function(model, data, suppress_alt = TRUE){
  model <- cna::noblanks(model)
  lhss <- lhs(model) #grab the left hand side
  rhs <- rhs(model) # and the right hand side
  disjuncts <- unlist(strsplit(lhss, "\\+")) # separate disjuncts
  facs <- lapply(disjuncts, cna:::tryparse) # and then the conjuncts
  facs <- lapply(facs, all.vars) # and then the conjuncts
  
  
  res <- vector("list", length(disjuncts))
  names(res) <- disjuncts
  
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
        
        if(length(cofactors) >= 1){
          cofactors_present <- paste0(cofactors, collapse = "*")
          # remove possible trailing "*":
          cofactors_present <- gsub("\\*$", "", cofactors_present)
          # wrap in parentheses to be extra careful
          cofactors_present <- paste0("(", cofactors_present, ")")  
        } else {
          cofactors_present <- NULL
        }
        
        # now create the expression to negate alternative disjuncts
        if(suppress_alt & length(alt_disjuncts) >= 1){
          alt_disjuncts_suppressed <- paste0(alt_disjuncts, collapse = "+")
          # remove possible trailing "+":
          alt_disjuncts_suppressed <- gsub("\\+$", "", alt_disjuncts_suppressed)
          alt_disjuncts_suppressed <- paste0("!(", alt_disjuncts_suppressed, ")")
        } else {
          alt_disjuncts_suppressed <- NULL
        }
        # now we can create the expression that goes to selectCases()
        if(is.null(alt_disjuncts_suppressed) & is.null(cofactors_present)){
          test_data <- data
        } else {
          ctrl_expression <- paste0(cofactors_present, 
                                    "*", 
                                    alt_disjuncts_suppressed)
          ctrl_expression <- gsub("^\\*|\\*$", "", ctrl_expression)
          test_data <- ct2df(selectCases(ctrl_expression, data))
        }
        
        # addition to make it work with single disjunct
        # if (cofactors_present == "()"){
        #   ctrl_expression <- alt_disjuncts_suppressed
        # }
        # 
        # if (alt_disjuncts_suppressed == "!()"){
        #   ctrl_expression <- cofactors_present
        # }
        # addition to make it work with single conjunct
        
        # now subset the data
        
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

generate_cand_models <- function(dat, 
                                 regex = "[Xx]", #regex to detect a factor
                                 outcome,
                                 cc_values = 
                                   expand.grid(seq(.95,.7,-.05),
                                               seq(.95,.7,-.05))
                                 ){
  if(nrow(cc_values) == 0){
    return(NULL)
  } else {
    #shuffle <- sample(1:nrow(cc_values), nrow(cc_values), replace = FALSE)
    #cc_values <- cc_values[shuffle,]
    rrow <- sample(1:nrow(cc_values), 1)
    cnares <- cna(dat, 
                  con = cc_values[rrow,1], 
                  cov = cc_values[rrow,2],
                  outcome = outcome)
    asfs <- asf(cnares)[,2]
    #regex <- paste0(toupper(factor_to_pick), "|", tolower(factor_to_pick))
    oks <- sapply(asfs, function(x) grepl(regex, x))
    if(all(!oks)){
      cc_values <- cc_values[-rrow,]
      generate_cand_models(dat = dat, 
                           #factor_to_pick = factor_to_pick,
                           regex = regex,
                           outcome = outcome,
                           cc_values = cc_values)
    } else {
      out <- asfs[oks]
      return(out)
    }
  }
}


