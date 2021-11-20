################################################################################
#  cScoreRand
################################################################################

#  FUNCTION DESCRIPTION: This function calculates the numerical distribution
#                        of c-scores using randomizations of species occurrence
#                        a pair of species in a collection of samples. The
#                        algorithm follows the SIM2 randomization model of
#                        Goteli (2000).

#  ARGUMENTS
# m: A matrix of two species (columns) by n samples (rows).
# standardize: Logical. If TRUE, c-score is standardize to number of potential
#              combinations between samples.
# rand: A number indicating the number of random collections of samples to 
#       generate.


cScoreRand <- function(m, ...) {
  
  extra.arg <- list(...)
  
  if( "standardize" %in% names(extra.arg) ) {
    
    standardize <- extra.arg[["standardize"]]
    
  } else {
    
    standardize <- FALSE
    
  }
  
  if( "rand" %in% names(extra.arg) ) {
    
    rand <- extra.arg[["rand"]]
    
  } else {
    
    rand <- 999
    
  }
  
  if( is.matrix(m) ) {
    
    if( ncol(m) == 2 ) {
      
      if( all.equal( sort( unique( c(m, 0, 1) ) ), c(0, 1) ) ){
        
        if( sum( rowSums(m) > 0 ) >= 1 ) {
          
          cScoreRand <- rep(NA, rand)
          
          for(i in 1:rand) {
            
            m.rand <- m.tmp
            
            m.rand[,2] <- m.rand[ sample( nrow(m.rand) ) ,2] 
            
            cScoreRand[i] <- cScore(m.rand, standardize = standardize)
            
          }
          
          return(cScoreRand)
          
        }        
        
      }
      
    }
    
  }
  
  stop("'m' is not a two species (columns) by n sites (rows) 
       presence-absence matrix")
}

# REFERENCE
# Stone. L. and A. Roberts. 1990. The checkerboard score and species 
# distributions. Oecologia 85: 74-79
# Gotelli, N.J.. 2000. Null model analysis of species co‐occurrence 
# patterns. Ecology, 81: 2606-2621