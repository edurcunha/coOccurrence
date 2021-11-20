################################################################################
#  cScore
################################################################################

#  FUNCTION DESCRIPTION: This function calculates the c-score index for a pair 
#                        of species in a collection of samples.

#  ARGUMENTS
# m: A matrix of two species (columns) by n samples (rows).
# standardize: Logical. If TRUE, c-score is standardize to number of potential
#              combinations between samples.


cScore <- function(m, ...) {

    
  extra.arg <- list(...)
  
  if( "standardize" %in% names(extra.arg) ) {
    
    standardize <- extra.arg[["standardize"]]
    
  } else {
    
    standardize <- FALSE
    
  }
  
  if( is.matrix(m) ) {
    
    if( ncol(m) == 2 ) {
      
      if( all.equal( sort( unique( c(m, 0, 1) ) ), c(0, 1) ) ){

        if( sum( rowSums(m) > 0 ) >= 1 ) {
          
          S_ij <- sum( ( m[,1] == 1 ) & (m[,2] == 1) )
          r_i <- sum( m[,1] )
          r_j <- sum( m[,2] )
          cScore <- ( r_i - S_ij ) * ( r_j - S_ij )
          
          if( standardize ) {
            cScore <- cScore / ( ( ( nrow(m)^ 2 ) / 2 ) - (nrow(m) / 2) ) 
          }
          
          return(cScore)
          
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