################################################################################
#  cScore
################################################################################

#  FUNCTION DESCRIPTION: This function calculates the c-score index for a pair 
#                        of species in a collection of samples

#  ARGUMENTS
# m: A matrix of two species (columns) by n samples (rows)


cScore <- function(m) {
  
  if( is.matrix(m) ) {
    
    if( ncol(m) == 2 ) {
      
      if( all.equal( unique( c(m) ), c(0, 1) ) ){

        if( sum( rowSums(m) > 0 ) >= 2 ) {
          
          S_ij <- sum( ( m[,1] == 1 ) & (m[,2] == 1) )
          r_i <- sum( ( m[,1] == 1 ) & (m[,2] == 0) )
          r_j <- sum( ( m[,1] == 0 ) & (m[,2] == 1) )
          cScore <- ( r_i - S_ij ) * ( r_j - S_ij )
          
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