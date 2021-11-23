################################################################################
#  cScore
################################################################################

#  FUNCTION DESCRIPTION: This function calculates the c-score index for a pair 
#                        of species in a collection of samples.

#  ARGUMENTS
# m: A matrix of two species (columns) by n samples (rows).
# standardize: Logical. If TRUE, c-score is standardized to the number 
#              of potential combinations between samples.
# null.model: Logical. If TRUE, the mean and sd of the c-score null
#                model is provided in the output.
# rand: A number indicating the number of random collections of samples to 
#       generate.

cScore <- function(m, ...) {
    
  extra.arg <- list(
    standardize = FALSE,
    null.model = FALSE,
    rand = 999
  )
  
  ellipsis <- list(...)
  
  arg.replace <- names(extra.arg) %in% names(ellipsis)
  
  extra.arg[arg.replace] <- ellipsis[names(extra.arg)[arg.replace] ] 

  
  if( is.matrix(m) ) {
    
    if( ncol(m) == 2 ) {
      
      if( all.equal( sort( unique( c(m, 0, 1) ) ), c(0, 1) ) ){

        if( sum( rowSums(m) > 0 ) >= 1 ) {
          
          S_ij <- sum( ( m[,1] == 1 ) & (m[,2] == 1) )
          r_i <- sum( m[,1] )
          r_j <- sum( m[,2] )
          cScore <- ( r_i - S_ij ) * ( r_j - S_ij )
          
          if( extra.arg[['standardize']] ) {
            
            cScore <- cScore / ( ( ( nrow(m)^ 2 ) / 2 ) - (nrow(m) / 2) )
            
          }
          
          if( extra.arg[['null.model']] ){
            
            null.model <- cScoreRand(m, standardize = extra.arg[['standardize']],
                       rand = extra.arg[['rand']] )
            
            cScore <- list(c.score = cScore,
                           null.dist = null.model)
            
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