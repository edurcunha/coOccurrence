################################################################################
#  cScorePairWise
################################################################################

#  FUNCTION DESCRIPTION: This function calculates the c-score index for species  
#                        pairs in a in a collection of samples.

#  ARGUMENTS
# m: A matrix of species by samples.

cScorePairWise <- function(m) {
  
  if( is.matrix(m) ) {
    
    if( ncol(m) == 2 ) {
       
      if( all.equal( unique( c(m) ), c(0, 1) ) ){
        
        if( sum( rowSums(m) > 0 ) >= 2 ) {
          
          n.spp <- ncol(m)
          
          spp.names <- colnames(m)
          
          cScoreMatrix <- matrix(NA, n.spp, n.spp, 
                                 dimnames = list(spp.names, spp.names))
          
          for(i in 1:(n.spp-1)) {
            
            for(j in (i + 1):n.spp){
              
              m.tmp <- m[ , c(i, j)]
              
              m.tmp <- m.tmp[rowSums(m.tmp) > 0, ]
              
              cScoreMatrix[j,i] <- cScore(m.tmp)
              
            } 
            
          }
          
          return(cScoreMatrix)
          
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




