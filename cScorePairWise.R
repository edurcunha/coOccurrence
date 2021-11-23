################################################################################
#  cScorePairWise
################################################################################

#  FUNCTION DESCRIPTION: This function calculates the c-score index for species  
#                        pairs in a in a collection of samples.

#  ARGUMENTS
# m: A matrix of species by samples.
# standardize: Logical. If TRUE, c-score is standardized to the number 
#              of potential combinations between samples.
# null.model: Logical. If TRUE, the mean and sd of the c-score null
#             model is provided in the output.
# rand: A number indicating the number of random collections of samples to 
#       generate.

cScorePairWise <- function(m, ...) {
  
  extra.arg <- list(
    standardize = FALSE,
    null.model = FALSE,
    rand = 999
  )
  
  ellipsis <- list(...)
  
  arg.replace <- names(extra.arg) %in% names(ellipsis)
  
  extra.arg[arg.replace] <- ellipsis[names(extra.arg)[arg.replace] ] 
  
  
  if( is.matrix(m) ) {
    
    if( ncol(m) >= 2 ) {
       
      if( all.equal( sort( unique( c(m, 0, 1) ) ), c(0, 1) ) ){
        
        if( sum( rowSums(m) > 0 ) >= 2 ) {
          
          n.spp <- ncol(m)
          
          spp.names <- colnames(m)
          
          cScoreMatrix <- matrix(NA, n.spp, n.spp, 
                                 dimnames = list(spp.names, spp.names))
          
          if( extra.arg[['null.model']] ) {
            
            NullMatrixMean <- matrix(NA, n.spp, n.spp, 
                                   dimnames = list(spp.names, spp.names))
            
            NullMatrixSD <- matrix(NA, n.spp, n.spp, 
                                   dimnames = list(spp.names, spp.names))
            
          }
          
          for(i in 1:(n.spp-1)) {
            
            for(j in (i + 1):n.spp){
              
              m.tmp <- m[ , c(i, j)]
              
              c.score.tmp <- cScore(m.tmp, 
                                    standardize = extra.arg[['standardize']],
                                    null.model = extra.arg[['null.model']],
                                    rand = extra.arg[['rand']])
              
              if( extra.arg[['null.model']] ) {
                
                cScoreMatrix[j,i] <- c.score.tmp$c.score
                
                NullMatrixMean[j,i] <- mean(c.score.tmp$null.dist)
                
                NullMatrixSD[j,i] <- sd(c.score.tmp$null.dist)
                
              } else {
                
                cScoreMatrix[j,i] <- c.score.tmp
                
              }
              
            } 
            
          }
          
          if( extra.arg[['null.model']] ) {
            
            cScoreMatrix <- list(c.scores = cScoreMatrix,
                                 null.dist.mean = NullMatrixMean,
                                 null.dist.sd = NullMatrixSD)
            
          }
          
          return(cScoreMatrix)
          
        }        
        
      }
      
    }
    
  }
  
  stop("'m' is not a species (columns) by sites (rows) presence-absence matrix. Please, check if values are only 1 and 0, and if the 'm' object is a matrix")
  
}

# REFERENCE
# Stone. L. and A. Roberts. 1990. The checkerboard score and species 
# distributions. Oecologia 85: 74-79




