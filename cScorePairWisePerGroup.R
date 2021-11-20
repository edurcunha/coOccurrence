################################################################################
#  cScorePairWisePerGroup
################################################################################

#  FUNCTION DESCRIPTION: This function calculates the c-score index for pairs 
#                        of species in sets of sample collections

#  ARGUMENTS
# m: A matrix of species by samples.
# time: A vector of labels indicating the time period of each sample.
# site: A vector of labels indicating the site of each sample.
# spp.group: A vector of labels indicating the category of each species.
# occurrence.tresh: A number indicating the threshold of occurrence to 
#                   include species in the c-score calculation.
# standardize: Logical. If TRUE, c-score is standardize to number of potential
#              combinations between samples.


cScorePairWisePerGroup <- function(m, time, site, spp.group, 
                                   occurrence.tresh, ...) {
  

  extra.arg <- list(...)
  
  if( "standardize" %in% names(extra.arg) ) {
    
    standardize <- extra.arg[["standardize"]]
    
  } else {
    
    standardize <- FALSE 
    
  }
  
  dataset <- matrix(0, 0, 5, dimnames = list( NULL, 
    c('Site', 'Time', 'Spp.Groups', 'Spp.names', 'c-score')) )
  
  site.lab <- levels( as.factor(site) )
  
  time.lab <- levels( as.factor(time) )
  
  spp.group.lab <- levels( as.factor(spp.group) )
  
  for(i in seq_along(site.lab) ) {
    
    for(j in seq_along(time.lab) ) {
      
      m.tmp <- m[site == site.lab[i] & time == time.lab[j], ]
      
      spp.include <- colSums(m.tmp) >= occurrence.tresh
      
      m.tmp <- m.tmp[ , spp.include, drop=F]      
      
      if( ncol(m.tmp) > 1 & nrow(m.tmp) > 1 ) {
        
        n.spp <- ncol(m.tmp)
        
        spp.names <- colnames(m.tmp)
        
        spp.group.tmp <- spp.group[spp.include]  
        
        c.score.pairwise.tmp <- c( cScorePairWise(m.tmp,
                                                  standardize = standardize) )
        
        index.pos.tmp <- seq_along(c.score.pairwise.tmp)
        
        valid.index.tmp <- !is.na(c.score.pairwise.tmp)
        
        c.score.vals <- c.score.pairwise.tmp[valid.index.tmp]
        
        sp.pos1 <- ( ( index.pos.tmp[valid.index.tmp] - 1 ) %/% n.spp ) + 1
        
        sp.pos2 <- ( ( index.pos.tmp[valid.index.tmp] - 1 ) %% n.spp ) + 1
        
        sp.pos <- cbind(sp.pos1, sp.pos2)
        
        tags.tmp <- sapply(1:nrow(sp.pos), function(x) {
          
          spp.order <- order( c( spp.group.tmp[sp.pos[x,1]], 
                             spp.group.tmp[sp.pos[x,2]] ) )
          
          tags.tmp <- c( 
            paste(spp.group.tmp[sp.pos[x,][spp.order]], 
                               collapse = ' vs '),
            paste(spp.names[sp.pos[x,][spp.order]], 
                  collapse = ' vs ') )
          
          return( tags.tmp )

        }, simplify = T )
        
        tags.tmp <- t(tags.tmp)
        
        dataset.tmp <- cbind(
              rep( site.lab[i], nrow(sp.pos) ), 
              rep( time.lab[j], nrow(sp.pos) ),
              tags.tmp,
              c.score.vals)
        
        dataset <- rbind(dataset, dataset.tmp)
        
        }
      
    } 
    
  }
  
  dataset <- as.data.frame(dataset)
  
  return(dataset)
  
}