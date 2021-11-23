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
# standardize: Logical. If TRUE, c-score is standardized to the number 
#              of potential combinations between samples.
# null.model: Logical. If TRUE, the mean and sd of the c-score null
#             model is provided in the output.
# rand: A number indicating the number of random collections of samples to 
#       generate.

cScorePairWisePerGroup <- function(m, time, site, spp.group, 
                                   occurrence.tresh, ...) {
  

  extra.arg <- list(
    standardize = FALSE,
    null.model = FALSE,
    rand = 999
  )
  
  ellipsis <- list(...)
  
  arg.replace <- names(extra.arg) %in% names(ellipsis)
  
  extra.arg[arg.replace] <- ellipsis[names(extra.arg)[arg.replace] ] 
  
  
  
  if( extra.arg[['null.model']] == FALSE ) {
    
    dataset <- matrix(0, 0, 5, dimnames = list( NULL, 
                  c('Site', 'Time', 'Sp.group', 'Sp.name', 'c.score')) )
    
  } else {
    
    dataset <- matrix(0, 0, 8, dimnames = list( NULL, 
                  c('Site', 'Time', 'Sp.group', 'Sp.name', 'c.score',
                    'null.c.score.mean', 'null.c.score.sd', 
                    'c.score.scaled')) )
    
  } 

  
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
        
        c.score.pairwise.tmp <- cScorePairWise(m.tmp,
                                    standardize = extra.arg[['standardize']], 
                                    null.model = extra.arg[['null.model']],
                                    rand = extra.arg[['rand']] )          
        
        if( extra.arg[['null.model']] ) {
          
          null.pairwise.mean.tmp <- c( c.score.pairwise.tmp$null.dist.mean )
          null.pairwise.sd.tmp <- c( c.score.pairwise.tmp$null.dist.sd )
          c.score.pairwise.tmp <- c( c.score.pairwise.tmp$c.score )
          
        } else {
          
          c.score.pairwise.tmp <- c( c.score.pairwise.tmp )
          
        }
        
        
        index.pos.tmp <- seq_along(c.score.pairwise.tmp)
        
        valid.index.tmp <- !is.na(c.score.pairwise.tmp)
        
        c.score.vals <- c.score.pairwise.tmp[valid.index.tmp]
        
        if( extra.arg[['null.model']] ) {
          
          null.mean.vals <- null.pairwise.mean.tmp[valid.index.tmp]
          
          null.sd.vals <- null.pairwise.sd.tmp[valid.index.tmp]
          
          scaled.c.score.vals <- ( c.score.vals - null.mean.vals ) / 
                                                          null.sd.vals
          
        } 
        
        sp.pos1 <- ( ( index.pos.tmp[valid.index.tmp] - 1 ) %/% n.spp ) + 1
        
        sp.pos2 <- ( ( index.pos.tmp[valid.index.tmp] - 1 ) %% n.spp ) + 1
        
        sp.pos <- cbind(sp.pos1, sp.pos2)
        
        tags.tmp <- sapply( 1:nrow(sp.pos), function(x) {
          
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
        
        if( extra.arg[['null.model']] ) {
          
          dataset.tmp <- cbind(
            rep( site.lab[i], nrow(sp.pos) ), 
            rep( time.lab[j], nrow(sp.pos) ),
            tags.tmp,
            c.score.vals,
            null.mean.vals,
            null.sd.vals,
            scaled.c.score.vals)
          
        } else {
          
          dataset.tmp <- cbind(
            rep( site.lab[i], nrow(sp.pos) ), 
            rep( time.lab[j], nrow(sp.pos) ),
            tags.tmp,
            c.score.vals)
          
        }
        
        dataset <- rbind(dataset, dataset.tmp)
        
        }
      
    } 
    
  }
  
  dataset <- as.data.frame(dataset)
  
  numeric.labs <- c('c.score', 'null.c.score.mean', 'null.c.score.sd', 
    'c.score.scaled')

  factor.labs <- c('Site', 'Time', 'Sp.group', 'Sp.name')
    
  for(i in seq_along(numeric.labs) ) {
    
    if(numeric.labs[i] %in% colnames(dataset) ) {
      
      dataset[ , numeric.labs[i] ] <- suppressWarnings( 
        as.numeric(dataset[ , numeric.labs[i] ] ) )
      
    }
    
  }
  
  for(i in seq_along(factor.labs) ) {
    
    if(factor.labs[i] %in% colnames(dataset) ) {
      
      dataset[ , factor.labs[i] ] <- as.factor(dataset[ , factor.labs[i] ] )
      
    }
    
  } 
  

  return(dataset)
  
}
