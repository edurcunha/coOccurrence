################################################################################
#  cScorePairWisePerGroup
################################################################################

#  FUNCTION DESCRIPTION: This function calculates the c-score index for pairs 
#                        of species in sets of sample collections

#  ARGUMENTS
# m: 
# time: 
# site: 
# spp.group: 
# occurrence.tresh: 

cScorePairWisePerGroup <- function(m, time, site, spp.group, occurrence.tresh) {
  
  site.lab <- levels( as.factor(site) )
  
  time.lab <- levels( as.factor(time) )
  
  spp.group.lab <- levels( as.factor(spp.group) )
  
  
  spp.group.contrasts <- list()
  
  for(i in seq_along(spp.group.lab)) {
    
    for(j in i:length(spp.group.lab)) {
      
      spp.group.contrasts <- c(spp.group.contrasts, 
                               list( sort( c(spp.group.lab[i], 
                                             spp.group.lab[j]) ) ) )
      
    } 
    
  }
  
  for(i in seq_along(site.lab) ) {
    
    for(j in seq_along(time.lab) ) {
      
      m.tmp <- m[site == site.lab[i] && time == time.lab[j], ]
      
      spp.include <- colSums(m.tmp) >= occurrence.tresh
      
      samples.include <- rowSums(m.tmp) > 0
      
      m.tmp <- m.tmp[samples.include, spp.include, drop=F]      
      
      if( ncol(m.tmp) > 1 && nrow(m.tmp) > 1 ) {
        
        n.spp <- ncol(m.tmp)
        
        spp.group.tmp <- spp.group[spp.include]  
        
        c.score.pairwise.tmp <- c( cScorePairWise(m.tmp) )
        
        index.pos.tmp <- seq_along(c.score.pairwise.tmp)
        
        valid.index.tmp <- !is.na(c.score.pairwise.tmp)
        
        c.score.vals <- c.score.pairwise.tmp[valid.index.tmp]
        
        sp.pos1 <- index.pos.tmp[valid.index.tmp] %% n.spp
        
        sp.pos2 <- index.pos.tmp[valid.index.tmp] %/% n.spp
      
        }
      


      
      

      
    } 
    
  }
  
}