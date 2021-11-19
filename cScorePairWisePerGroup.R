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
                               list(c(spp.group.lab[i], spp.group.lab[j])))
      
    } 
    
  }
  
  
  for(i in seq_along(site.lab) ) {
    
    for(j in seq_along(time.lab) ) {
      
      m.tmp <- m[site == site.lab[i] && time == time.lab[j], ]
      
      spp.include <- colSums(m.tmp) >= occurrence.tresh
      
      m.tmp <- m.tmp[ , spp.include]
      
      spp.group.tmp <- spp.group[spp.include]
      
      
      
      
      for(k in  seq_along(spp.group.contrasts) ) {
        
      }          
    }    
  }
  
  
  
  
  
  
}