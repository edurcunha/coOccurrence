################################################################################
#  coOccurrence_EXAMPLE
################################################################################

#  FUNCTION DESCRIPTION: This is a simple guided example to  calculates 
#                        the c-score index for pairs of species in 
#                        collections of samples 



# IMPORTANT NOTE:
# Copy the files of the GitHub repository 'coOccurrence' 
# ('https://github.com/edurcunha/coOccurrence') into your 
# R project folder.

# Loads the package
source('coOccurrence.R')

# Creates a presence=absence matrix of species by samples
m <- matrix( rbinom(200, 1, 0.5), 20, 10)   

# Assigns species names 
colnames(m) <- LETTERS[1:10]               

# Creates a vector indicating the time period of the samples 
time <- rep(c("After", "Before"), each = 10)

# Creates a vector indicating the site of the samples
site <- rep( rep(c("Parana", "Baia"), each = 5), 2 )

# Creates a vector indicating the group the species belong 
spp.group <- rep( c("native", "non-native"), each = 5)

# Assigns the minimum threshold of occurrence to include species in the c-score
# calculations
occurrence.tresh <- 1

# Builds the c-score dataset  
dataset <- cScorePairWisePerGroup(m, time, site, spp.group, occurrence.tresh, standardize = standardize)


