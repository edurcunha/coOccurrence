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

# Load the package
source('coOccurrence.R')

# Create a presence=absence matrix of species by samples
m <- matrix( rbinom(200, 1, 0.5), 20, 10)   

# Assign species names 
colnames(m) <- LETTERS[1:10]               

# Create a vector indicating the time period of the samples 
time <- rep( c("After", "Before"), each = 10)

# Create a vector indicating the site of the samples
site <- rep( rep( c("Parana", "Baia"), each = 5), 2 )

# Create a vector indicating the group the species belong 
spp.group <- rep( c("native", "non-native"), each = 5)

# Assign the minimum threshold of occurrence to include species in the c-score
# calculations
occurrence.tresh <- 1

# Assign the criterion for standardization of c-scores
standardize <- TRUE

# Assign the criterion for generating a null-model for the c-scores
null.model <- TRUE

# Assign the number of randomizations for creating the null-model 
rand = 9999

# Build the c-score dataset  
dataset <- cScorePairWisePerGroup(m = m, time = time, site = site, 
                    spp.group = spp.group, occurrence.tresh = occurrence.tresh, 
                    standardize = TRUE, null.model = TRUE, rand = 9999)

