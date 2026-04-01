# 1. Load the data --------------------------------------------------------------
pilot <- read.csv("G16.pilot.data.csv")
dim(pilot)

# Convert character colums to factors
factor_cls <- c("location", "ecotr_id", "tree_id", "pear_id", "species")
pilot[factor_cls] <- lapply(pilot[factor_cls], factor)


# 2. Pilot model -------------------------------------------------------------

# Fit
m2 <- lm(quality_index ~ location + species, data = pilot)
summary(m2)

# Extract model summary
intercept <- 54.822
species_effect <- -6.194
location_effect <- -6.610
climate_effect <- 10.0
residual_sd <- 20.33

# Simulation --------------------------------------------------------------

# Parameters
alpha = 0.05
test <- "greater"
n_grid <- seq(3, 100, by=5)

power_vec <- numeric(length(n_grid))
for (j in 1:length(n_grid)){
    N <- n_grid[j]
    total_N <- 8 * N  # 2 species × 2 climate × 2 location = 8 groups
    species <- rep(c(0, 1), each = 4 * N, times = 1)
    climate <- rep(c(0, 1), each = 2 * N, times = 2)
    location <- rep(c(0, 1), each = N, times = 4)
    
    numberSimulation <- 1000
    pval <- numeric(numberSimulation)
    tval <- numeric(numberSimulation)
    set.seed(54678)
    for (i in 1:numberSimulation){
        
        response <- intercept + 
            species_effect * species + 
            climate_effect * climate + 
            location_effect * location +
            rnorm(total_N, mean = 0, sd = residual_sd)
        
        simData <- data.frame(response = response,
                              species = species,
                              location = location,
                              climate = climate)
        
        
        lm_fit <- summary(lm(response ~ climate + species + location, data = simData))
        
        pval[i] <- lm_fit$coeff["climate", "Pr(>|t|)"]
        tval[i] <- lm_fit$coeff["climate", "t value"]
        
        if (test == "greater" & tval[i] > 0){
            pval[i] <- pval[i]/2
        }
        if (test == "greater" & tval[i] < 0){
            pval[i] <- 1 - (pval[i]/2)
        }
        if (test == "less" & tval[i] < 0){
            pval[i] <- pval[i]/2
        }
        if (test == "less" & tval[i] > 0){
            pval[i] <- 1 - (pval[i]/2)
        }
    }
    
    ## TODO: add multiplicity adjustment
    
    power_vec[j] = sum(pval < alpha) / numberSimulation
}

plot(n_grid, power_vec, type = "b", 
     pch = 19, lwd = 1,
     ylim = c(0, 1))
abline(h = 0.8, col = "red", lty = 2, lwd = 2)

################################################################
# NEXT STEPS
################################################################

# 1. estimate variance-convariance matrix for random effects from pilot data
# 2. simulate data with random effect
# 3. use lmer(response ~ climate + species + location + (1|ecotron) + (1|tree), data = simData)
# 4. add multiplicity adjustment: Bonferroni, Holm or BH (FDR) - less conservative (more powerful)

################################################################
# QUESTIONS
################################################################

# 1. how to determine if random effect is significant? By ICC threshold?
# 2. how to define hierarchical structure between random effects? What type of variance-covariance matrix to use
# Option 1 (independent)
# [ var econtr 0]
# [ 0 var tree]

# Option 2 (full)
# [ var econt covarXXX]
# [ covarXXX var tree]

# 3. simulation approach vs single formula approach: which one works better? 
# Answer: simulation is more flexiable and hence universal for any study setup

# 4. how to account for interaction effect between climate and species during simulation?
# Assume some range of values interaction = seq(0, 10, by = 0.2)???




