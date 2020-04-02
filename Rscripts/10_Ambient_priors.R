library(aster)

load("Routput/ambient_dataprep.RData")

#-------------------------------------------
# PRIORS: Fixed Effects and Random Effects
hyperparameters <- list()

# ------------------------------------- Fixed Effects - logistic priors We are
# using the fixed effect estimates from the original run 
load("./Routput/ambient_aster0.RData")  # load estimates

## mu - location parameter in logistic dist
hyperparameters$fixed$mu <- summary(aout)$coefficients[, 1]


## sigma - scale parameter in logistic dist
hyperparameters$fixed$sigma <- sqrt(3) * summary(aout)$coefficients[, 2]/pi


# ------------------------------------- Random Effects - square root exponential
# priors exponential priors on variance imply root exponential priors on std dev
agv <- 0.2  # from Mason's GC raster run, adjusted to what seems to be happening here.
hyperparameters$random <- list(lambda = 1/agv) # here we have NO IDEA, so just use this.

hyperparameters
rm(aout, agv)


# -------------------------------------- log prior for sigma (sigma^2 ~
# exp(lambda))
lp_sigma <- function(sigma, lambda) {
    # calculate the log prior for sigma f(sig^2) = lambda exp(-lambda sig^2) f(sig) =
    # 2 lambda sig exp(-lambda sig^2) l(sig) = log(2lambda) + log(sig) - lambda sig^2
    # drop constant log(2lambda) Watch out for overflow/underflow.  If sig < 0 return
    # -Inf If sig large, we want -Inf, but R gives Inf - Inf = NaN
    out <- ifelse(sigma < 0, -Inf, log(sigma) - lambda * sigma^2)
    out <- ifelse(is.finite(out), out, -Inf)
}
save.image("Routput/ambient_priors.RData")
ls()
