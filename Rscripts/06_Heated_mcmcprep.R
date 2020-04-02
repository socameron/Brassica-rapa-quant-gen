#library(individual)
library(aster)
library(compiler)
library(mcmc)
R.version.string
#packageVersion("individual")
packageVersion("aster")
packageVersion("mcmc")

set.seed(42)

load("Routput/heated_priors.RData")
## new minimal function included here to avoid the package stuff
## A tweaking of minimal to allow for a general pedigree

minimalnew <- function(random, alpha, sigma, pred, fam,
    y, root, modmat, fit, ind, deriv = 0,
    famlist = fam.default(), origin) {

    stopifnot(is.numeric(pred))
    stopifnot(length(pred) > 0)
    stopifnot(pred %in% 0:length(pred))
    stopifnot(pred < seq(along = pred))

    stopifnot(is.numeric(deriv))
    stopifnot(length(deriv) == 1)
    stopifnot(deriv %in% 0:2)

    stopifnot(is.list(famlist))
    stopifnot(sapply(famlist, function(x) inherits(x, "astfam")))

    stopifnot(is.numeric(fam))
    stopifnot(length(fam) == length(pred))
    stopifnot(fam %in% seq(along = famlist))

    stopifnot(is.numeric(y))
    stopifnot(is.finite(y))
    stopifnot(length(y) %% length(pred) == 0)

    stopifnot(is.numeric(root))
    stopifnot(is.finite(root))
    stopifnot(length(root) == length(y))

    if (! missing(origin)) {
        stopifnot(is.numeric(origin))
        stopifnot(is.finite(origin))
        stopifnot(length(origin) == length(y))
    }

    stopifnot(is.numeric(alpha))
    stopifnot(length(alpha) > 0)
    stopifnot(is.finite(alpha))

# this is the aster modmat
    stopifnot(is.numeric(modmat))
    stopifnot(is.finite(modmat))
    stopifnot(length(modmat) == length(y) * length(alpha))

    stopifnot(is.numeric(sigma))
#   stopifnot(length(sigma) == 3)  we want sigma to be a scalar here.
    stopifnot(length(sigma) == 1)
    stopifnot(is.finite(sigma))

    stopifnot(is.numeric(fit))
    stopifnot(is.finite(fit))
    stopifnot(length(fit) == length(y))
    stopifnot(fit %in% 0:1)

    stopifnot(is.numeric(ind))
    stopifnot(length(ind) == length(y))
    stopifnot(ind == round(ind))
    stopifnot(ind > 0)
    
# all the sire and dam stuff is here removed

    stopifnot(is.numeric(random))
    stopifnot(is.finite(random))
    stopifnot(length(random)==length(unique(ind)))

    # random    checked
    # alpha     checked
    # sigma     checked
    # pred      checked
    # fam       checked
    # y         checked
    # root      checked
    # modmat    checked
    # origin    checked
    # fit       checked
    # ind       checked
    # sire      checked
    # dam       checked
    # deriv     checked
    # famlist   checked

    nind <- max(ind)

    random.ind <- random[ind] * fit

    modmat <- cbind(modmat, random.ind) # just one extra column
    dimnames(modmat) <- NULL

    parm <- c(alpha, sigma)

    if (missing(origin)) {
        mout <- mlogl(parm, pred, fam, y, root, modmat, deriv,
            famlist = famlist)
    } else {
        mout <- mlogl(parm, pred, fam, y, root, modmat, deriv,
            famlist = famlist, origin = origin)
    }
    mout <- lapply(mout, function(x) (- x))
    return(mout)
}

minimalnew <- cmpfun(minimalnew)

## lupfun = log un-normalized posterior function
lupfun <- function(x) {
    stopifnot(is.numeric(x))
    stopifnot(is.finite(x))
    stopifnot(length(x) == nped + ncol(modmat) + 1)
    u <- x[1:nped] # the random effects for the whole pedigree
    psi <- x[-(1:nped)] # the prarmeters
    alpha <- psi[1:ncol(modmat)] # the fixed effect parameters
    sigma <- psi[-(1:ncol(modmat))] # sigma
    tau <- c(sigma/sqrt(2))  # ind, sire, dam components
    if (any(sigma < 0)) 
        return(-Inf)

### here is code that takes the N(0,1) vector u and makes it into a vector of breeding values
### first unnormalize u
    u <- u * sigma/sqrt(2)
    u[1:maxP_F2] <- u[1:maxP_F2]*sqrt(2)
### create vector of breeding values
    bv <- u
    for (a in (nf+1):maxP_F2) {
     bv[a] <- u[a] + (bv[pedigree[a,2]] + bv[pedigree[a,3]])/2
    }  ## takes care of intermediate unobserved generations
    bv[(maxP_F2+1):nped] <- bv[(maxP_F2+1):nped] + (bv[pedigree[(maxP_F2+1):nped,2]] + 
    bv[pedigree[(maxP_F2+1):nped,3]])/2
### this requires the individuals to be in generational order

# just the observed goes to minimalnew

    bv <- bv[-(1:maxP_F2)]

###

    Likelihood <- minimalnew(bv, alpha, tau, pred, fam, resp, root, modmat, fit, ind)$value
    Likelihood <- Likelihood - sum(x[1:nped]^2) / 2  # here random is N(0,1). Do this on outside
    # FE priors
    Priors <- dlogis(x = alpha, location = hyperparameters$fixed$mu, scale = hyperparameters$fixed$sigma, 
        log = T)
    # RE priors
    Priors <- c(Priors, lp_sigma(sigma = sigma, lambda = hyperparameters$random$lambda))
    
    Posterior <- Likelihood + sum(Priors)
    return(Posterior)
    print(Posterior)
}
lupfun <- cmpfun(lupfun)

# redefine the dimensions of the modmat once and for all

nnode <- length(pred)
modmat <- matrix(as.vector(modmat), nrow = nind * nnode, ncol = length(alpha))
nf <- sum(is.na(pedigree[,2]))  # number of founders

state.start <- c(rep(0, nind + maxP_F2), alpha, sigma) # maxP_F2 is the part of the pedigree above observed.

ntotal <- length(state.start)
ntotal
is.parm <- seq(1, ntotal) > nind + maxP_F2

mout <- metrop(lupfun, state.start, nspac = 2e4, nbatch = 1, scale = 0.0018, outfun = is.parm)
mout$accept

save.image("Routput/heated_mcmcprep01.RData")













#Archive

mout <- metrop(lupfun, state.start, blen = 2e4, nbatch = 1000, scale = 0.0018, outfun = is.parm)
mout$accept

save.image("Routput/heated_mcmcprep02.RData")

#For the new 03 run
load("Routput/heated_mcmcprep02.RData")

mout <- metrop(mout, nspac = 2e4, blen = 1, nbatch = 1000, scale = 0.0018)
mout$accept

# Started at around 9:30pm Eastern 2020-03-13, Completed 10:30pm 2020-03-18
save.image("Routput/heated_mcmcprep03.RData")
