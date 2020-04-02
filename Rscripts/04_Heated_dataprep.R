# ------------------------------------- Load aster package and aster object
library(aster)
R.version.string
packageVersion("aster")
load("./Routput/heated_aster0.RData")
# what have we got

ls()

# ------------------------------------- Invesitgate aster object
class(aout)  # aster object
names(aout)
aout$coefficients  # estimates of fixed effects
class(aout$modmat)  # model matrix for fixed effects (M)

class(aout$x)  # aster graph data
dim(aout$x)  # 5 aster graph nodes
nind <- dim(aout$x)[1]
nnode <- dim(aout$x)[2]
nfix <- length(aout$coefficients)

# now have to construct a random vector with individuals (observed) and parents
# (unobserved) comprising founders and two more generations (F1 and F2).
# F3 is the observed generation (ind).
ind <- rep(1:nind, times = nnode)  # 1,2, ... repeated 5 times

# create fitness
fit <- as.numeric(!seq(along = aout$pred) %in% aout$pred)  # 1 if terminal, 0 ow
fit <- rep(fit, each = nind)  # indicator corresponds to fitness var totalseeds
head(fit)
length(fit)

# read in the F1-F2 pedigree
fpeddata <- read.csv("Rdata/heatarrays_animal_pedigree.csv")
inddata <- read.csv("Rdata/heated_aster.csv")  ## this is for the H treatment only
# want a pedigree that corresponds to the ind vector, the individuals as they
# appear in the aster dataset.
maxP_F2 <- max(c(inddata$matID,inddata$patID)) # this is the max ind for the nonobserved.
indped <- inddata[,c(3,5,6)]
foundped <- fpeddata[1:maxP_F2,]
pedigree <- rbind(foundped,indped)
nped <- length(pedigree$animal) 
nobs <- nped - maxP_F2   
nobs
# this had better make sense!
# try it out

# get parameters from aster obj fixed effect parameters
alpha <- aout$coefficients
# random effect parameter
sigma <- 0.01  # as a start
names(sigma) <- c("additive genetic sd")
alpha
sigma

fam <- aout$fam  # ber ber pois ber pois
identical(aout$famlist, fam.default())

resp <- redata$resp
length(resp)
root <- aout$root
dim(root)
modmat <- aout$modmat
dim(modmat)
# make it into a matrix again
modmat <- matrix(as.vector(modmat),ncol=length(alpha))
dim(modmat)
is.null(aout$origin)

save.image("Routput/heated_dataprep.RData")

