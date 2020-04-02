
setwd("D:/Users/Cameron/Google Drive/University of Toronto/Graduate/Methods/Aster Modelling/Aster Thesis 2019")

pdf(file = "Routput/ambient_results03.pdf")

R.version.string

load("Routput/ambient_mcmcprep03.RData")

batch <- mout$batch
nparm <- ncol(batch)
names.psi <- c(names(alpha), "sigma_A")

for (i in 1:nparm) plot(batch[ , i], ylab = names.psi[i],
    main = paste("batch means for", names.psi[i]))

for (i in 1:nparm) acf(batch[ , i], lag.max = 200,
    main = paste("batch means for", names.psi[i]))

dev.off()

save.image("Routput/ambient_results03.RData")
