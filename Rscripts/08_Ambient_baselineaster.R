#This code is adapted from the revised data to estimate mean fitness and Va(w) in Grey Cloud Dunes
#read data
data.C<-read.csv("Rdata/ambient_aster.csv")

##########################################################################
##SECTION 1a: Aster to estimate mean expected fitness                   ##
##########################################################################

#load aster package
library(aster)

#####################################################
##################Aster Modeling#####################
#####################################################

#make response data a vector
vars<- c("germ","flower","flwr_clstr","seed","seed_pods")

#reshape data so that all response variables are located in a single vector in a new data
#set called "redata

redata <- reshape(data.C, varying = list(vars), direction = "long",timevar = "varb", times = as.factor(vars), v.names = "resp")

#Designation of fitness variable
#to estimate mean fitness, we need to identify (via a new variable) which node (seed_pods) will be used
#to calcualte fitness.  In this case, we want to identiy seed_pods as that variable.

fit <- grepl("seed_pods", as.character(redata$varb)) 

#make fit numeric 0 = false, 1 = true
fit<- as.numeric(fit)

#add fit to redata
redata$fit <- fit
#data checks
with(redata, sort(unique(as.character(varb)[fit == 0])))
with(redata, sort(unique(as.character(varb)[fit == 1])))

#add a variable "root" to redata, where value is 1
redata<- data.frame(redata, root=1)

#make plot a factor
redata$plot <- as.factor(redata$plot)

#set graphical mode and dist. for preds.
pred<- c(0,1,2,3,4)
fam<- c(1,1,2,1,2)

#describe dist. of preds.
sapply(fam.default(), as.character)[fam]

#basic first model with "block"
aout<- aster(resp~varb + fit:plot, pred, fam, varb, id, root, data=redata)

#save model aout as "aster0H.RData"
save.image("Routput/ambient_aster0.RData")

#show summary of model
summary(aout, show.graph=TRUE)


#generate MLE of saturated model mean value parameter vector: mu
pout<- predict.aster(aout, se.fit=TRUE)

