#Set working directory. This is for Mason's computer, you will need to change this to where you saved
# the required data/files

setwd("C:/Users/Mason Kulbaba/Dropbox/git/Miksha_bee_aster")

#load data - this is the revised data (Nov 10, 2019) from Ron

dat<- read.csv("data/BumbleBeesAsterSuccessDensity2.csv")


##############################################################################
#
#First some preliminary checks that will potentially save us some troubleshooting later.

#What is the class of each column/variable?

sapply(dat, class)

#OK, look fine except for class of "HB.Host" that describes if a honeybee hive was present (1) or not (0)
# Need to make this a 'factor' so that aster correctly treats this as a 2-level treatment and not a
# continuous variable

dat$HB.Host<- as.factor(dat$HB.Host)


#check to make sure HB.Host is now a 'factor' variable

sapply(dat, class)# good, is it a factor

#############################################################################

# Prepare data for aster analysis

names(dat)

#designate sequence of response variables -> these represent variables in graphical model

#4 node graphical model:

# 1) Fluff - any evidence of a visit? (0/1)

# 2) Pupae.1 - did nest produce at least 1 cocoon (0/1)?

# 3) Pupae.8 - did nest reach cell cup/honey stage (0/1)?

# 4) Success - number of cups/cocoons produced by nest (1-5, integer)


vars<- c("Fluff","Pupae.1","Pupae.8", "Success")


#reshape data so that all response variables are located in a single vector in a new data
#set called "redata"
redata <- reshape(dat, varying = list(vars), direction = "long",timevar = "varb", times = as.factor(vars), v.names = "resp")


#Designation of "fitness" variable: success
fit <- grepl("Success", as.character(redata$varb))
fit<- as.numeric(fit)

#add this designation to the reshaped data, 'redata'
redata$fit <- fit

#check
with(redata, sort(unique(as.character(varb)[fit == 0])))
with(redata, sort(unique(as.character(varb)[fit == 1]))) #all good!

#add a variable "root" to redata files, where value is 1. This represents the 'planting' of a 
#domocile. This represents the variable 'EmptyDomicile' in the original data 

redata<- data.frame(redata, root=1)

#load aster package
library(aster)

#set graphical mode and dist. for fitness nodes (preds)
pred<- c(0,1,2,3)
fam<- c(1,1, 1,2)

#describe dist. of preds.
sapply(fam.default(), as.character)[fam]

#fixed effect model
aout<- aster(resp~varb, pred, fam, varb, id, root, data=redata)

#show summary of first aster model...although we can't learn much from this
summary(aout, show.graph = T)

#add HB.Host treatment
aout.2<- aster(resp~varb + fit:(HBDensity), pred, fam, varb, id, root, data=redata)

summary(aout.2, show.graph = T)

#liklihood ratio test of two models

anova(aout, aout.2)

