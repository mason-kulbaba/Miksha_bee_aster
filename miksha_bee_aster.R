#Set working directory. This is for Mason's computer, you will need to change this to where you saved
# the required data/files

setwd("C:/Users/Mason Kulbaba/Dropbox/git/Miksha_bee_aster")

#load data - this is the revised data (Nov 10, 2019) from Ron

dat<- read.csv("data/BBSuccess_Aster.csv")


##############################################################################
#
#First some preliminary checks that will potentially save us some troubleshooting later.

#What is the class of each column/variable?

sapply(dat, class)

#OK, looks fine except for class of "HB.Host" that describes if a honeybee hive was present (1) or not (0)
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

# 2) Cells - did nest produce any cells (0/1)?

# 3) Cell.Count - how many cells were produced?



vars<- c("Fluff","Cells", "Cell.Count")


#reshape data so that all response variables are located in a single vector in a new data
#set called "redata"
redata <- reshape(dat, varying = list(vars), direction = "long",timevar = "varb", times = as.factor(vars), v.names = "resp")


#Designation of "fitness" variable: success
fit <- grepl("Cell.Count", as.character(redata$varb))
fit<- as.numeric(fit)

#add this designation to the reshaped data, 'redata'
redata$fit <- fit

#check
with(redata, sort(unique(as.character(varb)[fit == 0])))
with(redata, sort(unique(as.character(varb)[fit == 1]))) #all good!

#add a variable "root" to redata file, where value is 1. This represents the 'planting' of a 
#domocile. This represents the variable 'EmptyDomicile' in the original data 

redata<- data.frame(redata, root=1)


#The distribution for one node (Success) is the only debatable distribution (others are clearly Bernoulli)

library(MASS)

count<- dat$Cell.Count

fit.1<- fitdistr(count, "normal")
fit.2<- fitdistr(count, "negative binomial")#size: 0.042093032 
fit.3<- fitdistr(count, "poisson")

AIC(fit.1, fit.2, fit.3)

fit.2 #Just for now (we need to discuss this), set Success as neg.binomial

#Therefore, we need to set a custom family list of statistical distributions for aster
#This wilil be easy, as there are only two unique distribution

famlist<- list(fam.bernoulli(), fam.negative.binomial(0.042093032 ))


#load aster package
library(aster)

#set graphical mode and dist. for fitness nodes (preds)
pred<- c(0,1,2)
fam<- c(1,1,2)

#describe dist. of preds.
sapply(fam.default(), as.character)[fam] #this will still say bernoulli & Poisson, but the designation
                                          # of custom family list above will override this in aster

#fixed effect model
aout<- aster(resp~varb, pred, fam, varb, id, root, data=redata, famlist= famlist)

#show summary of first aster model...although we can't learn much from this
#NOTE: that success is now treated as neg.binomial with size 0.11105926
summary(aout, show.graph = T)

#add presence of hive/host
aout2<- aster(resp~varb + fit:HB.Host, pred, fam, varb, id, root, data=redata, famlist=famlist)

summary(aout2, show.graph = T)

#add bee density
aout3<- aster(resp~varb + fit:HB.Host + fit:HBDensity, pred, fam, varb, id, root, data=redata, famlist=famlist)

summary(aout3, show.graph = T)

#liklihood ratio test
anova(aout, aout2, aout3)

#only bee density
aout4<- aster(resp~varb + fit:HBDensity, pred, fam, varb, id, root, data=redata, famlist=famlist)

summary(aout4, show.graph = T)


anova(aout, aout4)


############################################
# For prelimiary resutls, lets generate estiamtes of "Success" with standar errors for the two levels
#of the treatment "HB.Host"

#First, make a design matrix to hold the appropriate estiamtes. Call this matrix 'fred' and make
# all variable equal to 1 for now, just so we can make an actual matrix (the value doesn't matter).
# However, make sure both levels of 'HB.Host' are represented.

fred <- data.frame(HB.Host=levels(dat$HB.Host), Fluff=1, Cells=1, Cell.Count=1, root = 1)

# reshape the design matrix just as we did for the actual data
renewdata <- reshape(fred, varying = list(vars),
                     direction = "long", timevar = "varb",
                     times = as.factor(vars), v.names = "resp")

# make character string from "varb" of renewdata, without actual values (i.e., the layers of varb in renewdata)
layer<- gsub("[0-9]", "", as.character(renewdata$varb))

# add this layer to renewdata
renewdata<- data.frame(renewdata, layer= layer)

# add 'Success' in new layer column of renewdata as numeric, called fit
fit<- as.numeric(layer=="Cell.Count")

# add fit to renewdata
renewdata<- data.frame(renewdata, fit = fit)


#Generate fintess estimates and standard errors for HB.Host treatment
nHost<- nrow(fred)#all data has same number of blocks so any file will do
nnode<- length(vars)
amat<- array(0, c(nHost, nnode, nHost))
dim(amat)# makes a 2 x 4 x 2 matrix (2 Host types and 4 nodes of graphicla model)

#only want prediction for k'th individual that contribute to expected
#fitness, and want to add only success entries

foo<- grepl("Cell.Count", vars)
for(k in 1:nHost)
  amat[k, foo, k]<- 1

#check
foo #yes, only last node is "true"; corresponds to "Success"

#generate predicted valuses using region-specific aout object, with renewdata, and amat format
pout.amat<- predict(aout2, newdata= renewdata, varvar= varb,
                        idvar= id, root = root, se.fit=TRUE, amat = amat)

#Combine estimates wiht standard errors
ests<- cbind(pout.amat$fit, pout.amat$se.fit)

#assign row names: 0= no colony, 1= presence of colony 
rownames(ests)<- as.character(fred$HB.Host)

#assign column names
colnames(ests)<- c("Success", "SE")

#show estimates
ests



################################################################################
# Below is later work that should be discussed as a group...that is to say that
# Mason doesn't know enough about the project to do the things that he is doing...but it 
# hasn't seemed to have stoped him. 

#add individual site as a random effect. Note: we have to use the 'reaster' function to perform
# an aster analsyis with random effects. However, we can still still test if variables explain 
# significantly more variation in models between fixed and random effects aster modesl the usual way.
aout.2<- reaster(resp~varb,list(Host=~0 + fit:Host), 
               + pred, fam, varb, id, root, data=redata)

summary(aout.2, show.graph = T)

#liklihood ratio test of two models

anova(aout, aout.2)# presence of 'Host' explains significantly more variation than the model without




