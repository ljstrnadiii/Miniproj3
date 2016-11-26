# Boston House Data set Model Selection
# Methods: AIC, BIC, CV, and test/train

library(MASS)
source("Mini3.R")

data(Boston)
df = Boston

X = df[,-14]
Y = df$medv

#standardize predictors
X = data.frame(scale(X))


### AIC ###

#ridge AIC
ridgeaic = ridgeAIC(X,Y, 0:300)
print("ridge: lambda, AIC")
print(ridgeaic$minlambda)

#lasso AIC
lassoaic = lassoAIC(X,Y,seq(.00001,1,.01))
print("lasso: t, AIC")
print(lassoaic$minlambda)


#pcr AIC
pcraic = pcrAIC(X,Y, ncol(X))
print("pcr: numcomps, AIC")
print(pcraic$mincomps)

#pls AIC
plsaic = plsAIC(X,Y, ncol(X))
print("pls: numcomps, AIC")
print(plsaic$mincomps)


### BIC ###

#ridge BIC
ridgebic = ridgeBIC(X,Y, 0:300)
print("ridge: lambda, BIC")
print(ridgebic$minlambda)

#lasso BIC
lassobic = lassoBIC(X,Y, seq(.0001,1,.01))
print("lasso: t, BIC")
print(lassobic$minlambda)

#pcr BIC
pcrbic = pcrBIC(X,Y,ncol(X))
print("pcr: numcomps, BIC")
print(pcrbic$mincomps)

#pls BIC
plsbic = plsBIC(X,Y, ncol(X))
print("pls: numcomps, BIC")
print(plsbic$mincomps)


### CV ### (need to multiply gcv by N to get the Err estimate)
#ridge CV 
ridgecv = ridgeCV(X,Y, 0:300)
print("ridgeCV: lambda, Err estimate")
print(ridgecv$minlam)

#lasso CV
lassocv = lassoCV(X,Y, seq(.01,1,.01))
print("lassoCV: t, Err estimate")
print(lassocv$mint)

#pcr CV (no mincomp returned, refer to plot output)
pcrcv = pcrCV(X,Y, ncol(X))

#plsr CV (no mincomp returned, refer to plot output)
plsrcv = plsrCV(X,Y,ncol(X))
