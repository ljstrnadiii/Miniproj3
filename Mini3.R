#Basic collection of functions:
# 1. train and test
# 2. standardize
# 3. ridgeAIC
# 4. 

# TODO: work with plotting to output a png file 

library(MASS)
library(lasso2)
library(pls)

# train and test set
tandt <- function(dat) {
    # samples and partition into 2/3 train and 1/3 test
    I = sample(1:dim(dat)[1], round((2/3)*dim(dat)[1]))
    not.I = setdiff(1:dim(dat)[1],I)
    dat.train = dat[I,]
    dat.test = dat[not.I,]

    return(list("train"=dat.train,"test"=dat.test))
}

standardize <- function(dat.train, dat.test, response,test.response){
    # Standardize Predictors and use predictors stats on test standardization
    # dat.train, dat.test are predictors only and response/test.response
    # is the response of the training and test.

    standardized.predictors = data.frame(scale(dat.train))

    # standardize
    standardized.train = cbind(response,standardized.predictors)
    standardized.train = as.data.frame(standardized.train)

    # must apply standardization of train mean and var to test.
    mean_train = apply(dat.train,2,mean)
    var_train = apply(dat.train,2,var)
    standardized.test.predictors = 
      cbind(as.data.frame(t((t(dat.test)-mean_train)/sqrt(var_train))))
    standardized.test = cbind(standardized.test.predictors,response = test.response)

    return(list("train" = standardized.train,"test"= standardized.test))
}


#################################
#AIC Methods
ridgeAIC <- function(dat, response, lambdalist){
    # Ride regression with parameter lambda chosen using AIC
    # dat is dat df of predictors, response is response of model
    
    lambda.and.aic <- c()
    df = dat
    dat = dat[complete.cases(response),]
    df[["response"]] <- response
    df = df[complete.cases(response),]

    N <- nrow(df)
	
    for(i in lambdalist){
		lambda <- i
		X <- as.matrix(dat)
		W <- X %*% solve(t(X)%*%X + lambda * diag(dim(X)[2])) %*% t(X)
		d <- sum(diag(W))  # the trace of W
		ans <- lm.ridge(df$response ~ ., lambda = i, data=df)
		predicted.y <- ans$ym + as.matrix(dat) %*% coef(ans)[2:(dim(X)[2]+1)] 
        
        pe <- mean((predicted.y - df$response)^2)
		se2 <- sum((predicted.y - df$response)^2)/(N-d-1)

        lambda.and.aic <- rbind(lambda.and.aic,
                                  c(lambda, pe + 2*d/N*se2))
    }

    minlambda = lambda.and.aic[lambda.and.aic[,2]==min(lambda.and.aic[,2]),]
    png("LambdaRidgeAIC.png")         
    plot(lambda.and.aic, type="l", xlab=expression(lambda), ylab="AIC")
    dev.off()

    ans <- lm.ridge(df$response ~ ., lambda = minlambda[1], data=df)
    norm = sum(coef(ans)[2:diag(X)[2]]^2)

    return(list("model"=ans, "norm" = norm, "minlambda" = minlambda))
}


lassoAIC <- function(dat, response, lambdalist){
    # Lasso Regression using AIC to choose lambda
    
    lambda.and.aic <- c()
    df = dat
    df[["response"]] <- response
    df = df[complete.cases(df),]
    N = nrow(df)
    
    for(i in lambdalist){
		ans <- l1ce(df$response ~ ., data=df, bound=i, absolute.t=FALSE)
		predicted.y <- predict(ans, df)
		pe <- mean((predicted.y - df$response)^2)
		d <- sum(coef(ans) != 0)
		se2 <- pe*N/(N-d-1)
		aic <- pe + 2*d/N*se2
		lambda.and.aic <- rbind(lambda.and.aic, c(i,aic))
	}

	minlambda = lambda.and.aic[lambda.and.aic[,2]==min(lambda.and.aic[,2]),]
	
    png("LambdaLassoAIC.png")
    plot(lambda.and.aic, type="l", xlab=expression(t), ylab="AIC")
    dev.off()
	ans <- l1ce(df$response ~ ., data=df, bound=minlambda[1], absolute.t=FALSE)
	norm = sum(abs(coef(ans))[2:ncol(df)])
	
    return(list("model" = ans, "norm" = norm, "minlambda"=minlambda))

}


pcrAIC <- function(dat, response, numcomp){
    # Pricipal Component Regression with AIC to choose num components
    
    df = dat 
    df[["response"]] <- response
    df = df[complete.cases(df),]
    
    ans <- pcr(df$response ~ ., data=df, ncomp=numcomp)
	predicted.y <- predict(ans, df)
    
    d.and.aic <- c()
	N <- nrow(df)
	for(i in 1:numcomp){
		pe <- mean((predicted.y[,,i] - df$response)^2)
		d <- i
		se2 <- pe*N/(N-d-1)
		aic <- pe + 2*d/N*se2
		d.and.aic <- rbind(d.and.aic, c(d,aic))
	}

	mincomps = d.and.aic[d.and.aic[,2]==min(d.and.aic[,2]),]
    
    ans <- pcr(df$response ~ ., data=df, ncomp=mincomps[1])
    
    png("pcrAIC.png")
    plot(d.and.aic, type="l", xlab="d", ylab="AIC")
    dev.off()
    
    return(list("model"=ans, "mincomps" = mincomps))
}

plsAIC <- function(dat, response, numcomp) {
    # partial least squares regression using AIC to choose num components
    # numcomps is max numcomps to check using AIC in case p is large
    # 
    df = dat
    df[["response"]] <- response
    df = df[complete.cases(df),]
    
    ans <- plsr(df$response ~ ., data=df, ncomp=numcomp)
	predicted.y <- predict(ans, df)

	d.and.aic <- c()
	N <- nrow(df)
	for(i in 1:numcomp){
		pe <- mean((predicted.y[,,i] - df$response)^2)
		d <- i
		se2 <- pe*N/(N-d-1)
		aic <- pe + 2*d/N*se2
		d.and.aic <- rbind(d.and.aic, c(d,aic))
	}

	mincomps = d.and.aic[d.and.aic[,2]==min(d.and.aic[,2]),]
	
    png("plsAIC.png")
    plot(d.and.aic, type="l", xlab="d", ylab="AIC")
    dev.off()

    ans <- plsr(response ~ ., data=df, ncomp=mincomps[1])

    return(list("model" = ans, "mincomps" = mincomps))
}

lmBackAIC <- function(dat, response){
   	df = dat 
    df[["response"]] <- response
    df = df[complete.cases(df),]
    
    step <- stepAIC(lm(df$response ~ . , data=df))   # using backward selection

	se2 <- summary(step)$sigma^2
	d <- length(step$coefficients)-1 
	N <- nrow(df)
	pe <- mean( (predict(step, df) - df$response)^2)
	aic <- pe + 2*d/N*se2
	

	# get AIC for full model

	ans <- lm(df$response ~ ., data=df)
	se2 <- summary(ans)$sigma^2
	d <- length(ans$coefficients)-1 
	N <- nrow(df)
	pe <- mean( (predict(ans, df) - df$response)^2)
	aic2 <- pe + 2*d/N*se2

    return(list("model"=step, "aic"=aic, "aicback"=aic2, "full_model" = ans ))

}



##########################
# BIC Methods
ridgeBIC <- function(dat, response, lambdalist){
    # function runs BIC ridge regression. Returns the model and
    # lambda with BIC vector 

    lambda.and.bic <- c()
    df = dat
    dat = dat[complete.cases(response),]
    df[["response"]] <- response
    df = df[complete.cases(response),]


    N <- nrow(df)
	
    for(i in lambdalist){
		lambda <- i
		X <- as.matrix(dat)
		W <- X %*% solve(t(X)%*%X + lambda * diag(dim(X)[2])) %*% t(X)
		d <- sum(diag(W))  # the trace of W
		
        ans <- lm.ridge(df$response ~ ., lambda = i, data=df)
		predicted.y <- ans$ym + as.matrix(dat) %*% coef(ans)[2:(dim(X)[2]+1)] 
        
        pe <- mean((predicted.y - df$response)^2)
		se2 <- sum((predicted.y - df$response)^2)/(N-d-1)
		lambda.and.bic <- rbind(lambda.and.bic,
                                  c(lambda, N/se2 * (pe + log(N)*d/N*se2)))
	}

	minlambda = lambda.and.bic[lambda.and.bic[,2]==min(lambda.and.bic[,2]),]

    mod = lm.ridge(df$response~. , lambda = minlambda[1], data = df)

    png("ridgeBIC.png")
    plot(lambda.and.bic, type="o", xlab=expression(lambda), ylab="BIC")
    dev.off()

    return(list("model" = mod, "minlambda" = minlambda))
}


lassoBIC <- function(dat, response, lambdalist) {
    # Lasso Regression using AIC to choose lambda
    
    lambda.and.bic <- c()
    df = dat
    df[["response"]] <- response
    df = df[complete.cases(df),]
    N = nrow(df)
    
    for(i in lambdalist){
		ans <- l1ce(response ~ ., data=df, bound=i, absolute.t=FALSE)
		predicted.y <- predict(ans, df)
		pe <- mean((predicted.y - df$response)^2)
		d <- sum(coef(ans) != 0)
		se2 <- pe*N/(N-d-1)
		bic <- N/se2 * (pe + log(N)*d/N*se2)		
        lambda.and.bic <- rbind(lambda.and.bic, c(i,bic))
	}

	minlambda = lambda.and.bic[lambda.and.bic[,2]==min(lambda.and.bic[,2]),]
    # get min lambda for min values that occurr more than once.
    minlambda = minlambda[minlambda[,1]==min(minlambda[,1]),]

    png("LambdaLassoBIC.png")
    plot(lambda.and.bic, type="o", xlab=expression(t), ylab="BIC")
    dev.off()
	ans <- l1ce(df$response ~ ., data=df, bound=minlambda[1], absolute.t=FALSE)
	norm = sum(abs(coef(ans))[2:ncol(df)])
	
    return(list("model" = ans, "norm" = norm, "minlambda"=minlambda))

}



pcrBIC <- function(dat, response, numcomp){
    # Pricipal Component Regression with AIC to choose num components
    
    df = dat 
    df[["response"]] <- response
    df = df[complete.cases(df),]
    
    ans <- pcr(df$response ~ ., data=df, ncomp=numcomp)
	predicted.y <- predict(ans, df)
    
    d.and.bic <- c()
	N <- nrow(df)
	for(i in 1:numcomp){
		pe <- mean((predicted.y[,,i] - df$response)^2)
		d <- i
		se2 <- pe*N/(N-d-1)
		bic <- N/se2 * (pe + log(N)*d/N*se2)
		d.and.bic <- rbind(d.and.bic, c(d,bic))
	}

	mincomps = d.and.bic[d.and.bic[,2]==min(d.and.bic[,2]),]
    
    ans <- pcr(df$response ~ ., data=df, ncomp=mincomps[1])
    png("pcrBIC.png")
    plot(d.and.bic, type="o", xlab="d", ylab="BIC")
    dev.off()
    
    return(list("model"=ans, "mincomps" = mincomps))
}



plsBIC <- function(dat, response, numcomp) {
    # partial least squares regression using AIC to choose num components
    # numcomps is max numcomps to check using AIC in case p is large
    # 
    df = dat
    df[["response"]] <- response
    df = df[complete.cases(df),]
    
    ans <- plsr(df$response ~ ., data=df, ncomp=numcomp)
	predicted.y <- predict(ans, df)

	d.and.bic <- c()
	N <- nrow(df)
	for(i in 1:numcomp){
		pe <- mean((predicted.y[,,i] - df$response)^2)
		d <- i
		se2 <- pe*N/(N-d-1)
		bic <- N/se2 * (pe + log(N)*d/N*se2)
		d.and.bic <- rbind(d.and.bic, c(d,bic))
	}

	mincomps = d.and.bic[d.and.bic[,2]==min(d.and.bic[,2]),]
	
    png("plsBIC.png")
    plot(d.and.bic, type="o", xlab="d", ylab="BIC")
    dev.off()
    ans <- plsr(response ~ ., data=df, ncomp=mincomps[1])

    return(list("model" = ans, "mincomps" = mincomps))
}


lmBackBIC <- function(dat, response){
   	df = dat 
    df[["response"]] <- response
    df = df[complete.cases(df),]
	N <- nrow(df)
    step <- stepAIC(lm(df$response ~ . , data=df), k=log(N))

	se2 <- summary(step)$sigma^2
	d <- length(step$coefficients)-1 
	pe <- mean( (predict(step, df) - df$response)^2)
	bic <- N/se2 * (pe + log(N)*d/N*se2)
	

	# get BIC for full model

	ans <- lm(df$response ~ ., data=df)
	se2 <- summary(ans)$sigma^2
	d <- length(ans$coefficients)-1 
	N <- nrow(df)
	pe <- mean( (predict(ans, df) - df$response)^2)
	bic2 <- pe + 2*d/N*se2

    return(list("model"=step, "bic"=bic, "bicfull"=bic2, "full_model" = ans ))


}

ridgeCV <- function(dat, response, lambdalist){
    lambda.and.bic <- c()
    df = dat
    dat = dat[complete.cases(response),]
    df[["response"]] <- response
    df = df[complete.cases(response),]
	
    
    ans <- lm.ridge(df$response ~ ., lambda = lambdalist, data=df)

    png("ridgeCV.png")
    plot(lambdalist, ans$GCV, xlab="Lambda", ylab="GCV", type="l")
    dev.off()
	
    minlambda = ans$GCV[ans$GCV == min(ans$GCV)]
    thing = which(ans$GCV == min(ans$GCV))
	minlam = c(c(lambdalist)[thing[[1]]], minlambda[[1]]*nrow(X))
    
    model = lm.ridge(df$response ~ ., lambda = minlam[1], data=df)

    return(list("model" = model, "minlambda"= minlam, "gcv" = ans$GCV))
}

lassoCV <- function(dat, response, lambdalist) {
    lambda.and.bic <- c()
    df = dat
    dat = dat[complete.cases(response),]
    df[["response"]] <- response
    df = df[complete.cases(response),]
	
    ans <- l1ce(df$response ~ ., data=df, bound=lambdalist, absolute.t=FALSE)
	GCV <- gcv(ans)
	
    png("lassoCV.png")
    plot(GCV[,1], GCV[,4], xlab="t/t_0 (relative t)", 
                            ylab="GCV", type="p")
    dev.off()

	minlambda = GCV[,4][GCV[,4] == min(GCV[,4])]
    thing = seq(.01,1,.01)[GCV[,4] == min(GCV[,4])]
	model = l1ce(response ~ ., data=df, bound=thing, absolute.t=FALSE)
    mint = c(thing,minlambda)
    return(list("model" = model, "mint" =mint , "GCV" = GCV))
}

pcrCV <- function(dat, response, numcomponents){
    df = dat
    dat = dat[complete.cases(response),]
    df[["response"]] <- response
    df = df[complete.cases(response),]

    png("pcrCV.png")
    par(mfrow = c(2,2))

    for(k in c(5,7,10,12)){
      save <- NULL  
      for(reps in 1:10){
        ans <- pcr(df$response ~ ., data=df, ncomp=numcomponents,
                   validation="CV", segments=k)
        cv.preds <- ans$validation$pred
        
        nc.and.err = c()
        for(j in 1:numcomponents){
 			EstErr <- mean( (df$response - cv.preds[,1,j])^2)
			nc.and.err <- rbind(nc.and.err, c(j,EstErr))
		}
		save <- rbind(save,nc.and.err[,2])
      }
      plot(1:numcomponents, save[1,], xlab="# Components in PCR",
                      ylab="Estimated Error", type="l")
	  for(reps in 2:10) lines(1:numcomponents, save[reps,], type="l")
	
	  ceil = ceiling(nrow(df)/k)*(k-1)
	  title(k, col.main="red", cex.main=.9)
    }
    dev.off()

    save <- NULL  
      for(reps in 1:10){
        ans <- pcr(df$response ~ ., data=df, ncomp=numcomponents,  
                   validation="LOO")
        cv.preds <- ans$validation$pred
        
        nc.and.err = c()
        for(j in 1:numcomponents){
 			EstErr <- mean( (df$response - cv.preds[,1,j])^2)
			nc.and.err <- rbind(nc.and.err, c(j,EstErr))
		}
		save <- rbind(save,nc.and.err[,2])
      }
      plot(1:numcomponents, save[1,], xlab="# Components in PCR",
                      ylab="Estimated Error", type="l")
	  for(reps in 2:10) lines(1:numcomponents, save[reps,], type="l")
	
	  ceil = ceiling(nrow(df)/k)*(k-1)
	  title(k, col.main="red", cex.main=.9)
    
    dev.off()

} 

plsrCV <- function(dat, response, numcomponents){
    df = dat
    dat = dat[complete.cases(response),]
    df[["response"]] <- response
    df = df[complete.cases(response),]

    png("plsrCV.png")
    par(mfrow = c(2,2))

    for(k in c(5,7,10,12)){
      save <- NULL  
      for(reps in 1:10){
        ans <- plsr(df$response ~ ., data=df, ncomp=numcomponents,
                   validation="CV", segments=k)
        cv.preds <- ans$validation$pred
        
        nc.and.err = c()
        for(j in 1:numcomponents){
 			EstErr <- mean( (df$response - cv.preds[,1,j])^2)
			nc.and.err <- rbind(nc.and.err, c(j,EstErr))
		}
		save <- rbind(save,nc.and.err[,2])
      }
      plot(1:numcomponents, save[1,], xlab="# Components in PLSR",
                      ylab="Estimated Error", type="l")
	  for(reps in 2:10) lines(1:numcomponents, save[reps,], type="l")
	
	  ceil = ceiling(nrow(df)/k)*(k-1)
	  title(k, col.main="red", cex.main=.9)
    }
    dev.off()
} 



