rm(list = ls())

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/UofT Admin and TA/STA302 S 2021/Lectures/Lecture 12")
## Install UsingR package ##
#  install.packages("UsingR")
## Scatter plot 3d ##
#install.packages("scatterplot3d")

library(glmnet)
library(rms)
library(MASS)

### Prostate Data ###
train <- read.table("https://gattonweb.uky.edu/sheather/book/docs/datasets/prostateTraining.txt", header = T)
head(train)
str(train)

test <- read.table("https://gattonweb.uky.edu/sheather/book/docs/datasets/prostateTest.txt", header = T)
    
### First fit a multiple linear regression ##
model.lm <- lm(lpsa ~ ., data = train[, -c(1)])
summary(model.lm)    

## Perform Prediction ##
pred.y <- predict(model.lm, newdata = test, type = "response")

## Prediction error ##
mean((test$lpsa - pred.y)^2)

## Fit a ridge penalty ##
model.ridge <- glmnet(x = as.matrix(train[,2:9]), y = train$lpsa, standardize = T, alpha = 0)

## Perform Prediction ##
pred.y.ridge <- predict(model.ridge, newx = as.matrix(test[,2:9]), type = "response")

## Prediction error ##
mean((test$lpsa - pred.y.ridge)^2)


## Fit a LASSO penalty ##
model.lasso <- glmnet(x = as.matrix(train[,2:9]), y = train$lpsa, standardize = T, alpha = 1)

## Perform Prediction ##
pred.y.lasso <- predict(model.lasso, newx = as.matrix(test[,2:9]), type = "response")

## Prediction error ##
mean((test$lpsa - pred.y.lasso)^2)

## Elastic net ##

model.EN <- glmnet(x = as.matrix(train[,2:9]), y = train$lpsa, standardize = T, alpha = 0.5)

## Perform Prediction ##
pred.y.EN <- predict(model.EN, newx = as.matrix(test[,2:9]), type = "response")

## Prediction error ##
mean((test$lpsa - pred.y.EN)^2)

####### Variable selection #######

## Step wise regression ###

sel.var.aic <- step(model.lm, trace = 0, k = 2) #log(nrow(dat1)))
select_var<-attr(terms(sel.var.aic), "term.labels")   
select_var

### LASSO selection ###

## Perform cross validation to choose lambda ##
set.seed(1002656486)
cv.out <- cv.glmnet(x = as.matrix(train[,2:9]), y = train$lpsa, standardize = T, alpha = 1)
plot(cv.out)
best.lambda <- cv.out$lambda.1se
best.lambda
co<-coef(cv.out, s = "lambda.1se")

#Selection of the significant features(predictors)

## threshold for variable selection ##

thresh <- 0.00
# select variables #
inds<-which(abs(co) > thresh )
variables<-row.names(co)[inds]
sel.var.lasso<-variables[!(variables %in% '(Intercept)')]
sel.var.lasso

## Step wise regression ###

### First fit a multiple linear regression ##
## Based on AIC ##
model.lm <- lm(lpsa ~ ., data = train[, -c(1)])
summary(model.lm)  
n <- nrow(train)
sel.var.aic <- step(model.lm, trace = 0, k = 2, direction = "both") 
sel.var.aic<-attr(terms(sel.var.aic), "term.labels")   
sel.var.aic

## Based on BIC ##
model.lm <- lm(lpsa ~ ., data = train[, -c(1)])
summary(model.lm)  
n <- nrow(train)
sel.var.bic <- step(model.lm, trace = 0, k = log(n), direction = "both") 
sel.var.bic<-attr(terms(sel.var.bic), "term.labels")   
sel.var.bic


### Cross Validation and prediction performance of AIC based selection ###
ols.aic <- ols(lpsa ~ ., data = train[,which(colnames(train) %in% c(sel.var.aic, "lpsa"))], 
               x=T, y=T, model = T)

## 10 fold cross validation ##    
aic.cross <- calibrate(ols.aic, method = "crossvalidation", B = 10)
## Calibration plot ##
pdf("aic_cross.pdf", height = 8, width = 16)
plot(aic.cross, las = 1, xlab = "Predicted LPSA", main = "Cross-Validation calibration with AIC")
dev.off()

## Test Error ##
pred.aic <- predict(ols.aic, newdata = test[,which(colnames(train) %in% c(sel.var.aic, "lpsa"))])
## Prediction error ##
pred.error.AIC <- mean((test$lpsa - pred.aic)^2)


### Cross Validation and prediction performance of BIC based selection ###
ols.bic <- ols(lpsa ~ ., data = train[,which(colnames(train) %in% c(sel.var.bic, "lpsa"))], 
               x=T, y=T, model = T)

## 10 fold cross validation ##    
bic.cross <- calibrate(ols.bic, method = "crossvalidation", B = 10)
## Calibration plot ##
pdf("bic_cross.pdf", height = 8, width = 16)
plot(bic.cross, las = 1, xlab = "Predicted LPSA", main = "Cross-Validation calibration with BIC")
dev.off()

## Test Error ##
pred.bic <- predict(ols.bic, newdata = test[,which(colnames(train) %in% c(sel.var.bic, "lpsa"))])
## Prediction error ##
pred.error.BIC <- mean((test$lpsa - pred.bic)^2)

### Cross Validation and prediction performance of lasso based selection ###
ols.lasso <- ols(lpsa ~ ., data = train[,which(colnames(train) %in% c(sel.var.lasso, "lpsa"))], 
                 x=T, y=T, model = T)

## 10 fold cross validation ##    
lasso.cross <- calibrate(ols.lasso, method = "crossvalidation", B = 10)
## Calibration plot ##
pdf("lasso_cross.pdf", height = 8, width = 16)
plot(lasso.cross, las = 1, xlab = "Predicted LPSA", main = "Cross-Validation calibration with LASSO")
dev.off()

## Test Error ##
pred.lasso <- predict(ols.lasso, newdata = test[,which(colnames(train) %in% c(sel.var.lasso, "lpsa"))])
## Prediction error ##
pred.error.lasso <- mean((test$lpsa - pred.lasso)^2)

print(c(pred.error.AIC, pred.error.BIC, pred.error.lasso))


#### GAM ###
gc()
rm(list = ls())
library(faraway)
library(gam)
library(mgcv)


### The Cherry Tree Data ###
data(trees)
## Height vs Volume ##
pdf("HvV.pdf", height = 8, width = 12)
par(family = 'serif', mfrow = c(1,2))
plot(trees$Height, trees$Volume, xlab = "Height", ylab = "Volume")
abline(lm(Volume ~ Height, data = trees))

plot(trees$Height, log(trees$Volume), xlab = "Height", ylab = "log(Volume)")
abline(lm(I(log(Volume)) ~ Height, data = trees))
dev.off()

## Girth vs Volume ##
pdf("GvV.pdf", height = 8, width = 12)
par(family = 'serif', mfrow = c(1,2))
plot(trees$Girth, trees$Volume, xlab = "Girth", ylab = "Volume")
abline(lm(Volume ~ Girth, data = trees))
plot(trees$Girth, log(trees$Volume), xlab = "Girth", ylab = "log(Volume)")
abline(lm(I(log(Volume)) ~ Girth, data = trees))
dev.off()

## First fit a gamma glm ##
ct.glm<-glm(Volume~Height+Girth,
            family=Gamma(link=log),data=trees)
summary(ct.glm)

## First fit a gamma gam ##
ct1.gam<-gam(Volume~s(Height)+s(Girth),
             family=Gamma(link=log),data=trees)
ct1.gam
par(mfrow = c(1,2))
plot(ct1.gam,residuals=T)


### The Ozone data ###
data(ozone)

## Fit a linear model ##
olm <- lm(O3 ~ temp + ibh  + ibt, data = ozone)
summary(olm)

## Fit a Additive model ##
ogam<-gam(O3~s(temp) + s(ibh) + s(ibt),
          data=ozone)
summary(ogam)
pdf("ogam.pdf", height = 8, width = 12)
par(family = 'serif', mfrow =c(2,2))
plot(ogam,residuals=T)
dev.off()

## The cubic splines ##
ogam.cs<-gam(O3~s(temp, bs = "cr") + s(ibh, bs = "cr") + s(ibt, bs = "cr"),
             data=ozone)
summary(ogam)
pdf("ogamcs.pdf", height = 8, width = 12)
par(family = 'serif', mfrow =c(2,2))
plot(ogam.cs,residuals=T)
dev.off()

## 20 basis functions ##
ogam.cs20<-gam(O3~s(temp, bs = "cr", k=20) + s(ibh, bs = "cr", k = 4) + s(ibt, bs = "cr", k = 20),
               data=ozone)
summary(ogam)
pdf("ogamcs20.pdf", height = 8, width = 12)
par(family = 'serif', mfrow =c(2,2))
plot(ogam.cs20,residuals=T)
dev.off()

## 4 basis functions ##
ogam.cs4<-gam(O3~s(temp, bs = "cr", k=4) + s(ibh, bs = "cr", k = 4) + s(ibt, bs = "cr", k = 4),
              data=ozone)
summary(ogam)
pdf("ogamcs4.pdf", height = 8, width = 12)
par(family = 'serif', mfrow =c(2,2))
plot(ogam.cs4,residuals=T)
dev.off()

## Anova or deviance tests ##
anova(ogam, ogam.cs, test = "F")
anova(ogam.cs20, ogam.cs4, test = "F")

## The interaction  ##
ogam.int<-gam(O3~s(temp,ibh) + s(ibt),
              data=ozone)
summary(ogam.int)
pdf("ogamint.pdf", height = 8, width = 12)
par(family = 'serif', mfrow =c(2,2))
plot(ogam.int,residuals=T)
vis.gam(ogam.int,theta=-45,color="heat")
dev.off()

## ## Anova or deviance tests ##
anova(ogam, ogam.int, test = "F")


