# Poster
análisis de consumo de combustible con parámetros de operación mecánico ambientales 
#PAQUETES
library(gamlss)
library(mctest)
library(car)
library(corrplot)
library(leaps)
library(rpart)
library(rpart.plot)
library(ggplot2)
library(xtable)
#Datos
datos <- read.csv2(file = 'Datos.csv', header = TRUE)
head(datos)
plot(datos[-c(6,7)])
datos$Modo<-as.factor(datos$Modo)
attach(datos)
class(Configuración)
par(mfrow=c(2,2))
plot(CC~RV, lty = 3, pch = 1)
lines(lowess(CC ~ RV),col='red', lty = 6, lwd = 2)
plot(CC~DR, lty = 3, pch = 1)
lines(lowess(CC ~ DR),col='red', lty = 6, lwd = 2)
plot(CC~X.CO, lty = 3, pch = 1)
lines(lowess(CC ~ X.CO),col='red', lty = 6, lwd = 2)
plot(CC~X.CH4, lty = 3, pch = 1)
lines(lowess(CC ~ X.CH4),col='red', lty = 6, lwd = 2)
mat <- cor(datos[,-c(6,7)])
par(mfrow=c(1,1))
corrplot.mixed(mat, lower.col = "black", number.cex = .7)
#Modelos con lm
#USANDO LEAPS
horizonte <- formula(CC~(RV+DR+X.CO+X.CH4)^2+I(RV^2)
                     +I(DR^2)+I(X.CH4^2)+I(RV^3))
m1 <- regsubsets(horizonte,data=datos,intercept = TRUE, nbest = 3, nvmax = 8)
summary(m1)
par(mfrow=c(1,1))
plot(m1,scale = "adjr2", main =expression(R[Adj]^2))
plot(m1,scale = "bic", main = 'BIC')
plot(m1, scale = "Cp", main = 'Cp')
#Mejores modelos
modr2 <- lm(CC ~ RV+DR+I(RV^2)+I(X.CH4^2)+I(RV^3)+RV*DR+RV*X.CO)
modbic <- lm(CC ~ I(RV^2)+I(X.CH4^2)+RV*DR)
modcp <- lm(CC~ RV+DR+I(RV^2)+I(X.CH4^2)+I(RV^3)+RV*DR)
nb <- c("R2","BIC","Cp")
ai1 <- AIC(modr2,k=3)
ai2 <- AIC(modbic,k=3)
ai3 <- AIC(modcp,k=3)
va <- c(ai1,ai2,ai3)
co1 <- cor(CC, modr2$fitted.values)
co2 <-cor(CC, modbic$fitted.values)
co3 <-cor(CC, modcp$fitted.values)
vc <- c(co1,co2,co3)
mat <- cbind(nb,va,vc)
colnames(mat)<- c("Críterio","AIC","Correlación")
xtable(mat)
#MEJOR MODELO
modrs <- lm(CC~ RV+DR+I(RV^2)+I(X.CH4^2)+I(RV^3)+RV*DR)
#Con Indicadoras
modrsi <- lm(CC~ RV+DR+I(RV^2)+I(X.CH4^2)+I(RV^3)+RV*DR+Modo+Configuración)
AIC(modrsi, k=3)
cor(CC,modrsi$fitted.values)
#EL MEJOR ES CON AMBAS INDICADORAS
vif(modrsi)
x <- datos[,c(2,3,5)]
y <- datos[,1]
omcdiag(x,y,detr = 0.001, red = 0.6, conf = 0.99,
        theil = 0.6, cn = 15)
imcdiag(x,y)
#MODELO BÁSICO
modb <- lm(CC~RV+DR+X.CH4+Modo+Configuración)
AIC(modb,k=3)
cor(CC,modb$fitted.values)


#CON STEPAIC
AIC(modb)
extractAIC(modb)
fullmod <- lm(CC~(RV+DR+X.CO+X.CH4)^2+I(RV^2)
              +I(DR^2)+I(X.CH4^2)+I(RV^3)+Modo + Configuración)
summary(fullmod)
fullmod2 <- lm(CC~(RV+DR+X.CH4)^2+I(RV^2)
              +I(DR^2)+I(X.CH4^2)+I(RV^3)+Modo + Configuración)

#Modelo Backwards
modback <- stepAIC(object = fullmod, direction = 'backward', k = 3)
modback2 <- stepAIC(object = fullmod2, direction = 'backward', k = 3)
AIC(modback)
  AIC(modback2)
extractAIC(fit=modback, k=2)
#Modelo Forward
hori <- formula(CC~(RV+DR+X.CO+X.CH4+Modo+Configuración)^2+I(RV^2)
                +I(DR^2)+I(X.CH4^2)+I(RV^3))
empty.mod <- lm(CC ~ 1)
modfor <- stepAIC(empty.mod, trace = FALSE, direction = 'forward',
                  scope = hori)
modfor$anova
AIC(modfor)

modboth <- stepAIC(empty.mod, trace = FALSE, direction = 'both',
                   scope = hori)
modboth$anova

summary(modback)$sigma
summary(modfor)$sigma
summary(modboth)$sigma
summary(modback2)$sigma
summary(modback)$adj.r.squared
summary(modfor)$adj.r.squared
summary(modboth)$adj.r.squared


par(mfrow=c(1, 2))
plot(modback, main="Backward", pch=19, cex=1, which=1)
plot(modfor, main="Forward", pch=19, cex=1, which=1)
plot(modback2, main="Backward", pch=19, cex=1, which=1)

qqnorm(modback$residuals, main="Backward")
qqline(modback$residuals, col = 'red')
qqnorm(modback2$residuals, main="Backward")
qqline(modback2$residuals, col = 'red')
qqnorm(modfor$residuals, main="Forward")
qqline(modfor$residuals, col = 'red')
shapiro.test(modback$residuals)
shapiro.test(modfor$residuals)
shapiro.test(modback2$residuals)
#MODELOS GAMLSS
fitt <- fitDist(CC, type = "realplus", k=3)
fitt$fits
#Step GAIC
mod1 <- gamlss(CC~(RV+DR+X.CO+X.CH4)^2+I(RV^2)
               +I(DR^2)+I(X.CH4^2)+I(RV^3)+Modo+Configuración,
               family = "WEI3", trace = FALSE)
summary(mod1)
vif(modback)
mod2 <- stepGAIC(mod1)
mod3 <- stepGAIC(mod1, scope=list(lower=~1,
                                  upper=~(RV+DR+X.CO+X.CH4)^2+I(RV^2)
                                  +I(DR^2)+I(X.CH4^2)+I(RV^3)+Modo+Configuración,
                                  k=2))
Rsq(mod3)
x1 <- cbind(x,X.CO)
od <- omcdiag(x1,y)
id <-imcdiag(x1,y)
tod <- od$odiags
  tid <- id$idiags
xtable(tod)
xtable(tid)
#STEPGAIC.ALL
m1 <- gamlss(CC~1,family = "WEI3", trace=FALSE)
summary(m1)
m2<-stepGAICAll.A(m1, scope=list(lower=~1, upper=~(RV+DR+X.CO+X.CH4)^2+I(RV^2)
                                 +I(DR^2)+I(X.CH4^2)+I(RV^3)+Modo+Configuración),
                  k=3)
m3 <- stepGAICAll.B(m1, scope=list(lower=~1, upper=~(RV+DR+X.CO+X.CH4)^2+I(RV^2)
                                   +I(DR^2)+I(X.CH4^2)+I(RV^3)+Modo+Configuración),
                    k=3)

#COMPARANDO
summary(m2)
AIC(m2, k = 3)
Rsq(m2)

#CON GAMMA

m4 <- gamlss(CC~1,family = "GA", trace=FALSE)
summary(m4)
m5<-stepGAICAll.A(m4, scope=list(lower=~1, upper=~(RV+DR+X.CO+X.CH4)^2+I(RV^2)
                                 +I(DR^2)+I(X.CH4^2)+I(RV^3)+Modo+Configuración),
                  k=3)
m6 <- stepGAICAll.B(m4, scope=list(lower=~1, upper=~(RV+DR+X.CO+X.CH4)^2+I(RV^2)
                                   +I(DR^2)+I(X.CH4^2)+I(RV^3)+Modo+Configuración),
                    k=3)
summary(m5)
AIC(m5, k=3)
Rsq(m5)
#GAMA GENERALIZADA
m7 <- gamlss(CC~1,family = "GG", trace=FALSE)
m8<-stepGAICAll.A(m7, scope=list(lower=~1, upper=~(RV+DR+X.CO+X.CH4)^2+I(RV^2)
                                 +I(DR^2)+I(X.CH4^2)+I(RV^3)+Modo+Configuración),
                  k=3)
m9 <- stepGAICAll.B(m7, scope=list(lower=~1, upper=~(RV+DR+X.CO+X.CH4)^2+I(RV^2)
                                   +I(DR^2)+I(X.CH4^2)+I(RV^3)+Modo+Configuración),
                    k=3)
summary(m8)
AIC(m8, k=3)
Rsq(m8)
#COMPARANDO
AIC1 <- AIC(m2, k = 3)
R21 <- Rsq(m2)
yf2 <- fitted(m2, what = "mu")
cor1 <- cor(CC, yf2)
AIC2 <- AIC(m5, k=3)
R22 <- Rsq(m5)
yf5 <- fitted(m5, what = "mu")
cor2 <-cor(CC, yf5)
AIC3 <- AIC(m8, k=3)
R23 <- Rsq(m8)
yf8 <- fitted(m8, what = "mu")
cor3 <-cor(CC, yf8)

enc <- c("AIC", "R2", "Correlación")
f1 <- c(AIC1,R21,cor1)
f2 <- c(AIC2,R22,cor2)
f3 <- c(AIC3,R23,cor3)
ML <- data.frame(enc,f1,f2,f3)
names(ML) <- c("Criterio", "WEI3", "GA", "GG")
ML
xtable(ML)
#ARBOLES DE REGRESIÓN
modar1 <- rpart(CC~., data = datos, method = "anova")
modar1
rpart.plot(modar1, type = 3, digits = 4, fallen.leaves = TRUE)
p1 <- predict(modar1,datos)                              
p1
CC

#BASE
try1 <- lm(CC~., data=datos)
summary(try1)
AIC(try1, k=3)
cor(CC, try1$fitted.values)

#LOESS

aux <- function(x){
  if (x>1|x<0)stop("Span invalido")
  modl <- loess(CC ~ RV + DR + X.CO + X.CH4,data=datos, degree = 2, span = x)
  cor(CC,modl$fitted)
}
aux <- Vectorize(aux)
sp <- seq(0.1,0.9, by=0.01)
cori<-aux(sp)
which.max(cori)
sp[40]

aux1 <- function(x){
  if (x>1|x<0)stop("Span invalido")
  modl <- loess(CC ~ RV + DR + X.CO + X.CH4,data=datos, degree = 1, span = x)
  cor(CC,modl$fitted)
}
aux1 <- Vectorize(aux1)
cori1<-aux1(sp)
which.max(cori1)
sp[16]
max(cori)
max(cori1)

modl1 <- loess(CC ~ RV + DR + X.CO + X.CH4,data=datos, degree = 1, span = 0.25)
modl2 <- loess(CC ~ RV + DR + X.CO + X.CH4,data=datos, degree = 2, span = 0.49)
AIC(modl1)
loess.aic <- function (x) {
  	if (!(inherits(x,"loess"))) stop("Error: argument must 
                                      > be a loess object")
  	# extract values from loess object
    	span <- x$pars$span
    	n <- x$n
    	traceL <- x$trace.hat
    	sigma2 <- sum( x$residuals^2 ) / (n-1)
    	delta1 <- x$one.delta
    	delta2 <- x$two.delta
    	enp <- x$enp
    
      	aicc <- log(sigma2) + 1 + 2* (2*(traceL+1)) / (n-traceL-2)
 
	aicc1<- n*log(sigma2) + n* ( 
  (delta1/delta2)*(n+enp)/(delta1^2/delta2)-2 )
	gcv  <- n*sigma2 / (n-traceL)^2
 	
   	result <- list(span=span, aicc=aicc, aicc1=aicc1, gcv=gcv)
   	return(result)
}
loess.aic(modl1)
loess.aic(modl2)
loess.aic(mlrv)
mlrv <- loess(CC~RV, degree = 1, span=0.1)
auxrv <- function(x){
  if (x>1|x<0)stop("Span invalido")
  modl <- loess(CC ~ RV,data=datos, degree = 1, span = x)
  cor(CC,modl$fitted)
}
auxrv <- Vectorize(auxrv)
corirv<-auxrv(sp)
which.max(corirv)
sp[1]
max(cori)
max(cori1)
cor(CC,modl1$fitted)
cor(CC,modl2$fitted)

#MEJORES MODELOS

modrsi <- lm(CC~ RV+DR+I(RV^2)+I(X.CH4^2)+I(RV^3)+RV*DR+Modo+Configuración)
res1 <- modrsi$residuals
res2 <- modfor$residuals
res3 <- m8$residuals
fit3 <- fitted(m8, what = "mu")

par(mfrow=c(2,3))
plot(modrsi, main="Leap", pch=19, cex=1, which = 1)
plot(modfor, main="stepAIC", pch=19, cex=1, which = 1)
plot(fit3, res3, main="stepGAIC.ALL", pch=19, cex=1)
lines(lowess(res3 ~ fit3),col='red', lty = 6, lwd = 1)
qqnorm(res1, main="Leap")
qqline(res1, col = 'red')
qqnorm(res2, main="stepAIC")
qqline(res2, col = 'red')
qqnorm(res3, main="stepGAIC.ALL")
qqline(res3, col = 'red')
s1 <- summary(modrsi)
s2 <- summary(modfor)
s3 <- summary(m8)

aic1 <- AIC(modrsi, k=3)
aic2 <- AIC(modfor, k=3)
aic3 <- AIC(m8, k=3)
vaic <- c(aic1,aic2,aic3)
cor1 <- cor(CC,modrsi$fitted.values)
cor2 <- cor(CC, modfor$fitted.values)
cor3 <- cor(CC, fit3)
vcor <- c(cor1,cor2,cor3)
ECM1 <- mean((CC-modrsi$fitted.values)^2)
ECM2 <- mean((CC-modfor$fitted.values)^2)
ECM3 <- mean((CC-fit3)^2)
vECM <- c(ECM1,ECM2,ECM3)
vmod <- c("Leaps", "stepAIC", "stepGAIC")
M <- cbind(vmod,vaic,vcor,vECM)
colnames(M) <- c("Modelo","AIC", "Correlación", "ECM")
M
xtable(M, digits = 4)
summary(modfor)
