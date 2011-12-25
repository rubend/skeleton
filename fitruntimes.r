#trying to fit the runtimes with a meaningful fit
en <- c(3,4,5)
p <- c(50,100,500,1000)
alpha <- 0.01
m1 <- c(323,331.5,166.1005,113.3078)
m2 <- c(324,324.75,166.2621,116.0109)
m3 <- c(326,323,159.3972,119.9674)
ratio <- c(m1,m2,m3)
data <- cbind(matrix(ratio,12,1),rep(p,3),c(rep(en[1],4),rep(en[2],4),rep(en[3],4)))
frame <- data.frame(data)
names(frame) <- c("r","p","en")
fit <- lm(r ~.,frame)
summary(fit) #this is without taking any logs
plot(frame$p,fitted(fit))
lines(frame$p[1:4],frame$r[1:4],type="l",col=3)
par(mfrow=c(1,2))
plot(fitted(fit),residuals(fit),col=3)# tukey-anscombe plot
abline(0,0)

qqnorm(residuals(fit))# checking if the residuals seem to be normally distributed
qqline(residuals(fit))

plot(residuals(fit),col=2)# serial correlation plot


frame$r <- log(frame$r)
frame$p <- log(frame$p)

logfit <- lm(r~.,frame)
sum <- summary(logfit) # looks good, p is significant.
par(mfrow=c(1,2))
plot(fitted(logfit),residuals(logfit),col=3)
abline(0,0)

qqnorm(residuals(logfit))
qqline(residuals(logfit))

par(mfrow=c(2,2))
plot(fitted(fit),residuals(fit))
plot(fitted(logfit),residuals(logfit))

qqnorm(residuals(fit))
qqline(residuals(fit))
qqnorm(residuals(logfit))
qqline(residuals(logfit))

#so the en variable seems to be unrelevant

newframe <- frame[,-3]

finalsum <- summary(lm(r~.,newframe))
finalsum$coefficients

## //NEW TEST\\ new measurements
## done for p= 50,100,250,500,750,1000
m4 <- c(325,318.75,216.6053,143.3178,131.9550,114.8082)#en = 3
m5 <- c(325,255.60,208.8974,160.5571,129.2362,116.1248)#en = 4
m6 <- c(325,318.75,210.1282,160.6588,129.4663,116.0660)#en = 5
ratio <- c(m4,m5,m6)
en <- c(3,4,5)
p2 <- c(50,100,250,500,750,1000)
data <- cbind(matrix(ratio,18,1),rep(p2,3),c(rep(en[1],6),rep(en[2],6),rep(en[3],6)))
frame <- data.frame(data)
names(frame) <- c("r","p","en")

logfit <- lm(r~.,frame)
sum <- summary(logfit) # looks good, p is significant.
par(mfrow=c(1,2))
plot(fitted(logfit),residuals(logfit),col=3)
abline(0,0)

qqnorm(residuals(logfit))
qqline(residuals(logfit))

frame$r <- log(frame$r)
frame$p <- log(frame$p)

logfit <- lm(r~.,frame)
sum <- summary(logfit) # looks good, p is significant.
par(mfrow=c(1,2))
plot(fitted(logfit),residuals(logfit),col=3)
abline(0,0)

qqnorm(residuals(logfit))
qqline(residuals(logfit))
#so the en variable seems to be unrelevant

newframe <- frame[,-3]

##model selection problem, need to use different data for the fit this time.
finalsum <- summary(lm(r~.,newframe))
finalsum$coefficients

## seems like en is really unrelevant and the log(p) coefficient is -0.36129
