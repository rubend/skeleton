#trying to fit the runtimes with a meaningful fit
en <- c(3,5)
p <- c(50,100,500,1000)
alpha <- 0.01

ratio <- c(323,331.5,166.1005,113.3078,326,323,159.3972,119.9674)
data <- cbind(matrix(ratio,8,1),rep(p,2),c(rep(en[1],4),rep(en[2],4)))
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
