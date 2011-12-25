timecpp <- c(0.001,0.005,0.04,0.22,0.621,1.34)
p <- c(50,100,250,500,750,1000)

fit <- lm(p~timecpp)

par(mfrow=c(1,2))
qqnorm(residuals(fit))
qqline(residuals(fit))
#TA-plot
plot(residuals(fit))
