x = -35:180
y = 0.2*log(x+50)

# corrected age: 
plot(y~c(x-35),type='l',lty=2,xlim=c(-35,180),xlab='chronological age')
lines(y~x,type='l',col='pink')
abline(v=-35)
