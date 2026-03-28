
# claim: given equal variances, distance between y1=y2,y3 is same
# as euclidean distance betwen y1, y3
a = matrix(c(1,0.99,0,0.99,1,0,0,0,1),nrow=3)

x1 = c(2,2,1.5)

sum(x1^2)
sum(2^2+1.5^2)
t(x1)%*%solve(a)%*%x1
