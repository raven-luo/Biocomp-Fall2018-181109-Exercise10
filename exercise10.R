#Task 1
setwd("C:/Users/raven/Desktop/biocomputing/exercise10/Biocomp-Fall2018-181109-Exercise10")
data=read.csv(file="data.txt",header=TRUE)
#Maximum likelihood
#Quadratic
nllike_quad=function(p,x,y){
  a=p[1]
  b=p[2]
  c=p[3]
  sigma=exp(p[4]) #Question: you still need an error term right?
  expected=a+b*x+c*(x^2)
  nll=-sum(dnorm(x=y,mean=expected,sd=sigma,log=TRUE))
  return(nll)
}
initialGuess_quad=c(1,1,1,1)
fit_quad=optim(par=initialGuess,fn=nllike_quad,x=data$x,y=data$y)
print(fit_quad)
quadNLL=nllike_quad(fit_quad$par,data$x,data$y) #Question: is this how you get quadNLL?
#Linear
nllike_lin=function(p,x,y){
  a=p[1]
  b=p[2]
  sigma=exp(p[3])
  expected=a+b*x
  nll=-sum(dnorm(x=y,mean=expected,sd=sigma,log=TRUE))
  return(nll)
}
initialGuess_lin=c(1,1,1)
fit_lin=optim(par=initialGuess_lin,fn=nllike_lin,x=data$x,y=data$y)
print(fit_lin)
linNLL=nllike_lin(fit_lin$par,data$x,data$y)
#Likelihood ratio test
testStatistic=2*(quadNLL-linNLL)
df=length(fit_quad$par)-length(fit_lin$par)
p_value=1-pchisq(q=testStatistic,df=df)
if (p_value<0.05) {print("The difference between quadratic and linear models are significant. The quadratic model is better.")} else {print("The difference between the models are not significant. The linear model is better")}
#Question: is this the way you compare which model is better? By the necessity/significance of the quadratic term?

#Task2