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
library(deSolve)
library(ggplot2)
ddSim=function(t,y,p){
  N1=y[1]
  N2=y[2]
  R1=p[1]
  R2=p[2]
  a11=p[3] #a12<a11, p[5]<p[3]
  a22=p[4]
  a12=p[5]
  a21=p[6] #a21<a22, p[6]<p[4]
  dN1dt=R1*(1-N1*a11-N2*a12)*N1
  dN2dt=R2*(1-N2*a22-N1*a21)*N2
  return(list(c(dN1dt,dN2dt)))
}
params1=c(0.5,0.5,0.6,0.4,0.4,0.2) #p[5]<p[3] & p[6]<p[4]
params2=c(0.5,0.5,0.5,10,10,2) #p[6]<p[4] only
params3=c(0.5,0.5,10,2,0.5,10) #p[5]<p[3] only
params4=c(0.5,0.5,0.5,2,10,10) #neither
N0=c(10,10) #Question: I just randomly generated this?
times=1:100
modelSim1=ode(y=N0,times=times,func=ddSim,parms=params1)
modelOutput1=data.frame(time=modelSim1[,1],N1=modelSim1[,2],N2=modelSim1[,3])
plot1=ggplot(data=modelOutput1)+geom_line(data=modelOutput1,aes(x=time,y=N1),color="blue")+geom_line(data=modelOutput1,aes(x=time,y=N2),color="red")+xlab("Time")+ylab("Population")
print(plot1) #Question: Is this the right way to try models and generate plots? Here there is coexistence, but it couldn't prove that other conditions canNOT generate coexistence?
modelSim2=ode(y=N0,times=times,func=ddSim,parms=params2)
modelOutput2=data.frame(time=modelSim2[,1],N1=modelSim2[,2],N2=modelSim2[,3])
plot2=ggplot(data=modelOutput2)+geom_line(data=modelOutput2,aes(x=time,y=N1),color="blue")+geom_line(data=modelOutput2,aes(x=time,y=N2),color="red")+xlab("Time")+ylab("Population")
print(plot2)
modelSim3=ode(y=N0,times=times,func=ddSim,parms=params3)
modelOutput3=data.frame(time=modelSim3[,1],N1=modelSim3[,2],N2=modelSim3[,3])
plot3=ggplot(data=modelOutput3)+geom_line(data=modelOutput3,aes(x=time,y=N1),color="blue")+geom_line(data=modelOutput3,aes(x=time,y=N2),color="red")+xlab("Time")+ylab("Population")
print(plot3)

#Question: to summarize into for loops? But how to return all the data frames and plots? rbind?
for (params in c(params1,params2,params3,params4)) {
  modelSim=ode(y=N0,times=times,func=ddSim,parms=params) #overwritten with each iteration
  modelOutput=data.frame(time=modelSim[,1],N1=modelSim[,2],N2=modelSim[,3]) #appended with each iteration, with a column indicating parameter set
  plot=ggplot(data=modelOutput)+geom_line(data=modelOutput,aes(x=time,y=N1),color="blue")+geom_line(data=modelOutput,aes(x=time,y=N2),color="red")+xlab("Time")+ylab("Population")
  print(plot)
  }
