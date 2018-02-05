source("sir_func.R")

npop = 499
I_0 = 1
R_0 = 0
S_0 = npop-I_0-R_0
tbegin = 0
tend   = 150
vt = seq(tbegin,tend,1)

gamma = c(0.1, 0.2, 0.3)
beta  = c(0.5, 0.5, 0.5)
vparameters1 = c(gamma=gamma[1],beta=beta[1])
vparameters2 = c(gamma=gamma[2],beta=beta[2])
vparameters3 = c(gamma=gamma[3],beta=beta[3])
inits = c(S=S_0,I=I_0,R=R_0)

solved_model1 = as.data.frame(lsoda(inits, vt, derivative_calc_func, vparameters1))
solved_model2 = as.data.frame(lsoda(inits, vt, derivative_calc_func, vparameters2))
solved_model3 = as.data.frame(lsoda(inits, vt, derivative_calc_func, vparameters3))

colnames(solved_model1) <- c("time", "S1", "I1", "R1")
colnames(solved_model2) <- c("time", "S2", "I2", "R2")
colnames(solved_model3) <- c("time", "S3", "I3", "R3")

plt1 <- cbind(solved_model1$time, solved_model1$S1, solved_model2$S2, solved_model3$S3)
plt2 <- cbind(solved_model1$time, solved_model1$I1, solved_model2$I2, solved_model3$I3)
plt3 <- cbind(solved_model1$time, solved_model1$R1, solved_model2$R2, solved_model3$R3)

vS = solved_model1$S1
vI = solved_model1$I1
vR = solved_model1$R1
vtime = solved_model1$time
vnpop = vS+vI+vR

#################################
#gamma plot                     #
#################################

mult.fig(4,main="SIR model with varied gamma, by S, I, and R separately")

ymax = max(plt1[,2])*1.2
ymin = min(plt1[,2])*0.8
plot(plt1[,1],plt1[,2],type="l",xlab="time",ylab="S(t)",ylim=c(ymin,ymax),lwd=3,col=4,main="Susceptible Population")
lines(plt1[,1], plt1[,3], type='l', lwd=3, col=2)
lines(plt1[,1], plt1[,4], type='l', lwd=3, col="green")

legend("topright",legend=c("gamma=0.1", "gamma=0.2", "gamma=0.3"),bty="n",lwd=3,col=c(4,2, "green"))

ymax = max(plt2[,2])*1.2
ymin = min(plt2[,2])*0.8
plot(plt2[,1],plt2[,2],type="l",xlab="time",ylab="I(t)",ylim=c(ymin,ymax),lwd=3,col=4,main="Infected Population")
lines(plt2[,1], plt2[,3], type='l', lwd=3, col=2)
lines(plt2[,1], plt2[,4], type='l', lwd=3, col="green")

legend("topright",legend=c("gamma=0.1", "gamma=0.2", "gamma=0.3"),bty="n",lwd=3,col=c(4,2, "green"))

ymax = max(plt3[,2])*1.2
ymin = min(plt3[,2])*0.8
plot(plt3[,1],plt3[,2],type="l",xlab="time",ylab="R(t)",ylim=c(ymin,ymax),lwd=3,col=4,main="Removed Population")
lines(plt3[,1], plt3[,3], type='l', lwd=3, col=2)
lines(plt3[,1], plt3[,4], type='l', lwd=3, col="green")

legend("bottomright",legend=c("gamma=0.1", "gamma=0.2", "gamma=0.3"),bty="n",lwd=3,col=c(4,2, "green"))

#################################
#beta plot                      #
#################################
npop = 499
I_0 = 1
R_0 = 0
S_0 = npop-I_0-R_0
tbegin = 0
tend   = 150
vt = seq(tbegin,tend,1)

gamma = c(0.1, 0.1, 0.1)
beta  = c(0.15, 0.2, 0.3)
vparameters1 = c(gamma=gamma[1],beta=beta[1])
vparameters2 = c(gamma=gamma[2],beta=beta[2])
vparameters3 = c(gamma=gamma[3],beta=beta[3])
inits = c(S=S_0,I=I_0,R=R_0)

solved_model1 = as.data.frame(lsoda(inits, vt, derivative_calc_func, vparameters1))
solved_model2 = as.data.frame(lsoda(inits, vt, derivative_calc_func, vparameters2))
solved_model3 = as.data.frame(lsoda(inits, vt, derivative_calc_func, vparameters3))

colnames(solved_model1) <- c("time", "S1", "I1", "R1")
colnames(solved_model2) <- c("time", "S2", "I2", "R2")
colnames(solved_model3) <- c("time", "S3", "I3", "R3")

plt1 <- cbind(solved_model1$time, solved_model1$S1, solved_model2$S2, solved_model3$S3)
plt2 <- cbind(solved_model1$time, solved_model1$I1, solved_model2$I2, solved_model3$I3)
plt3 <- cbind(solved_model1$time, solved_model1$R1, solved_model2$R2, solved_model3$R3)

vS = solved_model1$S1
vI = solved_model1$I1
vR = solved_model1$R1
vtime = solved_model1$time
vnpop = vS+vI+vR

mult.fig(4,main="SIR model with varied beta, by S, I, and R separately")

ymax = max(plt1[,2])*1.2
ymin = min(plt3[,2])*0.8
plot(plt1[,1],plt1[,2],type="l",xlab="time",ylab="S(t)",ylim=c(ymin,ymax),lwd=3,col=4,main="Susceptible Population")
lines(plt1[,1], plt1[,3], type='l', lwd=3, col=2)
lines(plt1[,1], plt1[,4], type='l', lwd=3, col="green")

legend("topright",legend=c("beta=0.15", "beta=0.2", "beta=0.3"),bty="n",lwd=3,col=c(4,2, "green"))

ymax = max(plt1[,3])*1.2
ymin = min(plt3[,3])*0.8
plot(plt2[,1],plt2[,2],type="l",xlab="time",ylab="I(t)",ylim=c(ymin,ymax),lwd=3,col=4,main="Infected Population")
lines(plt2[,1], plt2[,3], type='l', lwd=3, col=2)
lines(plt2[,1], plt2[,4], type='l', lwd=3, col="green")

legend("topright",legend=c("beta=0.15", "beta=0.2", "beta=0.3"),bty="n",lwd=3,col=c(4,2, "green"))

ymax = max(plt1[,4])*1.2
ymin = min(plt3[,4])*0.8
plot(plt3[,1],plt3[,2],type="l",xlab="time",ylab="R(t)",ylim=c(ymin,ymax),lwd=3,col=4,main="Removed Population")
lines(plt3[,1], plt3[,3], type='l', lwd=3, col=2)
lines(plt3[,1], plt3[,4], type='l', lwd=3, col="green")

legend("topleft",legend=c("beta=0.15", "beta=0.2", "beta=0.3"),bty="n",lwd=3,col=c(4,2, "green"))
