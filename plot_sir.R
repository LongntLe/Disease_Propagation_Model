require("sfsmisc")
source("sir_func.R")
#npop = 10000000
npop = 78000000
I_0 = 1
R_0 = 0
S_0 = npop-I_0-R_0
tbegin = 0
tend   = 150
vt = seq(tbegin,tend,1)  

gamma = 1/15        
#R0    = 1.50     
R0 = 3
beta  = R0*gamma
vparameters = c(gamma=gamma,beta=beta)
inits = c(S=S_0,I=I_0,R=R_0)

solved_model = as.data.frame(lsoda(inits, vt, derivative_calc_func, vparameters))

cat("The item names in the solved_model object are:",names(solved_model),"\n")

vS = solved_model$S
vI = solved_model$I
vR = solved_model$R
vtime = solved_model$time
vnpop = vS+vI+vR

mult.fig(4,main="SIR model of pandemic influenza with R0=1.5, N=300")

ymax = 1.4*max(vI/npop)
plot(vtime,vI/npop,type="l",xlab="time",ylab="fraction infected",ylim=c(0,ymax),lwd=3,col=4,main="Infected")

n=length(vtime)
lines(vtime[2:n],-diff(vS)/(diff(vtime)*vnpop[1:(n-1)]),type="l",lwd=3,col=2)

legend("topright",legend=c("total infected (prevalence)","newly infected/day (incidence)"),bty="n",lwd=3,col=c(4,2))

ymin = 0.9*min(vS/vnpop)
plot(vtime,vS/vnpop,type="l",xlab="time",ylab="fraction susceptible",ylim=c(ymin,1),lwd=3,main="Susceptible")

iind = which.min(abs(vS/vnpop-1/R0)) # find the index at which S/N is equal to 1/R0
lines(c(vtime[iind],vtime[iind]),c(-1000,1000),col=3,lwd=3)
legend("bottomleft",legend=c("time at which S=1/R0"),bty="n",lwd=3,col=c(3),cex=0.7)

plot(vtime,log(vI/vnpop),type="l",xlab="time",ylab="log(fraction infected)",lwd=3,col=4,main="log(Infected)")

text(40,-14,"Initial\n exponential\n rise",cex=0.7)

lines(c(vtime[iind],vtime[iind]),c(-1000,1000),col=3,lwd=3)
legend("topleft",legend=c("time at which S=1/R0","log(Infected)"),bty="n",lwd=3,col=c(3,4),cex=0.7)

epsilon = 0.00001
vsinf = seq(0,1-epsilon,epsilon)

LHS = -log(vsinf)+log(S_0/npop)
RHS = R0*(1-vsinf)
iind = which.min(abs(RHS-LHS))
sinf_predicted = vsinf[iind]

cat("The final fraction of susceptibles at the end of the epidemic from the model simulation is ",min(vS/vnpop),"\n")
cat("The final fraction of susceptibles at the end of the epidemic predicted by the final size relation is ",sinf_predicted,"\n")

