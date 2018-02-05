install.packages("markovchain")
require("markovchain")

mcSIR <- new("markovchain", states=c("S","I","R"),
             transitionMatrix=matrix(data=c(1-0.001,0.001,0,0,1-0.1,0.1,0,0,1),
             byrow=TRUE, nrow=3), name="SIR")
initialState <- c(99,1,0)

timesteps <- 300
sir.df <- data.frame( "timestep" = numeric(),
                      "S" = numeric(), "I" = numeric(),
                      "R" = numeric(), stringsAsFactors=FALSE)
for (i in 0:timesteps) {
  newrow <- as.list(c(i,round(as.numeric(initialState * mcSIR ^ i),0)))
  sir.df[nrow(sir.df) + 1, ] <- newrow
}

text_main = paste("Discrete Time Markov Chain Simulation, One Realization")
mult.fig(1,main=text_main,oma=c(1,2,4,1))
plot(sir.df$timestep,sir.df$S, xlab="time step", ylab="number of people")
points(sir.df$timestep,sir.df$I, col="red")
points(sir.df$timestep,sir.df$R, col="green")

legend("topright", legend=c("S", "I", "R"), bty="n",lwd=3,col=c("black", "red", "green"))

absorbingStates(mcSIR)
transientStates(mcSIR)
ab.state <- absorbingStates(mcSIR)
occurs.at <- min(which(sir.df[,ab.state]==max(sir.df[,ab.state])))

occurs.at

