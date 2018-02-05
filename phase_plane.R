require("phaseR")

sir <- function(t, y, parameters) {
  x <- y[1]
  y <- y[2]
  #n <- y[3]
  beta <- parameters[1]
  gamma <- parameters[2]
  dy <- numeric(2)
  dy[1] <-  - beta*x*y
  dy[2] <- beta*x*y - gamma*y
  list(dy)
}

flowField <- flowField(sir, x.lim = c(-1,100), y.lim = c(-1,100), 
                       parameters = c(0.3,0.1), points = 19, add = FALSE,
                       xlab = "S", ylab = "I", main = "Phase Plane with NullClines, beta=0.3, gamma=0.1, R>1")
nullclines <-
  nullclines(sir, x.lim = c(-1, 100), y.lim = c(-1, 100),
               parameters = c(1/5,1/15), points = 500)
y0 <- matrix(c(99, 1, 40, 80, 30, 70, 20,20), ncol = 2, nrow = 4, byrow = TRUE)
trajectory <-
  trajectory(sir, y0 = y0, t.end = 300,
               parameters = c(0.03,0.1), colour = rep("black", 3))

flowField <- flowField(sir, x.lim = c(-1,100), y.lim = c(-1,100), 
                       parameters = c(0.001,0.1), points = 19, add = FALSE,
                       xlab = "S", ylab = "I", main = "Phase Plane with NullClines, beta=0.001, gamma=0.1, R<1")
nullclines <-
  nullclines(sir, x.lim = c(-1, 100), y.lim = c(-1, 100),
             parameters = c(0.001,0.1), points = 500)
y0 <- matrix(c(99, 1, 40, 80, 30, 70, 20,20), ncol = 2, nrow = 4, byrow = TRUE)
trajectory <-
  trajectory(sir, y0 = y0, t.end = 300,
             parameters = c(0.001,0.1), colour = rep("black", 3))
