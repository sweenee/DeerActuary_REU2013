# ###############################################################
# R script to approximate the balance in an insurance
# fund with a focus on deer collisions.
#
# Authors: Lindong Zhou
#          Kelly Black
#

TESTING = TRUE

if(TESTING)
  {
    # The script is set up for testing purposes.
    # Set the seed using the value in the data file.
    load('randomSeed.Rdata')
    set.seed(theSeed)
  }

T <- 10        # Set the final time
N <- 3650      # Set the number of sample points
r <- 10			   # r tilde - The scaled growth factor
F <- 1000000	 # F tilde - The scaled carrying capacity
alpha <- 0.1   # Noise coeff. in the deer equation.

dt <- T/N
dw <- sqrt(dt) * rnorm(N)
w <- c(0, cumsum(dw))

x0 <- F
z0 <- F/x0
a <- alpha ^ 2 - r
b <- alpha
g0 <- 1 - r/(0.5 * b ^ 2 - a)
t  = seq(0,T,dt);
phi <- exp((a - 0.5 * b ^ 2) * seq(0, T, dt) - b * w)

integral <- vector(length = N + 1)
integral[1] <- 0
integral = c(0,cumsum(exp((0.5 * b ^ 2 - a) * t[1:N] + b * w[1:N]) * dw[1:N]));

Ztrue <- r/(0.5 * b ^ 2 - a) +
  g0 * phi - b * r/(0.5 * b ^ 2 - a) * phi * integral
Xtrue <- F/Ztrue
plot(seq(0, T, dt), Xtrue,
     type = "l", col = "orange", lwd = 2,
     ylim = c(900000, 1100000),
     xlab = "T", ylab = "X", main = "Deer Population")

L <- 730
R <- N/L
Dt <- R * dt
Xem <- Xmil <- vector(length = L + 1)
Xem[1] <- Xmil[1] <- x0
Xtemp1 <- Xtemp2 <- x0
for (i in 1:L) {
  Winc <- sum(dw[(R * (i - 1) + 1):(R * i)])
  Xtemp1 <- Xtemp1 + Dt * r * Xtemp1 * (1 - Xtemp1/F) + alpha * Xtemp1 * Winc
  Xem[i + 1] <- Xtemp1
  Xtemp2 <- Xtemp2 + Dt * r * Xtemp2 * (1 - Xtemp2/F) + alpha * Xtemp2 * Winc + 0.5 * alpha * alpha * Xtemp2 * (Winc * Winc - Dt)
  Xmil[i + 1] <- Xtemp2
}
points(seq(0, T, Dt), Xem, col = "green", pch = 16)
points(seq(0, T, Dt), Xmil, col = "red", pch = 4, lwd = 2)
legend("bottomleft", legend = c("Euler", "Milstein"),
       col = c("green", "red"), pch = c(16, 4), lwd = c(1, 2), cex = 0.8)


# Population Distribution
Z10 <- vector(length = 1000)
for (i in 1:1000) {
  dt <- T/N
  dw <- sqrt(dt) * rnorm(N)
  w <- c(0, cumsum(dw[1:(N - 1)]))
  x0 <- F
  z0 <- F/x0
  a <- alpha ^ 2 - r
  b <- alpha
  g0 <- 1 - r/(0.5 * b ^ 2 - a)
  phi <- exp((a - 0.5 * b ^ 2) * T - b * w[N])
  integral <- 0
  ds <- dt
  s = seq(0, T - ds, ds)
  integral = sum(exp((0.5 * b ^ 2 - a) * s + b * w) * dw)
  Ztrue <- r/(0.5 * b ^ 2 - a) + g0 * phi - b * r/(0.5 * b ^ 2 - a) * phi * integral
  Z10[i] <- Ztrue
}

mean(Z10)
sd(Z10)
hist(Z10, probability = TRUE, xlab = "Z", ylab = "Frequency density",
     main = "Distribution of Z", col = "grey")
qqnorm(Z10)

# Insurance Payout

I <- 0.0005
beta <- 0.1
P <- 10000
gamma <- 0.1
g <- 0.0001

m0 <- (beta * F - P)/(I - g)
Mem <- Mmil <- vector(length = L + 1)
Mem[1] <- Mmil[1] <- m0
Mtemp1 <- Mtemp2 <- m0
for (i in 1:L) {
  Winc <- sum(dw[(R * (i - 1) + 1):(R * i)])
  Mtemp1 <- Mtemp1 + Dt * (I * Mtemp1 - beta * Xtrue[i] + P) - Winc * gamma * Xtrue[i]
  Mem[i + 1] <- Mtemp1
  Mtemp2 <- Mtemp2 + Dt * (I * Mtemp2 - beta * Xtrue[i] + P) - Winc * gamma * Xtrue[i] + 0.5 * gamma * gamma * Xtrue[i] * (Winc * Winc - Dt)
  Mmil[i + 1] <- Mtemp2
}
plot(seq(0, T, Dt), Mem, type = "p", pch = 16,
     col = "light blue", lwd = 2, xlab = "T", ylab = "M",
     main = "Insurance Payout")
points(seq(0, T, Dt), Mmil, pch = 4, col = "purple", lwd = 2)
legend("bottomright", legend = c("Euler", "Milstein"),
       col = c("light blue", "purple"), pch = c(16, 4),
       lwd = c(1, 2), cex = 0.8)
