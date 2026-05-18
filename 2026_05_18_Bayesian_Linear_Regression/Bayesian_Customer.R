library(readr)
dat <- read_csv("supermarket_bayesian_regression_dataset.csv")

dat <- dat[dat$spending > 0, ]
dat$gender <- as.factor(dat$gender)
dat$status <- as.factor(dat$status)

lmod <- lm(spending ~ ., data = dat)
summary(lmod)

library(car)
vif(lmod)
durbinWatsonTest(lmod)
residualPlots(lmod)
plot(lmod$residuals)
shapiro.test(lmod$residuals)
hist(lmod$residuals)

# Predictions
pred_lm <- predict(lmod)

# Residuals
res_lm <- dat$spending - pred_lm

# Metrics
rmse_lm <- sqrt(mean(res_lm^2))

mae_lm <- mean(abs(res_lm))

r2_lm <- summary(lmod)$r.squared

adj_r2_lm <- summary(lmod)$adj.r.squared

cat("===== Classical LR =====\n")
cat("RMSE:", rmse_lm, "\n")
cat("MAE:", mae_lm, "\n")
cat("R2:", r2_lm, "\n")
cat("Adjusted R2:", adj_r2_lm, "\n")

###################################################################
# Observed vs fitted
plot(dat$spending, pred_lm,
     xlab = "Observed spending",
     ylab = "Fitted values",
     pch = 16, col = "blue",
     main = "Observed vs Fitted (LM)")

# 1:1 line
abline(0, 1, col = "red", lwd = 2)

legend("topleft",
       legend = c("Observed vs fitted", "1:1 line"),
       col = c("blue", "red"), pch = c(16, NA), lwd = c(NA, 2), bty = "n")



###################################################################
library(rjags)

mod1_string = "
model{
  for (i in 1:n){
    y[i] ~ dt(mu[i], prec, nu)
    mu[i] = b[1] +
            b[2]*age[i] +
            b[3]*gender[i] +
            b[4]*salary[i] +
            b[5]*status[i]
  }

  # Priors
  for (j in 1:5){
    b[j] ~ dnorm(0, 0.001)
  }
  prec ~ dgamma(0.001, 0.001)
  nu   ~ dexp(0.1)

  sig2 = 1/prec
  sig  = sqrt(sig2)
}
"

# Convert to numeric
set.seed(72)
dat$gender <- ifelse(dat$gender == "Male", 0,1)
dat$status  <- ifelse(dat$status == "Single", 0,1)

data_jags = list(
  y      = dat$spending,
  n      = nrow(dat),
  age    = dat$age,
  gender = dat$gender,
  salary = dat$salary,
  status = dat$status
)

params1 = c("b", "sig", "sig2", "nu")

inits1 = function(){
  list(
    "b[1]" = rnorm(1, 200, 10),
    "b[2]" = rnorm(1, 3, 5),
    "b[3]" = rnorm(1, 40, 10),
    "b[4]" = rnorm(1, 0, 1),
    "b[5]" = rnorm(1, 50, 10),
    prec   = rgamma(1, 100, 20),
    nu     = 10
  )
}
mod = jags.model(
  textConnection(mod1_string),
  data     = data_jags,
  inits    = inits1,
  n.chains = 3
)

update(mod, 1e3)

mod_sim = coda.samples(
  model          = mod,
  variable.names = params1,
  n.iter         = 5e3
)

mod_csim = as.mcmc(do.call(rbind, mod_sim))

#### Convergence Diagnosis
plot(mod_sim)
gelman.diag(mod_sim)
autocorr.diag(mod_sim)
autocorr.plot(mod_sim)
effectiveSize(mod_sim)
# =========================
# Bayesian Regression Metrics
# =========================

X <- model.matrix(
  ~ age + gender + salary + status,
  data = dat
)

# Posterior coefficients
beta_post <- mod_csim[, grep("^b\\[", colnames(mod_csim))]
beta_mean <- colMeans(beta_post)

# 95% credible intervals
beta_ci <- apply(beta_post, 2, quantile, probs = c(0.025, 0.975))

# Combine into a nice summary table
beta_summary <- cbind(
  Mean = beta_mean,
  t(beta_ci)   # transpose so intervals are columns
)

format(beta_summary, scientific = FALSE, digits = 4)
# Predictions
pred <- X %*% beta_mean

# Residuals
res <- dat$spending - pred

# RMSE and MAE
rmse_bayes <- sqrt(mean(res^2))
mae_bayes  <- mean(abs(res))

# R²
ss_res  <- sum(res^2)
ss_tot  <- sum((dat$spending - mean(dat$spending))^2)
r2      <- 1 - (ss_res / ss_tot)

# Adjusted R²
n      <- nrow(dat)
p      <- ncol(X) - 1
adj_r2 <- 1 - (1 - r2) * ((n - 1) / (n - p - 1))

# Posterior summary for nu
nu_post <- mod_csim[, "nu"]

cat("===== Bayesian Student-t Regression =====\n")
cat("RMSE               :", round(rmse_bayes,       4), "\n")
cat("MAE                :", round(mae_bayes,         4), "\n")
cat("R²                 :", round(r2,                4), "\n")
cat("Adjusted R²        :", round(adj_r2,            4), "\n")
cat("Posterior mean nu  :", round(mean(nu_post),     2),
    "  [95% CI:", round(quantile(nu_post, 0.025), 2),
    "-",         round(quantile(nu_post, 0.975), 2), "]\n")


dic = dic.samples(mod, n.iter = 1e3)

x1 = c(30, 1, 75000, 0)
x2 = c(35, 0, 70000, 1)

# Monte Carlo predictive distributions
pred1 <- mod_csim[, "b[1]"] + mod_csim[, c(2,3,4,5)] %*% x1
pred2 <- mod_csim[, "b[1]"] + mod_csim[, c(2,3,4,5)] %*% x2

# Densities
dens1 <- density(pred1)
dens2 <- density(pred2)


xrange <- range(c(dens1$x, dens2$x))
yrange <- range(c(dens1$y, dens2$y))

plot(dens1, col = "blue", lwd = 2,
     xlim = xrange,
     ylim = yrange,
     main = "Posterior predictive distributions",
     xlab = "Predicted spending", ylab = "Density")

# Overlay second density
lines(dens2, col = "red", lwd = 2)

legend("topright", legend = c("Person 1", "Person 2"),
       col = c("blue", "red"), lwd = 2, bty = "n")

mean(pred1 < pred2)


###########################################################

comparison <- data.frame(
  
  Model = c("Classical LR", "Bayesian LR"),
  
  RMSE = c(
    rmse_lm,
    rmse_bayes
  ),
  
  MAE = c(
    mae_lm,
    mae_bayes
  ),
  
  R2 = c(
    r2_lm,
    r2
  ),
  
  Adjusted_R2 = c(
    adj_r2_lm,
    adj_r2
  )
  
)

print(comparison)

###############################################################

# Posterior draws of coefficients
beta_post <- mod_csim[, grep("^b\\[", colnames(mod_csim))]

# Design matrix
X <- model.matrix(~ age + gender + salary + status, data = dat)

# Posterior fitted values (matrix: iterations x observations)
fitted_post <- beta_post %*% t(X)

# Posterior mean fitted values
fitted_mean <- colMeans(fitted_post)

# 95% credible intervals for fitted values
fitted_ci <- apply(fitted_post, 2, quantile, probs = c(0.025, 0.975))

# Observed response
y <- dat$spending

# Plot: observed vs fitted mean with 95% CI
plot(y, fitted_mean,
     xlab = "Observed spending",
     ylab = "Fitted (posterior mean)",
     pch = 16, col = "blue")

# Add 95% credible intervals as vertical lines
arrows(x0 = y, y0 = fitted_ci[1,],
       x1 = y, y1 = fitted_ci[2,],
       angle = 90, code = 3, length = 0.05, col = "darkgrey")

# Add 1:1 line for reference
abline(0, 1, col = "red", lwd = 2)

legend("topleft",
       legend = c("Observed vs fitted", "95% credible interval", "1:1 line"),
       col = c("blue", "darkgrey", "red"),
       pch = c(16, NA, NA), lwd = c(NA, 1, 2), bty = "n")
