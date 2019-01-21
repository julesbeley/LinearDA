rm(list = ls())
X2 <- data.frame(x = rnorm(2000, mean = 6, sd = 2), class = "blue")
X3 <- data.frame(x = rnorm(10000, mean = 15, sd = 2), class = "green")
X <- rbind(X2, X3)

lda1 <- function(value, data) {
  tab <- table(data$class)
  n_cla <- dim(tab)
  len <- length(value)
  nam <- names(tab)
  pi <- c()
  mu <- c()
  vark <- c()
  wei <- c()
  dis <- matrix(nrow = len, ncol = n_cla)
  exp <- matrix(nrow = len, ncol = n_cla)
  for (i in (1:n_cla)) {
    pi[i] <- table(data$class)[i] / dim(data)[1]
    mu[i] <- mean(data$x[data$class == nam[i]])
    vark[i] <- var(data$x[data$class == nam[i]])
    wei[i] <- ((tab[i] - 1) / (dim(data)[1] - n_cla)) * vark[i]
  }
  var <- sum(wei)
  for (i in (1:n_cla)) {
    for (j in (1:len)) {
      dis[j, i] <- value[j] * (mu[i] / var) - (mu[i] ^ 2) / (2 * var) + log(pi[i])
      exp[j, i] <- exp(dis[j, i])
    }
  }
  sum <- c()
  for (j in (1:len)) {
    sum[j] <- sum(exp[j, ])
  }
  prob <- matrix(nrow = len, ncol = n_cla)
  for (i in (1:n_cla)) {
    for (j in (1:len)) {
      prob[j, i] <- exp[j, i] / sum[j]
    }
  }
  colnames(dis) <- nam
  names(pi) <- nam
  exp <- cbind(value, exp)
  sum <- cbind(value, sum)
  class <- c()
  score <- c()
  proba <- c()
  for (j in (1:len)) {
    class[j] <- colnames(dis)[dis[j,] == max(dis[j,])]
    score[j] <- max(dis[j,])
    proba[j] <- max(prob[j,])
  }
  out <- as.data.frame(cbind(value, class, score, proba))
  for (j in (1:len)) {
    if (out$proba[j] == "NaN") {
      out$proba <- as.numeric(out$proba)
      out$proba[j] <- 1
      warning("NaN probability coerced to 1 for certain input values")
    }
  }
  min <- min(mu - 3 * sqrt(var))
  max <- max(mu + 3 * sqrt(var))
  y <- matrix(ncol = n_cla, nrow = 1000)
  subout <- rbind(pi, mu)
  subout <- subout[, order(subout[1,], decreasing = TRUE)]
  x <- seq(min, max, length = 1000)
  for (i in (1:n_cla)) {
    y[, i] <- subout[1, i] * dnorm(x, mean = subout[2, i], sd = sqrt(var))
  }
  colors <- heat.colors(n_cla)
  plot(x,
       y[, 1],
       type = "n",
       xlab = "",
       ylab = "")
  title("Class densities")
  lines(x, y[, 1], lwd = 2, col = colors[2])
  for (i in (2:n_cla)) {
    lines(x, y[, i], col = colors[i], lwd = 2)
  }
  abline(v = value[1])
  abline(v = value[length(value)])
  lines(density(data$x), col = "grey")
  for (i in (1:n_cla)) {
    text(subout[2, i], 0.9 * max(y[, i]), colnames(subout)[i])
  }
  subout <- subout[, order(subout[2,])]
  rownames(subout) <- c("priors", "averages")
  out <- list(subout, out, var)
  return(subout)
}
lda1(seq(7, 12, length=10), X)


