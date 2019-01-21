lda1 <- function(value, data) {
  if (is.numeric(data[, 1]) == TRUE) {
    names(data) <- c("x", "class")
  } else {
    names(data) <- c("class", "x")
  }
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
      dis[j, i] <-
        value[j] * (mu[i] / var) - (mu[i] ^ 2) / (2 * var) + log(pi[i])
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
    y[, i] <-
      subout[1, i] * dnorm(x, mean = subout[2, i], sd = sqrt(var))
  }
  colors <- heat.colors(n_cla, alpha = 0.7)
  plot(x,
       y[, 1],
       type = "n",
       xlab = "",
       ylab = "")
  title("Predicted class densities")
  lines(x, y[, 1], lwd = 2, col = colors[1])
  for (i in (2:n_cla)) {
    lines(x, y[, i], col = colors[i], lwd = 2)
  }
  lines(density(data$x), col = "grey", xlim = c(min, max))
  if (length(value) == 1) {
    abline(v = value, col = rgb(0.5, 0.5, 0, 0.15))
  } else {
    polygon(
      x = c(value[1], value[1], value[length(value)], value[length(value)]),
      y = c(0, 1, 1, 0),
      col = rgb(0.5, 0.5, 0, 0.15),
      border = NA
    )
  }
  for (i in (1:n_cla)) {
    text(subout[2, i], 0.9 * max(y[, i]), colnames(subout)[i])
  }
  subout <- subout[, order(subout[2,])]
  rownames(subout) <- c("priors", "averages")
  out <- list(subout, out)
  return(out)
}
lda1(1.2, X)
