lda1 <- function(value, data) {
# Make sure the data is composed of two variables and rename the variables 
    if (dim(data)[2] != 2) {
        warning("The data must include exactly two variables")
    } else {
        if (is.numeric(data[, 1])) {
            names(data) <- c("x", "class")
        } else {
            names(data) <- c("class", "x")
        }
# Retrieve class names, number of classes, and length of vector of values to classify 
        tab <- table(data$class)
        n_cla <- dim(tab)
        len <- length(value)
        nam <- names(tab)
# Create vectors and matrices used for the calculation of the discriminant score
        pi <- c() 
        mu <- c() 
        vark <- c()  
        wei <- c()
        dis <- matrix(nrow = len, ncol = n_cla)
        exp <- matrix(nrow = len, ncol = n_cla)
# Compute discriminant scores for each class and for each value of input vector 
        for (i in (1:n_cla)) {
            pi[i] <- tab[i] / dim(data)[1] # prior probabilities
            mu[i] <- mean(data$x[data$class == nam[i]]) # class means
            vark[i] <- var(data$x[data$class == nam[i]]) # class variances
            wei[i] <- ((tab[i] - 1) / (dim(data)[1] - n_cla)) * vark[i] # weighted class variances
        }
        var <- sum(wei) # sum of weighted class variances, i.e. global variance estimate
        for (i in (1:n_cla)) {
            for (j in (1:len)) {
                dis[j, i] <- value[j] * (mu[i] / var) - (mu[i] ^ 2) / (2 * var) + log(pi[i]) # discriminant score
                exp[j, i] <- exp(dis[j, i]) # exponentiated discriminant score (to calculate probabilities)
            }
        }
# Compute probabilities for each class
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
# Add column names to discrimant scores and priors, create vectors which will host output values
        colnames(dis) <- nam
        names(pi) <- nam
        class <- c()
        score <- c()
        proba <- c()
# Extract class name for which discrimant score is maximal, maximal discriminant score and maximal probability
        for (j in (1:len)) {
            class[j] <- colnames(dis)[dis[j,] == max(dis[j,])]
            score[j] <- max(dis[j,])
            proba[j] <- max(prob[j,])
        }
        out <- as.data.frame(cbind(value, class, score, proba))
# Coerce probability to 1 when numerator and denominator are too large for fraction to be computed (NaN), add warning
        for (j in (1:len)) {
            if (out$proba[j] == "NaN") {
                out$proba <- as.numeric(out$proba)
                out$proba[j] <- 1
                warning("NaN probability coerced to 1 for certain input values")
            }
        }
# Define minimum and maximum x's for each class, choose minimum among class minimums and maximum among class maxmimums (for plot)
        min <- min(mu - 3 * sqrt(var))
        max <- max(mu + 3 * sqrt(var))
# Create matrix which will host class distributions to be plotted, start loop with class with largest prior to get ylim right 
        y <- matrix(ncol = n_cla, nrow = 1000)
        x <- seq(min, max, length = 1000)
        subout <- rbind(pi, mu)
        subout <- subout[, order(subout[1,], decreasing = TRUE)]
        for (i in (1:n_cla)) {
            y[, i] <-
                subout[1, i] * dnorm(x, mean = subout[2, i], sd = sqrt(var))
        }
# Define color palette and plot class with largest prior; plot other classes 
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
# Add real density for the input variable in the background
        lines(density(data$x),
              col = "grey",
              xlim = c(min, max))
# Add vertical line or polygon showing value(s) to be classified by the function; add names of classes under the curves
        if (length(value) == 1) {
            abline(v = value,
                   col = rgb(0.5, 0.5, 0, 0.15))
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
# Reorder subout by ascending mean (mu)
        subout <- subout[, order(subout[2,])]
        rownames(subout) <- c("priors", "averages")
        out <- list(subout, out)
        return(out)
    }
}