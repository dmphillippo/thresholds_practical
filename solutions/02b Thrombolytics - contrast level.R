## 02b. Thrombolytics example - contrast level

library(nmathresh)


# Step 1: NMA and data prep -----------------------------------------------

# Again, we will use the saved posterior results
Thrombo.post.summary
Thrombo.post.cov

K <- 6   # Number of treatments

# Contrast design matrix is
X <- matrix(ncol = K-1, byrow = TRUE,
            c(1, 0, 0, 0, 0,
              0, 1, 0, 0, 0,
              0, 0, 1, 0, 0,
              0, 0, 0, 1, 0,
              0, -1, 1, 0, 0,
              0, -1, 0, 1, 0,
              0, -1, 0, 0, 1))


# Step 2: Reconstruct two-stage likelihood covariance matrix --------------

lik.cov <- recon_vcov(Thrombo.post.cov, prior.prec = .0001, X = X)

# The KL divergence is very small, which means that the reconstructed likelihood
# covariance matrix results in a posterior which is very close to that coming
# from the original NMA.


# Step 3: Derive thresholds -----------------------------------------------

thresh <- nma_thresh(mean.dk = Thrombo.post.summary$statistics[1:(K-1), "Mean"], 
                     lhood = lik.cov, 
                     post = Thrombo.post.cov, 
                     nmatype = "fixed", 
                     X = X, 
                     opt.max = FALSE)


# Step 4: Plot ------------------------------------------------------------

# Get treatment codes for the contrasts with data
d.a <- d.b <- vector(length = nrow(X))
for (i in 1:nrow(X)){
  d.a[i] <- ifelse(any(X[i, ] == -1), which(X[i, ] == -1), 0) + 1
  d.b[i] <- ifelse(any(X[i, ] == 1), which(X[i, ] == 1), 0) + 1
}

# Transform from d_ab style contrast references into d[i] style from the full
# set of contrasts for easy indexing in R
d.i <- d_ab2i(d.a, d.b, K = 6)

# Create plot data
plotdat <- data.frame(lab = paste0(d.b, " vs. ", d.a),
                      contr.mean = Thrombo.post.summary$statistics[d.i, "Mean"],
                      CI2.5 = Thrombo.post.summary$quantiles[d.i, "2.5%"],
                      CI97.5 = Thrombo.post.summary$quantiles[d.i, "97.5%"])

# Plot

# pdf("./figures/02b Thrombolytics - contrast level.pdf", width = 11, height = 3.5)

thresh_forest(thresh, contr.mean, CI2.5, CI97.5, label = lab, data = plotdat,
              label.title = "Contrast", xlab = "Log Odds Ratio", 
              CI.title = "95% Credible Interval",
              xlim = c(-.3, .3), refline = 0, digits = 2)

# dev.off()
