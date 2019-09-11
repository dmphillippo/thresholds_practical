## 02a. Thrombolytics example - study level

library(nmathresh)


# Step 1: NMA and data prep -----------------------------------------------


# The results of the NMA are available in the `nmathresh` package, as
# `Thrombo.post.summary` (for the posterior summaries) and `Thrombo.post.cov`
# (for the posterior covariance matrix).

# library(coda)   # Not required - but prints the summary in a nicer format
Thrombo.post.summary
Thrombo.post.cov


# To run a study level threshold analysis, we require the study data. This is
# available in `nmathresh` as a tab-delimited text file, and is read in like so:
dat.raww <- read.delim(system.file("extdata", "Thrombo_data.txt", package = "nmathresh"))

# print first few rows
head(dat.raww)

n <- nrow(dat.raww)   # number of studies


# Thresholds will be derived on the log odds ratio scale, so we derive log odds
# ratios against arm 1 as a reference.

# Log OR for two-arm trials, arm 2 vs. 1
dat.raww$lor.1 <- with(dat.raww, log(r.2 * (n.1 - r.1) / ((n.2 - r.2) * r.1)) )

dat.raww$k.1 <- dat.raww$t.2   # Arm 2 treatment
dat.raww$b.1 <- dat.raww$t.1   # Reference treatment

# Log OR for three-arm trials, arm 3 vs. 1
dat.raww$lor.2 <- with(dat.raww, log(r.3 * (n.1 - r.1) / ((n.3 - r.3) * r.1)) )

dat.raww$k.2 <- dat.raww$t.3   # Arm 3 treatment (NA if only 2 arms)
dat.raww$b.2 <- ifelse(is.na(dat.raww$k.2), NA, dat.raww$t.1)   # Reference treatment


# The likelihood covariance matrix is then constructed in a block diagonal
# manner, with the aid of `Matrix::bdiag`.

# LOR variances and covariances, likelihood covariance matrix V
V.diag <- as.list(rep(NA, n))
attach(dat.raww)
for (i in 1:n) {
  if (dat.raww$n_arms[i] == 2){
    V.diag[[i]] <- 1/r.1[i] + 1/r.2[i] + 1/(n.1[i]-r.1[i]) + 1/(n.2[i]-r.2[i])
  }
  else if (dat.raww$n_arms[i] == 3){
    v1 <- 1/r.1[i] + 1/r.2[i] + 1/(n.1[i]-r.1[i]) + 1/(n.2[i]-r.2[i])
    v2 <- 1/r.1[i] + 1/r.3[i] + 1/(n.1[i]-r.1[i]) + 1/(n.3[i]-r.3[i])
    # Covariance term
    c1 <- 1/r.1[i] + 1/(n.1[i] - r.1[i])
    V.diag[[i]] <- matrix(c(v1, c1, c1, v2), nrow = 2)
  }
}
detach(dat.raww)

library(Matrix)
V <- bdiag(V.diag)


# The raw data was imported in wide format, with one row per study. It is much
# easier to work with the data in long format, with one row per data point
# (contrast).

# Reshape the data to have one row per contrast per study
dat.rawl <- reshape(dat.raww, varying = c("lor.1", "b.1", "k.1", "lor.2", "b.2", "k.2"), 
                    timevar = "c", idvar = "studyID", direction = "long")

# Sort data by study and contrast, removing NA rows
dat.rawl <- dat.rawl[order(dat.rawl$studyID, dat.rawl$c, dat.rawl$b, na.last = NA), ]

N <- nrow(dat.rawl)   # number of data points
K <- length(unique(c(dat.rawl$b, dat.rawl$k)))   # Number of treatments


# Construct the design matrix, with a row for each contrast and K-1 columns (parameters)
X <- matrix(0, nrow = N, ncol = K-1)

for (i in 1:N){
  X[i, dat.rawl$k[i]-1] <- 1
  if (dat.rawl$b[i] != 1){
    X[i, dat.rawl$b[i]-1] <- -1
  }
}


# Step 2: Derive thresholds -----------------------------------------------

# We are now ready to perform a threshold analysis at the study level. The
# `nma_thresh` function takes the posterior means and covariance matrix of the
# treatment effect parameters ($d_k$), the likelihood covariance matrix, and the
# design matrix. We specify `nmatype = "fixed"` to derive thresholds for the FE
# model, and `opt.max = FALSE` since the optimal treatment is the one which
# minimises the log odds.

thresh <- nma_thresh(mean.dk = Thrombo.post.summary$statistics[1:(K-1), "Mean"], 
                     lhood = V, 
                     post = Thrombo.post.cov, 
                     nmatype = "fixed",
                     X = X,
                     opt.max = FALSE)

# The `nma_thresh` function prints some basic details, which can be used to
# verify that the input was as expected (the number of data points and
# treatments, and the base-case optimal treatment). These are also shown when
# the threshold object is printed:

thresh


# Step 3: Plot ------------------------------------------------------------

# Finally, we will use the function `thresh_forest` to display the thresholds on
# a forest plot. We sort the rows of the plot to display those with smallest
# thresholds first; this is achieved using the `orderby` option.

# Display using a forest plot, along with 95% confidence intervals for LORs
# Create row labels
dat.rawl$lab <- rep(NA, nrow(dat.rawl))
for (i in 1:nrow(dat.rawl)) {
  dat.rawl$lab[i] <- paste0(dat.rawl$studyID[i], " (", dat.rawl$k[i], " vs. ", dat.rawl$b[i], ")")
}

# Calculate 95% CIs
dat.rawl$CI2.5 <- dat.rawl$lor + qnorm(.025)*sqrt(diag(V))
dat.rawl$CI97.5 <- dat.rawl$lor + qnorm(.975)*sqrt(diag(V))


# We will sort the plot, smallest thresholds first
absmin <- function(x) min(abs(x))
dat.rawl$absmin.thresh <- apply(thresh$thresholds[, c("lo", "hi")], 1, absmin)

# Plot

# pdf("./figures/02a Thrombolytics - study level.pdf", width = 13, height = 5.5)

thresh_forest(thresh, 
              y = lor, CI.lo = CI2.5, CI.hi = CI97.5, 
              label = lab, orderby = absmin.thresh, data = dat.rawl,
              CI.title = "95% Confidence Interval", y.title = "Log OR", 
              label.title = "Study (Contrast)", xlab = "Log Odds Ratio", 
              xlim = c(-3, 2), refline = 0, digits = 2)

# dev.off()

