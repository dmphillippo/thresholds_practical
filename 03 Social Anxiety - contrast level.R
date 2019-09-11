## 03. Social Anxiety example - contrast level

library(nmathresh)


# Step 1: NMA and data prep -----------------------------------------------

# The results of the NMA are available in the `nmathresh` package, as
# `SocAnx.post.summary` (for the posterior summaries) and `SocAnx.post.cov` (for
# the posterior covariance matrix).

# SocAnx.post.summary
# SocAnx.post.cov


# For the contrast-level analysis, we need the contrast design matrix. Since
# there are many studies, we won't construct this by hand (though this is
# possible, just laborious); instead, we use the treatment details from the
# study data.

trt.dat <- read.delim(system.file("extdata", "SocAnx_data.txt", package = "nmathresh"))[, 1:6]

head(trt.dat)   # Print first few rows

K <- with(trt.dat, length(unique(c(t.1, t.2, t.3, t.4, t.5)))) - 1  # Number of treatments , -1 for NA

# work out which contrasts have data
contr.ab <- data.frame(a = c(), b = c())

for (i in 1:nrow(trt.dat)) {
   rowi <- trt.dat[i, which(!is.na(trt.dat[i, 2:6])) + 1] # non NA elements of ith row
   
   # get contrast from all combinations of treatments
   trtcomb <- combn(rowi, 2, function(x) sapply(x, as.numeric))
   
   a <- apply(trtcomb, 2, min)
   b <- apply(trtcomb, 2, max)
   
   # remove contrasts of treatments against themselves
   iseq <- a == b
   a <- a[!iseq]
   b <- b[!iseq]
   
   if (!all(iseq)) contr.ab <- rbind(contr.ab, cbind(a, b))
}

contr.ab <- unique(contr.ab[order(contr.ab$a, contr.ab$b), ])


# Contrast design matrix
X <- matrix(0, nrow = nrow(contr.ab), ncol = K - 1)
for (i in 1:nrow(X)) {
  if (contr.ab[i, "a"] > 1) X[i, contr.ab[i, "a"] - 1]  <- -1
  if (contr.ab[i, "b"] > 1)   X[i, contr.ab[i, "b"] - 1]    <- 1
}

# As the posterior summaries contain several variables, we pick out the indices
# of those which we need.

# Get indices of d, sd, diff in the CODA data
vnames <- sub("(.*)\\[.*","\\1", rownames(SocAnx.post.summary$statistics))
ind.d <- which(vnames == "d")
ind.sd <- which(vnames == "sd")
ind.diff <- which(vnames == "diff")
ind.delta <- which(vnames == "delta")


# Step 2: Reconstruct two-stage likelihood covariance matrix --------------

# Class model is used, so use the prior mean precision from the gamma
# distribution
prior.prec <- rep(3.9/0.35, 40)
# Other than for treatment 3
prior.prec[2] <- 0.0001

lik.cov <- recon_vcov(SocAnx.post.cov[1:(K-1), 1:(K-1)], prior.vcov = diag(1/prior.prec), X = X)

# The KL divergence indicates that the reconstructed likelihood covariance
# matrix results in a posterior that is reasonably close to the true posterior
# (values less than 1 indicate negligible differences, values greater than 3
# indicate considerable differences).


# Step 3: Derive thresholds -----------------------------------------------

thresh <- nma_thresh(mean.dk = SocAnx.post.summary$statistics[ind.d, "Mean"], 
                     lhood = lik.cov,
                     post = SocAnx.post.cov[1:(K-1), 1:(K-1)],
                     nmatype = "fixed", 
                     X = X, 
                     opt.max = FALSE)

# Q: How would you implement the following decision rules?
#      - A "do not do" decision, for the worst-ranked treatment
#      - A decision based on a MCID of 0.2 (a "small" effect)


# Step 4: Plot ------------------------------------------------------------

# Get indices of contrasts in likelihood
d.a <- d.b <- vector(length = nrow(X))
for (i in 1:nrow(X)) {
  d.a[i] <- ifelse(any(X[i, ] == -1), which(X[i, ] == -1), 0) + 1
  d.b[i] <- ifelse(any(X[i, ] == 1), which(X[i, ] == 1), 0) + 1
}

d.i <- d_ab2i(d.a, d.b, K = K)

# Contrast labels and credible intervals (from the posterior summary)
plotdat <- data.frame(
  contr = paste0(d.b, " vs. ", d.a),
  contr.mean = SocAnx.post.summary$statistics[ind.diff[d.i], "Mean"],
  CI2.5 = SocAnx.post.summary$quantiles[ind.diff[d.i], "2.5%"],
  CI97.5 = SocAnx.post.summary$quantiles[ind.diff[d.i], "97.5%"]
)

# We'll sort the plot to show contrasts with smallest thresholds first, and only
# show contrasts with thresholds < 1.5 SMD
cutoff <- 1.5
absmin <- function(x) min(abs(x))
plotdat$ord <- ifelse(thresh$thresholds$lo > -cutoff | thresh$thresholds$hi < cutoff,
                      apply(thresh$thresholds[, c("lo", "hi")], 1, absmin), NA)

# pdf("./figures/03 Social Anxiety - contrast level.pdf", width = 12, height = 8)

thresh_forest(thresh, contr.mean, CI2.5, CI97.5, 
              label = contr, orderby = list(ord, na.last = NA), data = plotdat,
              label.title = "Contrast", xlab = "SMD", xlim = c(-4, 3), 
              CI.title = "95% Credible Interval",
              refline = 0, digits = 2, calcdim = FALSE)

# dev.off()

# Q: What does the threshold analysis tell us? How would this affect the treatment recommendation?
# Hint: a SMD of 0.8 is considered "large"


# Step 5: More complex analyses -------------------------------------------

# More complex analyses are possible by manipulating the set of $u_{ak^*,m}$
# values, which are contained in the `thresh` object output by `nma_thresh` in
# the matrix `Ukstar`. For example, continuing with the contrast-level analysis
# of the social anxiety NMA, we could consider a common bias in all
# pharmacological treatments against inactive control, or all psychological
# treatments against inactive control.

# Firstly, consider a common pharmacological treatment bias. The overall
# influence on the contrast between treatments $a$ and $k^*$ of a common
# adjustment to all drug data points is found simply by summing over the
# individual influences of each drug data point (since a single common
# adjustment is to be made to all efficacy estimates). Therefore the point where
# treatment $a$ replaces $k^*$ as optimal is found by summing over the
# individual threshold solutions $u_{ak^*,m}$ for each drug data point $m$.

# Drug treatments (+ combi therapies)
drugtrts <- c(9:23, 37:41)

# Which data points compare these to an inactive trt?
drugdats <- which(contr.ab$a %in% 1:3 & contr.ab$b %in% drugtrts)

# Get U solutions by summing the indivual solutions of drug data points
U.drugs <- 1 / (rowSums(1 / thresh$Ukstar[,drugdats]))

# Which contrasts do the rows of Ukstar correspond to?
Ukstar.ab <- d_i2ab(1:(K*(K-1)/2), K)
Ukstar.ab <- Ukstar.ab[Ukstar.ab$a == thresh$kstar | Ukstar.ab$b == thresh$kstar, ]


# We then derive threshold values and new optimal treatments at each threshold
# by picking out the minimum positive value and maximum negative value from the
# overall drug threshold solution vector `U.drugs`. The function `get.int` makes
# this simple, with some additional handling of infinite thresholds behind the
# scenes.

# Thresholds are then
thresh.drugs <- get.int(U.drugs, thresh$kstar, 1:K, Ukstar.ab)

# Now we plot the invariant interval, along with the pharmacological treatment
# estimates.

## Function to plot the common invariant interval with the data
plotII <- function(thresh, contr.mean, CrI.lo, CrI.hi, rowlabs, xlim, xlab, ylab, ...){
  
  yaxs <- length(contr.mean):1
  
  # split plot in two
  layout(matrix(1:2,nrow=1), widths=c(.2,1))
  
  # plot row labels on left side
  gp <- par(mar=c(5,4,1,0))
  plot(rep(0,length(yaxs)), yaxs, pch=NA, ylim=c(.5,yaxs[1]+.5), ylab="",
       xlab="", yaxt="n", xaxt="n", bty="n")
  text(0, yaxs, labels=rowlabs,xpd=NA)
  title(ylab=ylab, line=2)
  
  # fake plot for setup of right side
  par(mar=c(5,1,1,2))
  plot(contr.mean, yaxs, pch=NA, yaxt="n", xlim=xlim,
       ylim=c(.5,yaxs[1]+.5), ylab="", xlab="",...)
  title(xlab=xlab, line=2)
  
  # reference line
  abline(v=0, lty=2, col="gray")
  
  # combined invariant region
  polygon(rep(c(contr.mean + thresh$lo, rev(contr.mean) + thresh$hi), each=2), 
          c(rep(yaxs,each=2) + c(.5,-.5), rep(rev(yaxs),each=2) + c(-.5,.5)),
          col=rgb(.7,.8,.9,.7),
          border=rgb(.6,.7,.8))
  
  # credible intervals
  segments(y0=yaxs, x0=CrI.lo, x1=CrI.hi, lend=1)
  
  # contrast means
  points(contr.mean, yaxs, pch=21, col="black", bg="white")
  
  # write new optimal treatments at either side
  text(xlim[1], round(length(yaxs)/2), 
       labels=as.expression(substitute(tilde(k)*"* = "*xx,list(xx=thresh$lo.newkstar))), 
       pos=4)
  text(xlim[2], round(length(yaxs)/2),
       labels=as.expression(substitute(tilde(k)*"* = "*xx,list(xx=thresh$hi.newkstar))),
       pos=2)
  
  # write invariant interval below plot
  with(thresh, title(xlab=paste0("Invariant interval about zero: ",lo.newkstar,
                                 " (",formatC(lo,digits=2,format="f"),", ",
                                 formatC(hi,digits=2,format="f"),") ",
                                 hi.newkstar), line=3))
  
  par(gp)
}


plotII(thresh.drugs, 
       SocAnx.post.summary$statistics[ind.diff[drugdats], "Mean"], 
       SocAnx.post.summary$quantiles[ind.diff[drugdats], "2.5%"],
       SocAnx.post.summary$quantiles[ind.diff[drugdats], "97.5%"],
       rowlabs = paste0(contr.ab[drugdats, "b"], " vs. ", contr.ab[drugdats, "a"]),
       xlim = c(-4, 1.5), ylab = "Drug vs. Inactive Contrasts", xlab = "SMD")


# Similarly for a common psychological treatment bias:

# Psych treatments
psychtrts <- c(4:8, 24:36)

# Which data points compare these to an inactive trt?
psychdats <- which(contr.ab$a %in% 1:3 & contr.ab$b %in% psychtrts)

# Get U solutions by summing the influences of drug data points
U.psych <- 1 / (rowSums(1 / thresh$Ukstar[,psychdats]))

# Thresholds are then
thresh.psych <- get.int(U.psych, thresh$kstar, 1:K, Ukstar.ab)

plotII(thresh.psych, 
       SocAnx.post.summary$statistics[ind.diff[psychdats], "Mean"], 
       SocAnx.post.summary$quantiles[ind.diff[psychdats], "2.5%"],
       SocAnx.post.summary$quantiles[ind.diff[psychdats], "97.5%"],
       rowlabs=paste0(contr.ab[psychdats,"b"]," vs. ",contr.ab[psychdats,"a"]),
       xlim=c(-3,2), ylab="Psych vs. Inactive Contrasts", xlab="SMD")

# The magnitude of the thresholds for both common pharmacological and common
# psychological treatment biases are large. Any such biases - if they exist -
# are likely to be much smaller than these thresholds.

# We may also assess the impact of adjusting for common biases of this nature
# simultaneously.

# Instead of deriving the threshold lines by hand using `U.drugs` and `U.psych`,
# we will assemble a bare-bones `thresh` object for input to `thresh_2d`, to
# take advantage of the automated plotting routines. The `Ukstar` matrix has a
# column for the threshold solutions for the combined drug data points
# `U.drugs`, and a column for threshold solutions for the combined psychological
# data points `U.psych`. `thresh_2d` also uses `kstar`, so we include that too.
# Internally, `thresh_2d` uses the `Ukstar` matrix to derive the threshold lines
# with gradient `- U.psych / U.drugs` and intercept `U.psych`.

thresh.drugpsych <- list(
  Ukstar = cbind(U.drugs, U.psych),
  kstar = thresh$kstar
)

thresh_2d(thresh.drugpsych, 1, 2,
          xlab = "SMD adjustment: Drug vs. Inactive",
          ylab = "SMD adjustment: Psych vs. Inactive",
          xlim = c(-6, 2), ylim = c(-6, 2),
          breaks = -6:2)
