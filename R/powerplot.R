#' @title Power Plot
#' @description Plots powers for the Kruskall-Wallis test, via Monte Carlo and two approximations.
#' @param numgrps Number of groups to compare
#' @param thetadagger Direction of effect
#' @param nnvec vector of numbers per group.
#' @param nmc Number of Monte Carlo trials
#' @param targetpower Target power for test
#' @param level level for test.
#' @export
powerplot <- function (numgrps = 3, thetadagger = NULL, nnvec = 5:30, nmc = 50000, 
    targetpower = 0.8, level = 0.05)
{
    if (is.null(thetadagger)) {
        thetadagger <- (seq(numgrps) - 1)/2
    }
    out <- array(NA, c(length(nnvec), 5))
    out[,5]<-nnvec
    for (jj in seq(length(nnvec))) {
        Delta <- kweffectsize(totsamp = numgrps * nnvec[jj], 
            shifts = thetadagger, distname = "normal", targetpower = 0.8, 
            proportions = rep(1, length(thetadagger))/length(thetadagger), 
            level = 0.05)
        out[jj,4]<- Delta
        out[jj, 1] <- kwpower(nreps = rep(nnvec[jj], numgrps), 
            shifts = Delta * thetadagger, distname = "normal", 
            level = 0.05, mc = 0, taylor = TRUE)$power
        out[jj, 2] <- kwpower(nreps = rep(nnvec[jj], numgrps), 
            shifts = Delta * thetadagger, distname = "normal", 
            level = 0.05, mc = 0)$power
        out[jj, 3] <- kwpower(nreps = rep(nnvec[jj], numgrps), 
            shifts = Delta * thetadagger, distname = "normal", 
            level = 0.05, mc = nmc)$power
    }
    plot(range(nnvec), range(out[,1:3]), type = "n", xlab = "Number per Group", 
          ylab = "Approximate Power", sub = paste(numgrps, "groups, level 0.05, power 0.8"), 
          main = "Approximate Powers for the Kruskall-Wallis Test")
    legend(median(nnvec)*.5, sum(range(out[, 3])*c(.85,.15)), lty = 1:3, 
         legend = c("Approximation with True Expectation", 
        "Approximation with Approximate Expectation", "Monte Carlo"))
    for (ii in 1:3) lines(nnvec, out[, ii], lty = ii)
    return(invisible(out))
}
