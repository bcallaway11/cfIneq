##-------------------------------------------------------##
##########################################################|
##                                                        |
################### COUNTERFACTUALS#######################|
##                                                        |
##########################################################|
##-------------------------------------------------------##

#' @title compute.counterfactual
#'
#' @description The internal code for computing counterfactual distributions
#'
#' @param outcome A vector of outcomes
#' @param wfun1 A function that takes in 'data' and returns a vector
#'  of the same length as 'outcome' that
#'  contains the first set of weights
#' @param wfun2 A function that takes in 'data' and returns a vector
#'  of the same length as 'outcome' that
#'  contains the first set of weights
#' @param data A data.frame for which the functions wf1(data) and
#'  wf2(data) will be called
#' @param tau The values at which to compute quantile effects, this should
#'  be a numeric vector with values strictly between 0 and 1
#' @param numyvals The number of values of the outcome to use in computing
#'  counterfactual distributions.  The default is to use all of them; it
#'  may speed up computation to set this equal to some smaller number (e.g.
#'  100)
#'  @param ... extra arguments to pass to the weighting functions
#'
#' @return A list that contains the difference in average outcomes using
#'  each set of weights, the distribution functions using each set of
#'  weights, the quantiles (and their difference) using each set of weights,
#'  tau, and the sequence of y values used to construct the distributions
#'
#' @keywords internal
#' @export
compute.counterfactual = function(outcome, wfun1, wfun2, data,
                                  tau=seq(.1,.9,.1), numyvals=NULL,
                                  ...) {

    weights1 <- wfun1(data)
    weights2 <- wfun2(data)

    ate <- BMisc::getWeightedMean(outcome, weights1) - BMisc::getWeightedMean(outcome, weights2)

    y.seq <- NULL ## this will eventually use all
    if (!is.null(numyvals)) {
        y.seq <- seq(quantile(outcome, .01, type="1"), quantile(outcome, .99, type="1"), length.out=numyvals)
    }

    Df1 <- BMisc::getWeightedDf(outcome, y.seq, weights1)
    Df2 <- BMisc::getWeightedDf(outcome, y.seq, weights2)

    q1 <- BMisc::getWeightedQuantiles(tau, outcome, weights=weights1)
    q2 <- BMisc::getWeightedQuantiles(tau, outcome, weights=weights2)

    out <- cfObj(ate=ate, Df1=Df1, Df2=Df2, q1=q1, q2=q2, qte=(q1-q2),
                tau=tau, y.seq=y.seq)
    return(out)

}

#' @title counterfactual
#'
#' @description a function for computing the counterfactual distribution
#'  of outcomes that would occur if untreated individuals received their
#'  current outcome but the distribution of outcomes for treated individuals
#'  is replaced by the counterfactual distribution of outcomes that
#'  they would have if they had the same conditional distribution of outcomes
#'  untreated individuals but retain their same characteristics.
#'
#' As an exammple, Callaway and Collins (2017) consider what a counterfactual
#'  economy with no unions in 1950 would look like.  Here the counterfactual
#'  distribution of earnings combines the observed distribution of earnings
#'  for non-union individuals with a counterfactual distribution of earnings
#'  for union individuals.  The latter counterfactual distribution replaces
#'  the conditional on covariates distribution of earnings for union individuals
#'  with the same one for non-union individuals, but holds constant the
#'  distribution of characteristics for union individuals.  Here,
#'  characteristics are things like age and education.
#'
#' @inheritParams compute.counterfactual
#'
#' @param se whether or not to compute standard errors
#' @param iters the number of bootstrap iterations to use to compute standard
#'  errors
#' @param pl whether or not to bootstrap in parallel
#' @param cores how many cores to uses when boostrapping in parallel
#'
#' @examples
#'
#' @return QTE object
#'
#' @export
counterfactual <- function(outcome, wfun1, wfun2, data,
                           tau=seq(.1,.9,.1), numyvals=NULL, se=TRUE,
                           iters=100, pl=FALSE, cores=1, ...) {

    ## qp <- qte::QTEparams(formla=formla, xformla=xformla, w=w,
    ##                 data=data, probs=probs, se=se, iters=iters,
    ##                 alp=alp, plot=plot, method=method,
    ##                 retEachIter=retEachIter, seedvec=seedvec,
    ##                 pl=pl, cores=cores, bootstrapiter=FALSE)

    ##first calculate the actual estimate
    counterfactual.res <- compute.counterfactual(outcome, wfun1,
                                                 wfun2, data, tau,
                                                 numyvals,
                                                 ...)

    if (se) {
        n <- length(outcome)
        #cat("Bootstrapping Standard Errors...\n")
        bootres <- pbapply::pblapply(1:iters, function(i) {
            booti <- sample(1:n,n,replace=TRUE)
            compute.counterfactual(outcome[booti], wfun1, wfun2,
                                   data[booti,], tau, numyvals)
        })

        ate.se <- sd(sapply(bootres, function(b) b$ate))
        qte.se <- apply(t(sapply(bootres, function(b) b$qte)), 2, sd)
        counterfactual.res$ate.se <- ate.se
        counterfactual.res$qte.se <- qte.se
        return(counterfactual.res)
    } else {
        return(counterfactual.res)
    }

}


#' @title cfObj
#'
#' @description Counterfactual objects
#'
#' @param ate Difference in average outcomes for 1st set of weights
#'  relative to second set of weights
#' @param Df1 Distribution function for first set of weights
#' @param Df2 Distribution function for second set of weights
#' @param q1 Quantiles for first set of weights
#' @param q2 Quantiles for second set of weights
#' @param qte q1-q2
#' @param tau Values at which quantiles are computed
#' @param y.seq Values at which the distribution functions are computed;
#'  if this is NULL, the default is to use all values of y
#' @param ate.se Standard error for ate
#' @param qte.se Standard errors for qte
#'
#' @return cfObj
#'
#' @export
cfObj <- function(ate=NULL, Df1=NULL, Df2=NULL, q1=NULL, q2=NULL,
                  qte=NULL, tau=NULL, y.seq=NULL, ate.se=NULL, qte.se=NULL) {

    out <- list(ate=ate, Df1=Df1, Df2=Df2, q1=q1, q2=q2, qte=(q1-q2),
                tau=tau, y.seq=y.seq, ate.se=ate.se, qte.se=qte.se)
    out
}

#' @title ggcf
#'
#' @description Plot counterfactual quantiles using ggplot2
#'
#' @import ggplot2
#'
#' @param cfObj A counterfactual object
#' @param plotate Whether or not to plot the ate too, default is FALSE
#'
#' @return ggplot object
#'
#' @export
ggcf <- function(cfObj, plotate=FALSE) {
    cmat <- data.frame(qte=cfObj$qte, qte.se=cfObj$qte.se, tau=cfObj$tau)

    p <- ggplot(data=cmat, mapping=aes(x=tau, y=qte)) +
        geom_line() +
        geom_point() +
        geom_line(aes(x=tau, y=qte+1.96*qte.se), linetype="dashed") +
        geom_line(aes(x=tau, y=qte-1.96*qte.se), linetype="dashed") +
        scale_x_continuous(breaks=cmat$tau) +
        theme_bw()
    if (plotate) {
        p <- p + geom_hline(yintercept=cfObj$ate, color="red")
        p <- p + geom_hline(yintercept=cfObj$ate + 1.96*cfObj$ate.se,
                            color="red", linetype="dashed")
        p <- p + geom_hline(yintercept=cfObj$ate - 1.96*cfObj$ate.se,
                            color="red", linetype="dashed")
    }
    p
}
