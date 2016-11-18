##-------------------------------------------------------##
##########################################################|
##                                                        |
################### COUNTERFACTUALS#######################|
##                                                        |
##########################################################|
##-------------------------------------------------------##

#### some code for producing counterfactuals
##counterfactual is "world without unions"
## the counterfactual will go in the slot for treated outcomes
## the observed world will go in the slot for untreated outcomess
#' @return QTE object
compute.counterfactual = function(qp) {

    setupData(qp)
    bootstrapiter <- qp$bootstrapiter
    
    ##no covariate att - will update if there are covariates
    ##should I make any changes to weights here
    ## not sure if these are right here, but should use weights
    ## anyway
    att <- getWeightedMean(untreated[,yname], untreated[,wname]) -
        getWeightedMean(data[,yname], data[,wname])

    qte <- getWeightedQuantiles(probs, untreated[,yname],
                                weights=untreated[,wname]) -
        getWeightedQuantiles(probs, data[,yname],
                                weights=data[,wname])
        
    
    n = nrow(data)

    
    ##set these up to access later
    observed.quantiles <- NULL
    counterfactual.quantiles <- NULL
    counterfactual.weights <- data[,wname]
    pscore.reg <- NULL
    if (!is.null(x)) {

        ##estimate the propensity score
        pscore.reg <- glm(data[,treat] ~ as.matrix(data[,x]),
                          family=binomial(link=method))
        pscore <- fitted(pscore.reg)
        p = rep(nrow(treated)/(nrow(treated) + nrow(untreated)), n)
        D <- data[,treat]
        y <- data[,yname]
        w <- data[,wname]
        n <- nrow(data)
        ##there are alternatives for how to compute the quantiles of 
        ##treated outcomes for the treated group:
        ##1) compute quantiles directly
        observed.quantiles = getWeightedQuantiles(probs, data[,yname],
                                                  weights=data[,wname])

      
        counterfactual.weights = w*(1-D)*(1/(1-pscore))

        counterfactual.quantiles = getWeightedQuantiles(probs, y,
                                                        weights=counterfactual.weights)
        
        qte <- counterfactual.quantiles - observed.quantiles
    }

    F.treated.t <- getWeightedDf(data[,yname], weights=counterfactual.weights)
    F.treated.t.cf <-getWeightedDf(untreated[,yname], weights=untreated[,wname])
        
    if (bootstrapiter) { ## do this to decrease size of each iteration
        out <- QTE(qte=qte, ate=att, probs=probs)
    } else {
        out <- QTE(F.treated.t=F.treated.t,
                   F.treated.t.cf=F.treated.t.cf,
                   qte=qte, pscore.reg=pscore.reg,  ate=att, probs=probs)
    }
    return(out)
    
}


#' @return QTE object
counterfactual <- function(formla, xformla=NULL, w=NULL, data,
                    probs=seq(0.05,0.95,0.05), se=TRUE,
                 iters=100, alp=0.05, plot=FALSE, method="logit",
                 retEachIter=FALSE, seedvec=NULL, pl=FALSE, cores=1) {

    qp <- QTEparams(formla=formla, xformla=xformla, w=w,
                    data=data, probs=probs, se=se, iters=iters,
                    alp=alp, plot=plot, method=method,
                    retEachIter=retEachIter, seedvec=seedvec,
                    pl=pl, cores=cores, bootstrapiter=FALSE)
    
    ##first calculate the actual estimate
    counterfactual.qte <- compute.counterfactual(qp)

    if (se) {

        qp$bootstrapiter <- TRUE

        ##bootstrap the standard errors
        SEobj <- bootstrap(qp, counterfactual.qte, compute.counterfactual)

        out <- QTE(qte=counterfactual.qte$qte, qte.upper=SEobj$qte.upper,
                    qte.lower=SEobj$qte.lower, ate=counterfactual.qte$ate,
                    ate.upper=SEobj$ate.upper, ate.lower=SEobj$ate.lower,
                    qte.se=SEobj$qte.se, ate.se=SEobj$ate.se,
                    pscore.reg=counterfactual.qte$pscore.reg,
                    F.treated.t=counterfactual.qte$F.treated.t,
                    F.treated.t.cf=counterfactual.qte$F.treated.t.cf,
                    eachIterList=eachIter,
                    probs=probs)
        return(out)
     
    } else {
        return(counterfactual.qte)
    }

}
