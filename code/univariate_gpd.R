# Class for doing a majority of the univariate threshold exceedance model
# Includes:
#   - Estimating the extremal index
#   - Declustering the observations
#   - Fitting de-clustered exceedances to the generalize Pareto distribution
#   - Calculating return levels

UnivariateGPD = setRefClass(
    Class = "UniGPD",

    fields = list(
        x = "numeric",                  # raw data
        threshold_quantile = "numeric", # threshold quantile used
        threshold = "numeric",          # threshold value used
        proc = "list"                   # processed data (for declustering)
        ),

    # proc is list with the following tags:
    #   exceed_index    index for exceedances (based on x)
    #   Tu              interexceedance times
    #   N               length of Tu
    #   m1              number of interexceedance times which equals 1
    #   probNotExceed   empirical probability of not exceeding threhold
    #   likelihood      either "ferro" or "suveges"
    #   method          either "classical" or "bayesian"
    #   classicEst      if method == "classical", the function for estimating
    #                   theta, based on the given likelihood
    #   theta           if method == "classical", then a length-1 numeric
    #                   representing the mean, else a vector of MCMC samples
    #                   of the posterior theta
    #   theta_vec       if method == "classical", a bootstrap sample of
    #                   estimated thetas
    #   T_C             list of clustered observations
    #   y               vector of independent, de-cluster exceedances
    

    methods = list(

        # Combine with MRL?

        # @method
        initialize = function(x, threshold_quantile, threshold){
            x <<- x

            if (!missing(threshold_quantile) && !missing(threshold))
                warning(
                    paste0("Both the quantile and the absolute value for ",
                        "threshold were provided when only 1 needs to be ",
                        "given. Be sure these quantities agree.")
                    )

            if (missing(threshold_quantile)){
                warning("Threshold quantile set to default of 0.9")
                threshold_quantile <<- 0.9
            } else {
                threshold_quantile <<- threshold_quantile
                }

            if (missing(threshold))
                threshold <<- quantile(x, threshold_quantile)

            # Option to just run all the functions on initialization?

            },


        # De-trending function?

        # Transforming to standard Frechet? (keep the indices if so)

        # Estimate the extremal index using either Ferro or Suveges likelihood.
        # Account for breaks in the time-series (e.g., in summer-only data, the
        # end of one summer and the beginning of next cannot have an
        # interexceedance time)
        # @method
        estimateTheta = function(
            likelihood = c("ferro", "suveges"),
            method = c("classical", "bayesian"),
            prior,
            ...){
            #' @param likelihood    character
            #' @param method        character
            #' @param prior         list, prior values for the bayesian method, the two
            #'                      likelihoods expect different items in the list
            #' @param ...           additional arguments passed to mwBASE::mcmc_sampler

            likelihood = match.arg(likelihood)
            method = match.arg(method)

            ### Set up
            # Index of exceedances
            proc$exceed_index <<- which(x > threshold)

            # Interexceedance times (has length one less than exceed)
            proc$Tu <<- diff(proc$exceed_index)
            proc$N <<- length(proc$exceed_index)

            # Number of interexceedance times equal to 1.
            proc$m1 <<- sum(proc$Tu == 1)

            # Probability of not exceeding
            proc$probNotExceed <<- mean(x <= threshold)

            proc$likelihood <<- likelihood
            proc$method <<- method

            ### Estimating theta (extremal index)
            # Classical estimation
            if (proc$method == "classical"){
                if (proc$likelihood == "ferro"){
                    proc$classicEst <<- function(Tu){
                        if (length(Tu) == 0)
                            return (1)
                        if (max(Tu) <= 2){
                            out = min(1, 2*(sum(Tu))^2 / ( length(Tu) * sum(Tu^2) ) )
                        } else {
                            out = min(1, 2*(sum(Tu - 1))^2 / (length(Tu)*sum((Tu-1)*(Tu-2))))
                            }
                        return (out)
                        }
                    }
                if (proc$likelihood == "suveges"){
                    proc$classicEst <<- function(Tu){
                        ((1-proc$probNotExceed)*sum(Tu - 1)+2*(proc$N-1)-proc$m1 -
                            sqrt(((1-proc$probNotExceed)*sum(Tu - 1) +
                                2*(proc$N-1)-proc$m1)^2 -
                                8*(proc$N-1-proc$m1)*(1-proc$probNotExceed)*sum(Tu-1))) / 
                            (2*(1-proc$probNotExceed)*sum(Tu-1))
                        }
                    }

                proc$theta <<- proc$classicEst(proc$Tu)
                }

            # Bayesian estimation
            if (proc$method == "bayesian"){

                # Gets the functions for doing the MCMC sampling
                require(mwBASE)

                # Default values for prior
                if (missing(prior)){
                    if (proc$likelihood == "ferro")
                        prior = list("theta_a" = 1, "theta_b" = 1/2,
                            "p_a" = proc$probNotExceed*100, "p_b" = (1-proc$probNotExceed)*100)
                    if (proc$likelihood == "suveges")
                        prior = list("theta_a" = 1, "theta_b" = 1/2)
                    }

                # Build the target function (i.e. the posterior) under each case.
                # Must be in accordance with mwBASE::mcmc_sampler
                if (proc$likelihood == "ferro"){
                    calcPosterior = function(dat, param){
                        theta = param[1]
                        p = param[2]
                        if (theta <= 0 || theta > 1)
                            return (-Inf)
                        if (p <= 0 || p >= 1)
                            return (-Inf)

                        # Likelihood (Ferro and Segers, Eq. 3)
                        out = proc$m1 * log(1 - theta*p^theta) +
                            (proc$N - 1 - proc$m1) * (log(theta) +
                            log(1-p^theta)) + theta*log(p)*sum(dat - 1)

                        # Priors
                        out = out + dbeta(theta, prior$theta_a, prior$theta_b, log = TRUE)
                        out = out + dbeta(p, prior$p_a, prior$p_b, log = TRUE)
                        return (out)
                        }
                    }
                if (proc$likelihood == "suveges"){
                    calcPosterior = function(dat, param){
                        theta = param[1]
                        if (theta <= 0 || theta > 1)
                            return (-Inf)

                        # Likelihood (Suveges, Eq. 1)
                        out = proc$m1 * log(1 - theta) +
                            2*(proc$N - 1 - proc$m1)*log(theta) -
                            theta * (1-proc$probNotExceed) * sum(dat - 1)

                        # Priors
                        out = out + dbeta(theta, prior$theta_a, prior$theta_b, log = TRUE)
                        return (out)
                        }
                    }

                # Do the sampling
                mcmc_out = mcmc_sampler(data = proc$Tu, target = calcPosterior,
                    nparam = ifelse(proc$likelihood == "ferro", 2, 1), ...)
                mcmc_out$param = as.matrix(mcmc_out$param)

                proc$theta <<- mcmc_out$param[,1]

                }

            },

        # estimateTheta() should be run first, return error otherwise
        # @method
        decluster = function(){
            if (is.null(proc$exceed_index))
                stop("Must run UnivariateGPD$estimateTheta() first")
            
            # C is the number of clusters
            C = floor(mean(proc$theta)*proc$N)+1
            C = min(C, proc$N-1)

            tmp = sort(proc$Tu, decreasing = TRUE)
            T_C = tmp[C]
            while (!(tmp[C-1] > T_C) && (C > 1)){
                C = C - 1
                T_C = tmp[C]
                }

            # The set of independent intercluster times
            inter.Clust = proc$Tu[proc$Tu > T_C]

            i_j = which(proc$Tu > T_C)
            i_j = c(0, i_j, proc$N)
            ind.seq = rep(list(NULL), C)
            intra.Clust = rep(list(NULL), C)    # The interexceedance times within each cluster
            nice.S = rep(list(NULL), C)         # List of independent clusters, marking when exceedances occur
            nice.C = rep(list(NULL), C)         # The observed value at the exceedance times

            for (k in 2:(C+1)){
                ind.seq[[k-1]] = seq(i_j[k-1]+1, i_j[k]-1)
                if (i_j[k-1]+1 == i_j[k]){
            #       nice.T[[j-1]] = NULL
                } else {
                    intra.Clust[[k-1]] = proc$Tu[seq(i_j[k-1]+1, i_j[k]-1)]
                    }
                nice.S[[k-1]] = proc$exceed_index[seq(i_j[k-1]+1, i_j[k])]
                nice.C[[k-1]] = x[nice.S[[k-1]]]
                }

            ### Bootstrap (for classical method)
            if (proc$method == "classical"){
                theta_vec = double(bsamp)
                for (i in 1:bsamp){

                    samp.inter = sample(C-1, replace = TRUE)
                    samp.intra = sample(C, replace = TRUE)

                    tmp = c(inter.Clust[samp.inter], unlist(intra.Clust[samp.intra]))
                    theta_vec[i] = proc$classicEst(tmp)
                    }

                proc$theta_vec <<- theta_vec
                }

            proc$T_C <<- T_C

            # Get the greatest value within each (independent) cluster
            proc$y <<- sapply(nice.C, max)


            },

        # @method
        fitGPD = function(){

            },

        # Requires fitGPD() to be run
        # @method
        doPredictions = function(){

            },

        # @method
        returnLevels = function(){
            

            }
        )
    )

#x = rnorm(10000)

n = 4000
x = double(n)
x[1] = rnorm(1)
for (i in 2:n)
    x[i] = 0.9 * x[i-1] + rnorm(1)

plot(x, type = 'l')


obj = UnivariateGPD$new(x, threshold_quantile = 0.9)

obj$estimateTheta(likelihood = "ferro", method = "bayesian")
obj$decluster()

hist(obj$proc$theta)
mean(obj$proc$theta)

obj$proc$y
