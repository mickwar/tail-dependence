# Class for doing a majority of the univariate threshold exceedance model
# Includes:
#   - Estimating the extremal index
#   - Declustering the observations
#   - Fitting de-clustered exceedances to the generalize Pareto distribution
#   - Calculating return levels

UnivariateGPD = setRefClass(
    Class = "UniGPD",

    fields = list(
        x = "numeric",
        theta = "numeric"
        ),

    methods = list(

        # Combine with MRL?

        # @method
        initialize = function(x, threshold){
            x <<- x
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
            threshold_quantile = 0.9,
            threshold,
            prior){

            likelihood = match.arg(likelihood)
            method = match.arg(method)

            # If not given, set threshold based on the given quantile
            if (missing(threshold))
                threshold = quantile(x, threshold_quantile)

            ### Set up
            # Index of exceedances
            exceed = which(x > threshold)

            # Interexceedance times (has length one less than exceed)
            Tu = diff(exceed)
            N = length(exceed)

            # Number of interexceedance times equal to 1.
            m1 = sum(Tu == 1)

            # Probability of not exceeding
            probNotExceed = mean(y <= threshold)

            ### Estimating theta (extremal index)
            # Classical estimation
            if (method == "classical"){
                if (likelihood == "ferro"){
                    classicEst = function(Tu){
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
                if (likelihood == "suveges"){
                    classicEst = function(Tu){
                        ((1-probNotExceed)*sum(Tu - 1)+2*(N-1)-m1 -
                            sqrt(((1-probNotExceed)*sum(Tu - 1)+2*(N-1)-m1)^2 - 8*(N-1-m1)*(1-probNotExceed)*sum(Tu-1))) / 
                            (2*(1-probNotExceed)*sum(Tu-1))
                        }
                    }

                theta <<- classicEst(Tu)
                }

            # Bayesian estimation
            if (method == "bayesian"){

                # Gets the functions for doing the MCMC sampling
                require(mwBASE)

                # Default values for prior
                if (missing(prior)){
                    if (likelihood == "ferro")
                        prior = list("theta_a" = 1, "theta_b" = 1/2,
                            "p_a" = probNotExceed*100, "p_b" = (1-probNotExceed)*100)
                    if (likelihood == "suveges")
                        prior = list("theta_a" = 1, "theta_b" = 1/2)
                    }

                # Build the target function (i.e. the posterior) under each case.
                # Must be in accordance with mwBASE::mcmc_sampler
                if (likelihood == "ferro"){
                    calcPosterior = function(dat, param){
                        theta = param[1]
                        p = param[2]
                        if (theta <= 0 || theta > 1)
                            return (-Inf)
                        if (p <= 0 || p >= 1)
                            return (-Inf)

                        # Likelihood (Ferro and Segers, Eq. 3)
                        out = m1 * log(1 - theta*p^theta) + (N - 1 - m1)*(log(theta) + log(1-p^theta)) +
                            theta*log(p)*sum(dat - 1)

                        # Priors
                        out = out + dbeta(theta, prior$theta_a, prior$theta_b, log = TRUE)
                        out = out + dbeta(p, prior$p_a, prior$p_b, log = TRUE)
                        return (out)
                        }
                    }
                if (likelihood == "suveges"){
                    calcPosterior = function(dat, param){
                        theta = param[1]
                        if (theta <= 0 || theta > 1)
                            return (-Inf)

                        # Likelihood (Suveges, Eq. 1)
                        out = m1 * log(1 - theta) + 2*(N - 1 - m1)*log(theta) - theta*(1-probNotExceed)*sum(dat - 1)

                        # Priors
                        out = out + dbeta(theta, prior$theta_a, prior$theta_b, log = TRUE)
                        return (out)
                        }
                    }

                # Do the sampling
                mcmc_out = mcmc_sampler(data = Tu, target = calcPosterior,
                    nparam = ifelse(likelihood == "ferro", 2, 1), ...)
                mcmc_out$param = as.matrix(mcmc_out$param)

                theta_hat = mean(mcmc_out$param[,1])

                }

            },

        # estimateTheta() should be run first, return error otherwise
        # @method
        decluster = function(){

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
            """

            """
            

            }
        )
    )

x = rnorm(100)
obj = UnivariateGPD$new(x, 0)

obj$returnLevels
