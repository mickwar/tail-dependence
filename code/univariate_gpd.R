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
        exceed_raw = "numeric",         # raw exceedances
        exceed_ind = "numeric",         # the indices
        threshold_quantile = "numeric", # threshold quantile used
        threshold = "numeric",          # threshold value used
        proc = "list",                  # processed data (for declustering)
        posterior = "list"
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
        initialize = function(x, u_q, u){
            "
            @param x   numeric vector, the observations
            @param u_q numeric, the threshold quantile, should be in (0, 1)
            @param u   numeric, the threshold
            "
            x <<- x

            if (missing(u_q) && missing(u))
                stop("Specify either threshold_quantile or threshold.")

            if (!missing(u_q) && !missing(u))
                warning(
                    paste0("Both the quantile and the absolute value for ",
                        "threshold were provided when only 1 needs to be ",
                        "given. Be sure these quantities agree.")
                    )

            if (missing(u_q)){
                warning("Threshold quantile set to default of 0.9")
                threshold_quantile <<- 0.9
            } else {
                threshold_quantile <<- u_q
                }

            if (missing(u)){
                threshold <<- quantile(x, threshold_quantile)
            } else {
                threshold <<- u
                }

            exceed_raw <<- x[x >= threshold] - threshold
            exceed_ind <<- which(x >= threshold)

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

            print("Estimating theta (extremal index).")

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

            print("Declustering exceedances based on mean estimate for theta.")
            
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
            proc$y <<- sapply(nice.C, max) - threshold

            # Do I need the index of where those declustered exceedances occur?
            # Need to remember how the declustering affect the bivariate analysis.

            },

        # @method
        # No zeta?
        fitGPD = function(m_ksi = 0, s_ksi = 10, nmcmc = 10000, nburn = 10000,
            window = 200, chainInit = c(1, 1e-6)){
            # Priors: ksi   ~ Normal(m_ksi, s_ksi^2)
            #         sigma ~ 1/sigma

            require(mwBASE)

            dat = list()

            # Check if decluster() was run. If so, fit the model on the
            # declustered exceedances. Otherwise, fit to the raw exceedances.
            if (is.null(proc$y)){
                print("Exceedances not declustered. Fitting model on raw exceedances.")
                dat$y = exceed_raw
                dat$n_c = length(dat$y)
            } else {
                print("Fitting model on declustered exceedances.")
                dat$y = proc$y
                dat$n_c = length(dat$y)
                }

            calcPosterior = function(dat, params){

                # dat is a list where dat$y is the exceedances, dat$n_c is the number of clusters
                sigma = params[1]
                ksi = params[2]

                # (lower) boundary check
                #if (any(1 + ksi*dat$y/sigma < 0))
                #    return (-Inf)

                # other boundary checks?
                if (ksi < 0 && max(dat$y) >= -sigma/ksi)
                    return (-Inf)
                if (sigma < 0)
                    return (-Inf)

                # Likelihood
                if (ksi != 0){
                    out = -dat$n_c*log(sigma) - (1 + 1/ksi)*sum(log(1+ksi*dat$y/sigma))
                } else {
                    out = -dat$n_c*log(sigma) - sum(dat$y)/sigma
                    }

                # Priors
                out = out - log(sigma)
                out = out + dnorm(ksi, m_ksi, s_ksi, log = TRUE)

                return (out)
                }

            posterior <<- mcmc_sampler(dat, calcPosterior, 2, nmcmc = nmcmc, nburn = nburn,
                nthin = 1, window = window, bounds = list("lower" = c(0, -Inf),
                "upper" = c(Inf, Inf)), chain_init = chainInit)

            posterior$sigma <<- as.numeric(posterior$params[,1])
            posterior$ksi <<- as.numeric(posterior$params[,2])

            },

        # Requires fitGPD() to be run
        # @method
        doPredictions = function(){

            },

        # @method
        returnLevels = function(){
            
            # Formula:
            # rl = u + sigma / ksi * [ (m * zeta * theta) ^ ksi - 1 ]
            #
            # rl is the quantity exceeded on average one every m observations.
            #
            # u is the threshold
            #
            # ksi is the shape for the GPD
            # sigma is the scale for the GPD
            #
            # zeta is the probability of exceeding the threshold
            # theta is the extremal index (probability of being in a cluster)

            }

        # Transform entire raw data vector to standard Frechet variates based
        # on the posterior parameters from the generalized Pareto distribution
        # fit.G
        # @method
        transformToFrechet = function(){
            if (is.null(posterior))
                stop("Must run UnivariateGPD$fitGPD() first.")

            # This looks a transformation to Pareto, not to Frechet. I got
            # this from the old code, but I don't understand it.
            #ksi = mean(posterior$ksi)
            #sigma = mean(posterior$sigma)
            #y = (1 + ksi/sigma * (x - threshold)) ^ (1/ksi)

            }

        )
    )

#x = rnorm(10000)

set.seed(1)
n = 40000
#x = rnorm(n)
x = double(n)
x[1] = rnorm(1, 0, 10)
for (i in 2:n)
    x[i] = 0.9 * x[i-1] + rnorm(1, 0, 10)




obj = UnivariateGPD(x, u_q = 0.95)

obj$estimateTheta(likelihood = "ferro", method = "bayesian")
obj$decluster()
obj$fitGPD()


mean(obj$proc$theta)
plot_hpd(obj$posterior$params[,1], main = "sigma", col1 = 'red')
plot_hpd(obj$posterior$params[,2], main = "ksi", col1 = 'blue')
