# Class for computing asymptotic tail dependence.
# Chi(q) = Pr(F_1(W_1) > q | F_2(W_2) > q) (W_1, W_2 are pareto processes)
# W(s) = Y * V(s),
#   Y is standard Pareto,
#   V(s) is stochastic process with sup V(s) = w_0 and E(V(s)) > 1
#
# lim q -> Inf Chi(q) = Chi = 

BivariateTail = setRefClass(
    "BivariateTail",

    # Using the list option for `fields` had problems with using the "BDPdensity"
    # class for `bdp`. Opted here to use the generic approach.
    fields = c(
        "cone",
        "angles",
        "chi_dat",    # chi calculated from data
        "chi_pp",      # chi calcualted from posterior predictive
        "bdp",
        "pred_angle",
        "pred_cone",
        "pred_y"
        ),
    #fields = list(
    #    cone = "matrix",
    #    angles = "numeric",
    #    chi_dat = "numeric",    # chi calculated from data
    #    bdp = "BDPdensity",
    #    pred_angle = "numeric",
    #    pred_cone = "matrix",
    #    chi_pp = "numeric"      # chi calcualted from posterior predictive
    #    ),

    methods = list(

        initialize = function(cone = NULL, angles = NULL){
            # cone   - n by 2 matrix with values in [0, 1] and the maximum
            #          between each pair (row) is 1.
            # angles - vector of angles from the origin to the cone in radians,
            #          should be in [0, pi/2]
            #
            # Able to handle raw values / magnitudes for the extreme values?

            if (is.null(cone) && is.null(angles))
                stop("One of `cone` or `angles` must be given.")

            # If the angles are given, calculate the their points on the cone
            if (is.null(cone)){
                angles <<- angles
                cone <<- angleToCone(angles)
                }

            # If the cone is given, calculate the scaled angles
            if (is.null(angles)){
                cone <<- cone
                angles <<- coneToAngle(cone)
                }

            # Compute chi just from the data
            chi_dat <<- calculateChi(cone)
            },

        # Convert angles on [0, pi/2] to the non-negative unit circle with
        # the supremum norm (i.e. {(x, y): max(x, y) == 1, x, y >= 0})
        angleToCone = function(angles){
            out = matrix(1, length(angles), 2)
            out[,1] = tan(pi/4 - abs(angles - pi/4))
            ind = which(angles <= pi/4)
            out[ind,] = out[ind, c(2,1)]
            return (out)
            },

        # The inverse of the above function
        coneToAngle = function(cone){
            out = atan2(cone[,2], cone[,1])
            return (out)
            },

        # Compute Chi
        calculateChi = function(cone){
            ind = (cone[,1] > 0 & cone[,2] > 0)
            z1 = cone[ind, 1] / mean(cone[, 1])
            z2 = cone[ind, 2] / mean(cone[, 2])
            z = cbind(z1, z2)
            m1 = apply(z, 1, min)
            return (mean(m1) * mean(ind))
            },

        # Fit Bernstein-Dirichlet prior to the angles
        doBDP = function(mcmc, prior, nsamples = 10000, scaled = FALSE){

            # DPpackage is deprecated
            require(DPpackage)

            # mcmc     - list that is passed to BDPdensity
            # prior    - list that is passed to BDPdensity
            # nsamples - number of posterior samples to use
            # scaled   - boolean, were the angles scaled to [0, 1]?
            if (missing(mcmc)){
                mcmc = list(nburn = 1000,
                    nsave = 1000,
                    nskip= 1,
                    ndisplay = 100)
                }
            if (missing(prior)){
                prior = list(aa0 = 1,
                    ab0 = 0.1,
                    kmax = 50,
                    a0 = 1,
                    b0 = 1)
                }

            if (scaled){
                scAngles = angles
            } else {
                scAngles = angles * 2/pi
                }

            # Don't fit the BDP with values too close to the edges
            ind_rm_0 = (scAngles <= 0.005)
            ind_rm_1 = (scAngles >= 0.995)
            ind_rm = (ind_rm_0 | ind_rm_1)

            # Fit Bernstein-Dirichlet prior
            # This can take a while
            cat(paste0("Fitting a Bernstein polynomial Dirichlet prior of order ", prior$kmax, "."))
            cat("\nThis may take a while (about 10 minutes).")
            bdp <<- BDPdensity(y = scAngles[!ind_rm],
                prior = prior, mcmc = mcmc,
                state = NULL, status = TRUE, support = 1)

            # Posterior predictive distribution
            pred_angle <<- sample(bdp$grid, nsamples, replace = TRUE, prob = bdp$fun)

            # Manually insert 0's and 1's
            pred_ind = sample(c(0,1,2), nsamples, replace = TRUE,
                prob = c(mean(ind_rm_0), mean(ind_rm_1), 1-mean(ind_rm)))
            pred_angle[pred_ind == 0] <<- 0
            pred_angle[pred_ind == 1] <<- 1

            # Convert posterior angle to posterior cone
            pred_cone <<- angleToCone(pred_angle * pi/2)

            # Calculate chi from the posterior predictive
            chi_pp <<- calculateChi(pred_cone)

            # The data should have previously been transformed to
            # a stadard Frechet, so now we can get posterior samples
            # of the exceedances
            # The support should be the positive quadrant without the
            # unit square [0, 1] * [0, 1].
            pred_v = 1/runif(nsamples)
            pred_y <<- pred_v * pred_cone

            # Compute chi as u approaches 1
            #u = seq(0.8, 0.999, length = 30)
            #calculateChi

            },


        # Plotting
        doPlot = function(){

            # TODO: The Y's (standard paretos) are the exceedances
            #denom1 = sapply(u, function(x) mean(Y * V[,1] > mean(V[,1]) / (1-x)) )
            #denom2 = sapply(u, function(x) mean(Y * V[,2] > mean(V[,2]) / (1-x)) )

            ## Plot marginal
            #plot(u, denom1, type = 'l', xlim = c(0, 1), ylim = c(0, 1), lwd = 2, col = 'blue')
            #lines(u, denom2, col = 'red', lty = 2, lwd = 2)
            #abline(1, -1)
            #abline(v = 1-mean(V[,1]), lty = 2, col = 'gray50')
            #abline(h = mean(V[,1]), lty = 2, col = 'gray50')
            #abline(v = 1-mean(V[,2]), lty = 2, col = 'gray50')
            #abline(h = mean(V[,2]), lty = 2, col = 'gray50')

            ## Compute joint probability
            #Vm1 = cbind(mean(V[,1]) / V[,1], mean(V[,2]) / V[,2])
            #Vm2 = cbind(V[,1] / mean(V[,1]), V[,2] / mean(V[,2]))
            #Vmax = apply(Vm1, 1, max)
            #Vmin = apply(Vm2, 1, min)
            #numer = sapply(u, function(x) mean(Y > Vmax / (1-x)) )
            }
        )
    )


angles = rbeta(1000, 2, 4) * pi/2
cone = matrix(1, length(angles), 2)
cone[,1] = tan(pi/4 - abs(angles - pi/4))
ind = which(angles <= pi/4)
cone[ind,] = cone[ind, c(2,1)]


obj = BivariateTail$new(cone, NULL)
#test2 = BivariateTail$new(NULL, angles)

head(obj$angles)
head(obj$cone)
obj$chi_dat

obj$doBDP()
