# Class for computing asymptotic tail dependence.
# Chi(q) = Pr(F_1(W_1) > q | F_2(W_2) > q) (W_1, W_2 are pareto processes)
# W(s) = Y * V(s),
#   Y is standard Pareto,
#   V(s) is stochastic process with sup V(s) = w_0 and E(V(s)) > 1
#
# lim q -> Inf Chi(q) = Chi = 

BivariateTail = setRefClass(
    "BivariateTail",

    fields = list(
        cone = "matrix",
        angles = "numeric",
        chi = "numeric"
        ),

    methods = list(

        initialize = function(x, iscone = TRUE){
            if (iscone){
                cone <<- x
                angles <<- coneToAngle(cone)
            } else {
                angles <<- x
                cone <<- angleToCone(angles)
                }
            calculateChi(cone)
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
            chi <<- mean(m1) * mean(ind)
            },

        # Plotting
        plot = function(){

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


test1 = BivariateTail$new(cone)
test2 = BivariateTail$new(angles, iscone = FALSE)

test1$chi

