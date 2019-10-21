# Mean residual life (MRL) helps to determine an appropriate exceedance threshold.
# When the plot begins to level out is about where the threshold should be.
#
# The MRL of a random variable X is defined as
#     m(t) = E(X - t | X > t),
# which is to the say the expected value of the exceedance (residual), given
# that the random variable is greater than a certain value (the threshold).
#
# In this class, we estimate the MRL based on a data vector. The data is
# sorted and each observation is then used as a `t` in the above equation.
# The mean and standard deviations are computed.

MRL = setRefClass(
    "MeanResidualLife",

    fields = list(
        x = "numeric",
        n = "numeric",
        sort_x = "numeric",
        mrl_est = "numeric",
        mrl_sig = "numeric"
        ),

    methods = list(

        # Initialize the data object
        initialize = function(x){
            x <<- x
            n <<- length(x)
            sort_x <<- sort(x)
            mrl_est <<- double(n-1)
            mrl_sig <<- double(n-1)
            calcMRL()
            },

        # There is probably a more optimized way of doing these
        # calculations, such as just working with a running mean
        # and sd.
        calcMRL = function(){
            for (i in 1:(n-1)){
                mrl_est[i] <<- mean(sort_x[(i+1):n] - sort_x[i])
                mrl_sig[i] <<- sd(sort_x[(i+1):n] - sort_x[i])
                }
            },

        # Function for plotting
        doPlot = function(){
            ylim = c(min(mrl_est - 1.96 * mrl_sig, na.rm = TRUE),
                     max(mrl_est + 1.96 * mrl_sig, na.rm = TRUE))
            plot(sort_x[-n], mrl_est, type = 'l', ylim = ylim)
            lines(sort_x[-n], mrl_est - 1.96 * mrl_sig, lty = 2)
            lines(sort_x[-n], mrl_est + 1.96 * mrl_sig, lty = 2)
            }
        )
    )


### Examples
x = rnorm(1000)
obj = MRL$new(x)
obj$doPlot()

x = rgamma(5000, 3, 1/2)
obj = MRL$new(x)
obj$doPlot()
