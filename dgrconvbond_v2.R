rm(list = ls(all.names = TRUE))
# Binomial Tree Test--------------------------
# "Tomaybedo" - restructure variables/inputs and rewrite the code as a function

# Timer
ptm <- proc.time()

# Variables
stock_0 <- 0.008
outstanding_shares <- 1657657946

bond_fv <- 1200000
bond_conprice <- 0.0160
conv_ratio <- bond_fv / bond_conprice

coupon_rate <- 0.00
yrly_coupons <- 4
coupon_value <- bond_fv * coupon_rate / yrly_coupons

time_maturity <- 2 # specification in years
n_steps <- 100 # number of steps in tree

# the following code identifies the nodes where coupon payments would occur
total_coupons <- yrly_coupons * time_maturity
coupon_times <- seq(0, n_steps, n_steps / total_coupons)
coupon_times <- round(coupon_times)
coupon_times <- coupon_times[1:total_coupons] # to remove the last coupon

coupon_values <- rep(0, n_steps)
coupon_values[coupon_times] <- coupon_value #gives the time vector of coupon payments

volatility <- 1.16
dt <- time_maturity / n_steps
u <- exp(volatility * sqrt(dt))
d <- 1 / u

risk_free <- 0.0162
div_yield <- 0.00
r_b <- risk_free - div_yield

# implementation of credit risk - with probability of default 
# at default, it is assumed that the value of the stock drops to 0

prob_default <- 0.507 #probability of default over the life of the conv bond
lambda_0 <- -log(1 - prob_default) / time_maturity
# prob_default_node <- 1 - exp(-lambda * dt) FROM CONSTANT INTENSITY MODEL
alpha <- -0.50

recovery_rate <- 0.00

# when stock prices fall below this lower bound, the probabilities of
# default, up and down movements will not belong to [0, 1]
lower_bound <- stock_0 * 
  (
    (1 / (lambda_0 * dt)) * 
     log((u - recovery_rate) / (exp(risk_free * dt) - recovery_rate))
  ) ^ (1 / alpha)

# r_blambda_dt <- exp((r_b + lambda) * dt) FROM CONSTANT INTENSITY MODEL
# q = (r_blambda_dt - d) / (u - d) FROM CONSTANT INTENSITY MODEL

rf_discountfactor <- exp(-risk_free * dt)

# creation of binomial trees --------------------------------------------------

binomial_tree_stock  <- data.frame() # initialise data frame for binomial tree
binomial_tree_stock[1, 1] <- stock_0 # set initial value for binomial tree

# the following binomial tree grows downwards (rows represent time) 
# and to the right (representing up movements)

for (tree_time in 2:n_steps) {
  binomial_tree_stock[tree_time, 1] <- binomial_tree_stock[tree_time -1, 1] * d
  for (nodes_at_time in 2:tree_time) {
    binomial_tree_stock[tree_time, nodes_at_time] <- binomial_tree_stock[tree_time - 1,
                                                                         nodes_at_time - 1] * u
  }
}

lower_bound_tree <- binomial_tree_stock < lower_bound
test <- colSums(lower_bound_tree, na.rm = TRUE)
test <- test[test > 0]
length(test) / n_steps

sum(rowSums(lower_bound_tree, na.rm = TRUE)) / sum(seq(1, n_steps))
# looks like with 100 steps, 19% of all stock prices fall below the lower bound.
# 500 steps, 31% of all stock prices fall below the lower bound.

# creates matrix of default intensities (lambdas) based on stock price
intensity_tree <- lambda_0 * (binomial_tree_stock / stock_0) ^ alpha

# creates matrix of probability of default (at each node) based on matrix of default
# intensities
default_tree <- 1 - exp(-intensity_tree * dt)

# creates matrix of probability of up movements (q) and renames it for brevity
probability_tree <- (exp((r_b + intensity_tree) * dt) - d) / (u - d)
q <- probability_tree
q[binomial_tree_stock < lower_bound] <- 0
# overrides values in the stock price tree where the stock price is below the lower
# bound. (sets values to 0)
# [binomial_tree_stock < lower_bound] <- 0

# overrides values in the probability of default tree where the stock price is
# below the lower bound. (sets prob of default to 1)
default_tree[binomial_tree_stock < lower_bound] <- 1

# creates adjusted stock price binomial tree with 0.83 VWAP factor
a_binomial_tree_stock <- binomial_tree_stock * 0.83
a_binomial_tree_stock[a_binomial_tree_stock > bond_conprice] <- bond_conprice

# creates a tree of new shares introduced by conversion using adjusted stock price tree
new_shares <- binomial_tree_stock 
new_shares[new_shares > 0 ] <- bond_fv
new_shares <- new_shares / a_binomial_tree_stock  

# creates a new tree of stock value prices which have been adjusted for the value of
# dilution and bond redemption
payoff_tree_stock <- (outstanding_shares * binomial_tree_stock + bond_fv) / 
                     (outstanding_shares + new_shares)

# overrides payoff tree stock to set values where the original share price was below
# the lower bound to 0
payoff_tree_stock[binomial_tree_stock < lower_bound] <- 0

# set convertible note value at final node
# assumed final node includes the final coupon payment + fv
binomial_tree_payoff <- binomial_tree_stock
binomial_tree_payoff[n_steps, ] <- pmax(bond_fv + coupon_value, 
                                        payoff_tree_stock[n_steps, ] * conv_ratio)

# overrides binomial tree payoff to the recovery amount at nodes below the lower bound.
# TODO - fix up this equation to actually use the recovery amount... will need to discount
binomial_tree_payoff[binomial_tree_stock < lower_bound] <- 0

continuation_tree <- binomial_tree_payoff # to initialise continuation values data frame 

# this loop calculates the payoffs for the convertible bond. It begins by calculating
# "continuation values" which are values of the convertible bond assuming it is not
# redeemed. The payoff is then calculated by comparing the continuation value with
# the value of the note if it is redeemed
for (tree_time in (n_steps - 1):1) {
  for (nodes_at_time in 1:tree_time) {
    continuation_tree[tree_time, nodes_at_time] <- 
      rf_discountfactor * (1 - default_tree[tree_time, nodes_at_time]) * 
       (
        q[tree_time, nodes_at_time] * 
         binomial_tree_payoff[tree_time + 1, nodes_at_time + 1] 
       + 
         (1 - q[tree_time, nodes_at_time]) * 
         binomial_tree_payoff[tree_time + 1, nodes_at_time]
       )
    +
      (1 - default_tree[tree_time, nodes_at_time]) * coupon_values[tree_time]
    +
      default_tree[tree_time, nodes_at_time] * rf_discountfactor * 
      recovery_rate * bond_fv * exp(-risk_free * (n_steps - tree_time))
    
    binomial_tree_payoff[tree_time, nodes_at_time] <-
      pmax(continuation_tree[tree_time, nodes_at_time],
           payoff_tree_stock[tree_time, nodes_at_time] * conv_ratio)
  }
}

binomial_tree_payoff[1, 1] #conv bond value at time 0

# Stop time
proc.time() - ptm

# checks and matrices --------------------------------------------------

# is the continuation value > conversion value??? -- matrix
conversion_test <- binomial_tree_payoff

for (tree_time in (n_steps):1) {
  for (nodes_at_time in 1:tree_time) {
    conversion_test[tree_time, nodes_at_time] <- 
      if (continuation_tree[tree_time, nodes_at_time] > 
          payoff_tree_stock[tree_time, nodes_at_time] * conv_ratio) 
      1
    else
      0
  }
}

# percentage of nodes where continuation value > conversion value
sum(rowSums(conversion_test, na.rm = TRUE)) / sum(seq(1, n_steps))

# difference between continuation value and conversion value
conversion_tree <- conv_ratio * payoff_tree_stock
difference_test <- continuation_tree - conversion_tree

