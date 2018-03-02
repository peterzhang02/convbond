#setwd(readClipboard())

rm(list = ls(all.names = TRUE))
# Binomial Tree Test--------------------------

# Timer
#ptm <- proc.time()

# Variables
stock_0 <- 100

bond_fv <- 100000
bond_conprice <- 125
conv_ratio <- bond_fv / bond_conprice

coupon_rate <- 0.05
yrly_coupons <- 4
coupon_value <- bond_fv * coupon_rate / yrly_coupons

time_maturity <- 2 # specification in years
n_steps <- 100 # number of steps in tree

# the following code identifies the nodes where coupon payments would occur
total_coupons <- yrly_coupons * time_maturity
coupon_times <- seq(0, n_steps, n_steps / total_coupons)
coupon_times <- round(coupon_times)
coupon_times <- coupon_times[1:total_coupons] # to remove the first and last coupons

coupon_values <- rep(0, n_steps)
coupon_values[coupon_times] <- coupon_value

volatility <- 0.25
dt <- time_maturity / n_steps
u <- exp(volatility * sqrt(dt))
d <- 1 / u

risk_free <- 0.05
div_yield <- 0.05
r_b <- risk_free - div_yield

# implementation of credit risk - with probability of default 
# at default, it is assumed that the value of the stock drops to 0

lambda <- 0.02
prob_default <- 1 - exp(-lambda * dt)
recovery_rate <- 0.50

r_blambda_dt <- exp((r_b + lambda) * dt)
q = (r_blambda_dt - d) / (u - d)
rf_discountfactor <- exp(-risk_free * dt)

# creation of binomial trees --------------------------------------------------

binomial_tree_stock  <- data.frame() # initialise data frame for binomial tree
binomial_tree_stock[1, 1] <- stock_0 # set initial value for binomial tree

# the following binomial tree grows downwards (rows represent time) and to the right (representing
# up movements)

for (tree_time in 2:n_steps) {
  binomial_tree_stock[tree_time, 1] <- binomial_tree_stock[tree_time -1, 1] * d
  for (nodes_at_time in 2:tree_time) {
    binomial_tree_stock[tree_time, nodes_at_time] <- binomial_tree_stock[tree_time - 1,
                                                                         nodes_at_time -1] * u
  }
}

# set convertible note value at final node

binomial_tree_payoff <- binomial_tree_stock
binomial_tree_payoff[n_steps, ] <- pmax(bond_fv + coupon_value, 
                                        binomial_tree_stock[n_steps, ] * conv_ratio)

continuation_tree <- binomial_tree_payoff # to initialise continuation values data frame 

for (tree_time in (n_steps-1):1) {
  for (nodes_at_time in 1:tree_time) {
    continuation_tree[tree_time, nodes_at_time] <- 
      rf_discountfactor * (1 - prob_default) * 
      (q * continuation_tree[tree_time + 1, nodes_at_time + 1] + 
      (1 - q) * continuation_tree[tree_time + 1, nodes_at_time])
    +
      (1 - prob_default) * coupon_values[tree_time]
    +
      prob_default * rf_discountfactor * 
      recovery_rate * bond_fv * exp(-risk_free * (n_steps - tree_time))
  }
}

for (tree_time in (n_steps-1):1) {
  for (nodes_at_time in 1:tree_time) {
    binomial_tree_payoff[tree_time, nodes_at_time] <- 
      pmax(continuation_tree[tree_time, nodes_at_time], 
           binomial_tree_stock[tree_time, nodes_at_time] * conv_ratio)
  }
}

binomial_tree_payoff[1, 1] #conv bond value at time 0

# Black Scholes

denom <- volatility * sqrt(time_maturity)

d1 <- (log(stock_0 / strike) + (risk_free - div_yield + (volatility * volatility) / 2) * 
         time_maturity) / denom

d2 <- (log(stock_0 / strike) + (risk_free - div_yield - (volatility * volatility) / 2) * 
         time_maturity) / denom

euro_call <- exp(-div_yield * time_maturity) * stock_0 * pnorm(d1) -
  exp(-risk_free * time_maturity) * strike * pnorm(d2)

euro_call
option_values[1, 1] - euro_call

# Stop time
proc.time() - ptm





u*q + d*(1-q) - exp((r_b)*dt)

q*(u^2) + (1-q)*(d^2) - exp(2*risk_free*dt)+(volatility^2)*dt
