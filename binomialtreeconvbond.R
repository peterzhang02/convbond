setwd(readClipboard())


rm(list = ls(all.names = TRUE))
# Binomial Tree Test--------------------------

# Timer
ptm <- proc.time()

# Variables
stock_0 <- 100

time_maturity <- 1 # specification in years
n_steps <- 100 # number of steps in tree
volatility <- 0.25
dt <- time_maturity / n_steps
u <- exp(volatility * sqrt(dt))
d <- 1 / u

risk_free <- 0.05
div_yield <- 0.05
risk_free_df <- exp(-risk_free * dt)
r_b <- risk_free - div_yield

r_b_dt <- exp(r_b * dt)
q = (r_b_dt - d) / (u - d)

strike <- 100

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

# set option price at end of tree (equal to payoff function)

binomial_tree_payoff <- binomial_tree_stock
binomial_tree_payoff[n_steps, ] <- binomial_tree_stock[n_steps, ] - strike
binomial_tree_payoff[binomial_tree_payoff < 0] <- 0

option_values <- binomial_tree_payoff # to initialise option values data frame 

for (tree_time in (n_steps-1):1) {
  for (nodes_at_time in 1:tree_time) {
    option_values[tree_time, nodes_at_time] <- (q * option_values[tree_time + 1,
                                                                  nodes_at_time + 1] +
                                                  (1 - q) * option_values[tree_time + 1,
                                                                          nodes_at_time]) *
      risk_free_df
  }
}

option_values[1, 1] #option value at time 0

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