##############################################################################
#
# allocate.r --- version 1.0.0
#
# Andrew P. Robinson
#
# 17/12/2010
#
# This file includes the 'inspect' function, which can be used to estimate
# the level of inspection effort required to provide a statistical guarantee
# that the future risk (as defined in the article) will be no higher than a
# specified cutoff with given probability.
#
# The code to recreate the results reported in Table 2 in the article 
# is appended.
# 
##############################################################################


inspect <- function(failed,             
                    inspected, 
                    approach,
                    risk.cutoff = 0.01, 
                    confidence = 0.95,
                    conservatism = 0.95,
                    verbose = FALSE) {
  if (inspected < 1 ||
      approach < 1 ||
      failed / inspected > risk.cutoff) {
    rate <- 1
  } else {
    n.tentative <- function(x) {       
      k.hat <- x * qbeta(conservatism,
                         shape1 = x * failed / inspected + 0.5, 
                         shape2 = x * (1 - failed / inspected) + 0.5)
      qbeta(confidence,
            shape1 = k.hat + 0.5, 
            shape2 = x - k.hat + 0.5) - risk.cutoff
    }
    if (n.tentative(approach) > 0) {
      rate <- 1
    } else {
      if (n.tentative(1) < 0 & n.tentative(approach) < 0) {
        rate <- 0
      } else {
        ideal.sample.size <- uniroot(f = n.tentative,
                                     interval = c(1, approach))$root
        rate <- min(ideal.sample.size/approach, 1)
      }
    }
  } 
  if (verbose) {
    sample.size <- rate * approach
    k.tilde <- sample.size *
      qbeta(conservatism,
            shape1 = sample.size * failed / inspected + 0.5, 
            shape2 = sample.size * (1 - failed / inspected) + 0.5)
    projected.risk <- qbeta(confidence,
                            shape1 = k.tilde + 0.5, 
                            shape2 = sample.size - k.tilde + 0.5)
    return(list(rate = rate,
                sample.size = sample.size,
                failed = failed,
                inspected = inspected,
                approach = approach,
                risk.cutoff = projected.risk,
                conservatism = conservatism,
                confidence = confidence,
                k.tilde = k.tilde))
  } else {
    return(rate)
  }
}

# ### South-east Queensland
# 
# inspect(58, 37743, 37743, risk.cutoff = 0.01) * 100
# 
# ### Far North Queensland
# 
# inspect(33, 2957, 2957, risk.cutoff = 0.01) * 100
# 
# ### New South Wales
# 
# inspect(137, 207764, 207764, risk.cutoff = 0.01) * 100
# 
# ### South Australia
# 
# inspect(59, 17510, 17510, risk.cutoff = 0.01) * 100
# 
# ### Victoria
# 
# inspect(24, 91491, 91491, risk.cutoff = 0.01) * 100
# 
# ### Western Australia
# 
# inspect(0, 14067, 14067, risk.cutoff = 0.01) * 100
# 
# ### National
# 
# inspect(311, 371532, 371532, risk.cutoff = 0.01) * 100

