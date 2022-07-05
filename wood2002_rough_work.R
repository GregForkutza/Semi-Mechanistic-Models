## recruitment function (Gaussian function)
R1 <- function( t, r_floor, r_height, r_ctr, r_scale) {
    r_floor+r_height*exp(-((t-r_ctr)/r_scale)^2)
}

n_stage <- 12

## stage duration
gamma <- rep(1, n_stage)

##death rate
delta<-function(t, d0) {
  d0
}

## recursion for recruitment function 
r <- function(i,...) {
  if (i ==1) return (R1(...))
  value<- v[i-1]*n[i-1]
  return(value)
}

##parms( what else gets included here)
parms0 <- c(r_floor=0, r_height=1, r_ctr=0, r_scale=1, size=1, prob =0.5, d0=1)
n_parms0 <- length(parms0)
parms1 <- c(parms0, gamma)

library(deSolve)
gradfun1 <- function(t, n, parms) {
    gamma <- parms[-(1:n_parms0)]
    pp <- as.list(parms1[1:n_parms0])  ## extract named parameters
    attach(pp)                         ## make them parameters available
    on.exit(detach(pp))  ## remove them from the search list when function ends
    grad <- rep(NA_real_, n_stage) 
    r <- R1(t, r_floor, r_height, r_ctr, r_scale)
    grad[1] <- r-gamma[1]*n[1]-d0*n[1]
    for (i in 2:n_stage) {
        r <- gamma[i-1]*n[i-1]
        grad[i] <- r-gamma[i]*n[i]-d0*n[i]
    }
    return(list(grad))  ## must return a list containing the gradient as the 1st element
}
gradfun1(t = 0, n = rep(1, n_stage), parms = parms1)

if (FALSE) {
    
    ## attach
    pp <- list(aaa = 1, bbb = 2)
    attach(pp)  ## if this is repeated
    aaa
    bbb
    search()
    detach(pp)

    ## example of with() being convenient:

    grad <- with(my_parm_and_state_vector_list, {
        incidence <- beta*S*I
        recovery <- gamma*I
        ## the last expression evaluated in the block is the gradient vector
        c(S = -incidence, I = incidence - recovery, R = recovery)
        ## this is the 'result' or the return value of with()
        ## so I assign it to 'grad'
    })
    ## in our case we needed
    grad <- with(..., {
        grad <- rep(NA, ...)
        for (...) {
            grad[i] <- ...
        }
        grad ## this is the last expression
    })
}
