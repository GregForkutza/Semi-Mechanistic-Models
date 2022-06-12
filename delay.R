## unpacking parameter and state vectors

## method 1, brute force
## tedious but clear.
gradfun <- function(t, x, parms) {
    S <- x[1]
    I <- x[2]
    R <- x[3]
    beta <- parms[1]
    gamma <- parms[2]
    ## etc ...
}

## named vectors
x <- c(S = 1, I = 0.001, R = 0)
p <- c(beta = 2, gamma = 1)

## 'atomic' vector: is an ordered set of values all of the same 'atomic' type
##  (numeric, integer, character, logical, complex, ...)
## can have names
## internally, stored as a contiguous block of memory

## list: is an ordered set of elements, not necessarily homogeneous, not necessarily
##   atomic types
## (most complex functions in R return their answers as a list)
m1 <- lm(mpg ~ hp, mtcars)
str(m1)

with( ## with(data, expr)
      ## evaluates the expression in an environment containing the elements of `data`
    as.list(## converts the combined, named *vector* to a *list*
            ## why? so we can pass it to with()
        c(  ## combine/catenate
            ## the only way this goes wrong is if I have states and parameters
            ##  with overlapping names!
            x, p
        )
    ),
    { ## curly brackets if I want to execute more than one statement
        incidence <- beta*S*I
        gradient <- incidence - gamma*I
    }
)

## using a pipe (|>)
(c(x, p)
    |> as.list()
    |> with({
        incidence <- beta*S*I
        gradient <- incidence - gamma*I
    })
)

## recruitment function
R1 <- function(t, floor, height, ctr, scale) {
    floor + height*exp(-((t-ctr)/scale)^2)
}
## recursion
ri <- function(t, i, ...) {
    if (i == 1) return(R1(...))
    ## index
    value <- ri(t-tau[i-1], i-1, ...)*S(t-tau[i-1], t, i-1, ...)
    return(value)
}
## (2) solve the recursion, write it out as a summation
##  sum (R1(t-tau)*exp(-sum(delta_i(tau)*tau_i))

##
## (3) return r, S as auxiliary variables
## and use lagvalues() to retrieve their values at previous times
## most of the time when using ode() etc. you return just list(grad) from
## the gradient function


## stages 1-6
delta1 <- function(t, d0, d1) {
    d0 * exp(-d1*t)
}
## stages 7-11
delta2 <- function(t, d0) {
    0
}
## visualize
curve(r(x, 1, 2, 1, 2), from = 0, to = 6, ylim = c(0, 3))
abline(h = c(1, 3), lty = 2)
abline(v = 1, lty = 2)
## log scale (set floor == 0 so we get a quadratic)
curve(r(x, floor = 0, 2, 1, 2), from = 0, to = 6, log = "y")

## parms <- c(r_floor = 1, r_height= 2, ...)
with(as.list(c(parms)), { ## don't include x here because it's not named
    ## because we want to use tau[i], pull out tau vector as well
    tau <- tail(parms, 12)
    grad <- rep(NA_real_, 12)  ## empty gradient vector (not numeric(12))
    r <- rep(NA_real_, 12)  ## empty recruitment vector
    r[1] <- R1(t, r_floor, r_height, r_ctr, r_scale)
    grad[1] <- r[1] - delta1(t, d0_1, d1_1)*x[1]
    for (i in 2:6) {
        ## r_1(t) = R1(t); r_{i+1}(t) = r_i(t-tau_i)*S_i(t-tau_i, t)
        ## → r_i(t) = r_{i-1}(t-tau_{i-1})*S_{i-1}(t-tau_{i-1}, t)
        ## n_i(t) = r_i(t) - r_i(t-tau_i)*S_i(t-tau_i, t) - d_i(t)*n_i(t)
        grad[i] <- r(t
                     r_floor, r_height, r_ctr, r_scale)
        - r(t-tau*(i-1),
                     r_floor, r_height, r_ctr, r_scale)
        ## S_i(a,b)
    }
    for (i in 7:12) {
        grad[i] <- ...
    }
})

tau_vec <- c(0.1, 0.5, 0.3, 0.5)
p <- c(r_floor = 1, r_height = 2, tau = tau_vec)

library(deSolve)
##
## lag < 0
## state vector: densities of 11 stages ?

sirgrad <- function(t, x, parms) {
    with(as.list(c(x,parms)),
    {
        incidence <- beta*S*I
        recovery <- gamma*I
        grad <- c(S = -incidence,
                  I = incidence - recovery,
                  R = recovery)
        return(list(grad, incidence = incidence, recovery = recovery))
    })
}
r1 <- ode(c(S=0.99, I = 0.01, R = 0),
          func = sirgrad,
          parms = c(beta = 2, gamma = 1),
          times = seq(0, 10, length.out = 101)
          )
head(r1)

ddesirgrad <- function(t, x, parms) {
    with(as.list(c(x,parms)),
    {
        incidence <- beta*S*I
        ## browser()
        recovery <- if (t<1/gamma) 0 else beta*lagvalue(t-1/gamma, 1)*lagvalue(t-1/gamma, 2)
        grad <- c(S = -incidence,
                  I = incidence - recovery,
                  R = recovery)
        return(list(grad, incidence = incidence))
    })
}
r2 <- dede(c(S=0.99, I = 0.01, R = 0),
          func = ddesirgrad,
          parms = c(beta = 2, gamma = 1),
          times = seq(0, 10, length.out = 101)
          )

matplot(r2[,1], r2[,2:4], type = "l")
          
