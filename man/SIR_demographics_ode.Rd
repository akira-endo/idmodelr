% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SIR.R
\name{SIR_demographics_ode}
\alias{SIR_demographics_ode}
\title{Susceptible-Infected-Recovered Model with Simple Demographics}
\usage{
SIR_demographics_ode(t, x, params)
}
\arguments{
\item{t}{The timestep over which to calculate derivatives}

\item{x}{A numeric vector of compartment populations.}

\item{params}{A named vector of parameter values.}
}
\value{
A vector of derivatives
}
\description{
Susceptible-Infected-Recovered Model with Simple Demographics
}
\examples{
##Model Input
S_0 <- 989
I_0 <- 1
R_0 <- 0
beta <- 3
tau <- 2
mu <- 1/81


parameters <- c(beta = beta, tau = tau, mu = mu)
inits <- c(S = S_0, I = I_0, R_0 = R_0)

SIR_demographics_ode(1, inits, parameters)
}
