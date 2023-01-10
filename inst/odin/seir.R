## Definition of the time-step and output as "time" 
dt <- user(1)
initial(time) <- 0
update(time) <- (step + 1) * dt

## Model parameters (default in parenthesis) 
beta <- user(0.2)
sigma <- user(0.1)
gamma <- user(0.1)

## Initial conditions 
initial(S) <- 1000 
initial(E) <- 10 
initial(I) <- 0 
initial(R) <- 0

## Core equations for transitions between compartments: 
update(S) <- S - n_SE
update(E) <- E + n_SE - n_EI
update(I) <- I + n_EI - n_IR
update(R) <- R + n_IR

## Individual probabilities of transition: 
N <- S + E + I + R # total population size
p_SE <- 1 - exp(-beta * I / N * dt) # S to I 
p_EI <- 1 - exp(-sigma * dt) # E to I
p_IR <- 1 - exp(-gamma * dt) # I to R

## Draws from binomial distributions for numbers changing between 
## compartments:
n_IR <- rbinom(I, p_IR)
n_EI <- rbinom(E, p_EI)
n_SE <- rbinom(S, p_SE)