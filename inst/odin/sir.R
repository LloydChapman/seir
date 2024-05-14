## Definition of the time-step and output as "time"
dt <- user(1)
initial(time) <- 0
update(time) <- (step + 1) * dt
steps_per_day <- 1 / dt

## Core equations for transitions between compartments:
update(S) <- S - n_SI
update(I) <- I + n_SI - n_IR
update(R) <- R + n_IR

## Individual probabilities of transition:
p_SI <- 1 - exp(-beta * I / N * dt) # S to I
p_IR <- 1 - exp(-gamma * dt) # I to R

## Draws from binomial distributions for numbers changing between
## compartments:
n_IR <- rbinom(I, p_IR)
n_SI <- rbinom(S, p_SI)

## Total population size
N <- S + I + R

## Initial states:
initial(S) <- S_ini
initial(I) <- I_ini
initial(R) <- 0

# Output incidence
update(I_inc) <- if (step %% steps_per_day == 0) n_SI else I_inc + n_SI
initial(I_inc) <- 0

## User defined parameters - default in parentheses:
S_ini <- user(1000)
I_ini <- user(10)
beta <- user(0.2)
gamma <- user(0.1)

# In-line compare function
cases <- data()
compare(cases) ~ poisson(I_inc)
