## Definition of the time-step and output as "time" 
dt <- user()
initial(time) <- 0
update(time) <- (step + 1) * dt
steps_per_day <- 1/dt

## Core equations for transitions between compartments: 
update(S_tot) <- S_tot - sum(n_SE)
update(E_tot) <- E_tot + sum(n_SE) - sum(n_EI)
update(I_tot) <- I_tot + sum(n_EI) - sum(n_IR)
update(R_tot) <- R_tot + sum(n_IR)
update(D_tot) <- D_tot + sum(n_ID)

## Equations for transitions between compartments by age group
update(S[]) <- S[i] - n_SE[i]
update(E[]) <- E[i] + n_SE[i] - n_EI[i]
update(I[]) <- I[i] + n_EI[i] - n_IR[i]
update(R[]) <- R[i] + n_IR[i]
update(D[]) <- D[i] + n_ID[i]

# ## To compute cumulative incidence
# update(cum_inc[]) <- cum_inc[i] + n_SE[i] # THIS MAY ONLY BE OK IF step IS NOT A MULTIPLE OF 1/dt -- CHECK!
# update(cum_inc_death[]) <- cum_inc_death[i] + n_ID[i]

## Individual probabilities of transition:
p_SE[] <- 1 - exp(-lambda[i] * dt) # S to E
p_EI <- 1 - exp(-sigma * dt) # E to I
p_Ix <- 1 - exp(-gamma * dt) # I to R or D
p_death[] <- user()

## Force of infection
m[, ] <- user() # age-structured contact matrix
s_ij[, ] <- m[i, j] * I[i]
lambda[] <- beta * sum(s_ij[, i])

## Draws from binomial distributions for numbers changing between 
## compartments:
n_SE[] <- rbinom(S[i], p_SE[i])
n_EI[] <- rbinom(E[i], p_EI)
n_Ix[] <- rbinom(I[i], p_Ix)
n_IR[] <- rbinom(n_Ix[i], 1 - p_death[i])
n_ID[] <- n_Ix[i] - n_IR[i]

## Initial conditions 
initial(S_tot) <- sum(S_ini) 
initial(E_tot) <- sum(E_ini)
initial(I_tot) <- sum(I_ini) 
initial(R_tot) <- 0
initial(D_tot) <- 0

initial(S[]) <- S_ini[i]
initial(E[]) <- E_ini[i]
initial(I[]) <- I_ini[i]
initial(R[]) <- 0
initial(D[]) <- 0
# initial(cum_inc[]) <- 0
# initial(cum_inc_death[]) <- 0

## Store infection and death incidence
initial(I_inc) <- 0
initial(D_inc) <- 0
update(I_inc) <- if (step %% steps_per_day == 0) sum(n_EI) else I_inc + sum(n_EI)
update(D_inc) <- if (step %% steps_per_day == 0) sum(n_ID) else D_inc + sum(n_ID)

## Store infection and death incidence by age
initial(I_inc_0_9) <- 0
update(I_inc_0_9) <- if (step %% steps_per_day == 0) n_EI[1] else I_inc_0_9 + n_EI[1]
initial(I_inc_10_19) <- 0
update(I_inc_10_19) <- if (step %% steps_per_day == 0) n_EI[2] else I_inc_10_19 + n_EI[2]
initial(I_inc_20_29) <- 0
update(I_inc_20_29) <- if (step %% steps_per_day == 0) n_EI[3] else I_inc_20_29 + n_EI[3]
initial(I_inc_30_39) <- 0
update(I_inc_30_39) <- if (step %% steps_per_day == 0) n_EI[4] else I_inc_30_39 + n_EI[4]
initial(I_inc_40_49) <- 0
update(I_inc_40_49) <- if (step %% steps_per_day == 0) n_EI[5] else I_inc_40_49 + n_EI[5]
initial(I_inc_50_59) <- 0
update(I_inc_50_59) <- if (step %% steps_per_day == 0) n_EI[6] else I_inc_50_59 + n_EI[6]
initial(I_inc_60_69) <- 0
update(I_inc_60_69) <- if (step %% steps_per_day == 0) n_EI[7] else I_inc_60_69 + n_EI[7]
initial(I_inc_70_plus) <- 0
update(I_inc_70_plus) <- if (step %% steps_per_day == 0) n_EI[8] else I_inc_70_plus + n_EI[8]

initial(D_inc_0_39) <- 0
update(D_inc_0_39) <- if (step %% steps_per_day == 0) sum(n_ID[1:4]) else D_inc_0_39 + sum(n_ID[1:4])
initial(D_inc_40_49) <- 0
update(D_inc_40_49) <- if (step %% steps_per_day == 0) n_ID[5] else D_inc_40_49 + n_ID[5]
initial(D_inc_50_59) <- 0
update(D_inc_50_59) <- if (step %% steps_per_day == 0) n_ID[6] else D_inc_50_59 + n_ID[6]
initial(D_inc_60_69) <- 0
update(D_inc_60_69) <- if (step %% steps_per_day == 0) n_ID[7] else D_inc_60_69 + n_ID[7]
initial(D_inc_70_plus) <- 0
update(D_inc_70_plus) <- if (step %% steps_per_day == 0) n_ID[8] else D_inc_70_plus + n_ID[8]

## Model parameters (default in parenthesis) 
S_ini[] <- user()
E_ini[] <- user()
I_ini[] <- user()
beta <- user(0.0165)
sigma <- user(0.1)
gamma <- user(0.1)

## Array dimensions
N_age <- user()
dim(S_ini) <- N_age
dim(E_ini) <- N_age
dim(I_ini) <- N_age
dim(S) <- N_age
dim(E) <- N_age
dim(I) <- N_age
dim(R) <- N_age
dim(D) <- N_age
# dim(cum_inc) <- N_age
# dim(cum_inc_death) <- N_age
dim(n_SE) <- N_age
dim(n_EI) <- N_age
dim(n_Ix) <- N_age
dim(n_IR) <- N_age
dim(n_ID) <- N_age
dim(p_SE) <- N_age
dim(lambda) <- N_age
dim(m) <- c(N_age, N_age)
dim(s_ij) <- c(N_age, N_age)
dim(p_death) <- N_age
