## Compartments E to R stratified by i and j with
## i for the age group
## j for the vaccine group

## Number of age and vaccination groups
n_age <- user()
n_vax <- user()

## Definition of the time-step and output as "time" 
dt <- user()
initial(time) <- 0
update(time) <- (step + 1) * dt
steps_per_day <- 1/dt

## Core equations for transitions between compartments:
update(S[, ]) <- S[i,j] - n_SE[i,j] - n_SV[i,j] + 
    (if (j == 1) n_SV[i,n_vax] else n_SV[i,j-1])
update(E[, ]) <- E[i,j] + n_SE[i,j] - n_EI[i,j] - n_EV[i,j] + 
    (if (j == 1) n_EV[i,n_vax] else n_EV[i,j-1])
update(I_P[, ]) <- I_P[i,j] + n_EI_P[i,j] - n_I_PI_C[i,j] - n_I_PV[i,j] + 
    (if (j == 1) n_I_PV[i,n_vax] else n_I_PV[i,j-1])
update(I_A[, ]) <- I_A[i,j] + n_EI_A[i,j] - n_I_AR[i,j] - n_I_AV[i,j] + 
    (if (j == 1) n_I_AV[i,n_vax] else n_I_AV[i,j-1])
update(I_C[, ]) <- I_C[i,j] + n_I_PI_C[i,j] - n_I_Cx[i,j]
update(R[, ]) <- R[i,j] + n_I_AR[i,j] + n_I_CR[i,j] + n_HR[i,j] - n_RV[i,j] + 
    (if (j == 1) n_RV[i,n_vax] else n_RV[i,j-1])
update(H[, ]) <- H[i,j] + n_I_CH[i,j] - n_HR[i,j] - n_HD[i,j]
update(G[, ]) <- G[i,j] + n_I_CG[i,j] - n_GD[i,j]
update(D[, ]) <- D[i,j] + n_HD[i,j] + n_GD[i,j]

## Individual probabilities of transition:
# infection progression:
p_SE[, ] <- 1 - exp(-lambda[i,j] * dt) # S to E
p_EI <- 1 - exp(-gamma_E * dt) # E to I
p_I_PI_C <- 1 - exp(-gamma_P * dt) # I_P to I_C
p_I_Cx <- 1 - exp(-gamma_C * dt) # I_C to x
p_I_AR <- 1 - exp(-gamma_A * dt) # I_A to R
p_Hx <- 1 - exp(-gamma_H * dt) # H to x
p_GD <- 1 - exp(-gamma_G * dt) # G to D

# vaccination:
p_SV[, ] <- p_V[i,j]
p_EV[, ] <- p_V[i,j]
p_I_AV[, ] <- p_V[i,j]
p_I_PV[, ] <- p_V[i,j]
p_RV[, ] <- p_V[i,j]

p_C[] <- user()
p_H[] <- user()
p_G[] <- user()
p_D[] <- user()

## Compute the force of infection
# Multiply numbers of infected individuals in each vaccine stratum j by their relative infectiousness
I_rel_inf[, ] <- rel_infectivity[i, j] * (theta_A*I_A[i,j] + I_P[i,j] + I_C[i,j])
# Calculate contact rate for individual in age group i with infectious individuals in age group j
s_ij[, ] <- m[i, j] * sum(I_rel_inf[j, ])
# Multiply by transmission rate and relative susceptibility to get force of infection on individual in age group i in vaccine stratum j
lambda[, ] <- beta * rel_susceptibility[i,j] * sum(s_ij[i, ])

## Draws from binomial distributions for numbers changing between 
## compartments:

# Flow out of S:
# infection
n_SE[, ] <- rbinom(S[i,j], p_SE[i,j])
# vaccination
n_SV[, ] <- rbinom(S[i,j] - n_SE[i,j], p_SV[i,j])

# Flow out of E:
# progression
n_EI[, ] <- rbinom(E[i,j], p_EI)
n_EI_P[, ] <- rbinom(n_EI[i,j], p_C[i])
n_EI_A[, ] <- n_EI[i,j] - n_EI_P[i,j]
# vaccination
n_EV[, ] <- rbinom(E[i,j] - n_EI[i,j], p_EV[i,j])

# Flow out of I_P:
# progression
n_I_PI_C[, ] <- rbinom(I_P[i,j], p_I_PI_C)
# vaccination
n_I_PV[, ] <- rbinom(I_P[i,j] - n_I_PI_C[i,j], p_I_PV[i,j])

# Flow out of I_A:
# progression
n_I_AR[, ] <- rbinom(I_A[i,j], p_I_AR)
# vaccination
n_I_AV[, ] <- rbinom(I_A[i,j] - n_I_AR[i,j], p_I_AV[i,j])

# Flow out of I_C:
# progression
n_I_Cx[, ] <- rbinom(I_C[i,j], p_I_Cx)
n_I_CR[, ] <- rbinom(n_I_Cx[i,j], 1 - p_H[i])
n_I_CHG[, ] <- n_I_Cx[i,j] - n_I_CR[i,j]
n_I_CH[, ] <- rbinom(n_I_CHG[i,j], 1 - p_G[i])
n_I_CG[, ] <- n_I_CHG[i,j] - n_I_CH[i,j]

# Flow out of H:
n_Hx[, ] <- rbinom(H[i,j], p_Hx)
n_HD[, ] <- rbinom(n_Hx[i,j], p_D[i])
n_HR[, ] <- n_Hx[i,j] - n_HD[i,j]

# Flow out of G:
n_GD[, ] <- rbinom(G[i,j], p_GD)

# Flow out of R:
n_RV[, ] <- rbinom(R[i,j],p_RV[i,j])

# Number vaccinated
n_V[, ] <- n_SV[i,j] + n_EV[i,j] + n_I_AV[i,j] + n_I_PV[i,j] + n_RV[i,j]

## Initial conditions 
initial(S[, ]) <- 0
initial(E[, ]) <- 0
initial(I_P[, ]) <- 0
initial(I_A[, ]) <- 0
initial(I_C[, ]) <- 0
initial(R[, ]) <- 0
initial(H[, ]) <- 0
initial(G[, ]) <- 0
initial(D[, ]) <- 0

## Seeding
seed_step_end <- seed_step_start + length(seed_value)
seed_rate <- if (step >= seed_step_start && step < seed_step_end)
    seed_value[as.integer(step - seed_step_start + 1)] else 0
seed <- rpois(seed_rate)

n_SE[seed_age,1] <- n_SE[seed_age,1] + min(S[seed_age,1] - n_SE[seed_age,1], seed)

## Store hospitalisation and death incidence
initial(H_inc) <- 0
initial(D_inc) <- 0
update(H_inc) <- if (step %% steps_per_day == 0) sum(n_I_CH) else H_inc + sum(n_I_CH)
update(D_inc) <- if (step %% steps_per_day == 0) sum(n_HD) else D_inc + sum(n_HD)

## Store hospitalisation and death incidence by age
initial(H_inc_0_39) <- 0
update(H_inc_0_39) <- if (step %% steps_per_day == 0) sum(n_I_CH[1:4, ]) else H_inc_0_39 + sum(n_I_CH[1:4, ])
initial(H_inc_40_49) <- 0
update(H_inc_40_49) <- if (step %% steps_per_day == 0) sum(n_I_CH[5, ]) else H_inc_40_49 + sum(n_I_CH[5, ])
initial(H_inc_50_59) <- 0
update(H_inc_50_59) <- if (step %% steps_per_day == 0) sum(n_I_CH[6, ]) else H_inc_50_59 + sum(n_I_CH[6, ])
initial(H_inc_60_69) <- 0
update(H_inc_60_69) <- if (step %% steps_per_day == 0) sum(n_I_CH[7, ]) else H_inc_60_69 + sum(n_I_CH[7, ])
initial(H_inc_70_plus) <- 0
update(H_inc_70_plus) <- if (step %% steps_per_day == 0) sum(n_I_CH[8, ]) else H_inc_70_plus + sum(n_I_CH[8, ])

initial(D_inc_0_39) <- 0
update(D_inc_0_39) <- if (step %% steps_per_day == 0) sum(n_HD[1:4, ]) else D_inc_0_39 + sum(n_HD[1:4, ])
initial(D_inc_40_49) <- 0
update(D_inc_40_49) <- if (step %% steps_per_day == 0) sum(n_HD[5, ]) else D_inc_40_49 + sum(n_HD[5, ])
initial(D_inc_50_59) <- 0
update(D_inc_50_59) <- if (step %% steps_per_day == 0) sum(n_HD[6, ]) else D_inc_50_59 + sum(n_HD[6, ])
initial(D_inc_60_69) <- 0
update(D_inc_60_69) <- if (step %% steps_per_day == 0) sum(n_HD[7, ]) else D_inc_60_69 + sum(n_HD[7, ])
initial(D_inc_70_plus) <- 0
update(D_inc_70_plus) <- if (step %% steps_per_day == 0) sum(n_HD[8, ]) else D_inc_70_plus + sum(n_HD[8, ])

## Model parameters (default in parenthesis) 
beta <- user(0.03)
gamma_E <- user(0.5)
gamma_P <- user(0.4)
gamma_C <- user(0.4)
gamma_A <- user(0.2)
gamma_H <- user(0.1)
gamma_G <- user(1/3)
theta_A <- user(0.5) 
m[, ] <- user() # age-structured contact matrix
rel_infectivity[, ] <- user() # average vaccine effectiveness against infection by age group and vaccine stratum
rel_susceptibility[, ] <- user() # infectiousness by age group and vaccine stratum
seed_age <- user(4) # 30-49 year band
seed_step_start <- user() # step at which to start seeding
seed_value[] <- user() # number of infections seeded into seed_age over length(seed_value) steps

# Vaccination
n_doses <- user()
index_dose[] <- user(integer = TRUE)
index_dose_inverse[] <- user(integer = TRUE)
vaccine_dose_step[, , ] <- user() # n_age x n_doses x n_time

# Number of vaccination candidates
vaccine_n_candidates[, ] <-
    S[i,index_dose[j]] +
    sum(E[i,index_dose[j]]) +
    sum(I_A[i,index_dose[j]]) +
    sum(I_P[i,index_dose[j]]) +
    sum(R[i,index_dose[j]])

# Work out the vaccination probability via doses, driven by the
# schedule
vaccine_probability_doses[, ] <- min(
    if (vaccine_n_candidates[i, j] > 0)
        vaccine_attempted_doses[i, j] / vaccine_n_candidates[i, j] else 0,
    as.numeric(1))

# Work out the total attempted doses
total_attempted_doses[, ] <- vaccine_missed_doses[i, j] + (
    if (as.integer(step) >= dim(vaccine_dose_step, 3)) 0
    else vaccine_dose_step[i, j, step + 1])


# No vaccine skip moves currently so attempted doses are just equal to total 
# attempted doses
vaccine_attempted_doses[, ] <- total_attempted_doses[i, j]

# Calculate catch-up vaccine doses, i.e. doses not delivered according to 
# schedule, e.g. because too many individuals are in the I or H states, that 
# are postponed till later
initial(vaccine_missed_doses[, ]) <- 0
update(vaccine_missed_doses[, ]) <-
    vaccine_catchup_fraction *
    max(total_attempted_doses[i, j] - n_V[i, index_dose[j]],
        as.numeric(0))
vaccine_catchup_fraction <- user(1)

# Calculate probability of vaccination with a certain dose based on progression
# at a constant rate (vaccine_progression_rate_base) or the supplied 
# time-varying probabilities (vaccine_probability_doses)
p_V[, ] <- if (index_dose_inverse[j] > 0)
        vaccine_probability_doses[i, index_dose_inverse[j]] else
        1 - exp(-vaccine_progression_rate_base[i, j] * dt)
vaccine_progression_rate_base[, ] <- user()

## Array dimensions
dim(S) <- c(n_age,n_vax)
dim(E) <- c(n_age,n_vax)
dim(I_P) <- c(n_age,n_vax)
dim(I_A) <- c(n_age,n_vax)
dim(I_C) <- c(n_age,n_vax)
dim(R) <- c(n_age,n_vax)
dim(H) <- c(n_age,n_vax)
dim(G) <- c(n_age,n_vax)
dim(D) <- c(n_age,n_vax)
dim(n_SE) <- c(n_age,n_vax)
dim(n_EI) <- c(n_age,n_vax)
dim(n_EI_P) <- c(n_age,n_vax)
dim(n_EI_A) <- c(n_age,n_vax)
dim(n_I_PI_C) <- c(n_age,n_vax)
dim(n_I_AR) <- c(n_age,n_vax)
dim(n_I_Cx) <- c(n_age,n_vax)
dim(n_I_CR) <- c(n_age,n_vax)
dim(n_I_CHG) <- c(n_age,n_vax)
dim(n_I_CH) <- c(n_age,n_vax)
dim(n_I_CG) <- c(n_age,n_vax)
dim(n_Hx) <- c(n_age,n_vax)
dim(n_HD) <- c(n_age,n_vax)
dim(n_HR) <- c(n_age,n_vax)
dim(n_GD) <- c(n_age,n_vax)
dim(n_SV) <- c(n_age,n_vax)
dim(n_EV) <- c(n_age,n_vax)
dim(n_I_PV) <- c(n_age,n_vax)
dim(n_I_AV) <- c(n_age,n_vax)
dim(n_RV) <- c(n_age,n_vax)
dim(n_V) <- c(n_age,n_vax)
dim(p_SE) <- c(n_age,n_vax)
dim(lambda) <- c(n_age,n_vax)
dim(m) <- c(n_age, n_age)
dim(s_ij) <- c(n_age, n_age)
dim(rel_infectivity) <- c(n_age,n_vax)
dim(I_rel_inf) <- c(n_age,n_vax)
dim(rel_susceptibility) <- c(n_age,n_vax)
dim(p_C) <- n_age
dim(p_H) <- n_age
dim(p_G) <- n_age
dim(p_D) <- n_age
dim(p_V) <- c(n_age,n_vax)
dim(vaccine_progression_rate_base) <- c(n_age, n_vax)
dim(p_SV) <- c(n_age,n_vax)
dim(p_EV) <- c(n_age,n_vax)
dim(p_I_PV) <- c(n_age,n_vax)
dim(p_I_AV) <- c(n_age,n_vax)
dim(p_RV) <- c(n_age,n_vax)
dim(seed_value) <- user()
dim(index_dose) <- n_doses
dim(index_dose_inverse) <- n_vax
dim(vaccine_dose_step) <- user()
dim(vaccine_n_candidates) <- c(n_age, n_doses)
dim(vaccine_attempted_doses) <- c(n_age, n_doses)
dim(vaccine_probability_doses) <- c(n_age, n_doses)
dim(total_attempted_doses) <- c(n_age, n_doses)
dim(vaccine_missed_doses) <- c(n_age, n_doses)