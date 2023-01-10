## Definition of the time-step and output as "time" 
dt <- user()
initial(time) <- 0
update(time) <- (step + 1) * dt
steps_per_day <- 1/dt

## Core equations for transitions between compartments:
# update(S_tot) <- S_tot - sum(n_SE)
# update(E_tot) <- E_tot + sum(n_SE) - sum(n_EI)
# update(I_tot) <- I_tot + sum(n_EI) - sum(n_IR)
# update(R_tot) <- R_tot + sum(n_IR)
# update(D_tot) <- D_tot + sum(n_ID)
# 
# # Equations for transitions between compartments by age group
update(S[]) <- S[i] - n_SE[i]
update(E[]) <- E[i] + n_SE[i] - n_EI[i]
update(I_P[]) <- I_P[i] + n_EI_P[i] - n_I_PI_C[i]
update(I_A[]) <- I_A[i] + n_EI_A[i] - n_I_AR[i]
update(I_C[]) <- I_C[i] + n_I_PI_C[i] - n_I_CR[i] - n_I_CH[i] - n_I_CG[i]
update(R[]) <- R[i] + n_I_AR[i] + n_I_CR[i] + n_HR[i]
update(H[]) <- H[i] + n_I_CH[i] - n_HR[i] - n_HD[i]
update(G[]) <- G[i] + n_I_CG[i] - n_GD[i]
update(D[]) <- D[i] + n_HD[i] + n_GD[i]

update(N_tot[]) <- S[i] + E[i] + I_P[i] + I_A[i] + I_C[i] + R[i] + H[i] + G[i] + D[i]

# ## To compute cumulative incidence
# update(cum_inc[]) <- cum_inc[i] + n_SE[i] # THIS MAY ONLY BE OK IF step IS NOT A MULTIPLE OF 1/dt -- CHECK!
# update(cum_inc_death[]) <- cum_inc_death[i] + n_ID[i]

## Individual probabilities of transition:
p_SE[] <- 1 - exp(-lambda[i] * dt) # S to E
p_EI <- 1 - exp(-sigma * dt) # E to I
p_I_PI_C <- 1 - exp(-gamma_P * dt) # I_P to I_C
p_I_Cx <- 1 - exp(-gamma_C * dt) # I_C to x
p_I_AR <- 1 - exp(-gamma_A * dt) # I_A to R
p_Hx <- 1 - exp(-gamma_H * dt) # H to x
p_GD <- 1 - exp(-gamma_G * dt) # G to D

p_C[] <- user()
p_H[] <- user()
p_G[] <- user()
p_D[] <- user()
# p_death[] <- user()

## Force of infection
theta_A <- user(0.5) 
m[, ] <- user() # age-structured contact matrix
s_ij[, ] <- m[i, j] * (theta_A*I_A[i] + I_P[i] + I_C[i])
lambda[] <- beta * sum(s_ij[, i])

## Draws from binomial distributions for numbers changing between 
## compartments:
n_SE[] <- rbinom(S[i], p_SE[i])
n_EI[] <- rbinom(E[i], p_EI)
n_EI_P[] <- rbinom(n_EI[i], p_C[i])
n_EI_A[] <- n_EI[i] - n_EI_P[i]
n_I_PI_C[] <- rbinom(I_P[i], p_I_PI_C)
n_I_Cx[] <- rbinom(I_C[i], p_I_Cx)
n_I_CR[] <- rbinom(n_I_Cx[i], 1 - p_H[i])
n_I_CHG[] <- n_I_Cx[i] - n_I_CR[i]
n_I_CH[] <- rbinom(n_I_CHG[i], 1 - p_G[i])
n_I_CG[] <- n_I_CHG[i] - n_I_CH[i]
n_I_AR[] <- rbinom(I_A[i], p_I_AR)
n_Hx[] <- rbinom(H[i], p_Hx)
n_HD[] <- rbinom(n_Hx[i], p_D[i])
n_HR[] <- n_Hx[i] - n_HD[i]
n_GD[] <- rbinom(G[i], p_GD)

# n_IR[] <- rbinom(n_Ix[i], 1 - p_death[i])
# n_ID[] <- n_Ix[i] - n_IR[i]

## Initial conditions 
# initial(S_tot) <- sum(S_ini) 
# initial(E_tot) <- sum(E_ini)
# initial(I_tot) <- sum(I_ini) 
# initial(R_tot) <- 0
# initial(D_tot) <- 0

initial(S[]) <- S_ini[i]
initial(E[]) <- E_ini[i]
initial(I_P[]) <- 0
initial(I_A[]) <- 0
initial(I_C[]) <- I_ini[i]
initial(R[]) <- 0
initial(H[]) <- 0
initial(G[]) <- 0
initial(D[]) <- 0
# initial(cum_inc[]) <- 0
# initial(cum_inc_death[]) <- 0

initial(N_tot[]) <- S_ini[i] + E_ini[i] + I_ini[i]

## Store hospitalisation and death incidence
# initial(I_inc) <- 0
initial(H_inc) <- 0
initial(D_inc) <- 0
# update(I_inc) <- if (step %% steps_per_day == 0) sum(n_EI) else I_inc + sum(n_EI)
update(H_inc) <- if (step %% steps_per_day == 0) sum(n_I_CH) else H_inc + sum(n_I_CH)
update(D_inc) <- if (step %% steps_per_day == 0) sum(n_HD) else D_inc + sum(n_HD)

## Store hospitalisation and death incidence by age
# initial(I_inc_0_9) <- 0
# update(I_inc_0_9) <- if (step %% steps_per_day == 0) n_EI[1] else I_inc_0_9 + n_EI[1]
# initial(I_inc_10_19) <- 0
# update(I_inc_10_19) <- if (step %% steps_per_day == 0) n_EI[2] else I_inc_10_19 + n_EI[2]
# initial(I_inc_20_29) <- 0
# update(I_inc_20_29) <- if (step %% steps_per_day == 0) n_EI[3] else I_inc_20_29 + n_EI[3]
# initial(I_inc_30_39) <- 0
# update(I_inc_30_39) <- if (step %% steps_per_day == 0) n_EI[4] else I_inc_30_39 + n_EI[4]
# initial(I_inc_40_49) <- 0
# update(I_inc_40_49) <- if (step %% steps_per_day == 0) n_EI[5] else I_inc_40_49 + n_EI[5]
# initial(I_inc_50_59) <- 0
# update(I_inc_50_59) <- if (step %% steps_per_day == 0) n_EI[6] else I_inc_50_59 + n_EI[6]
# initial(I_inc_60_69) <- 0
# update(I_inc_60_69) <- if (step %% steps_per_day == 0) n_EI[7] else I_inc_60_69 + n_EI[7]
# initial(I_inc_70_plus) <- 0
# update(I_inc_70_plus) <- if (step %% steps_per_day == 0) n_EI[8] else I_inc_70_plus + n_EI[8]

initial(H_inc_0_39) <- 0
update(H_inc_0_39) <- if (step %% steps_per_day == 0) sum(n_I_CH[1:4]) else H_inc_0_39 + sum(n_I_CH[1:4])
initial(H_inc_40_49) <- 0
update(H_inc_40_49) <- if (step %% steps_per_day == 0) n_I_CH[5] else H_inc_40_49 + n_I_CH[5]
initial(H_inc_50_59) <- 0
update(H_inc_50_59) <- if (step %% steps_per_day == 0) n_I_CH[6] else H_inc_50_59 + n_I_CH[6]
initial(H_inc_60_69) <- 0
update(H_inc_60_69) <- if (step %% steps_per_day == 0) n_I_CH[7] else H_inc_60_69 + n_I_CH[7]
initial(H_inc_70_plus) <- 0
update(H_inc_70_plus) <- if (step %% steps_per_day == 0) n_I_CH[8] else H_inc_70_plus + n_I_CH[8]

initial(D_inc_0_39) <- 0
update(D_inc_0_39) <- if (step %% steps_per_day == 0) sum(n_HD[1:4]) else D_inc_0_39 + sum(n_HD[1:4])
initial(D_inc_40_49) <- 0
update(D_inc_40_49) <- if (step %% steps_per_day == 0) n_HD[5] else D_inc_40_49 + n_HD[5]
initial(D_inc_50_59) <- 0
update(D_inc_50_59) <- if (step %% steps_per_day == 0) n_HD[6] else D_inc_50_59 + n_HD[6]
initial(D_inc_60_69) <- 0
update(D_inc_60_69) <- if (step %% steps_per_day == 0) n_HD[7] else D_inc_60_69 + n_HD[7]
initial(D_inc_70_plus) <- 0
update(D_inc_70_plus) <- if (step %% steps_per_day == 0) n_HD[8] else D_inc_70_plus + n_HD[8]

## Model parameters (default in parenthesis) 
S_ini[] <- user()
E_ini[] <- user()
I_ini[] <- user()
beta <- user(0.03)
sigma <- user(0.5)
gamma_P <- user(0.4)
gamma_C <- user(0.4)
gamma_A <- user(0.2)
gamma_H <- user(0.1)
gamma_G <- user(1/3)

## Array dimensions
n_age <- user()
dim(S_ini) <- n_age
dim(E_ini) <- n_age
dim(I_ini) <- n_age
dim(S) <- n_age
dim(E) <- n_age
dim(I_P) <- n_age
dim(I_A) <- n_age
dim(I_C) <- n_age
dim(R) <- n_age
dim(H) <- n_age
dim(G) <- n_age
dim(D) <- n_age
# dim(cum_inc) <- n_age
# dim(cum_inc_death) <- n_age
dim(N_tot) <- n_age
dim(n_SE) <- n_age
dim(n_EI) <- n_age
dim(n_EI_P) <- n_age
dim(n_EI_A) <- n_age
dim(n_I_PI_C) <- n_age
dim(n_I_AR) <- n_age
dim(n_I_Cx) <- n_age
dim(n_I_CR) <- n_age
dim(n_I_CHG) <- n_age
dim(n_I_CH) <- n_age
dim(n_I_CG) <- n_age
dim(n_Hx) <- n_age
dim(n_HD) <- n_age
dim(n_HR) <- n_age
dim(n_GD) <- n_age
# dim(n_IR) <- n_age
# dim(n_ID) <- n_age
dim(p_SE) <- n_age
dim(lambda) <- n_age
dim(m) <- c(n_age, n_age)
dim(s_ij) <- c(n_age, n_age)
# dim(p_death) <- n_age
dim(p_C) <- n_age
dim(p_H) <- n_age
dim(p_G) <- n_age
dim(p_D) <- n_age
