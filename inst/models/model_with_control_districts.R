## definition of time step
steps_per_day <- user(1)
steps_per_week <- steps_per_day * 7
steps_per_month <- steps_per_day * 30
dt <- 1 / steps_per_day
time <- step + 1
n_districts <- user()

# SMC Effect
SMC_effect[] <- decay[i, time] * eff_SMC * cov_SMC[i, time]
SMC_removal[] <- if (SMC[i, time] == 1) SMC_effect[i] else 0
dim(SMC_effect) <- n_districts
dim(SMC_removal) <- n_districts

# Children
update(SC[]) <- (SC[i] * r_C) * (1 - delta_a - delta_d - (1 - SMC_effect[i]) * mu_SE_C[i]) +
  delta_b * P[i] + mu_RS * RC[i] + mu_TS * TrC[i]
update(EC[]) <- (EC[i] * r_C) * (1 - delta_a - delta_d - mu_EI - SMC_removal[i]) +
  (1 - SMC_effect[i]) * mu_SE_C[i] * SC[i]
update(IC[]) <- (IC[i] * r_C) * (1 - delta_a - delta_d - mu_IR - SMC_removal[i]) +
  (1 - fT_C[i]) * mu_EI * EC[i]
update(TrC[]) <- (TrC[i] * r_C) * (1 - delta_a - delta_d - mu_TS) +
  fT_C[i] * mu_EI * EC[i] + SMC_removal[i] * (EC[i] + IC[i])
update(RC[]) <- (RC[i] * r_C) * (1 - delta_a - delta_d - mu_RS) + mu_IR * IC[i]

# Adults
update(SA[]) <- (SA[i] * r_A) * (1 - delta_d - mu_SE_A[i]) + delta_a * SC[i] + mu_RS * RA[i] + mu_TS * TrA[i]
update(EA[]) <- (EA[i] * r_A) * (1 - delta_d - mu_EI) + delta_a * EC[i] + mu_SE_A[i] * SA[i]
update(IA[]) <- (IA[i] * r_A) * (1 - delta_d - mu_IR) + delta_a * IC[i] + (1 - fT_A[i]) * mu_EI * EA[i]
update(TrA[]) <- (TrA[i] * r_A) * (1 - delta_d - mu_TS) + delta_a * TrC[i] + fT_A[i] * mu_EI * EA[i]
update(RA[]) <- (RA[i] * r_A) * (1 - delta_d - mu_RS) + delta_a * RC[i] + mu_IR * IA[i]


# Monthly incidence
initial(month_inc_C[]) <- 0
update(month_inc_C[]) <- if ((step) %% steps_per_month == 0) mu_EI * EC[i] * fT_C[i] else month_inc_C[i] + mu_EI * EC[i] * fT_C[i]
dim(month_inc_C) <- n_districts

initial(month_inc_A[]) <- 0
update(month_inc_A[]) <- if ((step) %% steps_per_month == 0) mu_EI * EA[i] * fT_A[i] else month_inc_A[i] + mu_EI * EA[i] * fT_A[i]
dim(month_inc_A) <- n_districts

# EIR
EIR[] <- alpha[i] * (X[i] / (b + X[i])) * temp_effect[i] * rain_effect[i]
temp_effect[] <- if (temp_current[i] <= T_opt) exp(-((temp_current[i] - T_opt)^2) / (2 * sigma_LT^2)) else exp(-((temp_current[i] - T_opt)^2) / (2 * sigma_RT^2))
rain_effect[] <- 1 / (1 + exp(-k1 * (c_R_D_current[i] - R_opt)))
dim(EIR) <- n_districts
dim(temp_effect) <- n_districts
dim(rain_effect) <- n_districts
X[] <- (IA[i] + IC[i] + qR * (RA[i] + RC[i])) / P[i] # Proportion of population infectious remains the same
dim(X) <- n_districts

# Temperature-entomological parameters
T_opt <- user()
sigma_LT <- user()
sigma_RT <- user()
k1 <- user()
R_opt <- user()
b <- user()
alpha[] <- user()
dim(alpha) <- n_districts

# Climate inputs
temp[,] <- user()     # Temperature data for each district and time
c_R_D[,] <- user()    # Rainfall data for each district and time
dim(temp) <- user()
dim(c_R_D) <- user()

# Current climate values
temp_current[] <- temp[i, time]
c_R_D_current[] <- c_R_D[i, time]
dim(temp_current) <- n_districts
dim(c_R_D_current) <- n_districts

# Growth rates
r_C <- user(1.000071) # daily growth rate u5 Chad
r_A <- user(1.000092) # daily growth rate o5 Chad

# Other parameters
fT_C[] <- user()
z[] <- user()
fT_A[] <- z[i] * fT_C[i]
dim(fT_C) <- n_districts
dim(z) <- n_districts
dim(fT_A) <- n_districts
qR <- user()
p_MH <- user()

# Transmission rates
mu_SE_C[] <- 1 - exp(-p_MH * EIR[i]) # Transmission rate for children
phi <- user()
mu_SE_A[] <- phi * mu_SE_C[i]   # Transmission rate for adults
mu_RS <- user()
mu_IR <- user()
mu_TS <- user()
mu_EI <- user()

dim(mu_SE_C) <- n_districts
dim(mu_SE_A) <- n_districts

# Rate parameters
delta_a <- user(1 / (5 * 365))
delta_b <- user(1 / (21 * 365))
delta_d <- user(1 / (21 * 365))

# Intervention parameters
eff_SMC <- user()      # SMC effectiveness per district
cov_SMC[,] <- user()     # SMC coverage over time for each district
SMC[,] <- user()         # SMC active status for each district and time
decay[,] <- user()       # Decay of SMC efficacy over time
dim(cov_SMC) <- user()
dim(SMC) <- user()
dim(decay) <- user()

# State variables for children
dim(SC) <- n_districts
dim(EC) <- n_districts
dim(IC) <- n_districts
dim(TrC) <- n_districts
dim(RC) <- n_districts

# State variables for adults
dim(SA) <- n_districts
dim(EA) <- n_districts
dim(IA) <- n_districts
dim(TrA) <- n_districts
dim(RA) <- n_districts

# Initial state variables
SC0 <- user(0.301)
EC0 <- user(0.071)
IC0 <- user(0.042)
TC0 <- user(0.054)
RC0 <- user(0.533)
SA0 <- user(0.281)
EA0 <- user(0.067)
IA0 <- user(0.040)
TA0 <- user(0.051)
RA0 <- user(0.562)

# Initial states for children
initial(SC[]) <- SC0 * N_C[i]
initial(EC[]) <- EC0 * N_C[i]
initial(IC[]) <- IC0 * N_C[i]
initial(TrC[]) <- TC0 * N_C[i]
initial(RC[]) <- RC0 * N_C[i]

# Initial states for adults
initial(SA[]) <- SA0 * N_A[i]
initial(EA[]) <- EA0 * N_A[i]
initial(IA[]) <- IA0 * N_A[i]
initial(TrA[]) <- TA0 * N_A[i]
initial(RA[]) <- RA0 * N_A[i]

# Initial population sizes
N[] <- user()
s[] <- user()
N_pop[] <- N[i] * s[i]
percAdult[] <- user()
percChild[] <- 1 - percAdult[i]

dim(N) <- n_districts
dim(N_pop) <- n_districts
dim(percAdult) <- n_districts
dim(percChild) <- n_districts
dim(s) <- n_districts

N_C[] <- N_pop[i] * percChild[i]
N_A[] <- N_pop[i] * percAdult[i]
dim(N_C) <- n_districts
dim(N_A) <- n_districts

# Population sizes
P_C[] <- SC[i] + EC[i] + IC[i] + TrC[i] + RC[i]
P_A[] <- SA[i] + EA[i] + IA[i] + TrA[i] + RA[i]
P[] <-  P_C[i] + P_A[i]
dim(P_C) <- n_districts
dim(P_A) <- n_districts
dim(P) <- n_districts
