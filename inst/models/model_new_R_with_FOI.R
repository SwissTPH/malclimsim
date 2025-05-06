## The model contained in the file 'age_model_equation_final' will be written
## using odin and inference in mcstate. Here, this is done with the
## simplification of not including the age structure or effect of SMC,
## which will be added later.

## rainfall in mL; temp in Celsius
## it is timeed an identical array at every time step

## definition of time step
steps_per_day <- user(1)
steps_per_week <- steps_per_day * 7
steps_per_month <- steps_per_day * 30
dt <- 1 / steps_per_day
time <- step + 1

#SMC_effect <- decay[time] * eff_SMC * cov_SMC[time]
lag_SMC <- user()
SMC_effect <- if (time > lag_SMC) decay[time - lag_SMC] * eff_SMC * cov_SMC[time - lag_SMC] else decay[time] * eff_SMC * cov_SMC[time]

initial(SMC_effect_2) <- decay[time] * eff_SMC * cov_SMC[time]
update(SMC_effect_2) <- decay[time] * eff_SMC * cov_SMC[time]
w1 <- user(0) # Gives some SMC to adults
w2 <- user(0) # Controls removal effect
SMC_effect_A <- w1 * SMC_effect
SMC_removal <- if (SMC[time] == 1) w2 * SMC_effect else 0

#SMC_removal <- user(0)

mu_SE_C <- (1 - SMC_effect) *(1 - exp(-p_MH_C * EIR))
initial(mu_SE_C_2) <-(1 - SMC_effect) * (1 - exp(-p_MH_C * EIR))
update(mu_SE_C_2) <- (1 - SMC_effect) *(1 - exp(-p_MH_C * EIR))

mu_SE_A <- (1 - SMC_effect_A) *(1 - exp(-rho * p_MH_C * EIR))
initial(mu_SE_A_2) <- (1 - SMC_effect_A) *(1 - exp(-rho * p_MH_C * EIR))
update(mu_SE_A_2) <- (1 - SMC_effect_A) *(1 - exp(-rho * p_MH_C * EIR))


# # Children

update(SC) <- SC -  mu_SE_C * SC + mu_RS_C * RC + mu_TS * TrC - (delta_d + delta_a) * SC + delta_b * P - SC * SMC_removal
update(EC) <- EC - mu_EI * EC +  mu_SE_C * SC - (delta_d + delta_a) * EC
update(IC) <- IC - mu_IR * IC + phi_1 * (1 - fT_C) * mu_EI * EC - (delta_d + delta_a) * IC - IC * SMC_removal
update(TrC) <- TrC - mu_TS * TrC + phi_1 * fT_C * mu_EI * EC - (delta_d + delta_a) * TrC + (IC + RC +SC) * SMC_removal
update(RC) <- RC - mu_RS_C * RC + mu_IR * IC + (1 - phi_1) * mu_EI * EC - (delta_d + delta_a) * RC  - RC * SMC_removal

# Adults
update(SA) <- SA -  mu_SE_A * SA + mu_RS_A * RA + mu_TS * TrA - delta_d * SA + delta_a * SC
update(EA) <- EA - mu_EI * EA +  mu_SE_A * SA - delta_d * EA + delta_a * EC
update(IA) <- IA - mu_IR * IA + phi_2 * (1 - fT_A) * mu_EI * EA - delta_d * IA + delta_a * IC
update(TrA) <- TrA - mu_TS * TrA + phi_2 * fT_A * mu_EI * EA - delta_d * TrA + delta_a * TrC
update(RA) <- RA - mu_RS_A * RA + mu_IR * IA + (1 - phi_2) * mu_EI * EA - delta_d * RA + delta_a * RC

# daily and weekly incidence
initial(day_inc_C) <- 0
update(day_inc_C) <- if ((step) %% steps_per_day == 0) mu_EI * EC * fT_C * phi_1 else day_inc_C + mu_EI * EC * fT_C * phi_1

initial(wk_inc_C) <- 0
update(wk_inc_C) <- if ((step) %% steps_per_week == 0) mu_EI * EC * fT_C * phi_1 else wk_inc_C + mu_EI * EC * fT_C * phi_1

initial(day_inc_A) <- 0
update(day_inc_A) <- if ((step) %% steps_per_day == 0) mu_EI * EA * fT_A * phi_2 else day_inc_A + mu_EI * EA * fT_A * phi_2

initial(wk_inc_A) <- 0
update(wk_inc_A) <- if ((step) %% steps_per_week == 0) mu_EI * EA * fT_A * phi_2 else wk_inc_A + mu_EI * EA * fT_A * phi_2

initial(day_inc_total) <- 0
update(day_inc_total) <- if ((step) %% steps_per_day == 0) mu_EI * (EC * fT_C * phi_1 + EA * fT_A * phi_2) else day_inc_total + mu_EI * (EC * fT_C * phi_1 + EA * fT_A * phi_2)

initial(wk_inc_total) <- 0
update(wk_inc_total) <-  if ((step) %% steps_per_week == 0) mu_EI * (EC * fT_C * phi_1 + EA * fT_A * phi_2) else wk_inc_total + mu_EI * (EC * fT_C * phi_1 + EA * fT_A * phi_2)

#beta_1 <- user()
#beta_2 <- user()
SMC_effective_coverage <- if (time > lag_SMC) decay[time - lag_SMC] * cov_SMC[time - lag_SMC] else decay[time] * cov_SMC[time]
#SMC_X <- beta_1 * SMC_effective_coverage + beta_2 * SMC_effective_coverage * c_R_D_shift
#SMC_X <- beta_1 * SMC_effective_coverage + beta_2 * c_R_D_shift * SMC_effective_coverage

# initial(month_inc_C) <- 0
# update(month_inc_C) <- if ((step) %% steps_per_month == 0) mu_EI * EC * fT_C * phi_1 * (1 - beta_1 * SMC_effective_coverage)  else month_inc_C + mu_EI * EC * fT_C * phi_1 * (1 - beta_1 * SMC_effective_coverage)
#
# initial(month_inc_A) <- 0
# update(month_inc_A) <- if ((step) %% steps_per_month == 0) mu_EI * EA * fT_A * phi_2 else month_inc_A + mu_EI * EA * fT_A * phi_2

initial(month_inc_C) <- 0
update(month_inc_C) <- if ((step) %% steps_per_month == 0) mu_EI * EC * fT_C * phi_1 else month_inc_C + mu_EI * EC * fT_C * phi_1

initial(month_inc_A) <- 0
update(month_inc_A) <- if ((step) %% steps_per_month == 0) mu_EI * EA * fT_A * phi_2 else month_inc_A + mu_EI * EA * fT_A * phi_2


initial(month_inc_total) <- 0
update(month_inc_total) <- if ((step) %% steps_per_month == 0) mu_EI * (EC * fT_C * phi_1 + EA * fT_A * phi_2) else month_inc_total + mu_EI * (EC * fT_C * phi_1 + EA * fT_A * phi_2)

# Calculating Prevalence
initial(prev_total_1) <- (IA + RA * s_2  + IC + RC * s_1) / N
update(prev_total_1) <- (IA + RA * s_2  + IC + RC * s_1) / N

initial(prev_C_1) <- (IC + RC * s_1) / N_C
update(prev_C_1) <- (IC + RC * s_1) / N_C

initial(prev_A_1) <- (IA + RA * s_2) / N_A
update(prev_A_1) <- (IA + RA * s_2) / N_A

initial(prev_total_2) <- (IA  + IC) / N
update(prev_total_2) <- (IA  + IC) / N

initial(prev_C_2) <- IC / N_C
update(prev_C_2) <- IC  / N_C

initial(prev_A_2) <- IA / N_A
update(prev_A_2) <- IA / N_A

size <- user(10)
size_1 <- user(10)
size_2 <- user(10)
kappa_C <- user(30)
kappa_A <- user(30)

# To be fitted for better reproducing prevalence
mu_RS_C <- user()
duration_infection <- user(200)
mu_IR <- 1 / (duration_infection -  1 / mu_RS_C) # I to R is a maximum of 200 days

mu_RS_A <- eta * mu_RS_C
eta <- user(1)

# Immunity related parameters
#s_1 <- user(1)
#s_2 <- c_s * s_1
#c_s <- prop_p_2
qR2 <- user()

#qR2 <- c_qR * qR1

qR1 <- qR2 + (1 - qR2) * c_qR

c_phi <- user() # alters
c_qR <- user()
c_s <- user()
prop_p_1 <- user() # prop of all malaria infections that are patent in <5
#prop_p_2 <- user() # prop of all malaria infections that are patent in >=5
c_p_1 <- user() # prop of patent infections that are symptomatic in <5
#c_p_2 <- user() # prop of patent infections that are symptomatic in >=5


# phi_1 <- user(0.5)
# phi_2 <- c_phi * phi_1

prop_p_2 <- c_s * prop_p_1
c_p_2 <- c_phi * c_p_1

# Proportion of population that is symptomatic
phi_1 <- c_p_1 * prop_p_1
phi_2 <- c_p_2 * prop_p_2


# Proprtion of asymptomatic population that is patent
s_1 <- 1 - ((1 - prop_p_1) / (1 - phi_1))
s_2 <- 1 - ((1 - prop_p_2) / (1 - phi_2))


# Others
mu_EI <- user()
mu_TS <- user()
#mu_IR <- user()
fT_C <- user()
fT_A <- z * fT_C
z <- user(1) # z < 1, reporting rate of adults vs. children
rho <- user(1)
p_MH_C <- user(0.50)
delta_b <- user()
delta_d <- user()
delta_a <- user()
alpha <- user(1)
T_opt <- user(28)
sigma_LT <- user(4)
sigma_RT <- user(4)
R_opt <- user(1)
k1 <- user(0.2)
b <- user()
lag_R <- user(0)  # Default lag = 0
lag_T <- user(0)

# defining SMC efficacy
eff_SMC <- user() # SMC effectiveness
cov_SMC[] <- user() # SMC coverage
SMC[] <- user()
decay[] <- user()
temp[] <- user()
c_R_D[] <- user() # this will be informed by data

initial(c_R_D_shift) <- c_R_D[time]  # Start with the first rainfall value
update(c_R_D_shift) <- if (time > lag_R) c_R_D[time - lag_R] else c_R_D[time]
initial(temp_shift) <- temp[time]
update(temp_shift) <- if (time > lag_T) temp[time - lag_T] else temp[time]

# Defining EIR
EIR <- alpha * (X / (b + X)) * temp_effect * rain_effect # Multiplicative effects
initial(EIR2) <- alpha * (X / (b + X)) * temp_effect * rain_effect # Multiplicative effects
update(EIR2) <- alpha * (X / (b + X)) * temp_effect * rain_effect # Multiplicative effects
#temp_effect <- exp(-((temp_shift - T_opt)^2) / (2 * sigma_LT^2)) # Gaussian term for temperature
temp_effect <- if (temp_shift <= T_opt) exp(-((temp_shift - T_opt)^2) / (2 * sigma_LT^2)) else exp(-((temp_shift - T_opt)^2) / (2 * sigma_RT^2))
rain_effect <- 1 / (1 + exp(-k1 * (c_R_D_shift - R_opt))) # Logistic term for rainfall

initial(temp_effect_2) <- if (temp_shift <= T_opt) exp(-((temp_shift - T_opt)^2) / (2 * sigma_LT^2)) else exp(-((temp_shift - T_opt)^2) / (2 * sigma_RT^2))
update(temp_effect_2) <- if (temp_shift <= T_opt) exp(-((temp_shift - T_opt)^2) / (2 * sigma_LT^2)) else exp(-((temp_shift - T_opt)^2) / (2 * sigma_RT^2))

initial(rain_effect_2) <- 1 / (1 + exp(-k1 * (c_R_D_shift - R_opt))) # Logistic term for rainfall
update(rain_effect_2) <-  1 / (1 + exp(-k1 * (c_R_D_shift - R_opt))) # Logistic term for rainfall

X <- (IC + IA + qR1 * (s_1 * RC + s_2 * RA) + qR2 * ((1 - s_1) * RC + (1 - s_2) * RA)) / P # Proportion of population infectious remains the same
# Egg-adult Development time due to temperature

# Dimensions of arrays
dim(c_R_D) <- user()
dim(temp) <- user()
dim(SMC) <- user()
dim(decay) <- user()
dim(cov_SMC) <- user()

## Initial conditions - user defined, defaults in parenthesis
N <- user()
s <- user()

N_pop <- N * s
percAdult <- user()
percChild <- 1 - percAdult
N_C <- N_pop * percChild
N_A <- N_pop * percAdult

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

initial(SC) <- SC0 * N_C
initial(EC) <- EC0 * N_C
initial(IC) <- IC0 * N_C
initial(TrC) <- TC0 * N_C
initial(RC) <- RC0 * N_C
initial(SA) <- SA0 * N_A
initial(EA) <- EA0 * N_A
initial(IA) <- IA0 * N_A
initial(TrA) <- TA0 * N_A
initial(RA) <- RA0 * N_A

## Total population size
initial(P_C) <- SC + EC + IC + TrC + RC
initial(P_A) <- SA + EA + IA + TrA + RA
update(P_C) <- SC + EC + IC + TrC + RC
update(P_A) <- SA + EA + IA + TrA + RA
P = (SC + EC + IC + TrC + RC + SA + EA + IA + TrA + RA)

## Relevant outputs to look at
initial(X2) <- (IC + IA + qR1 * (s_1 * RC + s_2 * RA) + qR2 * ((1 - s_1) * RC + (1 - s_2) * RA)) / P  # Proportion of population infectious remains the same
update(X2) <- (IC + IA + qR1 * (s_1 * RC + s_2 * RA) + qR2 * ((1 - s_1) * RC + (1 - s_2) * RA)) / P  # Proportion of population infectious remains the same
initial(X_I) <- (IC + IA) / (IC + IA + qR1 * (s_1 * RC + s_2 * RA) + qR2 * ((1 - s_1) * RC + (1 - s_2) * RA))
update(X_I) <- (IC + IA) / (IC + IA + qR1 * (s_1 * RC + s_2 * RA) + qR2 * ((1 - s_1) * RC + (1 - s_2) * RA))
initial(X_AP) <- (qR1 * (s_1 * RC + s_2 * RA)) / (IC + IA + qR1 * (s_1 * RC + s_2 * RA) + qR2 * ((1 - s_1) * RC + (1 - s_2) * RA))
update(X_AP) <- (qR1 * (s_1 * RC + s_2 * RA)) / (IC + IA + qR1 * (s_1 * RC + s_2 * RA) + qR2 * ((1 - s_1) * RC + (1 - s_2) * RA))
initial(X_ASP) <- (qR2 * ((1 - s_1) * RC + (1 - s_2) * RA)) / (IC + IA + qR1 * (s_1 * RC + s_2 * RA) + qR2 * ((1 - s_1) * RC + (1 - s_2) * RA))
update(X_ASP) <- (qR2 * ((1 - s_1) * RC + (1 - s_2) * RA)) / (IC + IA + qR1 * (s_1 * RC + s_2 * RA) + qR2 * ((1 - s_1) * RC + (1 - s_2) * RA))

