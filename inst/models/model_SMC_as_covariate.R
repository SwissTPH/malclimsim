## The model contained in the file 'age_model_equation_final' will be written
## using odin and inference in mcstate. Here, this is done with the
## simplification of not including the age structure or effect of SMC,
## which will be added later.

## rainfall in mL; temp in Celsius
## it is timeed an identical array at every time step

# Definition of time step
steps_per_day <- user(1)
steps_per_week <- steps_per_day * 7
steps_per_month <- steps_per_day * 30
dt <- 1 / steps_per_day
time <- step + 1

# Force of infection
mu_SE_C <- (1 - exp(-p_MH_C * EIR))
initial(mu_SE_C_2) <- (1 - exp(-p_MH_C * EIR))
update(mu_SE_C_2) <- (1 - exp(-p_MH_C * EIR))

mu_SE_A <- (1 - exp(-rho * p_MH_C * EIR))
initial(mu_SE_A_2) <- (1 - exp(-rho * p_MH_C * EIR))
update(mu_SE_A_2) <- (1 - exp(-rho * p_MH_C * EIR))

######################################
## State Equations                   #
## - split into adults and children  #
######################################
# Children - <5
update(SC) <- SC - mu_SE_C * SC + mu_RS_C * RC + mu_TS * TrC - (delta_d + delta_a) * SC + delta_b * P
update(EC) <- EC - mu_EI * EC + mu_SE_C * SC - (delta_d + delta_a) * EC
update(IC) <- IC - mu_IR * IC + phi_1 * (1 - fT_C) * mu_EI * EC - (delta_d + delta_a) * IC
update(TrC) <- TrC - mu_TS * TrC + phi_1 * fT_C * mu_EI * EC - (delta_d + delta_a) * TrC
update(RC) <- RC - mu_RS_C * RC + mu_IR * IC + (1 - phi_1) * mu_EI * EC - (delta_d + delta_a) * RC

# Adults - >=5
update(SA) <- SA - mu_SE_A * SA + mu_RS_A * RA + mu_TS * TrA - delta_d * SA + delta_a * SC
update(EA) <- EA - mu_EI * EA + mu_SE_A * SA - delta_d * EA + delta_a * EC
update(IA) <- IA - mu_IR * IA + phi_2 * (1 - fT_A) * mu_EI * EA - delta_d * IA + delta_a * IC
update(TrA) <- TrA - mu_TS * TrA + phi_2 * fT_A * mu_EI * EA - delta_d * TrA + delta_a * TrC
update(RA) <- RA - mu_RS_A * RA + mu_IR * IA + (1 - phi_2) * mu_EI * EA - delta_d * RA + delta_a * RC

#################################################
## Defining Incidence                           #
## - each variables represents number of cases  #
## in the given time period - day, week, month  #
#################################################
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

initial(month_inc_C) <- 0
update(month_inc_C) <- if ((step) %% steps_per_month == 0) mu_EI * EC * fT_C * phi_1 else month_inc_C + mu_EI * EC * fT_C * phi_1

initial(month_inc_A) <- 0
update(month_inc_A) <- if ((step) %% steps_per_month == 0) mu_EI * EA * fT_A * phi_2 else month_inc_A + mu_EI * EA * fT_A * phi_2

initial(month_inc_total) <- 0
update(month_inc_total) <- if ((step) %% steps_per_month == 0) mu_EI * (EC * fT_C * phi_1 + EA * fT_A * phi_2) else month_inc_total + mu_EI * (EC * fT_C * phi_1 + EA * fT_A * phi_2)

#########################
## Defining Prevalence  #
#########################
initial(prev_total_1) <- (IA + RA * s_2  + IC + RC * s_1) / N
update(prev_total_1) <- (IA + RA * s_2  + IC + RC * s_1) / N

initial(prev_C_1) <- (IC + RC * s_1) / N_C
update(prev_C_1) <- (IC + RC * s_1) / N_C

initial(prev_A_1) <- (IA + RA * s_2) / N_A
update(prev_A_1) <- (IA + RA * s_2) / N_A

########################################
## State Equation and EIR Parameters   #
########################################
# To be fitted for better reproducing prevalence
mu_RS_C <- user()
duration_infection <- user(200)
mu_IR <- 1 / (duration_infection -  1 / mu_RS_C) # I to R is a maximum of 200 days

mu_RS_A <- eta * mu_RS_C
eta <- user(1)

# Immunity related parameters
qR2 <- user()
qR1 <- qR2 + (1 - qR2) * c_qR
c_phi <- user() # alters
c_qR <- user()
c_s <- user()
prop_p_1 <- user() # prop of all malaria infections that are patent in <5
c_p_1 <- user() # prop of patent infections that are symptomatic in <5
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

##############################
## Defining climate vectors  #
##############################
temp[] <- user()
c_R_D[] <- user() # this will be informed by data
c_R_D_shift <- if (time > lag_R) c_R_D[time - lag_R] else c_R_D[time]
temp_shift <- if (time > lag_T) temp[time - lag_T] else temp[time]

# Defining EIR
EIR <- alpha * (X / (b + X)) * temp_effect * rain_effect # Multiplicative effects
temp_effect <- if (temp_shift <= T_opt) exp(-((temp_shift - T_opt)^2) / (2 * sigma_LT^2)) else exp(-((temp_shift - T_opt)^2) / (2 * sigma_RT^2))
rain_effect <- 1 / (1 + exp(-k1 * (c_R_D_shift - R_opt))) # Logistic term for rainfall

X <- (IC + IA + qR1 * (s_1 * RC + s_2 * RA) + qR2 * ((1 - s_1) * RC + (1 - s_2) * RA)) / P # Proportion of population infectious remains the same

# Dimensions of arrays
dim(c_R_D) <- user()
dim(temp) <- user()

###############################################
## Initial conditions; Population parameters  #
###############################################
# Default values in parenthesis
N <- user()
N_pop <- N
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

###################################################################
## Parameters that are used in observation function               #
## - some parameters are defined within the observation function  #
## that is defined with the `observation_functions.R` file        #
###################################################################
beta_1 <- user() # SMC effectiveness parameter
size_1 <- user(10) # Negative binomial dispersion parameter for <5
size_2 <- user(10) # Negative binomial dispersion parameter for >=5
kappa_C <- user(30) # Prevalence dispersion parameter for <5
kappa_A <- user(30) # Prevalence dispersion parameter for >=5

######################################
## Variables we are interested in    #
## viewing outside the model         #
######################################
# Entomological innoculation rate
initial(EIR2) <- alpha * (X / (b + X)) * temp_effect * rain_effect # Multiplicative effects
update(EIR2) <- alpha * (X / (b + X)) * temp_effect * rain_effect # Multiplicative effectsS

# Contribution of climate to EIR
initial(temp_effect_2) <- if (temp_shift <= T_opt) exp(-((temp_shift - T_opt)^2) / (2 * sigma_LT^2)) else exp(-((temp_shift - T_opt)^2) / (2 * sigma_RT^2))
update(temp_effect_2) <- if (temp_shift <= T_opt) exp(-((temp_shift - T_opt)^2) / (2 * sigma_LT^2)) else exp(-((temp_shift - T_opt)^2) / (2 * sigma_RT^2))
initial(rain_effect_2) <- 1 / (1 + exp(-k1 * (c_R_D_shift - R_opt))) # Logistic term for rainfall
update(rain_effect_2) <-  1 / (1 + exp(-k1 * (c_R_D_shift - R_opt))) # Logistic term for rainfall

# Contribution of infectivity to EIR
initial(X2) <- (IC + IA + qR1 * (s_1 * RC + s_2 * RA) + qR2 * ((1 - s_1) * RC + (1 - s_2) * RA)) / P  # Proportion of population infectious remains the same
update(X2) <- (IC + IA + qR1 * (s_1 * RC + s_2 * RA) + qR2 * ((1 - s_1) * RC + (1 - s_2) * RA)) / P  # Proportion of population infectious remains the same

# Contribution of different compartments to EIR
initial(X_I) <- (IC + IA) / (IC + IA + qR1 * (s_1 * RC + s_2 * RA) + qR2 * ((1 - s_1) * RC + (1 - s_2) * RA))
update(X_I) <- (IC + IA) / (IC + IA + qR1 * (s_1 * RC + s_2 * RA) + qR2 * ((1 - s_1) * RC + (1 - s_2) * RA))
initial(X_AP) <- (qR1 * (s_1 * RC + s_2 * RA)) / (IC + IA + qR1 * (s_1 * RC + s_2 * RA) + qR2 * ((1 - s_1) * RC + (1 - s_2) * RA))
update(X_AP) <- (qR1 * (s_1 * RC + s_2 * RA)) / (IC + IA + qR1 * (s_1 * RC + s_2 * RA) + qR2 * ((1 - s_1) * RC + (1 - s_2) * RA))
initial(X_ASP) <- (qR2 * ((1 - s_1) * RC + (1 - s_2) * RA)) / (IC + IA + qR1 * (s_1 * RC + s_2 * RA) + qR2 * ((1 - s_1) * RC + (1 - s_2) * RA))
update(X_ASP) <- (qR2 * ((1 - s_1) * RC + (1 - s_2) * RA)) / (IC + IA + qR1 * (s_1 * RC + s_2 * RA) + qR2 * ((1 - s_1) * RC + (1 - s_2) * RA))

# Misc
initial(c_R_D_shift_2) <- c_R_D[time]  # Start with the first rainfall value
update(c_R_D_shift_2) <- if (time > lag_R) c_R_D[time - lag_R] else c_R_D[time]
initial(temp_shift_2) <- temp[time]
update(temp_shift_2) <- if (time > lag_T) temp[time - lag_T] else temp[time]
