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

SMC_effect <- decay[time] * eff_SMC * cov_SMC[time]
#SMC_removal <- if (SMC[time] == 1) SMC_effect else 0
SMC_removal <- user(0)

# Children
update(SC) <- (SC * r_C) * (1 - delta_a - delta_d - (1 - SMC_effect) * mu_SE_C) + delta_b * P + mu_RS_C * RC + mu_TS * TrC
update(EC) <- (EC * r_C) * (1 - delta_a - delta_d - mu_EI - SMC_removal) + (1 - SMC_effect) * mu_SE_C * SC
update(IC) <- (IC * r_C) * (1 - delta_a - delta_d - mu_IR - SMC_removal) + (1 - fT_C) * mu_EI * EC
update(TrC) <- (TrC * r_C) * (1 - delta_a - delta_d - mu_TS) + fT_C * mu_EI * EC + SMC_removal * (EC + IC)
update(RC) <- (RC * r_C) * (1 - delta_a - delta_d - mu_RS_C) + mu_IR * IC


# Adults
update(SA) <- (SA * r_A) * (1 - delta_d - mu_SE_A) + delta_a * SC + mu_RS_A * RA + mu_TS * TrA
update(EA) <- (EA * r_A) * (1 - delta_d - mu_EI) + delta_a * EC + mu_SE_A * SA
update(IA) <- (IA * r_A) * (1 - delta_d - mu_IR) + delta_a * IC + (1 - fT_A) * mu_EI * EA
update(TrA) <- (TrA * r_A) * (1 - delta_d - mu_TS) + delta_a * TrC + fT_A * mu_EI * EA
update(RA) <- (RA * r_A) *  (1 - delta_d - mu_RS_A) + delta_a * RC + mu_IR * IA

# daily and weekly incidence
initial(day_inc_C) <- 0
update(day_inc_C) <- if ((step) %% steps_per_day == 0) mu_EI * EC * fT_C else day_inc_C + mu_EI * EC * fT_C

initial(wk_inc_C) <- 0
update(wk_inc_C) <- if ((step) %% steps_per_week == 0) mu_EI * EC * fT_C else wk_inc_C + mu_EI * EC * fT_C

initial(day_inc_A) <- 0
update(day_inc_A) <- if ((step) %% steps_per_day == 0) mu_EI * EA * fT_A else day_inc_A + mu_EI * EA * fT_A

initial(wk_inc_A) <- 0
update(wk_inc_A) <- if ((step) %% steps_per_week == 0) mu_EI * EA * fT_A else wk_inc_A + mu_EI * EA * fT_A

initial(day_inc_total) <- 0
update(day_inc_total) <- if ((step) %% steps_per_day == 0) mu_EI * (EC * fT_C + EA * fT_A) else day_inc_total + mu_EI * (EC * fT_C + EA * fT_A)

initial(wk_inc_total) <- 0
update(wk_inc_total) <-  if ((step) %% steps_per_week == 0) mu_EI * (EC * fT_C + EA * fT_A) else wk_inc_total + mu_EI * (EC * fT_C + EA * fT_A)

initial(month_inc_C) <- 0
update(month_inc_C) <- if ((step) %% steps_per_month == 0) mu_EI * EC * fT_C else month_inc_C + mu_EI * EC * fT_C

initial(month_inc_A) <- 0
update(month_inc_A) <- if ((step) %% steps_per_month == 0) mu_EI * EA * fT_A else month_inc_A + mu_EI * EA * fT_A

initial(month_inc_total) <- 0
update(month_inc_total) <- if ((step) %% steps_per_month == 0) mu_EI * (EC * fT_C + EA * fT_A) else month_inc_total + mu_EI * (EC * fT_C + EA * fT_A)

size <- user(10)

# User defined parameters
# Growth rates
#r_C <- user(1.000071) # daily growth rate u5 Chad
#r_A <- user(1.000092) # daily growth rate o5 Chad
r_C <- user(1.0000) # daily growth rate u5 Chad
r_A <- user(1.0000) # daily growth rate o5 Chad
phi <- user(1)
mu_SE_C <- 1 - exp(-p_MH_C * EIR)
mu_SE_A <- phi * (1 - exp(-rho * p_MH_C * EIR))
mu_RS_C <- user()
mu_RS_A <- eta * mu_RS_C
eta <- user(1)
mu_EI <- user()
mu_TS <- user()
mu_IR <- user()
fT_C <- user()
fT_A <- z * fT_C
z <- user(1) # z < 1, reporting rate of adults vs. children
rho <- user(1)
p_MH_C <- user()
delta_b <- user()
delta_d <- user()
delta_a <- user()
alpha <- user(1)
T_opt <- user(28)
sigma_LT <- user(4)
sigma_RT <- user(4)
R_opt <- user(1)
k1 <- user(0.2)
qR <- user()
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
qR2 <- user(1)
EIR <- alpha * (X / (b + X)) * temp_effect * rain_effect # Multiplicative effects
#temp_effect <- exp(-((temp_shift - T_opt)^2) / (2 * sigma_LT^2)) # Gaussian term for temperature
temp_effect <- if (temp_shift <= T_opt) exp(-((temp_shift - T_opt)^2) / (2 * sigma_LT^2)) else exp(-((temp_shift - T_opt)^2) / (2 * sigma_RT^2))
rain_effect <- 1 / (1 + exp(-k1 * (c_R_D_shift - R_opt))) # Logistic term for rainfall
X <- (qR2 * IA + IC + qR * (qR2 * RA + RC)) / P # Proportion of population infectious remains the same

# Egg-adult Development time due to temperature

# Dimensions of arrays
dim(c_R_D) <- user()
dim(temp) <- user()
dim(SMC) <- user()
dim(decay) <- user()
dim(cov_SMC) <- user()

## Initial conditions - user defined, defaults in parenthesis
N <- user()
N_pop <- s * N
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
s <- user(1)
initial(P_C) <- SC + EC + IC + TrC + RC
initial(P_A) <- SA + EA + IA + TrA + RA
update(P_C) <- SC + EC + IC + TrC + RC
update(P_A) <- SA + EA + IA + TrA + RA
P = (SC + EC + IC + TrC + RC + SA + EA + IA + TrA + RA)
