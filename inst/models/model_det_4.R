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
time_t <- step + shift_t
time_r <- step + shift_r

SMC_effect <- decay[time] * eff_SMC * cov_SMC[time]
SMC_removal <- if (SMC[time] == 1) SMC_effect else 0

# Children
update(SC) <- SC * (1 - delta_a - delta_d - (1 - SMC_effect) * mu_SE_C) + delta_b * P + mu_RS_C * RC + mu_TS * TrC
update(EC) <- EC * (1 - delta_a - delta_d - mu_EI - SMC_removal) + (1 - SMC_effect) * mu_SE_C * SC
update(IC) <- IC * (1 - delta_a - delta_d - mu_IR - SMC_removal) + (1 - fT_C) * mu_EI * EC
update(TrC) <- TrC * (1 - delta_a - delta_d - mu_TS) + fT_C * mu_EI * EC + SMC_removal * (EC + IC)
update(RC) <- RC * (1 - delta_a - delta_d - mu_RS_C) + mu_IR * IC


# Adults
update(SA) <- SA * (1 - delta_d - mu_SE_A) + delta_a * SC + mu_RS_A * RA + mu_TS * TrA
update(EA) <- EA * (1 - delta_d - mu_EI) + delta_a * EC + mu_SE_A * SA
update(IA) <- IA * (1 - delta_d - mu_IR) + delta_a * IC + (1 - fT_A) * mu_EI * EA
update(TrA) <- TrA * (1 - delta_d - mu_TS) + delta_a * TrC + fT_A * mu_EI * EA
update(RA) <- RA * (1 - delta_d - mu_RS_A) + delta_a * RC + mu_IR * IA

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

# Dispersion parameter for negative binomial distribution
size <- user()

# User defined parameters
phi <- user()
qR <- user()
p_HM <- user()
mu_RS_C <- user()
mu_RS_A <- eta * mu_RS_C
eta <- user()
mu_EI <- user()
mu_TS <- user()
mu_IR <- user()
fT_C <- user()
fT_A <- z * fT_C
z <- user() # z < 1, reporting rate of adults vs. children
rho <- user()
p_MH_C <- user()
delta_b <- user()
delta_d <- user()
delta_a <- user()
shift_r <- user(1)
shift_t <- user(1)

# defining SMC efficacy
eff_SMC <- user() # SMC effectiveness
cov_SMC[] <- user() # SMC coverage
SMC[] <- user()
decay[] <- user()

## Defining climate-driven entomological inoculation rate (EIR)
A <- m * a^2 * exp(-g * n)
B <- g
C <- a
EIR <- (A * p_HM * X / (B + (C * p_HM)))
X <- (IA + IC + qR * (RA + RC)) / P

# Defining Force of Infection
mu_SE_C <- 1 - exp(-p_MH_C * EIR)
mu_SE_A <- phi * (1 - exp(-rho * p_MH_C * EIR))

# Mosquito density
initial(m) <- (B_egg * p_EA) / (tau_EA * g)
update(m) <- (B_egg * p_EA) / (tau_EA * g)

# Lifetime number of eggs
B_egg <- e / (exp(GP * g)-1)
p_surv <- user(0.98)
g <- -log(p_surv^dt)  # mu_M in paper
e0 <- user(50)
e <- e0 * dt
GP <- 1 / a
temp[] <- user()

# Temperature relationships
k1 <- user(1)        # Coefficient for temperature update, default 1
c1 <- user(0)        # Constant for temperature update, default 0
temp_adj <- k1 * temp[time_t] + c1  # Adjusted temperature equation

# Weighted temperature relationship
k2 <- user(1)        # Coefficient for weighted temperature update, default 1
delta_temp <- user(2)  # Temperature shift, default 2 (based on original value in #2)
temp_w <- k2 * temp_adj + delta_temp  # Weighted temperature equation

# Mosquito biting rate
k3 <- user(0.017)    # Coefficient for mosquito biting rate, default 0.017
c3 <- user(0.165)    # Constant for mosquito biting rate, default 0.165
a <- (k3 * temp_adj - c3) * dt  # Updated biting rate equation, includes dt

# Egg-to-adult survival probability based on temperature
k4 <- user(-0.00924) # Coefficient for quadratic term, default -0.00924
k5 <- user(0.453)    # Coefficient for linear term, default 0.453
c4 <- user(4.77)     # Constant for p_EA_T, default 4.77

p_EA_T <- if (temp_w >= 33.3 || temp_w < 15.38) 1e-6 else (k4 * temp_w^2 + k5 * temp_w - c4)^dt  # Includes dt

# Egg-to-adult survival probability based on rainfall
a_R <- user()         # Rainfall response coefficient
b_R <- user()         # Rainfall threshold parameter
c_R_D[] <- user()     # Rainfall data array
p_EA_R <- (1 / (1 + exp(-a_R * (c_R_D[time_r] - b_R))))^dt  # Includes dt

# Total probability of egg-to-adult survival
p_EA <- p_EA_T * p_EA_R

# Egg-adult development time due to temperature
c5 <- user(14.7)      # Offset constant for bM, default 14.7
c6 <- user(34)        # Upper limit constant for bM, default 34

bM <- if (temp_w >= c6 || temp_w <= c5) 1e-6 else (k5 * temp_w * (temp_w - c5) * (c6 - temp_w)^(1/2)) * dt  # Includes dt

tau_EA <- 1 / bM  # Egg-to-adult development time

# Sporogony (development of sporozites in mosquitos)
k7 <- user(0.000112)   # Coefficient for sporogony equation, default 0.000112
c7 <- user(15.384)     # Lower temperature threshold, default 15.384
c8 <- user(35)         # Upper temperature threshold, default 35
c9 <- user(0.01)       # Constant sporogony rate, default 0.01

pT <- if (temp_adj >= c8) 1e-6 else dt * (k7 * temp_adj * (temp_adj - c7) * (c8 - temp_adj)^(1/2))

n <- 1 / pT  # Mosquito sporogony time


# Dimensions of arrays
dim(c_R_D) <- user()
dim(temp) <- user()
dim(SMC) <- user()
dim(decay) <- user()
dim(cov_SMC) <- user()
#dim(day_count) <- user()

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
