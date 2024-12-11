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


# SMC_effect <- decay[time] * eff_SMC * cov_SMC[time]
#
# # Children
# update(SC) <- SC * (1 - delta_a - delta_d - (1 - SMC_effect) * mu_SE_C) + delta_b * P + mu_RS_C * RC + mu_TS * TrC
# update(EC) <- EC * (1 - delta_a - delta_d - mu_EI) + (1 - SMC_effect) * mu_SE_C * SC
# update(IC) <- IC * (1 - delta_a - delta_d - mu_IR) + (1 - fT_C) * mu_EI * EC
# update(TrC) <- TrC * (1 - delta_a - delta_d - mu_TS) + fT_C * mu_EI * EC
# update(RC) <- RC * (1 - delta_a - delta_d - mu_RS_C) + mu_IR * IC
#
# # Adults
# update(SA) <- SA * (1 - delta_d - mu_SE_A) + delta_a * SC + mu_RS_A * RA + mu_TS * TrA
# update(EA) <- EA * (1 - delta_d - mu_EI) + delta_a * EC + mu_SE_A * SA
# update(IA) <- IA * (1 - delta_d - mu_IR) + delta_a * IC + (1 - fT_A) * mu_EI * EA
# update(TrA) <- TrA * (1 - delta_d - mu_TS) + delta_a * TrC + fT_A * mu_EI * EA
# update(RA) <- RA * (1 - delta_d - mu_RS_A) + delta_a * RC + mu_IR * IA



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
#dt <- user(1)
phi <- user(1)

mu_SE_C <- 1 - exp(-p_MH_C * EIR)
mu_SE_A <- phi * (1 - exp(-rho * p_MH_C * EIR))

initial(mu_SE_C_up) <- mu_SE_C
update(mu_SE_C_up) <- mu_SE_C

#mu_RS_C_0 <- user()
#mu_RS_C <- mu_RS_C_0 * dt
mu_RS_C <- user()
mu_RS_A <- eta * mu_RS_C
eta <- user(1)
#mu_EI_0 <- user()
#mu_EI <- mu_EI_0 * dt
mu_EI <- user()
#mu_TS_0 <- user()
#mu_TS <- mu_TS_0 * dt
mu_TS <- user()
#mu_IR_0 <- user()
#mu_IR <- mu_IR_0 * dt
mu_IR <- user()
#fT <- user()
fT_C <- user()
fT_A <- z * fT_C
z <- user(1) # z < 1, reporting rate of adults vs. children
rho <- user(1)
p_MH_C <- user()
#p_MH_A <- rho * p_MH_C
#delta_b_0 <- user()
#delta_b <- delta_b_0 * dt
delta_b <- user()
#delta_d_0 <- user()
#delta_d <- delta_d_0 * dt
delta_d <- user()
#delta_a_0 <- user()
#delta_a <- delta_a_0 * dt
delta_a <- user()
#day_count[] <- user()

# defining SMC efficacy
eff_SMC <- user() # SMC effectiveness
cov_SMC[] <- user() # SMC coverage
SMC[] <- user()
decay[] <- user()

## Defining climate-driven entomological inoculation rate (EIR)
A <- m * a^2 * exp(-g * n)
B <- g
C <- a
#EIR_C <- (A * p_HM_C * X / (B + ( C * p_HM_C * X)))
#EIR_A <- (A * p_HM_A * X / (B + ( C * p_HM_A * X)))
#p_HM_C <- user()
#p_HM_A <- theta * p_HM_C
# initial(EIR) <- (A * p_HM * X / (B + (C * p_HM * X)))
# update(EIR) <- (A * p_HM * X / (B + (C * p_HM * X)))
EIR <- (A * p_HM * X / (B + (C * p_HM * X)))
#update(EIR) <- (A * p_HM * X / (B + (C * p_HM * X)))
#EIR2 <- (A * p_HM * X / (B + (C * p_HM * X)))

# Define monthly EIR
initial(EIR_monthly) <- EIR
update(EIR_monthly) <- if ((step) %% steps_per_month == 0) EIR2 else EIR_monthly + EIR

# Define Daily EIR
initial(EIR2) <- EIR
update(EIR2) <- EIR

p_HM <- user()

qR <- user()
#qR2 <- user()
X <- (IA + IC + qR * (RA + RC)) / P
#X <- (qR2 * IA + IC + qR * (qR2 * RA + RC)) / P
#initial(X2) <- X
#update(X2) <- X

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

# Mosquito biting rate
a <- (0.017 * temp[time] - 0.165) * dt

# Sporogony (development of sporozites in mosquitos)
pT <- if (temp[time] >= 35) dt * 0.01 else (dt * (0.000112 * temp[time] * (temp[time] - 15.384) * (35 - temp[time])^(1/2)))
n <- 1 / pT

# Egg to adult survivorship
p_EA_T <- if (temp_w >= 33.3 || temp_w < 15.38) 0.01 else (-0.00924 * (temp_w)^2 + 0.453 * (temp_w) - 4.77)^dt
p_EA_R <- (1 / (1 + exp(-a_R * (c_R_D[time] - b_R))))^dt


a_R <- user()
b_R <- user()
c_R_D[] <- user() # this will be informed by data
p_EA <- p_EA_T * p_EA_R
shift <- user(1)

# Egg-adult Development time due to temperature
bM <- if (temp_w >= 34 || temp_w <= 14.8) 0.01 else (0.000111 * temp_w * (temp_w - 14.7) * (34-temp_w)^(1/2)) * dt
tau_EA <- 1 / bM

temp_w <- k_par * temp[time] + delta_temp
k_par <- user(1)
delta_temp <- user(2)

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
