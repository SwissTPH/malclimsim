########################
## Defining time steps #
########################
steps_per_day <- user(1)
steps_per_week <- steps_per_day * 7
steps_per_month <- steps_per_day * 30
dt <- 1 / steps_per_day
time <- step + 1

################################
## Defining SMC Vector Inputs  #
################################
cov_SMC[] <- user() # SMC coverage
SMC[] <- user()
decay[] <- user()

###########################
## Defining Effect of SMC #
###########################
eff_SMC <- user() # SMC effectiveness
SMC_effect <- decay[time] * eff_SMC * cov_SMC[time]

#######################
## Force of Infection #
#######################
mu_SE <- (1 - SMC_effect) *(1 - exp(-p_MH_C * EIR))
p_MH_C <- user(0.50)

######################################
## State Equations                   #
## - split into adults and children  #
######################################
# Children
update(SC) <- SC -  mu_SE * SC + mu_RS * RC + mu_TS * TrC - (delta_d + delta_a) * SC + delta_b * P
update(EC) <- EC - mu_EI * EC +  mu_SE * SC - (delta_d + delta_a) * EC
update(IC) <- IC - mu_IR * IC + pi_s_1 * (1 - fT_C) * mu_EI * EC - (delta_d + delta_a) * IC
update(TrC) <- TrC - mu_TS * TrC + pi_s_1 * fT_C * mu_EI * EC - (delta_d + delta_a) * TrC
update(RC) <- RC - mu_RS * RC + mu_IR * IC + (1 - pi_s_1) * mu_EI * EC - (delta_d + delta_a) * RC

# Adults
update(SA) <- SA -  mu_SE * SA + mu_RS * RA + mu_TS * TrA - delta_d * SA + delta_a * SC
update(EA) <- EA - mu_EI * EA +  mu_SE * SA - delta_d * EA + delta_a * EC
update(IA) <- IA - mu_IR * IA + pi_s_2 * (1 - fT_A) * mu_EI * EA - delta_d * IA + delta_a * IC
update(TrA) <- TrA - mu_TS * TrA + pi_s_2 * fT_A * mu_EI * EA - delta_d * TrA + delta_a * TrC
update(RA) <- RA - mu_RS * RA + mu_IR * IA + (1 - pi_s_2) * mu_EI * EA - delta_d * RA + delta_a * RC

#################################################
## Defining Incidence                           #
## - each variables represents number of cases  #
## in the given time period - day, week, month  #
#################################################
initial(day_inc_C) <- 0
update(day_inc_C) <- if ((step) %% steps_per_day == 0) mu_EI * EC * fT_C * pi_s_1 else day_inc_C + mu_EI * EC * fT_C * pi_s_1

initial(wk_inc_C) <- 0
update(wk_inc_C) <- if ((step) %% steps_per_week == 0) mu_EI * EC * fT_C * pi_s_1 else wk_inc_C + mu_EI * EC * fT_C * pi_s_1

initial(day_inc_A) <- 0
update(day_inc_A) <- if ((step) %% steps_per_day == 0) mu_EI * EA * fT_A * pi_s_2 else day_inc_A + mu_EI * EA * fT_A * pi_s_2

initial(wk_inc_A) <- 0
update(wk_inc_A) <- if ((step) %% steps_per_week == 0) mu_EI * EA * fT_A * pi_s_2 else wk_inc_A + mu_EI * EA * fT_A * pi_s_2

initial(day_inc_total) <- 0
update(day_inc_total) <- if ((step) %% steps_per_day == 0) mu_EI * (EC * fT_C * pi_s_1 + EA * fT_A * pi_s_2) else day_inc_total + mu_EI * (EC * fT_C * pi_s_1 + EA * fT_A * pi_s_2)

initial(wk_inc_total) <- 0
update(wk_inc_total) <-  if ((step) %% steps_per_week == 0) mu_EI * (EC * fT_C * pi_s_1 + EA * fT_A * pi_s_2) else wk_inc_total + mu_EI * (EC * fT_C * pi_s_1 + EA * fT_A * pi_s_2)

initial(month_inc_C) <- 0
update(month_inc_C) <- if ((step) %% steps_per_month == 0) mu_EI * EC * fT_C * pi_s_1 else month_inc_C + mu_EI * EC * fT_C * pi_s_1

initial(month_inc_A) <- 0
update(month_inc_A) <- if ((step) %% steps_per_month == 0) mu_EI * EA * fT_A * pi_s_2 else month_inc_A + mu_EI * EA * fT_A * pi_s_2

initial(month_inc_total) <- 0
update(month_inc_total) <- if ((step) %% steps_per_month == 0) mu_EI * (EC * fT_C * pi_s_1 + EA * fT_A * pi_s_2) else month_inc_total + mu_EI * (EC * fT_C * pi_s_1 + EA * fT_A * pi_s_2)

########################################
## Defining Transition Rate Parameters #
########################################
mu_RS <- user()
duration_infection <- user(200)
mu_IR <- 1 / (duration_infection -  1 / mu_RS) # I to R is a maximum of 200 days
mu_EI <- user()
mu_TS <- user()

#################################
## Defining Case Reporting Rate #
#################################
fT_C <- user()
fT_A <- z * fT_C
z <- user(1) # z < 1, reporting rate of adults vs. children

####################################################
## Defining Different Infectivity by Compartment   #
####################################################
qR <- user()

######################################################
## Defining Proportion Symptomatic and Asymptomatic  #
######################################################
# s => symptomatic
# c => some constant
# pi => proportion

# Proportion of population that is symptomatic
c_s <- user()
pi_s_1 <- user()
pi_s_2 <- c_s * pi_s_1 # >=5

##########################################
## Defining Demographic Rate Parameters  #
##########################################
delta_b <- user()
delta_d <- user()
delta_a <- user()

################################
## Defining Daily Growth Rates #
################################
r_C_0 <- user(1.000079)
r_A_0 <- user(1.000108)

initial(r_C) <- r_C_0
update(r_C) <- r_C * r_C_0

initial(r_A) <- r_A_0
update(r_A) <- r_A * r_A_0

###################################
## Defining Climate Vector Inputs #
###################################
temp[] <- user()
c_R_D[] <- user() # this will be informed by data

###########################
## Defining Climate Lags  #
###########################
lag_R <- user(0)  # Default lag = 0
lag_T <- user(0)
clim_SMC_lag <- user()

c_R_D_shift <- c_R_D[time + clim_SMC_lag - lag_R]
temp_shift <- temp[time + clim_SMC_lag - lag_T]

#initial(c_R_D_shift) <- c_R_D[time]  # Start with the first rainfall value
#update(c_R_D_shift) <- if (time > lag_R) c_R_D[time - lag_R] else c_R_D[time]
#initial(temp_shift) <- temp[time]
#update(temp_shift) <- if (time > lag_T) temp[time - lag_T] else temp[time]

##############################################
## Defining Entomological Inncoulation Rate  #
##############################################
EIR <- alpha * (X / (b + X)) * temp_effect * rain_effect # Multiplicative effects
temp_effect <- if (temp_shift <= T_opt) exp(-((temp_shift - T_opt)^2) / (2 * sigma_LT^2)) else exp(-((temp_shift - T_opt)^2) / (2 * sigma_RT^2))
rain_effect <- 1 / (1 + exp(-k1 * (c_R_D_shift - R_opt))) # Logistic term for rainfall

# Others
alpha <- user(1) # parameter to adjust magnitude of EIR
T_opt <- user(27.82) # temperature at which EIR is maximized
sigma_LT <- user(4) # standard deviation to left of T_opt
sigma_RT <- user(4) # standard deviation to right of T_opt
R_opt <- user(1) # rainfall at which EIR is maximized
k1 <- user(0.2) # controls the "rate" that EIR changes with rainfall
b <- user(0.65) # controls saturation rate of impact of population infectivity on EIR
#c_X <- user(0.35) # controls saturation rate of impact of population infectivity on EIR

############################################
## Defining Infectivity of the Population  #
############################################
X <- (IC + IA + qR * (RC + RA)) / P

########################################
## Dispersion Parameters               #
## - these are necessary for inference #
########################################
size <- user(10)
size_1 <- user(10)
size_2 <- user(10)

###############################################
## Initial conditions; Population parameters  #
###############################################
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

# Total population size
initial(P_C) <- SC + EC + IC + TrC + RC
initial(P_A) <- SA + EA + IA + TrA + RA
update(P_C) <- SC + EC + IC + TrC + RC
update(P_A) <- SA + EA + IA + TrA + RA
P = (SC + EC + IC + TrC + RC + SA + EA + IA + TrA + RA)

##################################################################
## Defining Vector Input Dimensions                              #
## - these need to be defined for technical reasons              #
## - they should "know" the correct dimensions from the inputs   #
##################################################################
dim(c_R_D) <- user()
dim(temp) <- user()
dim(SMC) <- user()
dim(decay) <- user()
dim(cov_SMC) <- user()

#################################
## Relevant outputs to look at  #
#################################
initial(SMC_effect_2) <- decay[time] * eff_SMC * cov_SMC[time]
update(SMC_effect_2) <- decay[time] * eff_SMC * cov_SMC[time]
initial(mu_SE_C_2) <-(1 - SMC_effect) * (1 - exp(-p_MH_C * EIR))
update(mu_SE_C_2) <- (1 - SMC_effect) *(1 - exp(-p_MH_C * EIR))
#initial(EIR2) <- alpha * (X / (b + X)) * temp_effect * rain_effect # Multiplicative effects
#update(EIR2) <- alpha * (X / (b + X)) * temp_effect * rain_effect # Multiplicative effects
initial(temp_effect_2) <- if (temp_shift <= T_opt) exp(-((temp_shift - T_opt)^2) / (2 * sigma_LT^2)) else exp(-((temp_shift - T_opt)^2) / (2 * sigma_RT^2))
update(temp_effect_2) <- if (temp_shift <= T_opt) exp(-((temp_shift - T_opt)^2) / (2 * sigma_LT^2)) else exp(-((temp_shift - T_opt)^2) / (2 * sigma_RT^2))
initial(rain_effect_2) <- 1 / (1 + exp(-k1 * (c_R_D_shift - R_opt))) # Logistic term for rainfall
update(rain_effect_2) <-  1 / (1 + exp(-k1 * (c_R_D_shift - R_opt))) # Logistic term for rainfall
