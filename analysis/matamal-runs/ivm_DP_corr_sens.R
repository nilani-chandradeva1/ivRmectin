#getting model in right parameter space
unlink("src/*.o")
unlink("src/*.so")
devtools::load_all()
require(tidyverse)

init_EIR_vec <- c(12,14,16)
Q0_vec <- c(0.20, 0.82, 0.93)
bites_Bed_vec <- c(0.84, 0.95, 1)

acc_refusal_rate <- c(1-0.05, 1-0.21)
ivm_cov_vec <- c(0.609, 0.682, 0.747)

#getting lowest possible coverage and highest possible coverage of ivm
ivm_cov_vec_ref_adj <- c(ivm_cov_vec[1]*acc_refusal_rate[2], ivm_cov_vec[2], ivm_cov_vec[3])

years <- seq(2012, 2022, by = 1)
distr_years <- c(years[3], years[6], years[9])
distr_campaign <- (6*30)/365 #distributions are in June of the distr years

distr_campaign_days <- distr_campaign*365

#how many days into times is the first itn campaign?
time_period <- 365*13

net_seq <- seq(distr_campaign_days+(365*3), time_period, by = 3*365)
itn_on <- net_seq[1]

#the prevalence survey measuring ~14% qPCR all age prevalence and ~7.8% all age microscopy is
nov <- 30*11
y_2018 <- 365*7 # 8 years into simulation
prev_survey_date <- nov+y_2018

#put nets in the model: 55% phenotypic resistance acc Moss et al 2023
pyr_nets_df <- read.csv("../../Glasgow/data/pyrethroid_only_nets.csv")
pyr_nets_res <- pyr_nets_df %>%
  filter(resistance == 0.55)

pyr_param_df <- data.frame(pyr_nets_res$dn0_med, pyr_nets_res$rn0_med, pyr_nets_res$gamman_med)
names(pyr_param_df) <- c("dn0", "rn0", "gamman")

mod_space_df <- expand.grid(init_EIR_vec, Q0_vec, bites_Bed_vec)
names(mod_space_df) <- c("init_EIR","Q0", "bites_Bed")

mod_space_df$dn0 <- pyr_param_df$dn0
mod_space_df$rn0 <- pyr_param_df$rn0
mod_space_df$gamman <- pyr_param_df$gamman
mod_space_df$itn_cov <- 0.75



mod_space_df <- mod_space_df %>%
  filter(bites_Bed == bites_Bed_vec[2] & Q0 == Q0_vec[2])

mod_space_list <- list()

for(i in seq_len(nrow(mod_space_df))){
  mod_space_list[[i]] <- as.numeric(mod_space_df[i,])
}


#getting the ivermectin information
ivm_haz <- read.table("IVM_derivation/ivermectin_hazards.txt", header=TRUE)
colnames(ivm_haz) = c("Day", "IVM_400_1_HS", "IVM_300_3_HS")
# Sourcing the extra functions required to generate the endectocide-specific parameters
source("R/mda_ivm_functions.R")



runfun <- function(mod_name){
  mod <- mod_name$generator$new(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period, tcrit = net_seq)
  op<- mod$transform_variables(modx)
  return(op)
}

#endectocide set up####
mda_int <- 30
eff_len <- 23
july <- 30*7
y_2021 <- 365*10
y_2022 <- 365*11
endec_y1_start <- july + y_2021 #start time of first MDA
endec_y2_start <- july+ y_2022

endec_start <- c(endec_y1_start, endec_y1_start+mda_int, endec_y1_start+mda_int+mda_int,
                 endec_y2_start, endec_y2_start+mda_int, endec_y2_start+mda_int+mda_int)


#set IVM params for Hannah's mosq model with hazards#
ivm_parms1 <- ivRmectin::ivm_fun(
  #IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = endec_start,
  time_period = time_period,
  hazard_profile = ivm_haz$IVM_300_3_HS[1:23],
  #hazard_profile = hazzy,
  ivm_coverage= 0.7, #gets updated later on!!
  ivm_min_age=5,
  ivm_max_age = 90)

ivm_parms2 <- ivRmectin::ivm_fun(
  #IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = endec_start,
  time_period = time_period,
  hazard_profile = rep(1,23),
  #hazard_profile = hazzy,
  ivm_coverage= 0.7, #gets updated later on!!
  ivm_min_age=5,
  ivm_max_age = 90)


eff_len <- 23

#ENDECTOCIDE MDA######

ivm_params_df <- mod_space_df %>%
  filter(init_EIR == 14) #may need to update this to init_EIR_vec[best_ref]]

DP_scaler_vec <- c(0.1, 0.5, 0.9, 1)

ivm_params_df <- ivm_params_df %>%
  crossing(DP_scaler = DP_scaler_vec)




######################
#ivm_cov_par <- ivm_parms1$ivm_cov_par*0.95 #if 5% refusal rate. or by ~0.8 if 20% refusal rate

#ivm_cov = ivm_cov_par*(exp(-ivm_parms1$ivm_min_age/21) - exp(-ivm_parms1$ivm_max_age/21))

#######################
#so 70% amongst eligibles goes to 54% in total population, or 51% w 5% refusal (for the median)


ivm_parms_list <- list()

for(i in seq_len(nrow(ivm_params_df))){
  ivm_parms_list[[i]] <- as.numeric(ivm_params_df[i,])
}




#then set up the endectocide model (model A) with the MDAs in 2021 and 2022


ivm_mda_mod <- function(itn_type_ivm_param){
  init_EIR_in <- itn_type_ivm_param[1]
  Q0_in  <- itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  d_ITN0_in <- itn_type_ivm_param[4]
  r_ITN0_in <- itn_type_ivm_param[5]
  itn_half_life_in <- itn_type_ivm_param[6]
  itn_cov_in <- itn_type_ivm_param[7]
  DP_scaler_in <- itn_type_ivm_param[8]
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide_DP_cor.R", package = "ivRmectin"),
    num_int = 2,
    #num_int = 2, # number of vector control (IRS and ITN) population groups
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR_in, # the Entomological Innoculation Rate
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = 0.682, # proportion of population receiving the endectocide
    ivm_min_age = 5, # youngest age group receiving endectocide
    ivm_max_age = 90, # oldest age group receiving endectocide
    IVRM_start = ivm_parms1$IVRM_start,
    Q0 = Q0_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in*365,
    itn_cov = itn_cov_in,
    ITN_IRS_on = itn_on,
    DP_scaler = DP_scaler_in
  )
  return(output)
}

my_sim_ivm_mda_mod <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(ivm_parms_list, ivm_mda_mod)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(ivm_params_df)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, Ivtot, DP_scaler))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "ivm_DP_int_arm"))}, simplify = F))
  return(out_df_antag)

}

my_sim_ivm_mda_out <- my_sim_ivm_mda_mod()

#then set up the control arm
ivm_mda_mod_control <- function(itn_type_ivm_param){
  init_EIR_in <- itn_type_ivm_param[1]
  Q0_in  <- itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  d_ITN0_in <- itn_type_ivm_param[4]
  r_ITN0_in <- itn_type_ivm_param[5]
  itn_half_life_in <- itn_type_ivm_param[6]
  itn_cov_in <- itn_type_ivm_param[7]
  DP_scaler_in <- itn_type_ivm_param[8]
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide_DP_cor.R", package = "ivRmectin"),
    num_int = 2,
    #num_int = 2, # number of vector control (IRS and ITN) population groups
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR_in, # the Entomological Innoculation Rate
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms2$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms2$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms2$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = 0.682, # control arm:falsly turn on, but HR is 1
    ivm_min_age = 5, # youngest age group receiving endectocide
    ivm_max_age = 90, # oldest age group receiving endectocide
    IVRM_start = ivm_parms2$IVRM_start,
    Q0 = Q0_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in*365,
    itn_cov = itn_cov_in,
    ITN_IRS_on = itn_on,
    DP_scaler = DP_scaler_in
  )
  return(output)
}

my_sim_ivm_mda_mod_control <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(ivm_parms_list, ivm_mda_mod_control)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(ivm_params_df)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, Ivtot, DP_scaler))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "ivm_DP_control_arm"))}, simplify = F))
  return(out_df_antag)

}

my_sim_ivm_mda_out_control <- my_sim_ivm_mda_mod_control()

DP_cor_sens_df <- rbind(my_sim_ivm_mda_out, my_sim_ivm_mda_out_control)
saveRDS(DP_cor_sens_df, file = "analysis/matamal-runs/DP_cor_sens_df.rds")

