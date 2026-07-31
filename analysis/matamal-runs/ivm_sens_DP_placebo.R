#dates of the MDAs
timestep_from_2012 <- function(year, month, day) {
  # Days in each month (non-leap year)
  month_lengths <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

  # Years since 2012 (each year = 365 days)
  years_since <- year - 2012
  days_from_years <- years_since * 365

  # Days from previous months in the same year
  if (month > 1) {
    days_from_months <- sum(month_lengths[1:(month - 1)])
  } else {
    days_from_months <- 0
  }

  # Days within the current month (start at day 1 → offset 0)
  days_from_days <- day - 1

  # Total timestep (timestep = 1 corresponds to 1 Jan 2012)
  timestep <- 1 + days_from_years + days_from_months + days_from_days
  return(timestep)
}

#dates of the prevalence surveys
prev_survey_date <- timestep_from_2012(2018, 11, 1)
prev_survey_date_2021 <- timestep_from_2012(2021, 11, 1)
prev_survey_date_2022 <- timestep_from_2012(2022, 11, 1)

DP_july_21 <- timestep_from_2012(2021, 7, 1)
DP_aug_21 <- timestep_from_2012(2021, 8, 1)
DP_sept_21 <- timestep_from_2012(2021, 9, 1)

DP_july_22 <- timestep_from_2012(2022, 7, 1)
DP_aug_22 <- timestep_from_2012(2022, 8, 1)
DP_sept_22 <- timestep_from_2012(2022, 9, 1)

#DP_mda_y_2021 <- c(y_2021+july, y_2021+aug, y_2021+sept)
DP_mda_y_2021 <- c(DP_july_21, DP_aug_21, DP_sept_21)
DP_mda_y_2022 <- c(DP_july_22, DP_aug_22, DP_sept_22)

mda_events <- c(DP_mda_y_2021, DP_mda_y_2022)



#dates of the MDAs
mda_int <- 30
eff_len <- 23
july <- 30*7
y_2021 <- 365*10
y_2022 <- 365*11
endec_y1_start <- july + y_2021 #start time of first MDA
endec_y2_start <- july+ y_2022

endec_start <- c(endec_y1_start, endec_y1_start+mda_int, endec_y1_start+mda_int+mda_int,
                 endec_y2_start, endec_y2_start+mda_int, endec_y2_start+mda_int+mda_int)

#dates of the prevalence surveys
prev_survey_date_2021 <- endec_start[3]+60
prev_survey_date_2022 <- endec_start[6]+60

nov <- 30*11
y_2018 <- 365*7 # 8 years into simulation
prev_survey_date <- nov+y_2018

mda_campaign_arrows <- data.frame(
  x = ((endec_start)/365)+2011,
  xend = ((endec_start)/365)+2011,
  y = rep(0.2*100, length(distr_years)),
  yend = rep(0.17*100, length(distr_years))
)


#prev survey 2021 in int arm

placebo_trial <- readRDS("analysis/matamal-runs/placebo_sens_df.rds")
DP_trial <- readRDS("analysis/matamal-runs/DP_cor_sens_df.rds")

ggplot(placebo_trial)+
  geom_line(aes(x = ((t)/365)+2011,
                y = slide_prev0to80*100,
                lty = as.factor(scenario)))+
  theme_bw()+
  coord_cartesian(xlim = c(2018, 2024), ylim = c(0,40)) +
  scale_x_continuous(breaks = 2018:2024)+
  xlab("Year")+
  ylab("All-age prevalence (%)")+
  #show MDA distributions
  geom_segment(
    data = mda_campaign_arrows,
    aes(x = x, xend = xend, y = y, yend = yend),
    arrow = arrow(length = unit(0.25, "cm")),
    size = 1, col = "orange"
  )+
  scale_linetype_manual(
    values = c(
      ivm_placebo_control_arm = "solid",
      ivm_placebo_int_arm = "dashed"
    ),
    labels = c(
      ivm_placebo_control_arm = "Control: placebo + standard of care",
      ivm_placebo_int_arm = "Intervention arm: ivermectin + standard of care"
    ),
    name = "Trial arm"
  )+
  theme(legend.position = c(0.5, 0.9))


DP_pals <- c('#e5f5f9','#99d8c9','#2ca25f', "black")



plot_dynamics_sens <- ggplot(DP_trial)+
  geom_line(aes(x = ((t-1)/365)+2011,
                y = slide_prev0to80*100,
                lty = as.factor(scenario),
            col = as.factor(DP_scaler)), size = 0.7)+
  theme_bw()+
  coord_cartesian(xlim = c(2018, 2023.5), ylim = c(0,40)) +
  scale_x_continuous(breaks = 2018:2023.5)+
  xlab("Year")+
  ylab("All-age prevalence (%)")+
  #show MDA distributions
  geom_segment(
    data = mda_campaign_arrows,
    aes(x = x, xend = xend, y = y, yend = yend),
    arrow = arrow(length = unit(0.25, "cm")),
    size = 1, col = "orange"
  )+
  scale_linetype_manual(
    values = c(
      ivm_DP_control_arm = "solid",
      ivm_DP_int_arm = "dashed"
    ),
    labels = c(
      ivm_DP_control_arm = "Control: DP + standard of care",
      ivm_DP_int_arm = "Intervention arm: ivermectin + standard of care"
    ),
    name = "Trial arm"
  )+
  theme(legend.position = c(0.3, 0.7))+
  scale_color_manual(labels = c("0.1" = "90%",
                                "0.5" = "50%",
                                "0.9" = "10%",
                                "1" = "No reduction"),
                     values = DP_pals,
                     name = "Assumed reduction in transmissability due to DP uptake amongst\n
                     individuals treated with an ivermectin-like endectocide")

##doesn't seem to change the control arm prevalence that much

eff_placebo_trial <- placebo_trial %>%
  filter(t == prev_survey_date_2022) %>%
  group_by(scenario) %>%
  select(slide_prev0to80) %>%
  pivot_wider(names_from = scenario,
              values_from = slide_prev0to80) %>%
  mutate(impact_trial_rel = ((ivm_placebo_control_arm - ivm_placebo_int_arm)/ivm_placebo_control_arm)*100)

k0_2021 <- 0.67
k1_2021 <- 0.37
k0_2022 <- 0.55
k1_2022 <- 0.54

eta <- 167.3134 #for MATAMAL

z_alpha <- 1.96#z_alpha/2
zb <- function(za = z_alpha, k0, k1, pi0, pi1, c, m){
  aux1 <- -za + sqrt(((c-1)*(pi0 - pi1)**2)/((pi0*(1-pi0)/m) +(pi1*(1-pi1)/m) + (k0**2)*(pi0**2) + (k1**2)*(pi1**2)))
  #aux2 <- -za - sqrt(((c-1)*(pi0 - pi1)**2)/((pi0*(1-pi0)/m) +(pi1*(1-pi1)/m) + (k0**2)*(pi0**2) + (k1**2)*(pi1**2)))
  return(aux1)#return(c(aux1, aux2))
}


eff_DP_trial <- DP_trial %>%
  filter(t %in% c(prev_survey_date_2021, prev_survey_date_2022)) %>%
  group_by(t,scenario, DP_scaler) %>%
  select(t, scenario, DP_scaler, slide_prev0to80) %>%
  pivot_wider(names_from = scenario,
               values_from = slide_prev0to80) %>%
  mutate(impact_trial_rel = ((ivm_DP_control_arm - ivm_DP_int_arm)/ivm_DP_control_arm)*100,
         survey_date = case_when(t == prev_survey_date_2021 ~ "survey_2021",
                                 t == prev_survey_date_2022 ~ "survey_2022")) %>%
  mutate(k0_obs = case_when(survey_date == "survey_2021" ~ k0_2021,
                            survey_date == "survey_2022" ~ k0_2022),
         k1_obs = case_when(survey_date == "survey_2021" ~ k1_2021,
                            survey_date == "survey_2022" ~ k1_2022),
         zb_k = zb(
           k0 = k0_obs,
           k1 = k1_obs,
           pi0 = ivm_DP_control_arm,
           pi1 = ivm_DP_int_arm,
           c = 12,
           m = eta
         ),
         power = pnorm(zb_k, mean = 0, sd = 1, lower.tail = T)*100)

ggplot(eff_DP_trial, aes(x = as.factor(survey_date), y = as.factor(DP_scaler), fill = impact_trial_rel))+
  geom_tile()
#need to make this a bit better, add the power as text (all are underpowered.) Check by also add in the 0.05 and 0.1 and see if recover 80% power





#checking step, can see that the eff_placebo_trial and eff_DP_trial are same for the DP_scaler of 1, so implementation works

1-((26.7-14.7)/26.7) #efficacy could be up to 55% lower because of this correlation

#then given the k, what would the power be?
