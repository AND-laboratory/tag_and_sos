### SOS (TAG) bio data final analyses
library(utils)
workdir <- 'S:/MNHS-Psych/ANDL-Lab-TAG-Study/SOS Study'
sos_hr <- read.csv(file.path(workdir,"/physio/???.csv", fsep=""))

sos_crp <- read.csv(file.path(workdir,"/saliva/Cleaned_datasets/SOS_CRP_cleaned.csv", fsep="")) 
sos_cytokines <- read.csv(file.path(workdir,"/saliva/Cleaned_datasets/SOS_cytokines_cleaned.csv", fsep=""))
sos_hormones <- read.csv(file.path(workdir,"/saliva/Cleaned_datasets/SOS_hormones_cleaned.csv", fsep=""))

### Heart rate

# Report heart rate in:
# Resting 1
# Resting 2
# Average resting state
# Average all task
# Self-connected condition
# Self-disconnected condition
# Other condition

# Check significant differences between:
# Average resting state vs. average all task
# Average resting state vs. self-connected
# Self-connected vs. self-disconnected
# Self-connected vs. other


### Saliva

# Report some basic plots taken from:
# S:\MNHS-Psych\ANDL-Lab-TAG-Study\SOS Study\saliva\Descriptives_plots
# and
# S:\MNHS-Psych\ANDL-Lab-TAG-Study\SOS Study\saliva\Biodata_scripts_output

## For each analyte, calculate each person's peak from T3 OR T4
# We do this because of the inconsistencies in the literature regarding when
# peak should happen - could be different by analyte and due to individual differences

library(panelr)
sos_crp_panel <- panel_data(sos_crp[ , c("SampleID", "Time","CRP_log")], 
                            id = "SampleID", wave = "Time")
sos_crp_wide <- widen_panel(sos_crp_panel)
sos_crp_wide$CRPlog_t3t4peak <- ifelse(sos_crp_wide$CRP_log_3 > sos_crp_wide$CRP_log_4, 
                                       sos_crp_wide$CRP_log_3, sos_crp_wide$CRP_log_4)

sos_cytokines_panel <- panel_data(sos_cytokines[ , c("SampleID", "Time","IL10_log","IL6_log","TNFalpha_log")], 
                                  id = "SampleID", wave = "Time")
sos_cytokines_wide <- widen_panel(sos_cytokines_panel)
sos_cytokines_wide$IL10log_t3t4peak <- ifelse(sos_cytokines_wide$IL10_log_3 > sos_cytokines_wide$IL10_log_4, 
                                              sos_cytokines_wide$IL10_log_3, sos_cytokines_wide$IL10_log_4)
sos_cytokines_wide$IL6log_t3t4peak <- ifelse(sos_cytokines_wide$IL6_log_3 > sos_cytokines_wide$IL6_log_4, 
                                              sos_cytokines_wide$IL6_log_3, sos_cytokines_wide$IL6_log_4)
sos_cytokines_wide$TNFalog_t3t4peak <- ifelse(sos_cytokines_wide$TNFalpha_log_3 > sos_cytokines_wide$TNFalpha_log_4, 
                                              sos_cytokines_wide$TNFalpha_log_3, sos_cytokines_wide$TNFalpha_log_4)

sos_hormones_panel <- panel_data(sos_hormones[ , c("SampleID", "Time","DHEA_log","Testosterone_log","Cortisol_log")], 
                                  id = "SampleID", wave = "Time")
sos_hormones_wide <- widen_panel(sos_hormones_panel)
sos_hormones_wide$DHEAlog_t3t4peak <- ifelse(sos_hormones_wide$DHEA_log_3 > sos_hormones_wide$DHEA_log_4, 
                                              sos_hormones_wide$DHEA_log_3, sos_hormones_wide$DHEA_log_4)
sos_hormones_wide$Testlog_t3t4peak <- ifelse(sos_hormones_wide$Testosterone_log_3 > sos_hormones_wide$Testosterone_log_4, 
                                             sos_hormones_wide$Testosterone_log_3, sos_hormones_wide$Testosterone_log_4)
sos_hormones_wide$Cortlog_t3t4peak <- ifelse(sos_hormones_wide$Cortisol_log_3 > sos_hormones_wide$Cortisol_log_4, 
                                              sos_hormones_wide$Cortisol_log_3, sos_hormones_wide$Cortisol_log_4)

## Check if that "peak" is significantly higher than T1 levels
mean(sos_crp_wide$CRPlog_t3t4peak, na.rm = TRUE) # -2.052981
mean(sos_crp_wide$CRP_log_1, na.rm = TRUE) # -2.543955
t.test(sos_crp_wide$CRPlog_t3t4peak, sos_crp_wide$CRP_log_1, alternative = "greater", paired = TRUE)

mean(sos_cytokines_wide$IL10log_t3t4peak, na.rm = TRUE) # 1.626901
mean(sos_cytokines_wide$IL10_log_1, na.rm = TRUE) # 1.089902
t.test(sos_cytokines_wide$IL10log_t3t4peak, sos_cytokines_wide$IL10_log_1, alternative = "greater", paired = TRUE)

mean(sos_cytokines_wide$IL6log_t3t4peak, na.rm = TRUE) # 1.01348
mean(sos_cytokines_wide$IL6_log_1, na.rm = TRUE) # 0.7275804
t.test(sos_cytokines_wide$IL6log_t3t4peak, sos_cytokines_wide$IL6_log_1, alternative = "greater", paired = TRUE)

mean(sos_cytokines_wide$TNFalog_t3t4peak, na.rm = TRUE) # 2.152535
mean(sos_cytokines_wide$TNFalpha_log_1, na.rm = TRUE) # 1.796177
t.test(sos_cytokines_wide$TNFalog_t3t4peak, sos_cytokines_wide$TNFalpha_log_1, alternative = "greater", paired = TRUE)

mean(sos_hormones_wide$DHEAlog_t3t4peak, na.rm = TRUE) # 4.621193
mean(sos_hormones_wide$DHEA_log_1, na.rm = TRUE) # 4.340595
t.test(sos_hormones_wide$DHEAlog_t3t4peak, sos_hormones_wide$DHEA_log_1, alternative = "greater", paired = TRUE)

mean(sos_hormones_wide$Testlog_t3t4peak, na.rm = TRUE) # 3.741303
mean(sos_hormones_wide$Testosterone_log_1, na.rm = TRUE) # 3.684813
t.test(sos_hormones_wide$Testlog_t3t4peak, sos_hormones_wide$Testosterone_log_1, alternative = "greater", paired = TRUE)

mean(sos_hormones_wide$Cortlog_t3t4peak, na.rm = TRUE) # -1.904597
mean(sos_hormones_wide$Cortisol_log_1, na.rm = TRUE) # -1.785707
t.test(sos_hormones_wide$Cortlog_t3t4peak, sos_hormones_wide$Cortisol_log_1, alternative = "greater", paired = TRUE)

## Acknowledge that x ppl (by analyte) had peak at T1, so they will have lower
# levels at T3/T4 than T1 - see if those people have lower subjective stress ratings
# This may be especially important for testosterone and cortisol

# Plot each analyte (and separate by immune markers vs. hormones)

# Some people had peak at T2 which is due to scanner environment, but T3/T4 
# may still be higher than T1 for those people, so it's ok. But may want to
# check significant differences between T2 and T1 for everyone and report

# May also want to report on fit of a quad model from T1 - T4.

# Then calculate each person's "reactivity" variable for use in further analyses
# This could be a simple difference score from T1 to T3/T4, or regress T1 on outcomes
# Could also calculate area under the curve