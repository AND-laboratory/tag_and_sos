################################################# 
####### SOS (TAG) bio data final analyses ####### 
################################################# 

library(utils)
library(tidyverse)
library(panelr)
library(readxl)
library(ggplot2)
library(viridis)
library(ggpubr)
library(rstatix)

workdir <- 'S:/MNHS-Psych/ANDL-Lab-TAG-Study/SOS Study'
HR_Resting1 <- read.csv(file.path(workdir,"/physio/Mean_HR_Resting1.csv", fsep=""))
HR_Resting2 <- read.csv(file.path(workdir,"/physio/Mean_HR_Resting2.csv", fsep=""))
HR_Task <- read.csv(file.path(workdir,"/physio/Mean_HR_Task.csv", fsep=""))

HR_conditions_raw <- read_excel(file.path(workdir,"/physio/Heart_Rate_Data_Task_Sequence_All_Slices_Conditions.xlsx", fsep=""))

sos_crp <- read.csv(file.path(workdir,"/saliva/Cleaned_datasets/SOS_CRP_cleaned.csv", fsep="")) 
sos_cytokines <- read.csv(file.path(workdir,"/saliva/Cleaned_datasets/SOS_cytokines_cleaned.csv", fsep=""))
sos_hormones <- read.csv(file.path(workdir,"/saliva/Cleaned_datasets/SOS_hormones_cleaned.csv", fsep=""))

####### Heart rate ####### 

colnames(HR_Resting1)[1] <- "ID"
colnames(HR_Resting2)[1] <- "ID"
colnames(HR_Task)[1] <- "ID"


# get rid of the repeat 313
HR_Resting1 <- filter(HR_Resting1, ID!="TAG313p2_20191117_151413")
HR_Resting2 <- filter(HR_Resting2, ID!="TAG313p2_20191117_151413")
HR_Task <- filter(HR_Task, ID!="TAG313p2_20191117_151413")

# There are two subjects (077 and 188) who had very low task HRs (high 30s) 
# despite having normal range resting state HRs, so they need to be removed
# as this was likely artifact during the task.

HR_Resting1 <- filter(HR_Resting1, ID!="TAG077_20190823_145422")
HR_Resting2 <- filter(HR_Resting2, ID!="TAG077_20190823_145422")
HR_Task <- filter(HR_Task, ID!="TAG077_20190823_145422")

HR_Resting1 <- filter(HR_Resting1, ID!="TAG188_20190731_145115")
HR_Resting2 <- filter(HR_Resting2, ID!="TAG188_20190731_145115")
HR_Task <- filter(HR_Task, ID!="TAG188_20190731_145115")

HR_Resting1$sequence <- rep(c("Resting 1"), each=50)
HR_Resting2$sequence <- rep(c("Resting 2"), each=49)
HR_Task$sequence <- rep(c("Task"), each=50)

HR_all_long <- bind_rows(HR_Resting1, HR_Resting2, HR_Task)

# average both resting states (for ID 080, use resting 1 as average because
# they don't have resting 2)

HR_all_wide <- spread(HR_all_long, sequence, Mean_HR)
HR_all_wide$RestingAvg <- ifelse(is.na(HR_all_wide$'Resting 2'), 
                                      HR_all_wide$'Resting 1', 
                                      rowMeans(HR_all_wide[,c('Resting 1', 'Resting 2')]))

# Add the individual task conditions:
HR_conditions <- pivot_longer(HR_conditions_raw, SOS068_20190303_150846:TAG223_20191116_145342, names_to = "ID", values_to = "HR")

HR_conditions$Task_condition <- as.factor(HR_conditions$Task_condition)

HR_conditions$ID <- as.factor(HR_conditions$ID)

HR_conditions_avg <- HR_conditions %>% group_by(ID, Task_condition) %>% summarize(mean_HR = mean(HR))

HR_conditions_clean <- HR_conditions_avg %>% filter(Task_condition!="NA")
HR_conditions_clean$Task_condition <- factor(HR_conditions_clean$Task_condition)
HR_conditions_clean_tbl <- as_tibble(HR_conditions_clean)
# Remove two outliers - their task HR data was impossibly low (probably artifact)
HR_conditions_clean_tbl <- filter(HR_conditions_clean_tbl, ID!="TAG077_20190823_145422")
HR_conditions_clean_tbl <- filter(HR_conditions_clean_tbl, ID!="TAG188_20190731_145115")

# Averages by condition
meanhr_resting1 <- mean(HR_all_wide$'Resting 1') # 73.439
meanhr_resting2 <- mean(HR_all_wide$'Resting 2', na.rm = TRUE) # 72.371
meanhr_restingavg <- mean(HR_all_wide$RestingAvg) # 72.927
meanhr_task <- mean(HR_all_wide$Task) # 76.942
meanhr_R <- mean(HR_conditions_clean_tbl[which(HR_conditions_clean_tbl$Task_condition=="R"), ]$mean_HR) # 76.913
meanhr_other <- mean(HR_conditions_clean_tbl[which(HR_conditions_clean_tbl$Task_condition=="1"), ]$mean_HR) # 75.981
meanhr_selfdis <- mean(HR_conditions_clean_tbl[which(HR_conditions_clean_tbl$Task_condition=="2"), ]$mean_HR) # 77.095
meanhr_selfcon <- mean(HR_conditions_clean_tbl[which(HR_conditions_clean_tbl$Task_condition=="3"), ]$mean_HR) # 77.493


# Check significant differences (Question: should this be paired t-tests or repeated measures ANOVAs?)

# paired t-tests
# Average all task vs. average resting state 
t.test(HR_all_wide$Task, HR_all_wide$RestingAvg, alternative = "greater", 
       paired = TRUE) # t = 7.4707, df = 49, p-value = 6.21e-10

# Self-connected vs. R condition (or should it be general resting state?) 
t.test(HR_conditions_clean_tbl[which(HR_conditions_clean_tbl$Task_condition=="3"), ]$mean_HR, 
       HR_conditions_clean_tbl[which(HR_conditions_clean_tbl$Task_condition=="R"), ]$mean_HR,
       alternative = "greater", paired = TRUE) # t = 2.9771, df = 50, p-value = 0.002238

# Self-connected vs. self-disconnected
t.test(HR_conditions_clean_tbl[which(HR_conditions_clean_tbl$Task_condition=="3"), ]$mean_HR, 
       HR_conditions_clean_tbl[which(HR_conditions_clean_tbl$Task_condition=="2"), ]$mean_HR,
       alternative = "greater", paired = TRUE) # t = 2.3487, df = 50, p-value = 0.01142

# Self-connected vs. Other
t.test(HR_conditions_clean_tbl[which(HR_conditions_clean_tbl$Task_condition=="3"), ]$mean_HR, 
       HR_conditions_clean_tbl[which(HR_conditions_clean_tbl$Task_condition=="1"), ]$mean_HR,
       alternative = "greater", paired = TRUE) # t = 6.2889, df = 50, p-value = 3.883e-08

# Plotting Heart rate

# Task vs. resting state: HR_all_wide
HR_avg_long <- HR_all_wide[,c("ID","Task","RestingAvg")]
HR_avg_long <- pivot_longer(HR_avg_long, cols = c("Task","RestingAvg"), 
                            names_to = "condition", 
                            values_to = "mean_HR")
HR_avg_long$condition <- factor(HR_avg_long$condition,
                                                 levels = c("Task","RestingAvg"),
                                                 labels = c("SOS Task","Resting State Scan"))

compare_means(mean_HR ~ condition, data = HR_avg_long, method = "t.test", paired = TRUE)

my_comparisons1 <- list( c("SOS Task","Resting State Scan"))

get_box_stats <- function(y, upper_limit = max(HR_avg_long$mean_HR) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Mean =", round(mean(y), 2), "\n",
      "Median =", round(median(y), 2), "\n"
    )
  ))
}

avgHR_plot <- 
  ggplot(HR_avg_long, aes(x = condition, y = mean_HR)) +
  geom_violin(aes(fill = condition), trim = FALSE) +
  geom_boxplot(width = 0.2) +
  scale_fill_viridis_d() +
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 9) +
  stat_compare_means(method = "t.test", paired = TRUE, 
                     comparisons = my_comparisons1, 
                     method.args = list(alternative = "greater")) + # one tailed based on my_comparisons order above
  ylab("Mean Heart Rate (bpm)") +
  xlab("Condition") +
  theme(legend.position = "none") 

# Within task: HR_conditions_clean_tbl
HR_conditions_clean_tbl$Task_condition <- factor(HR_conditions_clean_tbl$Task_condition,
                                                  levels = c("1","2","3","R"),
                                                  labels = c("Other", 
                                                             "Self-Disconnected",
                                                             "Self-Connected",
                                                             "Rest"))

compare_means(mean_HR ~ Task_condition, data = HR_conditions_clean_tbl, method = "t.test", paired = TRUE)

my_comparisons <- list( c("Self-Connected","Rest"),
                        c("Self-Connected","Self-Disconnected"),
                        c("Self-Connected","Other"),
                        c("Self-Disconnected","Other"),
                        c("Self-Disconnected","Rest"),
                        c("Other","Rest"))

get_box_stats <- function(y, upper_limit = max(HR_conditions_clean_tbl$mean_HR) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Mean =", round(mean(y), 2), "\n",
      "Median =", round(median(y), 2), "\n"
    )
  ))
}

taskHR_plot <- 
ggplot(HR_conditions_clean_tbl, aes(x = Task_condition, y = mean_HR)) +
  geom_violin(aes(fill = Task_condition), trim = FALSE) +
  geom_boxplot(width = 0.2) +
  scale_fill_viridis_d() +
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 7.2) +
  stat_compare_means(method = "t.test", paired = TRUE, 
                     comparisons = my_comparisons, 
                     method.args = list(alternative = "greater")) + # one tailed based on my_comparisons order above
  ylab("Mean Heart Rate (bpm)") +
  xlab("Condition") +
  theme(legend.position = "none") 



####### Salivary markers ####### 

# Report some basic plots taken from:
# S:\MNHS-Psych\ANDL-Lab-TAG-Study\SOS Study\saliva\Descriptives_plots
# and
# S:\MNHS-Psych\ANDL-Lab-TAG-Study\SOS Study\saliva\Biodata_scripts_output

# Also show plots for whole group with repeated measured (with standard errors)
# and growth curve modeling with quadratic curve (if it works)

## For each analyte, calculate each person's peak from T3 OR T4
# We do this because of the inconsistencies in the literature regarding when
# peak should happen - could be different by analyte and due to individual 
# differences. Doing this (picking T3 or T4) is basically a simpler version
# of that they do here (I think): https://pubmed.ncbi.nlm.nih.gov/24754834/
# Basically, what they did (they had heaps of time points) was manually select
# each person's peak, then they did a time-adjusted two-piece growth curve 
# model where the intercept was the person's peak. (They call this "landmark
# registration"*). Then the first part of the model was the
# activation slope (pre-peak slope) and the second part was the "regulation"
# slope (post-peak). I probably don't have enough time points for this to break
# into two parts. I think it's better to just select if T3 or T4 are that 
# person's peak. T2 was not designed to be the peak, so maybe for that we 
# just acknowledge some people had highest peak there due to MRI stress.
# Alternatively, could we include the T2 value as a covariate in further 
# analyses?
#
# * from the paper: "When alignment is based on a specific feature of the 
# curve, such as a response peak, the process is called landmark 
# registration (Kneip & Gasser, 1992)."

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
# OR - should these be repeated measures ANOVAs to do all time points at once?

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

# New long datasets with just T1 and the T3/T4 peak in it
sos_crp_long <- sos_crp_wide[,c("SampleID","CRP_log_1","CRPlog_t3t4peak")]
sos_crp_long <- pivot_longer(sos_crp_long, cols = c("CRP_log_1","CRPlog_t3t4peak"), 
                                            names_to = "time_point", 
                                            values_to = "CRP_log_conc")

sos_crp_long$time_point <- factor(sos_crp_long$time_point,
                                                 levels = c("CRP_log_1","CRPlog_t3t4peak"),
                                                 labels = c("T1","T3/T4 peak"))

compare_means(CRP_log_conc ~ time_point, data = sos_crp_long, method = "t.test", paired = TRUE)

my_comparisons_crp <- list( c("T3/T4 peak","T1"))

get_box_stats_crp <- function(y, upper_limit = max(sos_crp_long$CRP_log_conc) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Mean =", round(mean(y), 2), "\n",
      "Median =", round(median(y), 2), "\n"
    )
  ))
}

CRP_plot <- 
  ggplot(sos_crp_long, aes(x = time_point, y = CRP_log_conc)) +
  geom_violin(aes(fill = time_point), trim = FALSE) +
  geom_boxplot(width = 0.2) +
  scale_fill_viridis_d() +
  stat_summary(fun.data = get_box_stats_crp, geom = "text", hjust = 0.5, vjust = 7.2, na.rm = TRUE) + #printing the mean isn't working, maybe because there are some NAs
  stat_compare_means(method = "t.test", paired = TRUE, 
                     comparisons = my_comparisons_crp, 
                     method.args = list(alternative = "greater")) + # one tailed based on my_comparisons_crp order above
  ylab("Mean log(CRP)") +
  xlab("Time point") +
  theme(legend.position = "none") 

## Acknowledge that x ppl (by analyte) had peak at T1, so they will have lower
# levels at T3/T4 than T1 - see if those people have lower subjective stress ratings
# This may be especially important for testosterone and cortisol

# Some people had peak at T2 which is due to scanner environment, but T3/T4 
# may still be higher than T1 for those people, so it's ok. But may want to
# check significant differences between T2 and T1 for everyone and report

# Then calculate each person's "reactivity" variable for use in further analyses
# This could be a simple difference score from T1 to T3/T4, or regress T1 on outcomes
# Could also calculate area under the curve, or slope from growth curve, but
# growth curve would have to have intercept set at the individual's peak 
# (landmark registration) and the slope would be the slope before that intercept
# See: https://pubmed.ncbi.nlm.nih.gov/24754834/