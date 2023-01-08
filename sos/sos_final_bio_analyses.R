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


## for each analyte, calculate each person's peak from T3 OR T4

## Check if that "peak" is significantly higher than T1 levels

## Acknowledge that x ppl (by analyte) had peak at T1, so they will have lower
# levels at T3/T4 than T1 - see if those people have lower subjective stress ratings

# Some people had peak at T2 which is due to scanner environment, but T3/T4 
# may still be higher than T1 for those people, so it's ok. But may want to
# check significant differences between T2 and T1 for everyone and report

# May also want to report on fit of a quad model from T1 - T4.

# Then calculate each person's "reactivity" variable for use in further analyses
# This could be a simple difference score from T1 to T3/T4, or regress T1 on outcomes
# Could also calculate area under the curve