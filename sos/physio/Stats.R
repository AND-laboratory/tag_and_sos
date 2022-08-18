setwd('/Users/clairekelly/Documents/Work/Michelle/Run_PhysIO_All/Stats/')

library(tidyverse)

# Import/organise/clean data

HR_Resting1 <- read_csv('Mean_HR_Resting1.csv')
HR_Resting2 <- read_csv('Mean_HR_Resting2.csv')
HR_Task <- read_csv('Mean_HR_Task.csv')

colnames(HR_Resting1)[1] <- "ID"
colnames(HR_Resting2)[1] <- "ID"
colnames(HR_Task)[1] <- "ID"

# We couldn't get HR data for Resting2 for TAG 080. To make things easier, delete that subject across all sequences:
HR_Resting1 <- filter(HR_Resting1, ID!="TAG080_20190329_150631")
HR_Task <- filter(HR_Task, ID!="TAG080_20190329_150631")

HR_Resting1$sequence <- rep(c("Resting 1"), each=52)
HR_Resting2$sequence <- rep(c("Resting 2"), each=52)
HR_Task$sequence <- rep(c("Task"), each=52)

HR_all_sequences <- bind_rows(HR_Resting1, HR_Resting2, HR_Task)


# plotting 

library(ggplot2)

p <- HR_all_sequences %>% 
  group_by(sequence) %>% 
  summarize(m = mean(Mean_HR),
            s = sd(Mean_HR)) %>% 
  ggplot(aes(x = sequence, y = m)) + 
  geom_point(data=HR_all_sequences, aes(x=sequence,y=Mean_HR), alpha=0.5, size=2, fill='black', color='black')+
  geom_bar(stat = "identity", width=1, color='black', fill='#8A3324', alpha=0.5, lwd=0.2)+
  geom_errorbar(aes(ymin = m-s, ymax = m+s), width=0.2) +
  ylab("Mean Heart Rate")+xlab("")+
  theme_classic(base_size=16) 

ggsave('plot.png', p, width=120, height=100, units=c("mm"), dpi=300)


# check means

mean(HR_Resting1$Mean_HR) #73.20782
mean(HR_Resting2$Mean_HR) #72.10011
mean(HR_Task$Mean_HR) #75.5534


# paired t-tests

### resting1 vs resting2
HR_resting_only <- filter(HR_all_sequences, sequence!="Task")
t.test(Mean_HR ~ sequence, data = HR_resting_only, paired = TRUE)
#t = 3.8657, df = 51, p-value = 0.0003148


### resting2 vs task
HR_resting2_and_task <- filter(HR_all_sequences, sequence!="Resting 1")
t.test(Mean_HR ~ sequence, data = HR_resting2_and_task, paired = TRUE)
#t = -3.3042, df = 51, p-value = 0.001747

### resting1 vs task
HR_resting1_and_task <- filter(HR_all_sequences, sequence!="Resting 2")
t.test(Mean_HR ~ sequence, data = HR_resting1_and_task, paired = TRUE)
#t = -2.2252, df = 51, p-value = 0.03052





