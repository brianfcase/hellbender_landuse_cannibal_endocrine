setwd("~/R_parent")

library(survival)
library(lmerTest)
library(lme4)
library(ggplot2)
library(ggeffects)
library(rgl)
library(car)
library(MASS)
library(plot3D)
library(dplyr)
library(plotly)
library(NADA)
library(NADA2)
library(psych)
library(performance)
library(emmeans)
library(DHARMa)
library(AICcmodavg)
library(MuMIn)

parent <- read.csv("lcms_landuse_final.csv")

#fixing dataset categories -----

parent$s_f <- as.factor(parent$s_f)
parent$run_year <- as.factor(parent$run_year)
parent$forest_binary <- as.factor (parent$forest_binary)
parent$outcome <- as.factor(parent$outcome)
parent$period <- as.factor(parent$period)
parent$type <- as.factor(parent$type)
parent$year <- as.factor(parent$year)
parent$pit <- as.factor(parent$pit)
parent$pit_letter <- as.factor(parent$pit_letter)
parent$site <- as.factor(parent$site)
parent$leeches_binomial <- as.factor(parent$leeches_binomial)
parent$testosterone <- as.numeric(parent$testosterone)
parent$run <- as.factor(parent$run)
parent$river <- as.factor(parent$river)
parent$cortisol_nada <- as.factor(parent$cortisol_nada)
parent$cortisol_limit <- as.factor(parent$cortisol_limit)
parent$cortisol_detected <- as.factor(parent$cortisol_detected)

#subsetting dataset -----

parent_t0 <- subset(parent, type =="t0")
parent_t60 <- subset(parent, type =="t60")

init_t0 <- subset(parent_t0, period =="initiation")
init_t60 <- subset(parent_t60, period =="initiation")
one_month_t0 <- subset(parent_t0, period =="one month")
one_month_t60 <- subset(parent_t60, period =="one month")
one_month_t60_big <- subset(one_month_t60, tol >=38)
subset_data <- function(data) {
  # Ensure that the outcome column is treated as a factor or character
  data <- data %>% mutate(outcome = as.character(outcome))
  
  # Subset the data to remove points with forest cover <60 and outcome 'Success'
  subset_data <- data %>% filter(!(forest <57 & outcome == 'success'))
  
  return(subset_data)
}

filtered_data <- subset_data(one_month_t60)

#Method Statistics----

#Obtaining Samples-----
#(1) Days elapsed from detection to sampling, by sampling period
describeBy(x = parent_t0$dect_days, list(parent_t0$period))

#(2) Time elapsed from capture to sampling, T0
describeBy(x = parent$sample_elapsed_min, list(parent$type))

#(3) Time of day when samples collected
describeBy(x = parent$sample_time_hours, list(parent$type))
#min = 9.62 hours (0937) and max = 17.79 hours (1747)

#Table 1 - Blood Plasma Sample Sizes by Period, Nest Fate, Sample Type -----

parentdt <- data.table::as.data.table(parent)
parentdt[, .N, by = c("type","period")] #number of parental t0s and t60s by period
parentdt[, .N, by = c("type","period","outcome")] #number of parental t0s and t60s by period and outcome

#Table 2 - Hormone Detection Counts & Freq by Period, and Overall -----

#testosterone

parent_t0 %>% #by period
  group_by (period) %>%
  summarise(dect_num = sum(testosterone_nada == "FALSE"), total = n()) %>%
  mutate(dect_per = (dect_num/total*100))

parent_t0 %>% #overall
  summarise(dect_num = sum(testosterone_nada == "FALSE"), total = n()) %>%
  mutate(dect_per = (dect_num/total*100))

#dht

parent_t0 %>%
  group_by (period) %>%
  summarise(dect_num = sum(dht_nada == "FALSE"), total = n()) %>%
  mutate(dect_per = (dect_num/total*100))

parent_t0 %>% #overall
  summarise(dect_num = sum(dht_nada == "FALSE"), total = n()) %>%
  mutate(dect_per = (dect_num/total*100))

#cortisol

parent %>%
  group_by (period, type) %>% #by period
  summarise(dect_num = sum(cortisol_nada == "FALSE"), total = n()) %>%
  mutate(dect_per = (dect_num/total*100))

parent %>%
  group_by (type) %>% #overall
  summarise(dect_num = sum(cortisol_nada == "FALSE"), total = n()) %>%
  mutate(dect_per = (dect_num/total*100))

#corticosterone

parent %>%
  group_by (period, type) %>% #by period
  summarise(dect_num = sum(corticosterone_nada == "FALSE"), total = n()) %>%
  mutate(dect_per = (dect_num/total*100))

parent %>% 
  group_by (type) %>% #overall
  summarise(dect_num = sum(corticosterone_nada == "FALSE"), total = n()) %>%
  mutate(dect_per = (dect_num/total*100))


#2.4.2 - Potential Covariates of Hormone Concentrations-----

#(1) androgens and sample time of day

#testosterone and sample time of day, ATS regression
with(parent_t0, ATS(testosterone, testosterone_nada,sample_time_hours))

#dht and sample time of day, ATS regression
with(parent_t0, ATS(dht, dht_nada,sample_time_hours))

#(2) leeches and T60 cortisol, forest cover used as proxy for known variation 
#in leeches across watersheds

#survival terms
t60_cortisol_surv_init <- Surv(init_t60$cortisol, init_t60$cortisol_limit == 0, type = "left")
t60_cortisol_surv_om <- Surv(one_month_t60$cortisol, one_month_t60$cortisol_limit == 0, type = "left")

#stress cortisol and leech infection (leech magnitude)

#at initiation

cortisol_leech_mag_init <- survreg(t60_cortisol_surv_init ~ leeches * forest, 
                                   cluster = pit, dist = 'lognormal', data = init_t60)
summary(cortisol_leech_mag_init)

#at one month

cortisol_leech_mag_om <- survreg(t60_cortisol_surv_om ~ leeches * forest , 
                                 cluster = pit, dist = 'lognormal', data = one_month_t60)
summary(cortisol_leech_mag_om)

#(3) androgens, baseline cortisol detection (1/0), and t60 cortisol vs. SMI at nest initiation

with(init_t0, ATS(testosterone, testosterone_nada,smi))

with(init_t0, ATS(dht, dht_nada,smi))

t0_smi_init <- glmer(cortisol_detected ~ smi 
                     + (1|pit), data = init_t0, family = binomial)
summary(t0_smi_init)

with(init_t60, ATS(cortisol, cortisol_nada,smi))

#(4) androgens (initiation), baseline cortisol detection (1/0), 
#and t60 cortisol (initiation & one month) vs. total length

with(init_t0, ATS(testosterone, testosterone_nada,tol))

with(init_t0, ATS(dht, dht_nada,tol))

t0_tol_init <- glmer(cortisol_detected ~ tol 
                     + (1|pit), data = init_t0, family = binomial)
summary(t0_tol_init)

t0_tol_om <- glmer(cortisol_detected ~ tol 
                   + (1|pit), data = one_month_t0, family = binomial)
summary(t0_tol_om)

with(init_t60, ATS(cortisol, cortisol_nada,tol))

with(one_month_t60, ATS(cortisol, cortisol_nada,tol))

#__________________________________________________________________------
#Results Statistics ----
#Survey Summary ----

#(1) number of nests
length(unique(parent$nest))

#(2) number of unique males
length(unique(parent$pit))

#(3) number of t0 samples
nrow(parent_t0)

#(4) number of t60 samples 
nrow(parent_t60)

#(5) nest sample sizes & unique males by by nest outcome
parentdtsuccess <- subset(parentdt, outcome =="success")
parentdtsuccess[, .N, by = c("nest")] #number of nests successful
parentdtsuccess[, .N, by = c("pit")] #number of unique males successful

parentdtcannibal <- subset(parentdt, outcome =="cannibal")
parentdtcannibal[, .N, by = c("nest")] #number of nests cannibals
parentdtcannibal[, .N, by = c("pit")] #number of unique males cannibals

#(6) males that succeeded and cannibalized over course of study

# unique male IDs for each outcome
success_males <- unique(parentdtsuccess$pit)
cannibal_males <- unique(parentdtcannibal$pit)

# males that did both
both_males <- intersect(success_males, cannibal_males)

# count of unique males that both succeeded and cannibalized
length(both_males) #6 individuals

#(7) correlation between nest fate and forest cover

unique_values <- unique(parent$nest)

corresponding_values <- sapply(unique_values, function(x) parent$s_f[parent$nest == x][1])
forest <- sapply(unique_values, function(x) parent$forest[parent$nest == x][1])
pit <- sapply(unique_values, function(x) parent$pit[parent$nest == x][1])

result_df <- data.frame(column1 = unique_values, column2 = corresponding_values, column3 = forest, column4 = pit)
result_df$nest <- result_df$column1
result_df$s_f <- result_df$column2
result_df$forest <- result_df$column3
result_df$pit <- result_df$column4

result_df$forest <- scale(result_df$forest)

forest_can <- glmer(s_f ~ forest + (1|pit), data = result_df, family = binomial)
summary(forest_can)

#3.1 - Morphometrics of Nesting Males ----

#(1) total length at nest initiation
describeBy(x = init_t0$tol)

#(2) mass at nest initiation
describeBy(x = init_t0$mass)

#(3) summary statistics nesting males SMI
describeBy(x = init_t0$smi)

#(4) number of sampling events where SMI fell below 
#population's lower 10th percentile (SMI: 357)
sum(parent_t0$smi < 357, na.rm = TRUE)
9/163

#(5) body condition by nest outcome

hist(init_t0$smi)

#scale predictor variables
init_t0$smi_scaled <- scale(init_t0$smi)
init_t0$forest_scaled <- scale(init_t0$forest)

#glmer with s_f response (binomial) and smi and forest as predictor (scaled, eliminates errors)
smi_init_scaled <- glmer(s_f ~ smi_scaled * forest_scaled + (1|pit), data = init_t0, family = binomial)
summary(smi_init_scaled)

#smi not related to nest outcome, and relationship 
#between SMI and nest outcome does not depend on forest cover

#(6) total length by nest outcome

hist(init_t0$tol)

#scale predictor variables
init_t0$tol_scaled <- scale(init_t0$tol)
init_t0$forest_scaled <- scale(init_t0$forest)

#glmer with s_f response (binomial) and tol as predictor
tol_init <- glmer(s_f ~ tol_scaled * forest_scaled + (1|pit), data = init_t0, family = binomial)
summary(tol_init)

#3.2 - Steroids Across Nest Fate and Forest Cover----
###section (1) androgens------

#maximum successful testosterone
max(init_t0$testosterone[init_t0$s_f == 1], na.rm = TRUE)

#maximum cannibal testosterone
max(init_t0$testosterone[init_t0$s_f == 0], na.rm = TRUE)

#maximum successful dht
max(init_t0$dht[init_t0$s_f == 1], na.rm = TRUE)

#maximum successful dht
max(init_t0$dht[init_t0$s_f == 0], na.rm = TRUE)

#testosterone and days elapsed from nest detection
with(init_t0, ATS(testosterone, testosterone_nada,dect_days))

#dht and days elapsed from nest detection
with(init_t0, ATS(dht, dht_nada,dect_days))

#correlation between testosterone and DHT

#pulling samples where both compounds detected
t_dht_ratio <- init_t0[(init_t0$testosterone > 0.2) & (init_t0$dht > 1.0),]

#plotting above selection, dht by testosterone
plot(t_dht_ratio$dht,t_dht_ratio$testosterone)

#linear model of dht and testosterone
t_dht_ratio_lm <- lm(testosterone ~ dht, data = t_dht_ratio)
summary(t_dht_ratio_lm)

#testosterone models at initiation and AICC selection / comparison

t_1<- lmer(log(testosterone) ~ outcome*forest+ (1|pit), data = init_t0, REML = FALSE)
summary(t_1)

t_2 <- lmer(log(testosterone) ~ outcome+forest + (1|pit), data = init_t0, REML = FALSE)
summary(t_2)

t_3<- lmer(log(testosterone) ~ forest + (1|pit), data = init_t0, REML = FALSE)
summary(t_3)

t_4<- lmer(log(testosterone) ~ outcome + (1|pit), data = init_t0, REML = FALSE)
summary(t_4)

t_null<- lmer(log(testosterone) ~ (1|pit), data = init_t0, REML = FALSE)
summary(t_null)

testosterone_models <- list(
  t_1 = t_1,
  t_2 = t_2,
  t_3 = t_3,
  t_4 = t_4,
  t_null = t_null
)

t_aicc_table <- model.sel(testosterone_models, rank = AICc)

# Convert to data frame
t_aicc_table <- as.data.frame(t_aicc_table)

# Extract K (number of estimated parameters = df column in model.sel output)
t_aicc_table$K <- t_aicc_table$df

# Add cumulative Akaike weights
t_aicc_table$cum_weight <- cumsum(t_aicc_table$weight)

# Reorder columns to show K and cumulative weight clearly
t_aicc_table <- t_aicc_table[, c("df", "K", "logLik", "AICc", "delta", "weight", "cum_weight")]

# View testosterone AICc table
print(t_aicc_table)

#chi squared test between highest performing model and those including
#forest cover with outcome
anova(t_4, t_2, test = "Chisq")
anova(t_4, t_1, test = "Chisq")

#dht models at initiation and AICC selection / comparision

dht_1<- lmer(log(dht) ~ outcome*forest+ (1|pit), data = init_t0, REML = FALSE)
summary(dht_1)

dht_2 <- lmer(log(dht) ~ outcome+forest + (1|pit), data = init_t0, REML = FALSE)
summary(dht_2)

dht_3<- lmer(log(dht) ~ forest + (1|pit), data = init_t0, REML = FALSE)
summary(dht_3)

dht_4 <- lmer(log(dht) ~ outcome + (1|pit), data = init_t0, REML = FALSE)
summary(dht_4)

dht_null<- lmer(log(dht) ~ (1|pit), data = init_t0, REML = FALSE)
summary(dht_null)

dht_models <- list(
  dht_1 = dht_1,
  dht_2 = dht_2,
  dht_3 = dht_3,
  dht_4 = dht_4,
  dht_null = dht_null
)

dht_aicc_table <- model.sel(dht_models, rank = AICc)

# Convert to data frame
dht_aicc_table <- as.data.frame(dht_aicc_table)

# Extract K (number of estimated parameters = df column in model.sel output)
dht_aicc_table$K <- dht_aicc_table$df

# Add cumulative Akaike weights
dht_aicc_table$cum_weight <- cumsum(dht_aicc_table$weight)

# Reorder columns to show K and cumulative weight clearly
dht_aicc_table <- dht_aicc_table[, c("df", "K", "logLik", "AICc", "delta", "weight", "cum_weight")]

# View the DHT AICc table
print(dht_aicc_table)

#chi squared test between highest performing model and those including
#forest cover with outcome
anova(dht_4, dht_2, test = "Chisq")
anova(dht_4, dht_1, test = "Chisq")

#comparison between mean testosterone at initiation, successful vs. cannibal males
describeBy(x = init_t0$testosterone, group = init_t0$outcome)
21.5 - 14.0
#(successful males mean T = 21.5) - (cannibal males' mean T = 14.0) = 7.5
7.5/14
#increase from cannibals to successful males = 7.5/14.0, or 53.6% 

##comparison between mean dht at initiation, successful vs. cannibal males
describeBy(x = init_t0$dht, group = init_t0$outcome)
9.4 - 6.4
#(successful males mean DHT = 9.4) - (cannibal males' mean DHT = 6.4) = 3
3/6.4
#increase from cannibals to successful males = 3/6.4, or 47%

###section (2) baseline cortisol-----

#initiation
max(init_t0$cortisol[init_t0$s_f == 1], na.rm = TRUE)

sum(init_t0$cortisol[init_t0$s_f == 1] > 0.8) #8 above 0.8
sum(init_t0$s_f == 1) #44 successful samples
mean(init_t0$cortisol[init_t0$s_f == 1] > 0.8) * 100

max(init_t0$cortisol[init_t0$s_f == 0], na.rm = TRUE)
sum(init_t0$cortisol[init_t0$s_f == 0] > 0.8) #3 above 0.8
sum(init_t0$s_f == 0) #38 cannibal samples
mean(init_t0$cortisol[init_t0$s_f == 0] > 0.8) * 100

prop.test(x = c(8, 3), n = c(44, 38))

fisher.test(matrix(c(8, 36, 3, 35), nrow = 2))

#percentage of samples at initiation > 2.0 ng/mL

sum(init_t0$cortisol > 2) #1 above 2 ng/mL, max = 4.29
sum(init_t0$cortisol > 2) / nrow(init_t0) #1.2% greater than 2
high_cortisol_sf_init <- init_t0[init_t0$cortisol > 2, "s_f"]
high_cortisol_sf_init #sample from 1 successful male

#one month
max(one_month_t0$cortisol[one_month_t0$s_f == 1], na.rm = TRUE)
sum(one_month_t0$cortisol[one_month_t0$s_f == 1] > 0.8) #9 above 0.8
sum(one_month_t0$s_f == 1) #41 successful samples
mean(one_month_t0$cortisol[one_month_t0$s_f == 1] > 0.8) * 100

max(one_month_t0$cortisol[one_month_t0$s_f == 0], na.rm = TRUE)
sum(one_month_t0$cortisol[one_month_t0$s_f == 0] > 0.8) #11 above 0.8
sum(one_month_t0$s_f == 0) #40 cannibal samples
mean(one_month_t0$cortisol[one_month_t0$s_f == 0] > 0.8) * 100

prop.test(x = c(9, 11), n = c(41, 40))

fisher.test(matrix(c(9, 41, 11, 40), nrow = 2))

#percentage of samples at one month > 2.0 ng/mL

sum(one_month_t0$cortisol > 2) #7 above 2 ng/mL
sum(one_month_t0$cortisol > 2) / nrow(one_month_t0) #8.6% above 2ng/mL
high_cortisol_sf_om <- one_month_t0[one_month_t0$cortisol > 2, "s_f"]
high_cortisol_sf_om #samples from 2 successful males, 5 cannibals

#binomial model of detected or not detected

#initiation (random effect removed)

init_t0$forest_scaled <- scale(init_t0$forest)

t0_init_1 <- glm(cortisol_detected ~ outcome * forest_scaled,
                 data = init_t0, family = binomial)
summary(t0_init_1)

t0_init_2 <- glm(cortisol_detected ~ outcome + forest_scaled, 
                 data = init_t0, family = binomial)
summary(t0_init_2)

t0_init_3 <- glm(cortisol_detected ~ forest_scaled, 
                 data = init_t0, family = binomial)
summary(t0_init_3)

t0_init_4 <- glm(cortisol_detected ~ outcome, 
                 data = init_t0, family = binomial)
summary(t0_init_4)

t0_init_null <- glm(cortisol_detected ~ 1, 
                    data = init_t0, family = binomial)
summary(t0_init_null)

t0_init_models <- list(
  t0_init_1 = t0_init_1,
  t0_init_2 = t0_init_2,
  t0_init_3 = t0_init_3,
  t0_init_4 = t0_init_4,
  t0_init_null = t0_init_null
)

t0_init_aicc_table <- model.sel(t0_init_models, rank = AICc)

# Convert to data frame
t0_init_aicc_table  <- as.data.frame(t0_init_aicc_table)

# Extract K (number of estimated parameters = df column in model.sel output)
t0_init_aicc_table$K <- t0_init_aicc_table$df

# Add cumulative Akaike weights
t0_init_aicc_table$cum_weight <- cumsum(t0_init_aicc_table$weight)

# Reorder columns to show K and cumulative weight clearly
t0_init_aicc_table <- t0_init_aicc_table[, c("df", "K", "logLik", "AICc", "delta", "weight", "cum_weight")]

# View the cortisol detection at inititation AICc table
print(t0_init_aicc_table)

#one month

one_month_t0$forest_scaled <- scale(one_month_t0$forest)

t0_om_1 <- glmer(cortisol_detected ~ outcome * forest_scaled 
                 + (1|pit), data = one_month_t0, family = binomial)
summary(t0_om_1)

t0_om_2 <- glmer(cortisol_detected ~ outcome + forest_scaled 
                 + (1|pit), data = one_month_t0, family = binomial)
summary(t0_om_2)

t0_om_3 <- glmer(cortisol_detected ~ forest_scaled 
                 + (1|pit), data = one_month_t0, family = binomial)
summary(t0_om_3)

t0_om_4 <- glmer(cortisol_detected ~ outcome 
                 + (1|pit), data = one_month_t0, family = binomial)
summary(t0_om_4)

t0_om_null <- glmer(cortisol_detected ~ 1 + (1|pit), 
                    data = one_month_t0, family = binomial)
summary(t0_om_null)

t0_om_models <- list(
  t0_om_1 = t0_om_1,
  t0_om_2 = t0_om_2,
  t0_om_3 = t0_om_3,
  t0_om_4 = t0_om_4,
  t0_om_null = t0_om_null
)

t0_om_aicc_table <- model.sel(t0_om_models, rank = AICc)

# Convert to data frame
t0_om_aicc_table  <- as.data.frame(t0_om_aicc_table)

# Extract K (number of estimated parameters = df column in model.sel output)
t0_om_aicc_table$K <- t0_om_aicc_table$df

# Add cumulative Akaike weights
t0_om_aicc_table$cum_weight <- cumsum(t0_om_aicc_table$weight)

# Reorder columns to show K and cumulative weight clearly
t0_om_aicc_table<- t0_om_aicc_table[, c("df", "K", "logLik", "AICc", "delta", "weight", "cum_weight")]

# View the cortisol detection at inititation AICc table
print(t0_om_aicc_table)

###section (3) t60 cortisol----

#summary overall

#number of samples below 2 ng/mL
sum(parent_t60$cortisol < 2)
sum(parent_t60$cortisol < 2) / nrow(parent_t60)

#of detected samples (n = 143)
parent_t60_detected_cortisol <- parent_t60[parent_t60$cortisol > 0.8, ]
min(parent_t60_detected_cortisol$cortisol) #minimum was 1.5ng/mL

#initiation

#number of samples below 2 ng/mL
sum(init_t60$cortisol < 2)
sum(init_t60$cortisol < 2) / nrow(init_t60)

#summary statistics of stress-induced cortisol at initiation
#pulling all samples > 0.8 (detected)
init_t60_detected_cortisol <- init_t60[init_t60$cortisol > 0.8, ]

#summary statistics of t60 cortisol at initiation among detected samples
describeBy(x = init_t60_detected_cortisol$cortisol)

#scaling forest cover and total length
init_t60$forest_scaled <- scale(init_t60$forest)
init_t60$smi_scaled <- scale(init_t60$smi)

#model selection
f_init_t60_null <- lm(cortisol ~ 1, data = init_t60)
summary(f_init_t60_null)

f_init_t60_1 <- lm(cortisol ~ outcome + forest_scaled + smi_scaled + 
                     outcome:forest_scaled + outcome:smi_scaled + forest_scaled:smi_scaled +
                     outcome:forest_scaled:smi_scaled, data = init_t60)
summary(f_init_t60_1)

f_init_t60_2 <- lm(cortisol ~ outcome + forest_scaled + smi_scaled + 
                     outcome:forest_scaled + outcome:smi_scaled + forest_scaled:smi_scaled, 
                   data = init_t60)
summary(f_init_t60_2)

f_init_t60_3 <- lm(cortisol ~ outcome + forest_scaled + smi_scaled + 
                     outcome:forest_scaled + forest_scaled:smi_scaled,
                   data = init_t60)
summary(f_init_t60_3)

f_init_t60_4 <- lm(cortisol ~ outcome + forest_scaled + smi_scaled + 
                     outcome:forest_scaled + outcome:smi_scaled, 
                   data = init_t60)
summary(f_init_t60_4)

f_init_t60_5 <- lm(cortisol ~ outcome + forest_scaled + smi_scaled + 
                     forest_scaled:smi_scaled + outcome:smi_scaled, 
                   data = init_t60)
summary(f_init_t60_5)

f_init_t60_6 <- lm(cortisol ~ outcome + forest_scaled + smi_scaled + 
                     forest_scaled:smi_scaled, data = init_t60)
summary(f_init_t60_6)

f_init_t60_7 <- lm(cortisol ~ outcome + forest_scaled + smi_scaled + 
                     outcome:smi_scaled, data = init_t60)
summary(f_init_t60_7)

f_init_t60_8 <- lm(cortisol ~ outcome + forest_scaled + smi_scaled + 
                     outcome:forest_scaled, data = init_t60)
summary(f_init_t60_8)

f_init_t60_9 <- lm(cortisol ~ outcome + forest_scaled + smi_scaled, 
                   data = init_t60)
summary(f_init_t60_9)

f_init_t60_10 <- lm(cortisol ~ outcome + forest_scaled + outcome:forest_scaled, 
                    data = init_t60)
summary(f_init_t60_10)

f_init_t60_11 <- lm(cortisol ~ outcome + smi_scaled + outcome:smi_scaled, 
                    data = init_t60)
summary(f_init_t60_11)

f_init_t60_12 <- lm(cortisol ~ forest_scaled + smi_scaled + forest_scaled:smi_scaled, 
                    data = init_t60)
summary(f_init_t60_12)

f_init_t60_13 <- lm(cortisol ~ outcome + forest_scaled, 
                    data = init_t60)
summary(f_init_t60_13)

f_init_t60_14 <- lm(cortisol ~ outcome + smi_scaled, 
                    data = init_t60)
summary(f_init_t60_14)

f_init_t60_15 <- lm(cortisol ~ forest_scaled + smi_scaled, 
                    data = init_t60)
summary(f_init_t60_15)

f_init_t60_16 <- lm(cortisol ~ outcome, 
                    data = init_t60)
summary(f_init_t60_16)

f_init_t60_17 <- lm(cortisol ~ forest_scaled, 
                    data = init_t60)
summary(f_init_t60_17)

f_init_t60_18 <- lm(cortisol ~ smi_scaled, 
                    data = init_t60)
summary(f_init_t60_18)

t60_init_models <- list(
  f_init_t60_1 = f_init_t60_1,
  f_init_t60_2 = f_init_t60_2,
  f_init_t60_3 = f_init_t60_3,
  f_init_t60_4 = f_init_t60_4,
  f_init_t60_5 = f_init_t60_5,
  f_init_t60_6 = f_init_t60_6,
  f_init_t60_7 = f_init_t60_7,
  f_init_t60_8 = f_init_t60_8,
  f_init_t60_9 = f_init_t60_9,
  f_init_t60_10 = f_init_t60_10,
  f_init_t60_11 = f_init_t60_11,
  f_init_t60_12 = f_init_t60_12,
  f_init_t60_13 = f_init_t60_13,
  f_init_t60_14 = f_init_t60_14,
  f_init_t60_15 = f_init_t60_15,
  f_init_t60_16 = f_init_t60_16,
  f_init_t60_17 = f_init_t60_17,
  f_init_t60_18 = f_init_t60_18,
  f_init_t60_null = f_init_t60_null
)

t60_init_aicc_table <- model.sel(t60_init_models, rank = AICc)

# Convert to data frame
t60_init_aicc_table  <- as.data.frame(t60_init_aicc_table)

# Extract K (number of estimated parameters = df column in model.sel output)
t60_init_aicc_table$K <- t60_init_aicc_table$df

# Add cumulative Akaike weights
t60_init_aicc_table$cum_weight <- cumsum(t60_init_aicc_table$weight)

# Reorder columns to show K and cumulative weight clearly
t60_init_aicc_table<- t60_init_aicc_table[, c("df", "K", "logLik", "AICc", "delta", "weight", "cum_weight")]

# View the cortisol detection at initiation AICc table
print(t60_init_aicc_table)

#one month

#number of samples below 2 ng/mL
sum(one_month_t60$cortisol < 2)
sum(one_month_t60$cortisol < 2) / nrow(one_month_t60)

#summary statistics of stress-induced cortisol at initiation
#pulling all samples > 0.8 (detected)
one_month_t60_detected_cortisol <- one_month_t60[one_month_t60$cortisol > 0.8, ]

#summary statistics of t60 cortisol at one month among detected samples
describeBy(x = one_month_t60_detected_cortisol$cortisol)

#scaling forest cover and total length
one_month_t60$forest_scaled <- scale(one_month_t60$forest)
one_month_t60$tol_scaled <- scale(one_month_t60$tol)

#making sure variables are the correct type
one_month_t60$forest_scaled <- as.numeric(one_month_t60$forest_scaled)
one_month_t60$tol_scaled <- as.numeric(one_month_t60$tol_scaled)
one_month_t60$outcome <- as.factor(one_month_t60$outcome)

#model selection_________________without random effects

f_om_t60_null <- lm(cortisol ~ 1, 
                    data = one_month_t60)
summary(f_om_t60_null)

f_om_t60_1 <- lm(cortisol ~ outcome + forest_scaled + tol_scaled + 
                   outcome:forest_scaled + outcome:tol_scaled + forest_scaled:tol_scaled +
                   outcome:forest_scaled:tol_scaled, data = one_month_t60)
summary(f_om_t60_1)

f_om_t60_2 <- lm(cortisol ~ outcome + forest_scaled + tol_scaled + 
                   outcome:forest_scaled + outcome:tol_scaled + forest_scaled:tol_scaled, 
                 data = one_month_t60)
summary(f_om_t60_2)

f_om_t60_3 <- lm(cortisol ~ outcome + forest_scaled + tol_scaled + 
                   outcome:forest_scaled + forest_scaled:tol_scaled,
                 data = one_month_t60)
summary(f_om_t60_3)

f_om_t60_4 <- lm(cortisol ~ outcome + forest_scaled + tol_scaled + 
                   outcome:forest_scaled + outcome:tol_scaled, 
                 data = one_month_t60)

f_om_t60_5 <- lm(cortisol ~ outcome + forest_scaled + tol_scaled + 
                   forest_scaled:tol_scaled + outcome:tol_scaled, 
                 data = one_month_t60)
summary(f_om_t60_5)

f_om_t60_6 <- lm(cortisol ~ outcome + forest_scaled + tol_scaled + 
                   forest_scaled:tol_scaled, data = one_month_t60)
summary(f_om_t60_6)

f_om_t60_7 <- lm(cortisol ~ outcome + forest_scaled + tol_scaled + 
                   outcome:tol_scaled, data = one_month_t60)
summary(f_om_t60_7)

f_om_t60_8 <- lm(cortisol ~ outcome + forest_scaled + tol_scaled + 
                   outcome:forest_scaled, data = one_month_t60)
summary(f_om_t60_8)

f_om_t60_9 <- lm(cortisol ~ outcome + forest_scaled + tol_scaled, 
                 data = one_month_t60)
summary(f_om_t60_9)

f_om_t60_10 <- lm(cortisol ~ outcome + forest_scaled + outcome:forest_scaled, 
                  data = one_month_t60)
summary(f_om_t60_10)

f_om_t60_11 <- lm(cortisol ~ outcome + tol_scaled + outcome:tol_scaled, 
                  data = one_month_t60)
summary(f_om_t60_11)

f_om_t60_12 <- lm(cortisol ~ forest_scaled + tol_scaled + forest_scaled:tol_scaled, 
                  data = one_month_t60)
summary(f_om_t60_12)

f_om_t60_13 <- lm(cortisol ~ outcome + forest_scaled, 
                  data = one_month_t60)
summary(f_om_t60_13)

f_om_t60_14 <- lm(cortisol ~ outcome + tol_scaled, 
                  data = one_month_t60)
summary(f_om_t60_14)

f_om_t60_15 <- lm(cortisol ~ forest_scaled + tol_scaled, 
                  data = one_month_t60)
summary(f_om_t60_15)

f_om_t60_16 <- lm(cortisol ~ outcome, 
                  data = one_month_t60)
summary(f_om_t60_16)

f_om_t60_17 <- lm(cortisol ~ forest_scaled, 
                  data = one_month_t60)
summary(f_om_t60_17)

f_om_t60_18 <- lm(cortisol ~ tol_scaled, 
                  data = one_month_t60)
summary(f_om_t60_18)

t60_om_models <- list(
  f_om_t60_1 = f_om_t60_1,
  f_om_t60_2 = f_om_t60_2,
  f_om_t60_3 = f_om_t60_3,
  f_om_t60_4 = f_om_t60_4,
  f_om_t60_5 = f_om_t60_5,
  f_om_t60_6 = f_om_t60_6,
  f_om_t60_7 = f_om_t60_7,
  f_om_t60_8 = f_om_t60_8,
  f_om_t60_9 = f_om_t60_9,
  f_om_t60_10 = f_om_t60_10,
  f_om_t60_11 = f_om_t60_11,
  f_om_t60_12 = f_om_t60_12,
  f_om_t60_13 = f_om_t60_13,
  f_om_t60_14 = f_om_t60_14,
  f_om_t60_15 = f_om_t60_15,
  f_om_t60_16 = f_om_t60_16,
  f_om_t60_17 = f_om_t60_17,
  f_om_t60_18 = f_om_t60_18,
  f_om_t60_null = f_om_t60_null
)

t60_om_aicc_table <- model.sel(t60_om_models, rank = AICc)

# Convert to data frame
t60_om_aicc_table  <- as.data.frame(t60_om_aicc_table)

# Extract K (number of estimated parameters = df column in model.sel output)
t60_om_aicc_table$K <- t60_om_aicc_table$df

# Add cumulative Akaike weights
t60_om_aicc_table$cum_weight <- cumsum(t60_om_aicc_table$weight)

# Reorder columns to show K and cumulative weight clearly
t60_om_aicc_table<- t60_om_aicc_table[, c("df", "K", "logLik", "AICc", "delta", "weight", "cum_weight")]

# View the cortisol detection at initiation AICc table
print(t60_om_aicc_table)

anova(f_om_t60_17, f_om_t60_3, test = "Chisq")
anova(f_om_t60_null, f_om_t60_3, test = "Chisq")

###section (3b) t60 cortisol predictions & inquiries----

#(1) interaction between forest cover and nest fate

t60_fate_low_forest <- ggpredict(f_om_t60_3, terms = c("forest_scaled [-1.9348200]", "outcome"))
t60_fate_high_forest <- ggpredict(f_om_t60_3, terms = c("forest_scaled [1.0185190]", "outcome"))

#creating 200 evenly spaced raw forest cover values
forest_raw_vals <- seq(53.8, 68.3, length.out = 200)
forest_scaled_vals <- (forest_raw_vals - 63.1) / 4.67

#build dataframe grid with both nest fates
om_t60_forest_pred <- expand.grid(
  forest_scaled = forest_scaled_vals,
  outcome = levels(one_month_t60$outcome),
  tol_scaled = 0  # hold tol_scaled at mean
)

#using the f_om_t60_3 model to predict t60 cortisol with SEs (for Figure 2)
pred <- predict(f_om_t60_3, newdata = om_t60_forest_pred, se.fit = TRUE)
om_t60_forest_pred$predicted <- pred$fit
om_t60_forest_pred$se <- pred$se.fit

#back-converting scaled and centered forest values to raw %s
om_t60_forest_pred$forest <- rep(forest_raw_vals, 
                                 times = length(levels(one_month_t60$outcome)))

#for reviwer inquiry removing 2 low successful males_____________
filtered_data$forest_scaled <- scale(filtered_data$forest)
filtered_data$tol_scaled <- scale(filtered_data$tol)

filtered_data$forest_scaled <- as.numeric(filtered_data$forest_scaled)
filtered_data$tol_scaled <- as.numeric(filtered_data$tol_scaled)
filtered_data$outcome <- as.factor(filtered_data$outcome)

f_om_t60_3_filter <- lm(cortisol ~ outcome + forest_scaled + tol_scaled + 
                          outcome:forest_scaled + forest_scaled:tol_scaled,
                        data = filtered_data)
summary(f_om_t60_3_filter)

#(2) interaction between forest cover and total length


#re-run of model for for males > 38cm__________________

one_month_t60_big$forest_scaled <- scale(one_month_t60_big$forest)
one_month_t60_big$tol_scaled <- scale(one_month_t60_big$tol)

f_om_t60_3_big <- lm(cortisol ~ outcome + forest_scaled + tol_scaled + 
                       outcome:forest_scaled + forest_scaled:tol_scaled, 
                     data = one_month_t60_big)
summary(f_om_t60_3_big)

#generating predictions for figure 3___________________

# 1. generating evenly spaced values of forest and tol
forest_raw_vals <- seq(53.7, 68.6, length.out = 20)
tol_raw_vals <- seq(38.0, 56.3, length.out = 20) #bottoming out at 38 cm minimum

# 2. scaling values using the means and SDs used in the model
mean(one_month_t60$forest)
sd(one_month_t60$forest)

forest_scaled_vals <- (forest_raw_vals - 63.14079) / 4.672677
tol_scaled_vals <- (tol_raw_vals - 45.56053) / 4.71582

# 3. prediction grid (across both nest fates)
om_t60_tol_pred <- expand.grid(
  forest_scaled = forest_scaled_vals,
  tol_scaled = tol_scaled_vals,
  outcome = levels(one_month_t60$outcome)
)

# 4. predict the prediction grid using the model
om_t60_tol_pred$cortisol <- predict(f_om_t60_3, newdata = om_t60_tol_pred)

# 5. converting scaled values to actual values for plotting
om_t60_tol_pred$forest <- om_t60_tol_pred$forest_scaled * 4.672677 + 63.14079
om_t60_tol_pred$tol <- om_t60_tol_pred$tol_scaled *  4.71582 + 45.56053


#generating predictions for figure 4___________________

#54% forest cover

om_t60_tol_54 <- ggpredict(f_om_t60_3, 
                           terms = c("tol_scaled [all]", "forest_scaled [-1.9348200]"))

tol_mean <- mean(one_month_t60$tol, na.rm = TRUE)
tol_sd <- sd(one_month_t60$tol, na.rm = TRUE)

om_t60_tol_54$tol <- om_t60_tol_54$x * tol_sd + tol_mean

om_t60_tol_54 <- subset(om_t60_tol_54, tol > 38) #cropping to values >38cm

om_t60_tol_54$cortisol <- om_t60_tol_54$predicted
om_t60_tol_54$selower <- om_t60_tol_54$cortisol - om_t60_tol_54$std.error
om_t60_tol_54$sehigher <- om_t60_tol_54$cortisol + om_t60_tol_54$std.error

om_t60_tol_54$ymin = om_t60_tol_54$cortisol - (1.96 * om_t60_tol_54$std.error)
om_t60_tol_54$ymax = om_t60_tol_54$cortisol + (1.96 * om_t60_tol_54$std.error)

#68% forest cover

om_t60_tol_68 <- ggpredict(f_om_t60_3, 
                           terms = c("tol_scaled [all]", "forest_scaled [1.0185190]"))

tol_mean <- mean(one_month_t60$tol, na.rm = TRUE)
tol_sd <- sd(one_month_t60$tol, na.rm = TRUE)

om_t60_tol_68$tol <- om_t60_tol_68$x * tol_sd + tol_mean

om_t60_tol_68 <- subset(om_t60_tol_68, tol > 38)  #cropping to values >38cm

om_t60_tol_68$cortisol <- om_t60_tol_68$predicted
om_t60_tol_68$selower <- om_t60_tol_68$cortisol - om_t60_tol_68$std.error
om_t60_tol_68$sehigher <- om_t60_tol_68$cortisol + om_t60_tol_68$std.error

om_t60_tol_68$ymin = om_t60_tol_68$cortisol - (1.96 * om_t60_tol_68$std.error)
om_t60_tol_68$ymax = om_t60_tol_68$cortisol + (1.96 * om_t60_tol_68$std.error)

###section (4) VIF analysis of t60 cortisol model----

f_om_t60_3_VIF <- vif(f_om_t60_3, type = "predictor")
f_om_t60_3_VIF

#PLOTS------

#figure_1 A_______________________________________

init_t0$outcome_ordered <- ordered(init_t0$outcome, c("success","cannibal"))

t_init <- ggplot(init_t0, aes(outcome_ordered,testosterone,fill = outcome_ordered))+
  geom_bar(position = "dodge", stat = 'summary', fun.y ='mean', alpha = 0.5, width = 0.9)+
  scale_fill_manual(values=c("black","blue"), labels=c('Success', 'Cannibal'))+
  geom_point(aes(x = outcome_ordered), shape = 21, size = 5, alpha = 0.5, position = position_jitterdodge(
    jitter.width = 0.9, jitter.height=0.0, dodge.width = 0.9))+
  stat_summary(fun.data = mean_se, position =position_dodge(0.9), width = 0.25, size = 2, geom = "errorbar")+
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black", size = 2, linetype = "solid"),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(0.3,"cm"),
        axis.title.x = element_text(size = 30, color = "black", margin=margin(12,0,0,0)),
        axis.title.y = element_text(size = 40, color = "black", margin=margin(0,12,0,0)),
        axis.text.x = element_text(size =30, color = "black"),
        axis.text.y = element_text(size =40, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size =20, color = "black"),
        legend.position = "none")+
  labs(x = "", y = "Testosterone (ng/mL)")+
  scale_x_discrete(labels=c("Success","Cannibal","Other Fail"))+
  expand_limits(y = c(0,100))+
  scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "black", size = 1.5)

line(x=c(0,1),y=c(20,20))

t_init <- t_init +
  geom_segment(aes(x=1, y=100, xend=2, yend=100), size = 2)+
  geom_text(aes(x = 1.5, y = 108), label = c("p = 0.06"), size = 10)

t_init

ggsave("figure_1a.tiff", plot = t_init,
       dpi = 300,
       width = 9, height = 10,
       units = "in")

#figure_1 B__________________________________________

init_t0$outcome_ordered <- ordered(init_t0$outcome, c("success","cannibal"))

dht_init <- ggplot(init_t0, aes(outcome_ordered,dht,fill = outcome_ordered))+
  geom_bar(position = "dodge", stat = 'summary', fun.y ='mean', alpha = 0.5, width = 0.9)+
  scale_fill_manual(values=c("black","blue"), labels=c('Success', 'Cannibal'))+
  geom_point(aes(x = outcome_ordered), shape = 21, size = 5, alpha = 0.5, position = position_jitterdodge(
    jitter.width = 0.9, jitter.height=0.0, dodge.width = 0.9))+
  stat_summary(fun.data = mean_se, position =position_dodge(0.9), width = 0.25, size = 2, geom = "errorbar")+
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black", size = 2, linetype = "solid"),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(0.3,"cm"),
        axis.title.x = element_text(size = 30, color = "black", margin=margin(12,0,0,0)),
        axis.title.y = element_text(size = 40, color = "black", margin=margin(0,12,0,0)),
        axis.text.x = element_text(size =30, color = "black"),
        axis.text.y = element_text(size =40, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size =20, color = "black"),
        legend.position = "none")+
  labs(x = "", y = "DHT (ng/mL)")+
  scale_x_discrete(labels=c("Success","Cannibal"))+
  expand_limits(y = c(0,31))+
  scale_y_continuous(breaks = c(0,5,10,15,20,25,30))+
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "black", size = 1.5)


line(x=c(0,1),y=c(20,20))

dht_init <- dht_init +
  geom_segment(aes(x=1, y=30, xend=2, yend=30), size = 2)+
  geom_text(aes(x = 1.5, y = 32), label = c("p = 0.03*"), size = 10)

dht_init

ggsave("figure_1b.tiff", plot = dht_init,
       dpi = 300,
       width = 9, height = 10,
       units = "in")

#figure 2___________________________________________

om_t60_forest_pred$outcome <- ordered(om_t60_forest_pred$outcome, c("success","cannibal"))
om_t60_forest_pred$cortisol <- om_t60_forest_pred$predicted
om_t60_forest_pred$selower <- (om_t60_forest_pred$cortisol - om_t60_forest_pred$se)
om_t60_forest_pred$sehigher <- (om_t60_forest_pred$cortisol + om_t60_forest_pred$se)

predict_om_t60_forest <- ggplot(om_t60_forest_pred, aes(x = forest, y = cortisol, fill = outcome, color = outcome))+
  geom_line(size = 1.8, show.legend = FALSE)+
  scale_color_manual(values=c("black", "blue"), labels=c('Success','Cannibal'))+
  scale_fill_manual(values=c("black", "blue"), labels=c('Success','Cannibal'))+
  geom_ribbon(aes(ymin=selower, ymax=sehigher), alpha = 0.3)+
  geom_jitter(data = one_month_t60, aes(x = forest, y = cortisol, fill = outcome), size = 5, alpha = 0.5, 
              width = 0.5, show.legend = FALSE)+
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black", size = 2, linetype = "solid"),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(0.3,"cm"),
        axis.title.x = element_text(size = 40, color = "black", margin=margin(12,0,0,0)),
        axis.title.y = element_text(size = 42, color = "black", margin=margin(0,12,0,0)),
        axis.text.x = element_text(size =45, color = "black"),
        axis.text.y = element_text(size =45, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size =20, color = "black"),
        legend.position = c(0.85,0.93),
        legend.title = element_blank(),
        legend.text = element_text(size = 35))+
  labs(x = "Forest Cover (%)", y = "Stress-Induced Cortisol (ng/mL)")+
  expand_limits(y = c(0,14.0))+
  scale_x_continuous(breaks = c(54,56,58,60,62,64,66,68))+
  scale_y_continuous(breaks = c(0,2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,20.0,22.0,24.0,26.0))+
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "black", size = 2)
predict_om_t60_forest 

ggsave("figure_2_test.tiff", plot = predict_om_t60_forest,
       dpi = 300,
       width = 10, height = 9,
       units = "in")

#figure 3_________________________________________

library(lattice)
library(latticeExtra)
library(viridisLite)

modifyList(myParSettings, list(axis.line = list(lwd = 2)))

my_breaks <- c(0, 1, 2, 3, 4, 6, 8, 10, 12, 14, 16)

figure3 <- levelplot(cortisol ~ forest * tol, om_t60_tol_pred,
                     panel = panel.levelplot.points, cex = 0,
                     scales = list(x = list(at = c(52, 54, 56, 58, 60, 62, 64, 66, 68, 70), cex = 2.2),
                                   y = list(at = c(38, 40, 42, 44, 46, 48, 50, 52, 54, 56), cex = 2.2)),
                     at = seq(2.15, 15.6, length.out = 20),
                     col.regions = viridis(17, begin = 0.01, end = .99, alpha = 0.8),
                     xlim = c(53.3,68.9), ylim = c(37.5,56.8),
                     par.settings = list(axis.line = list(lwd = 3)),          
                     xlab = list(label = "Forest Cover %", cex = 2.6),         
                     ylab = list(label = "Total Length (cm)", cex = 2.6),
                     colorkey = list(
                       at = seq(2, 16, length.out = 20),  
                       col = viridis(20, begin = 0.01, end = .99, alpha = 0.75),
                       labels = list(at = seq(0, 20, by = 2), cex = 2.2),
                       width = 2.7, height = 0.9,
                       title = list(label = "Stress-Induced Cortisol (ng/mL)", cex = 2.2)))+
  layer_(panel.2dsmoother(..., n = 200)) +
  layer(panel.points(jitter(one_month_t60$forest, factor = 15), 
                     jitter(one_month_t60$tol, factor = 5),
                     pch = 19, cex = 2.0, col = "black", alpha = 0.7))

print(figure3)


tiff("figure_3_remake.tiff", width = 8, height = 6.5, units = "in", res = 300)
print(figure3)
dev.off()


#figure 4______________________________________

figure4a <- ggplot(om_t60_tol_54, aes(x = tol, y = cortisol))+
  geom_line(size = 1.8, show.legend = FALSE)+
  geom_ribbon(aes(ymin=ymin, ymax=ymax), alpha = 0.2)+
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black", size = 2, linetype = "solid"),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(0.3,"cm"),
        axis.title.x = element_text(size = 40, color = "black", margin=margin(12,0,0,0)),
        axis.title.y = element_text(size = 38, color = "black", margin=margin(0,12,0,0)),
        axis.text.x = element_text(size =40, color = "black"),
        axis.text.y = element_text(size =45, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size =20, color = "black"))+
  labs(x = "Total Length (cm)", y = "Stress-Induced Cortisol (ng/mL)")+
  expand_limits(y = c(0,21.0))+
  scale_x_continuous(breaks = c(38,40,42,44,46,48,50,52,54,56))+
  scale_y_continuous(breaks = c(0,3,6,9,12,15,18,21))+
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "black", size = 2)+
  annotate("text", x = 56.5, y = 20.8, label = "Forest Cover = 54%", 
           size = 12, hjust = 1, vjust = 1)
figure4a 

ggsave("figure_4a.tiff", plot = figure4a,
       dpi = 300,
       width = 10, height = 9,
       units = "in")

figure4b <- ggplot(om_t60_tol_68, aes(x = tol, y = cortisol))+
  geom_line(size = 1.8, show.legend = FALSE)+
  geom_ribbon(aes(ymin=ymin, ymax=ymax), alpha = 0.2)+
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black", size = 2, linetype = "solid"),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(0.3,"cm"),
        axis.title.x = element_text(size = 40, color = "black", margin=margin(12,0,0,0)),
        axis.title.y = element_text(size = 38, color = "black", margin=margin(0,12,0,0)),
        axis.text.x = element_text(size =40, color = "black"),
        axis.text.y = element_text(size =45, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size =20, color = "black"))+
  labs(x = "Total Length (cm)", y = "Stress-Induced Cortisol (ng/mL)")+
  expand_limits(y = c(0,21.0))+
  scale_x_continuous(breaks = c(38,40,42,44,46,48,50,52,54,56))+
  scale_y_continuous(breaks = c(0,3,6,9,12,15,18,21))+
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "black", size = 2)+
  annotate("text", x = 56.5, y = 20.8, label = "Forest Cover = 68%", 
           size = 12, hjust = 1, vjust = 1)
figure4b 

ggsave("figure_4b.tiff", plot = figure4b,
       dpi = 300,
       width = 10, height = 9,
       units = "in")