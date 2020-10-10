setwd("U:/Private/IHF")
install.packages(c("fhidata","fhi","fhiplot"), repos = c("https://folkehelseinstituttet.github.io/drat", "https://cran.rstudio.com"))

#Task 1
################################################

###### Q5 - Import the file and create indlev as the dataset
indlev <- readRDS("U:/Private/IHF/individual_level_data.RDS")

#Check the characteristics of the dataset
head(indlev, n = 10)
tail(indlev, n = 10)
str(indlev)
unique(indlev$location_code)
table(indlev$location_code, indlev$value)

#Add a new var called 'num_case'
indlev$num_case <- 1

require(dplyr)
require(tidyr)

#Summarize cases by location and date
ind_sum <- summarize(group_by(indlev, location_code, date),sum(num_case))

#Extract location_code
location_code <- c(unique(ind_sum$location_code))

#Make a date list from 2010 - 2020
dayseq <- seq(as.Date("2010-01-01"), as.Date("2020-12-31"), by="day")

seqdate <- data.frame(date = dayseq)
seqdate$location_code <- list(location_code)

#A full dataset of location and date
seqdate <- unnest(seqdate, cols = c(location_code))
ind_sum <- as.data.frame(ind_sum)
seqdate <- as.data.frame(seqdate)

compdf <- full_join(ind_sum, seqdate)
names(compdf )[3] <- "num_case" 

###### Q6 & 7 Get the complete dataset, replace NA with 0 and sort the dataset
compdf <- compdf %>%
          mutate(num_case = if_else(is.na(num_case), 0, num_case))

compdf <- compdf[order(compdf$location_code, as.Date(compdf$date)),]

###### Q8
require(surveillance)
require(tidyverse)

iso_list <- isoWeekYear(compdf$date)
iso_frame <- as.data.frame(iso_list)
iso_frame$iso_vec <- paste(iso_list$ISOYear,iso_list$ISOWeek, sep="/")

compdf_1 <- cbind(compdf, iso_frame)
compdf_lo_iso <- summarize(group_by(compdf_1, location_code, iso_vec),sum(num_case))
compdf_lo_iso_1 <- separate(compdf_lo_iso, iso_vec, into = c("year", "week"), remove = FALSE)

compdf_lo_iso_1$year <- as.numeric(compdf_lo_iso_1$year)
compdf_lo_iso_1$week <- as.numeric(compdf_lo_iso_1$week)
compdf_lo_iso_2 <- compdf_lo_iso_1[order(compdf_lo_iso_1$location_code, compdf_lo_iso_1$year, compdf_lo_iso_1$week),]

###### Q9
norloc <- fhidata::norway_locations_b2020
norloc <- norloc[order(norloc$municip_code),]
names(norloc)[1] <- "location_code"

comp_fyk <- merge(x=compdf_lo_iso_2,y=norloc,by="location_code",all.x=TRUE)
comp_fyk$nation <- "norge"
names(comp_fyk)[5] <- "num_case"
 
comp_fyk_sum <- aggregate(num_case ~ county_name, comp_fyk, FUN = sum)
comp_nat_sum <- aggregate(num_case ~ nation, comp_fyk, FUN = sum)

###### Q10
comp <- subset(comp_fyk, select = -c(region_code, region_name, faregion_name, faregion_code))
# Add index to each municip, country and nation 
comp$ID_muni <- cumsum(!duplicated(comp[1]))
comp$ID_cnty <- cumsum(!duplicated(comp[7]))
comp$ID_nat <- cumsum(!duplicated(comp[9]))

################## LOOP for 
################## each of the 356 municipalities (municip*), 
################## 11 counties (county*), and the national level ("norge")

sub_comp <- subset(comp, ID_muni == 1)

training <- subset(sub_comp,sub_comp$year >= 2010 & sub_comp$year <= 2019)
testing <- subset(sub_comp,sub_comp$year >= 2020)

###### Q11: Both codes and the following plots indicate seasonal variance
#Create a time-series object
Xdisprog <- create.disProg(week = training$week, observed = training$num_case,
        state = training$location_code, 
        freq = 52, start = c(2010, 1))

#Surveilance time series
Xsts <- disProg2sts(Xdisprog)
dim(Xsts)

f_S1 <- addSeason2formula(f = ~ 1, S = 1, period = 52)

# fit the Poisson model
result0 <- hhh4(Xsts, control = list(end = list(f = f_S1), family = "Poisson"))
summary(result0)
#Negative binomial model
result1 <- update(result0, family = "NegBin1")
AIC(result0, result1)
# fit an autoregressive model
result2 <- update(result1, ar = list(f = ~ 1))

#exp(ar.1) is the autoregressive coefficient and can be interpreted as the epidemic proportion of disease incidence
coef(result2, se = TRUE, amplitudeShift = TRUE, idx2Exp = TRUE) 

#Clear sign of seasonal variation
plot(Xsts, type = observed ~ time | unit, same.scale = FALSE, col = "grey")
plot(result2)

#Generate formula for temporal and seasonal trends
model_basic <- list(end = list(f = addSeason2formula(~1 + t, period = 52),
             offset = population(Xsts)),
             ar = list(f = ~1),
             end = list(f =f.end),
             family = "NegBin1")

# Estimate model
fit_basic <- hhh4(stsObj = Xsts,control = model_basic)
summary(fit_basic, idx2Exp = TRUE, amplitudeShift = TRUE, maxEv = TRUE)

#The resulting multiplicative effect of seasonality is shown 
plot(fit_basic, type = "season", components = "end", main = "")

#Check if negative binomial distribution is better than Poisson
AIC(fit_basic, update(fit_basic, family = "Poisson"))
confint(fit_basic, parm = "overdisp")

#Predictive model: do sequential one-step-ahead predictions for all time points
tp <- c(2, 520)
Xpred1 <- oneStepAhead(fit_basic, tp = tp, type = "final")
sc1 <- scores(Xpred1 , reverse = FALSE)

######Q12 95% CI
predci <- confint(Xpred1)

######Q13 Find potential outbreaks
predci1 <- as.data.frame(predci)
predci1$seq <- (1:nrow(predci1) + 2)
names(predci1)[1] <- "pred_CIL"
names(predci1)[2] <- "pred_CIU"

training$seq <- 1:nrow(training)
training_re <- merge(x=training,y=predci1,by="seq",all.y=TRUE)
training_re$ob <- ifelse(training_re$num_case > training_re$pred_CIU, 1, 0)

######Q14 Exclude the potential outbreaks
training_NoOb <- subset(training_re, ob == 0)

######Q15 Refit the model using new dataset
Xdisprog_re <- create.disProg(week = training_NoOb$week, observed = training_NoOb$num_case,
        state = training_NoOb$location_code, 
        freq = 50, start = c(2010, 3))

#Surveilance time series
Xsts_re <- disProg2sts(Xdisprog_re)
dim(Xsts)

#Generate formula for temporal and seasonal trends
model2_basic <- list(end = list(f = addSeason2formula(~1 + t, period = 50),
             offset = population(Xsts_re)),
             ar = list(f = ~1),
             end = list(f =f.end),
             family = "NegBin1")

# Estimate model
fit_basic_re <- hhh4(stsObj = Xsts_re,control = model2_basic)
summary(fit_basic_re, idx2Exp = TRUE, amplitudeShift = TRUE, maxEv = TRUE)

#Check if negative binomial distribution is better than Poisson
AIC(fit_basic_re, update(fit_basic_re, family = "Poisson"))
confint(fit_basic_re, parm = "overdisp")

#Predictive model: do sequential one-step-ahead predictions for all time points
tp <- c(2, 507)
Xpred3 <- oneStepAhead(fit_basic_re, tp = tp, type = "final")
sc3 <- scores(Xpred3 , reverse = FALSE)


######Q16 95% CI
predci_re <- confint(Xpred3)

######Q13 Find potential outbreaks
predci2 <- as.data.frame(predci_re)
predci2$seq <- (1:nrow(predci2) + 2)
names(predci2)[1] <- "pred2_CIL"
names(predci2)[2] <- "pred2_CIU"

training$seq <- 1:nrow(training)
training_re2 <- merge(x=training_re,y=predci2,by="seq",all.y=TRUE)
training_re2$ob <- ifelse(training_re2$num_case > training_re2$pred2_CIU, 1, 0)

######Q17 Identify the potential outbreaks
training_re2_Ob <- subset(training_re2, ob == 1)

















