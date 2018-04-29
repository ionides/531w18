library(plyr)
library(reshape2)
library(readxl)
library(data.table)
library(seaborn)
library(ggplot2)
library(depmixS4)
library(ggplot2)
library(gridExtra)
library(kableExtra)
library(dplyr)
library(data.table)
library(MASS)
library(zoo)
library(xts)

PSID <- read_excel("C:/Users/Sandra/Desktop/U of Michigan/Work/PSID Excel/J241976.xlsx", #Change reference of file
                   col_types = c("skip", "skip", "skip", 
                                 "numeric", "numeric", "numeric", 
                                 "numeric", "numeric", "skip", 
                                 "numeric", "numeric", "numeric", 
                                 "numeric", "numeric", "numeric", 
                                 "numeric", "skip", "numeric", 
                                 "numeric", "numeric", "numeric", 
                                 "numeric", "numeric", "numeric", 
                                 "skip", "numeric", "numeric", 
                                 "numeric", "numeric", "numeric", 
                                 "numeric", "numeric", "skip", 
                                 "numeric", "numeric", "numeric", 
                                 "numeric", "numeric", "numeric", 
                                 "numeric", "skip", "numeric", 
                                 "numeric", "numeric", "numeric", 
                                 "numeric", "numeric", "numeric", 
                                 "skip", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric"))

PSID <- `colnames<-`(PSID, c("ReleaseNumber.2005", "TotalIncome.2005", "MaritalStatus.2005", "InterviewNum.2005", "SequenceNum.2005", "TypeofRecord.2005", "WhyNonresponsive.2005",
                             "ReleaseNumber.2007", "TotalIncome.2007", "MaritalStatus.2007", "InterviewNum.2007", "SequenceNum.2007", "TypeofRecord.2007", "WhyNonresponsive.2007",
                             "ReleaseNumber.2009", "TotalIncome.2009", "MaritalStatus.2009", "InterviewNum.2009", "SequenceNum.2009", "TypeofRecord.2009", "WhyNonresponsive.2009",
                             "ReleaseNumber.2011", "TotalIncome.2011", "MaritalStatus.2011", "InterviewNum.2011", "SequenceNum.2011", "TypeofRecord.2011", "WhyNonresponsive.2011",
                             "ReleaseNumber.2013", "TotalIncome.2013", "MaritalStatus.2013", "InterviewNum.2013", "SequenceNum.2013", "TypeofRecord.2013", "WhyNonresponsive.2013",
                             "ReleaseNumber.2015", "TotalIncome.2015", "MaritalStatus.2015", "InterviewNum.2015", "SequenceNum.2015", "TypeofRecord.2015", "WhyNonresponsive.2015",
                             "EmploymentStatus.2005", "EmploymentStatus.2007", "EmploymentStatus.2009", "EmploymentStatus.2011", "EmploymentStatus.2013", "EmploymentStatus.2015"))

PSID2 <- read_excel("C:/Users/Sandra/Desktop/U of Michigan/Work/PSID Excel/J243911.xlsx") #Change reference

PSID <- cbind(PSID, PSID2)

#Convert to long format -- Move sex to beginning of dataframe first
#x <- "Sex"
PSID <- PSID[c(x, setdiff(names(PSID), x))]
#PSID <- PSID[,2:79]
PSID <- reshape(PSID[,2:79], dir = "long", varying = 1:78, sep = ".")
#Remove all variables except covariates, ID and response
PSID.pomp <- subset(PSID[,c("TotalIncome", "EmploymentStatus", "OwnOrRent", "Race",
                            "Age", "YearsCompletedEducation", "id", "time", "InterviewNum", "SequenceNum")])

PSID$Race <- ifelse(PSID$WifeRace1 == 1 & PSID$HeadRace1 == 1, 1,
                    ifelse(PSID$WifeRace1 == 1 & PSID$HeadRace1 == 0, 1,
                           ifelse(PSID$WifeRace1 == 0 & PSID$HeadRace1 == 1, 1,
                                  ifelse(PSID$WifeRace1 == 2 & PSID$WifeRace1 == 2, 2,
                                         ifelse(PSID$WifeRace1 == 2 & PSID$HeadRace1 == 0, 2,
                                                ifelse(PSID$WifeRace1 == 0 & PSID$HeadRace1 == 2, 2,
                                                       ifelse(PSID$WifeRace1 == 3 & PSID$HeadRace1 == 3, 3, 
                                                              ifelse(PSID$WifeRace1 == 3 & PSID$HeadRace1 == 0, 3,
                                                                     ifelse(PSID$WifeRace1 == 0 & PSID$HeadRace1 == 3, 3,
                                                                            ifelse(PSID$WifeRace1 == 4 & PSID$HeadRace1 == 4, 4,
                                                                                   ifelse(PSID$WifeRace1 == 4 & PSID$HeadRace1 == 0, 4,
                                                                                          ifelse(PSID$WifeRace1 == 0 & PSID$HeadRace1 == 4, 4, 
                                                                                                 ifelse(PSID$WifeRace1 == 5 & PSID$HeadRace1 == 5, 5,
                                                                                                        ifelse(PSID$WifeRace1 == 5 & PSID$HeadRace1 == 0, 5,
                                                                                                               ifelse(PSID$WifeRace1 == 0 & PSID$HeadRace1 == 5, 5, 0)))))))))))))))



summary(PSID.pomp$TotalIncome)

#Collapse by interview ID, year and income level. Since income level is the same across all R's in a 
#household I am only using the unique value. Marital status is also the same. 
PSID.pomp <- as.data.table(PSID.pomp)
PSID.pomp <- PSID.pomp[, InterviewPersonNum := paste(InterviewNum,SequenceNum,sep = ".")]
PSID.pomp <- PSID.pomp[, SequenceNum := NULL]

#To observe latent states within the poor, I will subset the data using only households with incomes < median
PSID.pomp <- subset(PSID.pomp, PSID.pomp$TotalIncome < median(PSID.pomp$TotalIncome, na.rm = TRUE))
PSID.pomp <- na.omit(PSID.pomp)

write.csv(PSID.pomp, "PSID.pomp.csv")

matplot(t(PSID.pomp$TotalIncome),
        ylab="Total Income",xlab="6-month intervals", 
        type="l",xaxp=c(1,4,3))

tapply(PSID.pomp$TotalIncome, PSID.pomp$time, summary)

#Remove everyone with income loss. 
PSID.pomp <- subset(PSID.pomp, PSID.pomp$TotalIncome >= 0)
#Code no income as  0.1 instead of 0 for computations purposes
PSID.pomp$TotalIncome <- ifelse(PSID.pomp$TotalIncome == 0, 0.1, 
                                ifelse(PSID.pomp$TotalIncome <= 0, 0.1, PSID.pomp$TotalIncome))
PSID.pomp$EmploymentStatus <- ifelse(PSID.pomp$EmploymentStatus == 1, 1, 
                                     ifelse(PSID.pomp$EmploymentStatus == 2 | PSID.pomp$EmploymentStatus == 3 |
                                              PSID.pomp$EmploymentStatus == 6, 2,
                                            ifelse(PSID.pomp$EmploymentStatus == 4, 3, 0)))

PSID.pomp$Race2 <- ifelse(PSID.pomp$Race == 1, 1, 0) #1 is white, non white 



install.packages("depmixS4")
#Remove missing data
PSID <- na.omit(PSID)
set.seed(1234)
#Scale income
PSID.pom$TotalIncome <- scale(PSID.pom$TotalIncome)
  d <- density(PSID.pomp$TotalIncome)
  plot(d, type="n", main = "Density Plot Total Income")
  polygon(d, col="grey", border="gray")
summary(d)

set.seed(1)
x <- rlnorm(100) # random values from a log-normal distribution

# fit distributions to decide which parameters to use in model
library(MASS)
lognormal <- fitdistr(PSID.pomp$TotalIncome, "lognormal")
gamma <- fitdistr(PSID.pomp$TotalIncome, "gamma")
gaussian <- fitdistr(PSID.pomp$TotalIncome, "normal")
weibull < fitdistr(PSID.pomp$TotalIncome, "weibull")
# compare AICs
AIC(lognormal)
AIC(gaussian)
#Normal distribution fits the model better. Weibull optimization failed, 
#gamma algorithm hit a non-finite finite-difference value
#CCF
library(zoo)
library(xts)
PSID.pomp$time <- as.Date(as.yearmon(PSID.pomp$time))
PSID.pomp$time <- as.Date(PSID.pomp$time, "%y%m%d")
psid.ts <- xts(PSID.pomp, order.by=as.Date(PSID.pomp$time))
inc <- xts(PSID.pomp$TotalIncome, order.by=as.Date(PSID.pomp$time))
inc <- as.vector(inc)
plot(inc)

educ <- xts(PSID.pomp$YearsCompletedEducation, order.by=as.Date(PSID.pomp$time))
educ <- as.vector(educ)
ccf(educ, inc)

race <- xts(PSID.pomp$Race2, order.by=as.Date(PSID.pomp$time))
race <- as.vector(race)
ccf(race, inc)

age <- xts(PSID.pomp$Age, order.by=as.Date(PSID.pomp$time))
age <- as.vector(age)
ccf(age, inc)


##GLM Regression. 
my.mod <- glm(log(TotalIncome) ~ Age + as.factor(EmploymentStatus) + as.factor(Race) + YearsCompletedEducation, data = PSIDdat)

summary(my.mod)
plot(my.mod$residuals)
qplot(my.mod$residuals, type = "l")
wald.test(b=coef(object=my.mod), Sigma=vcov(object=my.mod), Terms=6)

##HMM Model
PSIDdat <- PSID.pomp
PSIDdat$EmploymentStatus <- factor(PSIDdat$EmploymentStatus, 
                                 levels = c(0,1,2,3),
                                 labels = c("Other", "Employed", "Unemployed", "Retired"))
PSIDdat$Race <- factor(PSIDdat$Race,
                       levels = c(0, 1, 2, 3, 4, 5),
                       labels = c("Other", "White", "Black/African-American", "American Indian", "Asian", "Native Hawaiian/Islander"))
PSIDdat$OwnOrRent <- factor(PSIDdat$OwnOrRent,
                         levels = c(1, 5, 8, 9),
                         labels = c("Owns Home", "Rents", "Neither", "Other"))
PSIDdat$Age <- ifelse(PSIDdat$Age == 999, NA, PSIDdat$Age)
PSIDdat$YearsCompletedEducation <- ifelse(PSIDdat$YearsCompletedEducation == 99, NA, PSIDdat$YearsCompletedEducation)
PSIDdat <- na.omit(PSIDdat)
PSIDdat$time <- substring(PSIDdat$time,1, 4)
PSIDdat$time <- as.numeric(PSIDdat$time)

library(depmixS4)
set.seed(1234)
trst <- c(0.9, 0.1, 0, 0, 0.1, 0.9, 0, 0)
psid.ghmm <- depmix(TotalIncome ~ 1, 
                   data = PSIDdat, nstates = 2, instart = runif(2), 
                   ntimes = c(17033, 15767, 14324, 14185, 12857, 11078),
                   family = gaussian(),
                   transition = ~ Race + Age + YearsCompletedEducation)

fm <- fit(psid.ghmm, verbose = FALSE, emc=em.control(maxit = 500))
summary(fm, which = "response")
summary(fm, which = "transition")





low <- tmat[,2]
PSIDdat$low <- low
PSIDdat$time <- as.character(PSIDdat$time)
PSIDdat$time <- substring(PSIDdat$time,1, 4)
PSIDdat$time <- as.numeric(PSIDdat$time)
pdat <- data.frame(cbind(PSIDdat[,1], PSIDdat[,8], PSIDdat[,11]))
colnames(pdat) <- c("TotalIncome", "time", "low")
pdat <- ts(pdat)
pdat$time <- as.character(pdat$time)
as.Date()


#Check other hidden state values
  fmx <- fit(depmix(TotalIncome ~ 1, 
                   data = PSIDdat, nstates = 3, instart = runif(3), 
                   ntimes = c(17033, 15767, 14324, 14185, 12857, 11078),
                   family = gaussian(),
                   transition = ~ Race + Age + YearsCompletedEducation), verbose = FALSE, 
          emc=em.control(maxit = 500))
  
  fmxx <- fit(depmix(TotalIncome ~ 1, 
                    data = PSIDdat, nstates = 4, instart = runif(4), 
                    ntimes = c(17033, 15767, 14324, 14185, 12857, 11078),
                    family = gaussian(),
                    transition = ~ Race + Age + YearsCompletedEducation), 
              
              verbose = FALSE, 
             emc=em.control(maxit = 500))
  

  f <- fit(depmix(list(TotalIncome ~ 1, EmploymentStatus ~ 1),
                     data = PSIDdat, nstates = 2, instart = runif(2), 
                     ntimes = c(17033, 15767, 14324, 14185, 12857, 11078),
                     family = list(gaussian(), multinomial("identity"))), 
            verbose = FALSE, 
            emc=em.control(maxit = 500))
  
  ff <- fit(depmix(TotalIncome ~ 1, 
                   data = PSIDdat, nstates = 2, instart = runif(2), 
                   ntimes = c(17033, 15767, 14324, 14185, 12857, 11078),
                   family = gaussian(),
                   transition = ~ Race + Age + YearsCompletedEducation + EmploymentStatus), 
           verbose = FALSE, 
            emc=em.control(maxit = 500))
  
  
  fr <- fit(depmix(TotalIncome ~ 1,
                   data = PSIDdat, nstates = 3, instart = runif(3), 
                   ntimes = c(17033, 15767, 14324, 14185, 12857, 11078),
                   family = gaussian()), 
            verbose = FALSE, 
            emc=em.control(maxit = 500))
  
  fe <- depmix(TotalIncome ~ 1, 
               data = PSIDdat, nstates = 2, instart = runif(2), 
               ntimes = c(17033, 15767, 14324, 14185, 12857, 11078),
               family = gaussian(),
               transition = ~ Race + Age)
  
  fe.fit <- fit(fe, verbose = FALSE, emc=em.control(maxit = 500))
  
  fee <- fit(depmix(TotalIncome ~ 1, 
                   data = PSIDdat, nstates = 3, instart = runif(3), 
                   ntimes = c(17033, 15767, 14324, 14185, 12857, 11078),
                   family = gaussian(),
                   transition = ~ Race + Age), 
            verbose = FALSE, 
            emc=em.control(maxit = 500))
  
probs <- posterior(fmx)             # Compute probability of being in each state
head(probs)

# Lets change the name
colnames(probs)[2:6] <- paste("P",1:5, sep="-")
# Create dta.frame
dfu <- cbind(PSIDdat[,8], probs[,2:6])
# to Long format
dfu <- melt(dfu,id="time")

# Get the states values
stts<-getpars(fmx)[seq(43, by=2, length=5)]
# Change the names
names(stts) <- paste("St", (1:5), sep="-")
# Plot the data along with the time series of probabilities
qplot(time,value,data=dfu,geom="line",
      main = paste("States", paste(names(stts), stts, collapse=": "), collapse="; "),
      ylab = "State Probabilities") + 
  facet_grid(variable ~ ., scales="free_y") + theme_bw()


#Forward Backward Algorithm to estimate MLE vs iterated filtering
d <- forwardbackward(fm, return.all=TRUE)
d.sca <- d$sca
KL.dist(d.sca)
pam(d.sca)

#Results and graphs
c <- data.frame(States = c(2,3,4,2,2,3,2,3), 
                `Log Lik` = c(-922865.3, -934894, -934894, -1028984, -934894, -934894, -934894, -934894), 
                AIC = c(1845773, 1869900, 1870002, 2057994, 1869842, 1869816, 1869826, 1869888), 
                BIC = c(1845969, 1870424, 1871003, 2058115, 1870095, 1869947, 1870004, 1870356), 
                df = c(21, 56, 107, 13, 27, 14, 19, 50), 
                Response = c("Total Income", "Total Income", "Total Income", "Total Income and Employment Status", "Total Income", "Total Income",
                             "Total Income", "Total Income"),
                Covariates = c("Race, Age, Years of Education","Race, Age, Years of Education", "Race, Age, Years of Education", 
                                "None",  "Race, Age, Years of Education, Years of Employment", "None", "Race and Age", "Race and Age"))
income <- PSIDdat$TotalIncome
tim <- PSIDdat$time
post <- posterior(fm)
VitPath <- data.frame(post[,1])
d <- forwardbackward(fm, return.all=FALSE)

plo <- ggplot(PSIDdat, aes(TotalIncome)) + geom_histogram(bins = 30)
plo + facet_grid(time ~. )

PSIDdat$state <- post$state
state <- ggplot(post[10000:15000, ], aes(state, row.names(post)))
state

plot(VitPath[10000:15000, ], type='s', main='Implied States', xlab='', ylab='State', ylim = c(0, 3))


qplot(income,geom="line",
      main = "plot",
      ylab = "State Probabilities") + 
  facet_grid(time ~ ., scales="free_y") + theme_bw()





layout(1:3)
plot(income, tim)
geom_histogram()

fm.mod <- data.frame(Covariates = c("(Intercept)", "RaceWhite", "RaceBlack/African-American", "RaceAmerican Indian ", "RaceAsian", "RaceNative Hawaiian/Islander", 
                                    "Age", "YearsCompletedEducation", "Probalities at zero values of the covariates"),
                     `State 1.1` = c(0,0,0,0,0,0,0,0,0.8449153), 
                     `State 1.2` = c(-1.695265217, 0.345279128, -0.016206558, 0.428342941, -0.331647880, 0.095527697, 0.009451043, 0.005767248, 0.1550847),
                     `State 2.1` = c(0,0,0,0,0,0,0,0, 0.160579),
                     `State 2.2` = c(1.653926256, -0.161628929, -0.806133339, -0.658742499,0.399870158, -0.488206547, -0.006501024, -0.004946809, 0.839421))



results <- cbind(PSIDdat, post[,1])
setnames(results, "V2", "state")
results$state <- factor(results$state, 
                           levels = c(1,2),
                        labels = c("State 1", "State 2"))
results

res2 <- results
res2$state <- factor(res2$state, 
                        levels = c(1,2),
                        labels = c("State 1", "State 2"))

  
  
income.lev <- ggplot(res2, aes(TotalIncome)) + scale_fill_brewer(palette = "Spectral")
  
income.lev <- income.lev + geom_histogram(aes(fill = state),
                    bins=5, 
                    col="black", 
                    size=.1) +   # change number of bins
                    labs(title="State classification by income levels", 
                           subtitle="State 1 - Transitionary Poverty; State 2 - Permanent Poverty") 
  
 
age.lev <- ggplot(res2, aes(Age)) + scale_fill_brewer(palette = "Spectral")

age.lev <- age.lev + geom_histogram(aes(fill = state),
                                          bins=5, 
                                          col="black", 
                                          size=.1) +   # change number of bins
  labs(title="State classification by age levels", 
       subtitle="State 1 - Transitionary Poverty; State 2 - Permanent Poverty") 

age.lev

 ##RACE
race1 <- ggplot(subset(results, state == 1), aes(x = "", fill = factor(Race))) + 
  geom_bar(width = 1) +
  theme(axis.line = element_blank(), 
        plot.title = element_text(hjust=0.5)) + 
  labs(fill="Race", 
       x=NULL, 
       y=NULL, 
       title="State 1 - Permanent Poverty")

race1 <- race1 + coord_polar(theta = "y", start=0)


race2 <- ggplot(subset(results, state == 2), aes(x = "", fill = factor(Race))) + 
  geom_bar(width = 1) +
  theme(axis.line = element_blank(), 
        plot.title = element_text(hjust=0.5)) + 
  labs(fill="Race", 
       x=NULL, 
       y=NULL, 
       title="State 2 - Transitionary Poverty", 
       caption="Chart 1 - Race")

race2 <- race2 + coord_polar(theta = "y", start=0)
library(gridExtra)
grid.arrange(race1, race2)


##EMPLOYMENT
emp1 <- ggplot(subset(results, state == 1), aes(x = "", fill = factor(EmploymentStatus))) + 
  geom_bar(width = 1) +
  theme(axis.line = element_blank(), 
        plot.title = element_text(hjust=0.5)) + 
  labs(fill="Race", 
       x=NULL, 
       y=NULL, 
       title="State 1 - Permanent Poverty")

emp1 <- emp1 + coord_polar(theta = "y", start=0)


emp2 <- ggplot(subset(results, state == 2), aes(x = "", fill = factor(EmploymentStatus))) + 
  geom_bar(width = 1) +
  theme(axis.line = element_blank(), 
        plot.title = element_text(hjust=0.5)) + 
  labs(fill="Race", 
       x=NULL, 
       y=NULL, 
       title="State 2 - Transitionary Poverty", 
       caption="Chart 2 - Race")

emp2 <- emp2 + coord_polar(theta = "y", start=0)

grid.arrange(emp1, emp2)
