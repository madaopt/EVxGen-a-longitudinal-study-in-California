#####################################################################
# EV and gentrification, housing and rental price, longitudinal data 
# Futu Chen
# Last update: 2/13/2026
# Update notes:  Main analysis code, cleaned for publication; 
# Data available in the URL in citation
#
#####################################################################
# I LOVE PRESSURE #

library(lubridate)
library(dplyr)
library(tidyverse)
library(readxl)
library(stringr)
library(ggplot2)
library(ggtext)
library(lme4)
library(splines)
library(ggeffects)
library(gridExtra)
library(splines)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(ggpubr)
library(reshape2)
library(data.table)
library(lmerTest)
"%!in%" <- Negate("%in%")

# ======================
# Load Base Data
# ======================

dat <- evdat

# ======================
# 1. Gentrification Data
# ======================

sf   <- read.csv("displacement-typologies/data/downloads_for_public/sanfrancisco.csv")
scag <- read.csv("displacement-typologies/data/downloads_for_public/scag.csv")

gen <- bind_rows(sf, scag) %>%  mutate( Typology.f = factor(Typology),
                                        # Remove High Student Population & Unavailable or Unreliable Data
                                        Typology.f2 = case_when( Typology.f %in% c("High Student Population", "Unavailable or Unreliable Data") ~ NA_character_ ,
                                                                 TRUE ~ as.character(Typology.f)),
                                        Typology.f2 = factor( Typology.f2,
                                                              levels = c(
                                                                "Low-Income/Susceptible to Displacement",
                                                                "Ongoing Displacement",
                                                                "At Risk of Gentrification",
                                                                "Early/Ongoing Gentrification",
                                                                "Advanced Gentrification",
                                                                "Stable Moderate/Mixed Income",
                                                                "At Risk of Becoming Exclusive",
                                                                "Becoming Exclusive",
                                                                "Stable/Advanced Exclusive"
                                                              ), ordered = TRUE),
                                        
                                        # Try 3 categories
                                        Typology.f3 = case_when( Typology.f2 %in% c( "Low-Income/Susceptible to Displacement", "Ongoing Displacement", "At Risk of Gentrification", "Early/Ongoing Gentrification") ~ "Ongoing",
                                                                 Typology.f2 %in% c( "Advanced Gentrification", "Stable Moderate/Mixed Income" ) ~ "Stable",
                                                                 Typology.f2 %in% c( "At Risk of Becoming Exclusive","Becoming Exclusive","Stable/Advanced Exclusive" ) ~ "Exclusive",
                                                                 TRUE ~ NA_character_),
                                        # use my label
                                        Typology.f3 = factor( Typology.f3, levels = c("Ongoing", "Stable", "Exclusive"),labels = c("None/early", "Advanced/stable", "Exclusive")
                                        )
)

# ======================
# 2. Merge EV + Gentrification
# ======================

dat2 <- dat %>%
  left_join(gen, by = c("GEOID_CT" = "GEOID"))

# Subset to more than 300 pop #
dat.time.full <- dat2[ dat2$tpop >= 300,]
nrow(dat.time.full)/nrow(dat2)*100 #99.05

# ======================
# 3. Add DAC
# ======================

ap <- read.csv("CalEnviroScreen 4.0/calenviroscreen4.csv") 
ap <- ap[,c("Census.Tract", "PM2.5","CES.4.0.Percentile" )]
colnames(ap)[1] <- "GEOID_CT"
dat.time.full <- left_join(dat.time.full, ap, by="GEOID_CT")
dat.time.full$CES.bin <- factor(ifelse(dat.time.full$CES.4.0.Percentile >= 75, 1, 0), levels = c(0,1), labels = c("Non DAC","DAC"))
table(dat.time.full$CES.bin)

# ======================
# 4. Prepare Variables for Analysis
# ======================

dat.time.full <- dat.time.full %>%
  mutate(
    Year_2015 = Year - 2015, # center year
    bev_1000pop  = BEV  / tpop * 1000,
    phev_1000pop = PHEV / tpop * 1000,
    ice_1000pop  = ICE  / tpop * 1000,
    Year.f = factor(Year) # for plots
  )

# =====================================================
# Analysis 1
# =====================================================

library(splines)
library(lmerTest)

# Model Selection: Optimal Degrees of Freedom
# ----------------------------------
# BEV: Select optimal spline df
# ----------------------------------

temp <- NULL

for (i in c(1:8) ) {
  md <- lmer(bev_1000pop ~    ns(Year_2015, df=i) * CES.bin * Typology.f3  +  (1 | GEOID_CT), data = dat.time.full, REML=F)
  tb <- c(BIC(md), AIC(md), i)
  temp <- as.data.frame(rbind(temp,tb))
}

colnames(temp) <- c("BIC","AIC","Df")
temp[temp$BIC == min(temp$BIC),] #df=7
temp[temp$AIC == min(temp$AIC),] #df=8

# Decision: Use df=7 (df=8 likely overfitting)

# ----------------------------------
# PHEV: Select optimal spline df
# ----------------------------------

temp <- NULL
for (i in c(1:8) ) {
  md <- lmer(phev_1000pop ~    ns(Year_2015, df=i) * CES.bin * Typology.f3  +  (1 | GEOID_CT), data = dat.time.full, REML=F)
  tb <- c(BIC(md), AIC(md), i)
  temp <- as.data.frame(rbind(temp,tb))
}
colnames(temp) <- c("BIC","AIC","Df")
temp[temp$BIC == min(temp$BIC),] #df=7
temp[temp$AIC == min(temp$AIC),] #df=7

# Decision: Use df=7

# =====================================================
# Descriptive Statistics for Writing
# =====================================================

# BEV
summary(dat.time.full[dat.time.full$Year == 2015, ]$bev_1000pop)
summary(dat.time.full[dat.time.full$Year == 2023, ]$bev_1000pop)

# PHEV
summary(dat.time.full[dat.time.full$Year == 2015, ]$phev_1000pop)
summary(dat.time.full[dat.time.full$Year == 2023, ]$phev_1000pop)

# =====================================================
# Supplemental Figure: ZEV Trends
# =====================================================
md.b2  <- lmer(bev_1000pop ~    ns(Year_2015, df=7)   +  (1 | GEOID_CT), data = dat.time.full) 
pr.b2 <-  predict_response(md.b2 , terms =c("Year_2015 [n=50]"), type = "fixed")
pr.b2$x <- pr.b2$x+2015
pr.b2$group <- "BEV"

md.p2  <- lmer(phev_1000pop ~    ns(Year_2015, df=7) + (1 | GEOID_CT), data = dat.time.full)
pr.p2 <- predict_response(md.p2 , c("Year_2015 [n=50]"), type = "fixed")
pr.p2$x <- pr.p2$x+2015
pr.p2$group <- "PHEV"

plotdat <- rbind(pr.b2,pr.p2)

plot(plotdat) + scale_y_continuous(limits = c(0,25)) + scale_x_continuous(breaks=c(2015:2023))+
  ylab("Predicted nZEV/1000 pop") + ggtitle("") + scale_color_discrete(name = "ZEV type")+
  scale_color_manual(values=c("#1A85FF","#D41159"),name = "ZEV type")+
  scale_fill_manual(values=c("#1A85FF","#D41159"),name = "ZEV type")+
  xlab("Year") +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.position = "bottom",
        axis.text=element_text(size=14,angle = 45,  vjust = 0.5, hjust=1),
        axis.title.y.left=element_text(size=14, face="bold", vjust = 2),
        axis.title.x=element_text(size=14, face="bold", vjust = 0),
        plot.title = element_text(face="bold", size = 10),
        strip.background = element_rect(fill="gray"),
        strip.text.x = element_text(size = 14, color = "black", face = "bold"))

# =====================================================
# Analysis 1: Fit Full Interaction Models
# =====================================================

# ----------------------------------
# BEV Model
# ----------------------------------

# BEV #
md.b  <- lmer(bev_1000pop ~    ns(Year_2015, df=7) * CES.bin * Typology.f3  +  (1 | GEOID_CT), data = dat.time.full, REML=T)
pr.b <-  predict_response(md.b , terms =c("Year_2015 [n=50]","CES.bin","Typology.f3"), type = "fixed")
colnames(pr.b) <- c("x" , "predicted" ,"std.error", "conf.low",  "conf.high",  "facet","group"  )
pr.b$x <- pr.b$x+2015
p1 <- plot(pr.b) + scale_y_continuous(limits = c(-3,45)) + scale_x_continuous(breaks=c(2015:2023))+
  ylab("Estimated nBEV/1000 pop") + ggtitle("") + scale_color_discrete(name = "Gentrification status")+
  scale_color_manual(values=c("#648FFF","#DC267F","#FFB000"),name = "Gentrification status")+
  scale_fill_manual(values=c("#648FFF","#DC267F","#FFB000"),name = "Gentrification status")+
  xlab("Year") +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.position = "bottom",
        axis.text=element_text(size=14,angle = 45,  vjust = 0.5, hjust=1),
        axis.title.y.left=element_text(size=14, face="bold", vjust = 2),
        axis.title.x=element_text(size=14, face="bold", vjust = 0),
        plot.title = element_text(face="bold", size = 10),
        strip.background = element_rect(fill="gray"),
        strip.text.x = element_text(size = 14, color = "black", face = "bold")) 

# PHEV #
md.p  <- lmer(phev_1000pop ~    ns(Year_2015, df=7) * CES.bin * Typology.f3  + (1 | GEOID_CT), data = dat.time.full)
pr.p <- predict_response(md.p , c("Year_2015 [n=50]","CES.bin","Typology.f3"), type = "fixed")
colnames(pr.p) <- c("x" , "predicted" ,"std.error", "conf.low",  "conf.high",  "facet","group"  )
pr.p$x <- pr.p$x+2015
p2 <-  plot(pr.p) + scale_y_continuous(limits = c(-3,45)) + scale_x_continuous(breaks=c(2015:2023))+
  ylab("Estimated nPHEV/1000 pop") + ggtitle("") + scale_color_discrete(name = "Gentrification status")+
  scale_color_manual(values=c("#648FFF","#DC267F","#FFB000"),name = "Gentrification status")+
  scale_fill_manual(values=c("#648FFF","#DC267F","#FFB000"),name = "Gentrification status")+
  xlab("Year") +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.position = "bottom",
        axis.text=element_text(size=14,angle = 45,  vjust = 0.5, hjust=1),
        axis.title.y.left=element_text(size=14, face="bold", vjust = 2),
        axis.title.x=element_text(size=14, face="bold", vjust = 0),
        plot.title = element_text(face="bold", size = 10),
        strip.background = element_rect(fill="gray"),
        strip.text.x = element_text(size = 14, color = "black", face = "bold")) 

plot <- ggarrange(p1, p2,nrow=1, common.legend = TRUE, legend="bottom")

# =====================================================
# Analysis 2: Housing Analysis
# =====================================================

# Subset to none/early neighborhoods

dat.non <- dat.time.full[dat.time.full $Typology.f3 %in% "None/early",] #17763
nrow(gen[gen$Typology.f3 %in% "None/early",]) / nrow(gen) #1975, 29.65%

# ----------------------------------
# Visualize housing trend by change in EV, Figure 2
# ----------------------------------

# Among non/early gentrified #
housing <- dat.non[,c("GEOID_CT","Year","median_value","median_rent")] # data has been crosswalked from 2020 geo to 2010 geo, make sure housing data is crosswalked! 

# Making categorical ZEV predictors with 4 levels #
delta19_15 <- dat.non[,c("GEOID_CT" , "Year", "bev_1000pop","phev_1000pop","ice_1000pop")] %>% 
  filter(Year == 2015 | Year == 2019) %>% 
  setDT(.) %>% 
  data.table::dcast(GEOID_CT ~ Year, value.var = c("bev_1000pop","phev_1000pop","ice_1000pop")) %>%
  mutate( delta_bev1915= bev_1000pop_2019 - bev_1000pop_2015,
          delta_phev1915= phev_1000pop_2019 - phev_1000pop_2015,
          delta_ice1915= ice_1000pop_2019 - ice_1000pop_2015)%>%
  dplyr::select(c("GEOID_CT", "delta_bev1915","delta_phev1915","delta_ice1915")) %>% 
  mutate(
    delta_bev1915_cat= cut(delta_bev1915, unique(quantile(delta_bev1915,seq(0,1,.25), na.rm=T,include.lowest=TRUE, label = F))),
    delta_phev1915_cat= cut(delta_phev1915, unique(quantile(delta_phev1915,seq(0,1,.25), na.rm=T,include.lowest=TRUE, label = F))),
    delta_ice1915_cat= cut(delta_ice1915, unique(quantile(delta_ice1915,seq(0,1,.25), na.rm=T,include.lowest=TRUE, label = F))),
    
    delta_bev1915_cat.f = factor(delta_bev1915_cat, labels = c("Q1","Q2","Q3","Q4")),
    delta_phev1915_cat.f = factor(delta_phev1915_cat, labels = c("Q1","Q2","Q3","Q4")),
    delta_ice1915_cat.f = factor(delta_ice1915_cat, labels = c("Q1","Q2","Q3","Q4")),
    
  )
table(delta19_15$delta_bev1915_cat)
prop.table(table(delta19_15$delta_bev1915_cat.f))*100
table(delta19_15$delta_phev1915_cat)
prop.table(table(delta19_15$delta_phev1915_cat.f))*100

names(delta19_15)

housing <- dat.non[,c("GEOID_CT","Year","median_value","median_rent")]
housing <- left_join(housing,delta19_15, by="GEOID_CT")

housing.wide <- housing[,c("GEOID_CT" ,"Year" , "median_value"  , "median_rent" ,"delta_bev1915_cat.f" , "delta_phev1915_cat.f" ,"delta_ice1915_cat.f")]
housing.long <- gather(housing.wide, evtype, quntile, delta_bev1915_cat.f:delta_ice1915_cat.f, factor_key=TRUE)
levels(housing.long$evtype) <- c("BEV", "PHEV", "ICE")
housing.long$quntile <- as.factor(housing.long$quntile)
levels(housing.long$quntile) <- c("Q1","Q2","Q3","Q4")

# set up color scheme 
mycolor = c("#E5A81A","#E51ABD","#1A57E5","#009E73") 
# highlight the time frame 
rect_df <- data.frame(
  xmin = "2020",
  xmax = "2023"
)

p1 <- housing.long %>% filter( !is.na(quntile)) %>% filter( evtype %!in% "ICE") %>% 
  ggplot( aes(x=as.factor(Year), y=median_value, color = quntile))+
  geom_smooth(aes(group=quntile), method = "loess")+
  scale_color_manual(values=mycolor,name = "Delta nZEV/1000pop \n 2019-2015")+
  geom_rect(data = rect_df,aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            color = "gray60",  linetype = "dashed", fill = NA, linewidth = 1,  inherit.aes = FALSE ) +
  ylab("Median home value")+
  xlab("Year")+
  facet_grid(~evtype)+
  theme_bw()+
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        # legend.position = "bottom",
        axis.text=element_text(size=14,angle = 45,  vjust = 0.5, hjust=1),
        axis.title.y.left=element_text(size=14, face="bold", vjust = 2),
        axis.title.x=element_text(size=14, face="bold", vjust = 0),
        plot.title = element_text(face="bold", size = 10),
        strip.background = element_rect(fill="gray"),
        strip.text.x = element_text(size = 14, color = "black", face = "bold")) 


p2 <- housing.long %>% filter( !is.na(quntile)) %>%  filter( evtype %!in% "ICE") %>% 
  ggplot( aes(x=as.factor(Year), y=median_rent, color = quntile))+
  geom_smooth(aes(group=quntile), method = "loess")+
  scale_color_manual(values=mycolor,name = "Delta nZEV/1000pop \n 2019-2015")+
  geom_rect(data = rect_df,aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            color = "gray60",  linetype = "dashed", fill = NA, linewidth = 1,  inherit.aes = FALSE ) +
  ylab("Median rental price")+
  xlab("Year")+
  facet_grid(~evtype)+
  theme_bw()+
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        # legend.position = "bottom",
        axis.text=element_text(size=14,angle =45,  vjust = 0.5, hjust=1),
        axis.title.y.left=element_text(size=14, face="bold", vjust = 2),
        axis.title.x=element_text(size=14, face="bold", vjust = 0),
        plot.title = element_text(face="bold", size = 10),
        strip.background = element_rect(fill="gray"),
        strip.text.x = element_text(size = 14, color = "black", face = "bold")) 

ggarrange(p1,p2,nrow=2, common.legend = T, legend = "bottom")

# ----------------------------------
# Model housing trend, main analysis 2
# ----------------------------------
# Get cov needed for model
cov <- dat.non[,c("GEOID_CT","Year","pct_Bachelor","median_income","Typology.f2","tpop")]
housing <- dat.non[,c("GEOID_CT","Year","median_value","median_rent","tpop")]# none/early
summary(housing) # 969 NA median value, 86NA median rent 

housing <- left_join(housing,delta19_15, by="GEOID_CT")

# Subset to 2020 and beyond #
housing.2020 <- housing[housing$Year >=2020,]
# Center year
housing.2020$year_2020 <-  housing.2020$Year - 2020
housing.2020 <- left_join(housing.2020, cov, by= c("GEOID_CT","Year"))

# Consider linear based on previous visual and limited number of years

# log outcome
housing.2020$log_median_value=log(housing.2020$median_value)
housing.2020$log_median_rent=log(housing.2020$median_rent)

hist(housing.2020$log_median_value)
hist(housing.2020$log_median_rent)

# Function to summarize the interaction results #
# avg annual change in  Y for cat veh 
mdint.summary.fun <- function(md){
  temp <- as.data.frame(coef(summary(md)))
  #est avg annual change in  Y for cat bev, beta hat 
  ind.q2=grep("year_2020",rownames(temp))[2]
  ind.q3=grep("year_2020",rownames(temp))[3]
  ind.q4=grep("year_2020",rownames(temp))[4]
  est.q1 <- coef(summary(md))[2,1]
  est.q2 <- coef(summary(md))[2,1] + coef(summary(md))[ind.q2,1]
  est.q3 <- coef(summary(md))[2,1] + coef(summary(md))[ind.q3,1]
  est.q4 <- coef(summary(md))[2,1]  + coef(summary(md))[ind.q4,1]
  #var hat
  varhat.q1 <- vcov(summary(md))[2,2]
  varhat.q2 <- vcov(summary(md))[2,2] + vcov(summary(md))[ind.q2,ind.q2] + 2*vcov(md)[ind.q2,2]
  varhat.q3 <- vcov(summary(md))[2,2] + vcov(summary(md))[ind.q3,ind.q3] + 2*vcov(md)[ind.q3,2]
  varhat.q4 <- vcov(summary(md))[2,2] + vcov(summary(md))[ind.q4, ind.q4] + 2*vcov(md)[ind.q4,2]
  # CI
  estq1_ci_low <- est.q1-1.96*sqrt(varhat.q1)
  estq1_ci_up <- est.q1+1.96*sqrt(varhat.q1)
  estq2_ci_low <- est.q2-1.96*sqrt(varhat.q2)
  estq2_ci_up <- est.q2+1.96*sqrt(varhat.q2)
  estq3_ci_low <- est.q3-1.96*sqrt(varhat.q3)
  estq3_ci_up <- est.q3+1.96*sqrt(varhat.q3)
  estq4_ci_low <- est.q4-1.96*sqrt(varhat.q4)
  estq4_ci_up <- est.q4+1.96*sqrt(varhat.q4)
  # trans est
  est.q1_t <- (exp( est.q1)-1)*100 #% change 
  est.q2_t <- (exp( est.q2)-1)*100 #% change 
  est.q3_t <- (exp( est.q3)-1)*100 #% change 
  est.q4_t <- (exp( est.q4)-1)*100 #% change 
  
  # trans ci
  estq1_ci_low_t <- (exp( estq1_ci_low)-1)*100
  estq1_ci_up_t <-  (exp( estq1_ci_up)-1)*100
  estq2_ci_low_t <- (exp( estq2_ci_low)-1)*100
  estq2_ci_up_t <-  (exp( estq2_ci_up)-1)*100
  estq3_ci_low_t <- (exp( estq3_ci_low)-1)*100
  estq3_ci_up_t <-  (exp( estq3_ci_up)-1)*100
  estq4_ci_low_t <- (exp( estq4_ci_low)-1)*100
  estq4_ci_up_t <-  (exp( estq4_ci_up)-1)*100
  
  # others
  pval.q2 <- summary(md)$coefficients[ind.q2,"Pr(>|t|)"]
  pval.q3 <- summary(md)$coefficients[ind.q3,"Pr(>|t|)"]
  pval.q4 <- summary(md)$coefficients[ind.q4,"Pr(>|t|)"]
  nCT <- as.numeric(summary(md)$ngrps[1])
  nobs <- nobs(md)
  
  out <- c(nCT=nCT, nobs=nobs,est.q1=est.q1_t, est.q2=est.q2_t,est.q3=est.q3_t,est.q4=est.q4_t,
           estq1_ci_low=round(estq1_ci_low_t,2), estq1_ci_up= round(estq1_ci_up_t,2), 
           estq2_ci_low=round(estq2_ci_low_t,2), estq2_ci_up= round(estq2_ci_up_t,2), 
           estq3_ci_low=round(estq3_ci_low_t,2), estq3_ci_up= round(estq3_ci_up_t,2), 
           estq4_ci_low=round(estq4_ci_low_t,2), estq4_ci_up= round(estq4_ci_up_t,2), 
           p.q2=round(pval.q2,6),p.q3=round(pval.q3,6),p.q4=round(pval.q4,6))
  out
}

library(ggeffects)

# Prep inputs

ev.list <- c("delta_bev1915_cat.f","delta_phev1915_cat.f","delta_ice1915_cat.f")
var.list <- c("value","rent")
random.slope <- c("","+year_2020")
cov.list <- c("", "+scale(tpop.x)","+scale(pct_Bachelor)","+scale(median_income)","+scale(pct_Bachelor)+scale(median_income)")

# ----------------------------------
# Analysis 2: results on the abs scale
# ----------------------------------
# ----------------------------------
# Run models 
# ----------------------------------

### Annual average change ###
md.holder <- NULL
md.name <- NULL
dat.list <- housing.2020
# for (d in 1:length(dat.list)) {
for (i in 1:length(ev.list)) {
  for (j in 1:length(var.list)) {
    for (k in 1:length(random.slope)) {
      for (l in 1:length(cov.list)) {
        md <-  lmer(as.formula(paste0("log_median_",var.list[j], "~", "year_2020*",ev.list[i], cov.list[l],"+ (1",random.slope[k],"| GEOID_CT)")), dat.list) 
        md.holder <- c(md.holder,list(md))
        name <- paste0("log_median_",var.list[j], "~", "year_2020*",ev.list[i], cov.list[l],"+ (1",random.slope[k],"| GEOID_CT)","dat=")
        md.name <- c(md.name, name)
      } #l, cov
    } #k, random slope 
  }# j, var y 
} #i ev
# }#d in datlist 
length(md.holder)

# Get estimation: % increase 
result.holder.int <- NULL
for (i in 1:length(md.holder)) {
  temp <- lapply(md.holder[i], mdint.summary.fun) %>% bind_rows()
  result.holder.int  <- rbind(result.holder.int ,temp)
}

result.holder.int$cov <- rep(c("no cov","pop","edu","income","edu+income"),12)
result.holder.int$random <- rep(rep( c("no","year"), each=5),6)
result.holder.int$outcome <- rep( rep(c("value","rent"), each=10),3)
result.holder.int$ev <- rep(c("bev","phev","ice"), each=20)

# Visual #
result.holder.int$cov.f <- factor(result.holder.int$cov , levels = c("no cov", "pop","edu","income","edu+income"), labels = c("no cov","pop","edu","income","edu+income"))
result.holder.int$ev.f <- factor(result.holder.int$ev , levels = c("bev","phev","ice"), labels = c("BEV","PHEV","ICE"))
result.holder.int$outcome.f <- factor(result.holder.int$outcome , levels = c("value","rent"), labels = c("Home value","Rental price"))
result.holder.int$random.f <- factor(result.holder.int$random, levels = c("no", "year"), labels = c("no random slope", "random slope year"))

# Code significance level for ggplot 
result.holder.int$sig.p.q2 <- ifelse(result.holder.int$p.q2 < 0.001, "***", ifelse(result.holder.int$p.q2 < 0.01, "**",ifelse(result.holder.int$p.q2 < 0.05, "*",NA)))
result.holder.int$sig.p.q3 <- ifelse(result.holder.int$p.q3 < 0.001, "***", ifelse(result.holder.int$p.q3 < 0.01, "**",ifelse(result.holder.int$p.q3 < 0.05, "*",NA)))
result.holder.int$sig.p.q4 <- ifelse(result.holder.int$p.q4 < 0.001, "***", ifelse(result.holder.int$p.q4 < 0.01, "**",ifelse(result.holder.int$p.q4 < 0.05, "*",NA)))

# --------------------------------------------------
# Convert  results from wide to long; can be more efficient lol
# --------------------------------------------------


tempq1 <- result.holder.int[,c("est.q1", "estq1_ci_low" ,"estq1_ci_up","ev.f","outcome.f", "random.f" ,"cov.f")]
tempq1$sig.p <- NA
tempq2 <- result.holder.int[,c("est.q2", "estq2_ci_low" ,"estq2_ci_up","ev.f","outcome.f", "random.f" ,"cov.f","sig.p.q2")]
tempq3 <- result.holder.int[,c("est.q3", "estq3_ci_low" ,"estq3_ci_up","ev.f","outcome.f", "random.f" ,"cov.f","sig.p.q3")]
tempq4 <- result.holder.int[,c("est.q4", "estq4_ci_low" ,"estq4_ci_up","ev.f","outcome.f", "random.f" ,"cov.f","sig.p.q4")]

colnames(tempq1 ) <- c("est", "est_ci_low" ,"est_ci_up","ev.f" ,"outcome.f", "random.f" ,"cov.f","sig.p")
colnames(tempq2 ) <- c("est", "est_ci_low" ,"est_ci_up", "ev.f" ,"outcome.f", "random.f" ,"cov.f"  ,"sig.p")
colnames(tempq3 ) <- c("est", "est_ci_low" ,"est_ci_up", "ev.f" ,"outcome.f", "random.f" ,"cov.f"  , "sig.p")
colnames(tempq4 ) <- c("est", "est_ci_low" ,"est_ci_up", "ev.f" ,"outcome.f", "random.f" ,"cov.f" ,"sig.p")

tempq1$quartile <- "Q1"
tempq2$quartile <- "Q2"
tempq3$quartile <- "Q3"
tempq4$quartile <- "Q4"

result.holder.int2 <- rbind(tempq1,tempq2, tempq3, tempq4)
result.holder.int2$quartile.f <- factor(result.holder.int2$quartile)

# --------------------------------------------------
# Not get baseline $ value for each result 
# --------------------------------------------------

# function to summarize 
baseline.sum.fun <- function(md, ev){
  # Predicted log values at year = 0 for each category
  baseline_preds <- ggpredict(md, terms = c(ev, "year_2020 [0]"))
  # Calculate smearing factor
  resid_vals <- resid(md)
  smearing_factor <- mean(exp(resid_vals))
  # back trans to get oroginal
  baseline_preds$predicted_dollar <- exp(baseline_preds$predicted) * smearing_factor
  baseline_preds$ev <- ev
  baseline_preds
}

# get baseline estimation

ev.list <-  rep(c("delta_bev1915_cat.f","delta_phev1915_cat.f","delta_ice1915_cat.f"), each=20)
cov.list <-  rep(c("no cov","pop","edu","income","edu+income"),12)
random.list <- rep(rep( c("no","year"), each=5),6)
outcome.list <-  rep( rep(c("value","rent"), each=10),3)
evname.list <- rep(c("bev","phev","ice"), each=20)

result.holder.baseline <- NULL
for (i in 1:length(md.holder)) {
  md <- md.holder[[i]]
  baseline_preds <- ggpredict(md, terms = c(ev.list[i], "year_2020 [0]"))
  # Calculate smearing factor
  resid_vals <- resid(md)
  smearing_factor <- mean(exp(resid_vals))
  # back trans to get original
  baseline_preds$predicted_dollar <- exp(baseline_preds$predicted) * smearing_factor
  baseline_preds <-  baseline_preds[,c("x", "predicted_dollar")]
  #add labels
  baseline_preds$cov <- cov.list[i]
  baseline_preds$random <- random.list[i]
  baseline_preds$outcome <- outcome.list[i]
  baseline_preds$ev <- evname.list[i]
  # bind 
  result.holder.baseline <- rbind(result.holder.baseline,baseline_preds)
}
result.holder.baseline <- as.data.frame(result.holder.baseline)

result.holder.baseline$cov.f <- factor(result.holder.baseline$cov , levels = c("no cov","pop","edu","income","edu+income"), labels = c("no cov","pop","edu","income","edu+income"))
result.holder.baseline$ev.f <- factor(result.holder.baseline$ev , levels = c("bev","phev","ice"), labels = c("BEV","PHEV","ICE"))
result.holder.baseline$outcome.f <- factor(result.holder.baseline$outcome , levels = c("value","rent"), labels = c("Home value","Rental price"))
result.holder.baseline$random.f <- factor(result.holder.baseline$random, levels = c("no", "year"), labels = c("no random slope", "random slope year"))

# join to % increase
result.holder.int3 <- left_join(result.holder.int2,result.holder.baseline, by=c("quartile", "ev.f"  ,     "outcome.f" , "random.f" ,  "cov.f" ) )

# --------------------------------------------------
# Calc annual avg change in abs value  
# --------------------------------------------------

n_years <- 3  # from year 0 to year 3 → 3 years
result.holder.int3 <- result.holder.int3 %>%
  mutate(
    abs_change_mean  = predicted_dollar * ((1 + est/100)^n_years - 1) / n_years,
    abs_change_lower = predicted_dollar * ((1 + est_ci_low/100)^n_years - 1) / n_years,
    abs_change_upper = predicted_dollar * ((1 + est_ci_up/100)^n_years - 1) / n_years
  )

# --------------------------------------------------
# Figure 3
# --------------------------------------------------
result.holder.int3 <- result.holder.int3 %>%
  mutate(sig.p = ifelse(is.na(sig.p), "", sig.p))

mycolor = c("#E5A81A","#E51ABD","#1A57E5","#009E73")

result.holder.int3$cov.f2 <- result.holder.int3$cov.f # use new labels in response to reviewers
levels(result.holder.int3$cov.f2) <- c("Unadjusted","Population adjusted","Education adjusted","Income adjusted","Edu + Income adjusted")

result.holder.int3 %>%
  filter(random.f %in% "random slope year") %>%
  filter(ev.f %!in% "ICE") %>%
  ggplot(aes(x = abs_change_mean, y = cov.f2,color=quartile.f)) + ##
  #geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = abs_change_upper, xmin = abs_change_lower , group=quartile.f), size = .5, height = 
                   .2, color = "gray50", position=position_dodge(width=0.5)) +
  geom_point(aes(shape=quartile.f),size = 3.5,position=position_dodge(width=0.5)) +
  geom_point(aes(x=abs_change_mean, shape= quartile.f, group=quartile.f),size = 3.5,position=position_dodge(width=0.5)) +
  theme_bw()+
  theme(panel.grid.minor = element_blank()) +
  ylab("") +
  xlab("Annual change in housing price, US Dollar $") +
  coord_flip()+
  scale_color_discrete(name = "Delta change in veh/1000pop \n2015-2019")+
  scale_shape_discrete(name = "Delta change in veh/1000pop \n2015-2019")+
  geom_text(aes(label =  sig.p, group=quartile.f,color=quartile.f), vjust = 0.4, hjust = -3, position=position_dodge(width=0.6),show.legend = F)+
  # geom_text(aes(label= paste0(round(abs_change_mean,0)), group=as.factor(year)), vjust = 0, hjust = 1.5, color="black", position=position_dodge(width=0.5))+
  #facet_wrap(outcome.f~ev.f,scales="free_x")+ #random.f
  facet_grid(outcome.f ~ ev.f, scales = "free_y")+
  #facet_grid(ev.f~outcome.f, scales = "free_x")+
  # ggtitle("Annual average dollar change for each category of EV") +
  scale_color_manual(name = "Quartiles of Change in nZEV/1000pop (2019-2015)", values = mycolor)+
  scale_shape_discrete(name = "Quartiles of Change in nZEV/1000pop (2019-2015)")+
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.position = "bottom",
        axis.text=element_text(size=14,angle = 0,  vjust = 0.5, hjust=1),
        axis.title.y.left=element_text(size=14, face="bold", vjust = 2),
        axis.title.x=element_text(size=14, face="bold", vjust = 0),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(face="bold", size = 10),
        strip.background = element_rect(fill="gray"),
        strip.text.x = element_text(size = 14, color = "black", face = "bold"),
        strip.text.y = element_text(size = 14, color = "black", face = "bold"))


# ----------------------------------
# Analysis 2: results on the relative scale
# Reported as a sup table
# ----------------------------------

# list of data to loop on
ev.list <- c("delta_bev1915_cat.f","delta_phev1915_cat.f","delta_ice1915_cat.f")
var.list <- c("value","rent")
random.slope <- c("","+year_2020")
cov.list <- c("", "+scale(tpop.x)","+scale(pct_Bachelor)","+scale(median_income)","+scale(pct_Bachelor)+scale(median_income)")


md.holder <- NULL
md.name <- NULL
dat.list <- housing.2020


for (i in 1:length(ev.list)) {
  for (j in 1:length(var.list)) {
    for (k in 1:length(random.slope)) {
      for (l in 1:length(cov.list)) {
        md <-  lmer(as.formula(paste0("log_median_",var.list[j], "~", "year_2020*",ev.list[i], cov.list[l],"+ (1",random.slope[k],"| GEOID_CT)")), dat.list) 
        md.holder <- c(md.holder,list(md))
        name <- paste0("log_median_",var.list[j], "~", "year_2020*",ev.list[i], cov.list[l],"+ (1",random.slope[k],"| GEOID_CT)","dat=",dat.list)
        md.name <- c(md.name, name)
      } #l, cov
    } #k, random slope 
  }# j, var y 
} #i ev

length(md.holder)

# extract #
# 60 models # 
result.holder.int <- NULL
for (i in 1:length(md.holder)) {
  temp <- lapply(md.holder[i], mdint.summary.fun) %>% bind_rows()
  result.holder.int  <- rbind(result.holder.int ,temp)
}

result.holder.int$cov <- rep(c("no cov","pop","edu","income","edu+income"),12)
result.holder.int$random <- rep(rep( c("no","year"), each=5),6)
result.holder.int$outcome <- rep( rep(c("value","rent"), each=10),3)
result.holder.int$ev <- rep(c("bev","phev","ice"), each=20)


result.holder.int[16,c("cov","random","ev","outcome")]
summary(md.holder[[16]])
result.holder.int[29,c("cov","random","ev","outcome")]
summary(md.holder[[29]])

# visual #
result.holder.int$cov.f <- factor(result.holder.int$cov , levels = c("no cov", "pop","edu","income","edu+income"), labels = c("no cov","pop","edu","income","edu+income"))
result.holder.int$ev.f <- factor(result.holder.int$ev , levels = c("bev","phev","ice"), labels = c("BEV","PHEV","ICE"))
result.holder.int$outcome.f <- factor(result.holder.int$outcome , levels = c("value","rent"), labels = c("Home value","Rental price"))
result.holder.int$random.f <- factor(result.holder.int$random, levels = c("no", "year"), labels = c("no random slope", "random slope year"))

result.holder.int$sig.p.q2 <- ifelse(result.holder.int$p.q2 < 0.001, "**", ifelse(result.holder.int$p.q2 < 0.05,"*",NA))
result.holder.int$sig.p.q3 <- ifelse(result.holder.int$p.q3 < 0.001, "**", ifelse(result.holder.int$p.q3 < 0.05,"*",NA))
result.holder.int$sig.p.q4 <- ifelse(result.holder.int$p.q4 < 0.001, "**", ifelse(result.holder.int$p.q4 < 0.05,"*",NA))
#wide to long; dump idea, no gather()
tempq1 <- result.holder.int[,c("est.q1", "estq1_ci_low" ,"estq1_ci_up","ev.f","outcome.f", "random.f" ,"cov.f")]
tempq1$sig.p <- NA
tempq2 <- result.holder.int[,c("est.q2", "estq2_ci_low" ,"estq2_ci_up","ev.f","outcome.f", "random.f" ,"cov.f","sig.p.q2")]
tempq3 <- result.holder.int[,c("est.q3", "estq3_ci_low" ,"estq3_ci_up","ev.f","outcome.f", "random.f" ,"cov.f","sig.p.q3")]
tempq4 <- result.holder.int[,c("est.q4", "estq4_ci_low" ,"estq4_ci_up","ev.f","outcome.f", "random.f" ,"cov.f","sig.p.q4")]

colnames(tempq1 ) <- c("est", "est_ci_low" ,"est_ci_up","ev.f" ,"outcome.f", "random.f" ,"cov.f","sig.p")
colnames(tempq2 ) <- c("est", "est_ci_low" ,"est_ci_up", "ev.f" ,"outcome.f", "random.f" ,"cov.f"  ,"sig.p")
colnames(tempq3 ) <- c("est", "est_ci_low" ,"est_ci_up", "ev.f" ,"outcome.f", "random.f" ,"cov.f"  , "sig.p")
colnames(tempq4 ) <- c("est", "est_ci_low" ,"est_ci_up", "ev.f" ,"outcome.f", "random.f" ,"cov.f" ,"sig.p")

tempq1$quartile <- "Q1"
tempq2$quartile <- "Q2"
tempq3$quartile <- "Q3"
tempq4$quartile <- "Q4"

result.holder.int2 <- rbind(tempq1,tempq2, tempq3, tempq4)
result.holder.int2$quartile.f <- factor(result.holder.int2$quartile)

