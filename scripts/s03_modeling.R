## Created by: pharmacometrician
## Date: dd/mm/yyyy
## fit pre-specified model and make predictions

#Clean environment
rm(list=ls())

library(dplyr)
library(cowplot)
library(ggplot2)
library(lme4)  
library(lmerTest)
library(pbkrtest)
library(lsmeans)
library(broom)
library(broom.mixed)
library(tidyverse)
library(ggrepel)
library(flextable)
library(officer)
library(stringr)
library(ggthemes)

# add caption
time<-as.character(Sys.time())
script_name <- "s03_modeling.R"
caption  <-paste0("Datasource: r-script:", script_name,", ", time, sep="")

# ggplot themes
theme_set(theme_bw(base_size = 12))

## Settings for plots and labels
conc.label <- "Drug concentration (ng/mL)"
delta.QTcF.label <- expression(paste(Delta,"QTcF (msec)"))
deltadelta.QTcF.label <- expression(paste(Delta,Delta,"QTcF (msec)"))

## Import datasetS
d0.dof <- read.csv("data/derived/qtpk_dofetilide.csv") 
d0.ver <- read.csv("data/derived/qtpk_verapamil.csv") 

############# dofetilide #############

## Data preparations
exp.resp.dof <- d0.dof %>%
  select(USUBJID, TIME, QTcF.CFB, CONC, QTcF.mB, QTcF.B, TREAT, ACTIVE) %>%
  ## Time as factor 
  mutate(TIME=as.factor(TIME)) %>% 
  filter(!is.na(CONC)) %>% 
  filter(!is.na(QTcF.CFB)) %>%
  ## Ensure a dichotomous variable
  mutate(ACTIVE=as.factor(ACTIVE))

## Fitting prespecified model

Pre.spec.model.dof <- lmerTest::lmer(QTcF.CFB ~ TIME+
                                     ACTIVE+
                                     I(QTcF.B-QTcF.mB)+
                                     CONC+
                                     (CONC||USUBJID),
                                     data = exp.resp.dof)

summary(Pre.spec.model.dof)

## in case of convergence issues, one can try simplifying the model by eliminating random effect on slope, see code below:

Pre.spec.model.dof.v2 <- lmerTest::lmer(QTcF.CFB ~ TIME+
                                       ACTIVE+
                                       I(QTcF.B-QTcF.mB)+
                                       CONC+
                                       (1|USUBJID),
                                     data = exp.resp.dof)

## Extraction of parameter estimates 
Par.tab.temp.dof<-as.data.frame(coef(summary(Pre.spec.model.dof, 
                                         ddf="Kenward-Roger"))) 

## Formating the table and calculating CI
Par.tab.dof <-Par.tab.temp.dof %>%
  rownames_to_column(var="Fixed effect parameter") %>%
  rename(p  = `Pr(>|t|)`, 
         se = `Std. Error`) %>%
  mutate(`Relative standard error (%)`=se/Estimate*100, 
         `Lower 95% CI` = Estimate + qt(0.025, df=df) * se,
         `Upper 95% CI` = Estimate + qt(0.975, df=df) * se) %>%
  # rounding to 3 sign digits
  mutate(p=ifelse(p<0.001, "<0.001", signif(p,3)),
         Estimate = round(Estimate,3),
         `Relative standard error (%)`= round(`Relative standard error (%)`,3), 
         `Lower 95% CI` = round(`Lower 95% CI`,3),
         `Upper 95% CI` = round(`Upper 95% CI`,3)
         ) %>%
  select(`Fixed effect parameter`, Estimate, `Lower 95% CI`,
                `Upper 95% CI`, `Relative standard error (%)`, `p`)

Par.tab.dof

## Create table using package flextable
tab1 <- Par.tab.dof %>% flextable() %>%  
  ## Remove borders
  border_remove() %>% 
  ## Set borders for header and bottom
  border(i = NULL,
         j = NULL,
         border = NULL,
         border.top = fp_border(color="black", width=2),
         border.bottom = fp_border(color="black", width=2),
         border.left = NULL,
         border.right = NULL,
         part = "header") %>%
  hline_bottom(j = NULL, border = fp_border(color="black", width=2), part = "body") %>%
  ## Set header style
  bold(bold = TRUE, part = "header") %>%
  align(align="left", part = "all") %>% 
  align(align="right", j=2:6, part="body") %>%
  ## Set font and size
  font(fontname = "Times new roman", part = "all") %>% 
  fontsize(size = 10, part = "all") %>%
  fontsize(size = 8, part="footer") %>%
  autofit()
tab1

## Goodnes of fit plots------------------------------------------------------------------------------------------

GOF.dat<-augment(Pre.spec.model.dof) %>%
  mutate(Sd.rediduals=residuals(Pre.spec.model.dof, scaled = T))

IPRED.DV<-ggplot(data=GOF.dat, aes(x=.fitted, y=QTcF.CFB))+
  geom_point(alpha=0.2)+
  geom_abline(intercept=0, slope=1)+
  geom_smooth(aes(x=.fitted, y=QTcF.CFB),
              method="loess",se=T)+
  labs(x="Model-predicted value (msec)",
       y="Observed Value (msec)")+
  #coord_cartesian(ylim = c(, 140), xlim = c(-40, 140))+
  theme(legend.position="none")

QQ.plot<-ggplot(GOF.dat, aes(sample = Sd.rediduals)) +
  stat_qq(alpha=0.2) +
  #coord_cartesian(ylim = c(-5, 5), xlim = c(-5, 5))+
  geom_abline(slope=1)+
  labs(x="Theoretical Quantiles",
       y="Standardized Residuals") +
  scale_color_gdocs()+
  theme(legend.position="none")

RES.CONC<-ggplot(GOF.dat, aes(x=CONC, y= Sd.rediduals)) +
  geom_point(alpha=0.2) +
  geom_abline(slope=0) +
  labs(x="Drug concentration (ng/mL)",
       y="Standardized Residuals") +
  geom_smooth()+
  geom_abline(intercept=1.96,slope=0,linetype="dashed") +
  geom_abline(intercept=-1.96,slope=0,linetype="dashed") +
  scale_color_gdocs() +
  theme(legend.position="none")

RES.Base.QTcF<-ggplot(GOF.dat, aes(x=`I(QTcF.B - QTcF.mB)`, y= Sd.rediduals)) + 
  geom_point(alpha=0.2) +
  geom_abline(slope=0) +
  labs(x="Centered Baseline QTcF (msec)",
       y="Standardized Residuals") +
  geom_smooth() +
  geom_abline(intercept=1.96,slope=0,linetype="dashed") +
  geom_abline(intercept=-1.96,slope=0,linetype="dashed") +
  scale_color_gdocs() +
  theme(legend.position="none")

RES.TIME<-ggplot(GOF.dat, aes(x=TIME, y= Sd.rediduals)) +
  geom_boxplot(notch = TRUE, outlier.color="white") +
  geom_abline(slope=0, col="#f44336") +
  labs(x="Nominal Time (hours)",
       y="Standardized Residuals") +
  geom_abline(intercept=1.96,slope=0,linetype="dashed") +
  geom_abline(intercept=-1.96,slope=0,linetype="dashed") +
  scale_color_gdocs() +
  theme(legend.position="none")

# add ACTIVE label to indicate Placebo and Drug on the plot
GOF.dat$ACTIVE_label <- " "
GOF.dat$ACTIVE_label[GOF.dat$ACTIVE==0] <- "Placebo\n(ACTIVE=0)"
GOF.dat$ACTIVE_label[GOF.dat$ACTIVE==1] <- "Drug\n(ACTIVE=1)"

RES.ACTIVE<-ggplot(GOF.dat, aes(x=ACTIVE_label, y= Sd.rediduals)) + 
  geom_boxplot(notch = TRUE, outlier.color="white") +
  geom_abline(slope=0, col="#f44336") +
  labs(x="Active treatment",
       y="Standardized Residuals") +
  geom_abline(intercept=1.96,slope=0,linetype="dashed") +
  geom_abline(intercept=-1.96,slope=0,linetype="dashed") +
  scale_color_gdocs() +
  theme(legend.position="none")

allplot <- plot_grid(IPRED.DV, QQ.plot, RES.CONC, RES.Base.QTcF, RES.TIME, RES.ACTIVE, 
                     ncol=3, align = "hv", labels="auto", hjust = 0.1)

allplot

## Specify the concentrations of interest - this can be geo mean Cmax in this example

Cmax.dof <- d0.dof %>% 
  # Remove missing concentrations 
  filter(!is.na(CONC)) %>%
  # exclude placebo
  filter(TREAT!="Placebo") %>%                    
  group_by(USUBJID) %>% 
  ## Get Cmax for each ID
  mutate(Cmax.i = max(CONC)) %>%   
  ungroup() %>%   
  # Only select the first observation per ID 
  distinct(USUBJID, .keep_all=T) %>%  
  # Summary statistics 
  group_by(TREAT) %>%                 
  summarise(GMEAN    = exp(mean(log(Cmax.i))),
            GCV      = 100 * sqrt(exp(var(log(Cmax.i))) - 1),
            n.id     = n_distinct(USUBJID))

## Specify the concentrations of interest (here, it is geometric mean of active arm)
Cmax.ref.grid.dof <- ref.grid(Pre.spec.model.dof, 
                          at=list(CONC=c(0,Cmax.dof$GMEAN)), 
                          ACTIVE=c(0,1)) 

# for user-defined concentration of interest (for example, derived from a popPK model, rather than using geometric mean Cmax), 
# one can use the following example code (later in the code, Cmax.ref.grid.dof.user.coi should be used instead of Cmax.ref.grid.dof):

user.coi <- 3.1 # define concentration of interest here
Cmax.ref.grid.dof.user.coi <- ref.grid(Pre.spec.model.dof, 
                              at=list(CONC=c(0,user.coi)),  # concentration of interest is used here
                              ACTIVE=c(0,1)) 

## Generate LS-Means estimated of baseline adjusted QTcF for active and placebo
## by CONC, predictions are averaged over the levels of TIME. 
Cmax.lsm.dof <-lsmeans::lsmeans(Cmax.ref.grid.dof, c("CONC", "ACTIVE"))

## Get DQTCF estimates
Cmax.dQTcF.dof<-summary(Cmax.lsm.dof, level=0.9) %>% filter(ACTIVE==1, CONC!=0)

## Get  DDQTCF estimates
Cmax.ddQTcF.tmp.dof <-summary(contrast(Cmax.lsm.dof, method = "trt.vs.ctrl1"),
                         infer = c(TRUE, FALSE), 
                         level = .90, adjust = "none")
## Select only relevant contrasts
Cmax.ddQTcF.tmp2.dof <-Cmax.ddQTcF.tmp.dof[(nrow(Cmax.ddQTcF.tmp.dof)/2):nrow(Cmax.ddQTcF.tmp.dof)+1,]

## Formating of ddQTcF output table
Cmax.ddQTcF.dof<- Cmax.ddQTcF.tmp2.dof %>% 
  separate(contrast, into = "CONC", sep = "\\s",extra="drop") %>%
  mutate(CONC = sub('....', '', CONC)) %>%
  mutate(CONC = as.numeric(as.character(CONC))) %>%
  filter(CONC!=0) %>%  #remove zero conc
  dplyr::select(-SE, -df) %>%
  round(.,2)

## exposure-response plot with concentration bins

Bin.qtpk.dof <- exp.resp.dof %>%
  filter(ACTIVE==1) %>%
  group_by(Decile=ntile(CONC, 10), ACTIVE) %>% 
  mutate(CONC_min=min(CONC,na.rm=T),
         CONC_median=median(CONC,na.rm=T),
         DQTCF_UCL=mean(QTcF.CFB,na.rm=T)+qnorm(0.95)*(sd(QTcF.CFB,na.rm=T)/sqrt(sum(!is.na(QTcF.CFB)))),
         DQTCF_LCL= mean(QTcF.CFB,na.rm=T)+qnorm(0.05)*(sd(QTcF.CFB,na.rm=T)/sqrt(sum(!is.na(QTcF.CFB)))),
         DQTCF_mean=mean(QTcF.CFB,na.rm=T)
  ) %>%
  mutate(ACTIVE=as.factor(ACTIVE))

## Create data for plotting bin width
plot.bins.dof<-expand.grid(cutpoints=c(Bin.qtpk.dof$CONC_min, max(exp.resp.dof$CONC, na.rm=T)),
                       y.plot=min(exp.resp.dof$QTcF.CFB, na.rm=T)*1.2)

# extract model estimates for whole concentration range

tmp1.dof <- ref.grid(Pre.spec.model.dof,                                                   # Model 
                 at=list(CONC=seq(from=0,to=max(exp.resp.dof$CONC),length.out = 200)),     # Concentration to make prediction
                 ACTIVE=c(0,1))                                                             # Placebo and active

tmp2.dof<-lsmeans::lsmeans(tmp1.dof, c("CONC", "ACTIVE"))
tmp3.dof<-summary(tmp2.dof, level=0.9) %>% filter(ACTIVE==1, CONC!=0)
tmp4.dof<-summary(contrast(tmp2.dof, method = "trt.vs.ctrl1"),
              infer = c(TRUE, FALSE), 
              level = .90, adjust = "none") # Do not adjust for multiplicity

## Select only relevant contrasts
tmp5.dof<-tmp4.dof[(nrow(tmp4.dof)/2):nrow(tmp4.dof)+1,]
tmp6.dof<- tmp5.dof %>% 
  separate(contrast, into = "CONC", sep = "\\s",extra="drop") %>%
  mutate(CONC = sub('....', '', CONC)) %>%
  mutate(CONC = as.numeric(as.character(CONC))) %>%
  filter(CONC!=0) %>%  #remove zero conc
  dplyr::select(-SE, -df) %>%
  round(.,2)

# exp-resp plot, dQTc
exp.resp.dof$ACTIVE <- factor(exp.resp.dof$ACTIVE, levels = c(1,0))

Pred.plot.dQTcF.bin.dof <-ggplot()+
  # scatter of observed data
  geom_point(data=exp.resp.dof, aes(x=CONC, y=QTcF.CFB, col=as.factor(ACTIVE)), alpha=0.2) +
  geom_smooth(data=exp.resp.dof, aes(x=CONC, y=QTcF.CFB), method = "loess", color = "red", se = FALSE, linetype = "dashed") +
  # model-predictions
  geom_ribbon(data=tmp3.dof, 
              aes(x=CONC, ymin= lower.CL, ymax=upper.CL), fill="black", alpha=0.1)+
  geom_line(data=tmp3.dof,  aes(x=CONC, y=lsmean), size=1, col="black")+
  labs(title="Dofetilide",x=conc.label, y=delta.QTcF.label)+
  # add bins of observed data
  geom_pointrange(data=Bin.qtpk.dof, aes(x=CONC_median,
                                     ymax=DQTCF_UCL,
                                     ymin= DQTCF_LCL,
                                     y=DQTCF_mean,
                                     col=ACTIVE)) +
  geom_point(data=plot.bins.dof, aes(y=y.plot, x=cutpoints), shape="|", size=2 ) +
  geom_line(data=plot.bins.dof, aes(y=y.plot, x=cutpoints), size=0.25) +
  scale_color_gdocs()+
  theme(legend.position="none")
Pred.plot.dQTcF.bin.dof

# exp-resp plot, ddQTc

Pred.plot.ddQTcF.bin.dof <-ggplot()+
  # model-prediction
  geom_ribbon(data=tmp6.dof, 
              aes(x=CONC, ymin= lower.CL, ymax=upper.CL), fill="black", alpha=0.1)+
  geom_line(data=tmp6.dof,  aes(x=CONC, y=estimate), size=1, col="black")+
  geom_hline(aes(yintercept=10), linetype=2)+
  labs(title="Dofetilide",x=conc.label, y=deltadelta.QTcF.label) +
  geom_segment(data=Cmax.ddQTcF.dof[1,], 
               aes(x=CONC, xend=CONC,
                   y=-1, yend=upper.CL), col= "#f44336", size=1)+
  geom_segment(data=Cmax.ddQTcF.dof[1,], 
               aes(x=CONC, xend=0,
                   y=upper.CL, yend=upper.CL), col= "#f44336", size=1,
               arrow = arrow(length = unit(0.25,"cm")))+
  geom_label(data=Cmax.ddQTcF.dof[1,], 
             aes(x=0, y=upper.CL), 
             label=paste0("Upper 90% CI at geomean Cmax: ",Cmax.ddQTcF.dof$upper.CL, " (ms)"),
             hjust=0, vjust=-1, fill="#f44336", 
             col="white")
Pred.plot.ddQTcF.bin.dof

############# verapamil #############
## Data preparations
exp.resp.ver <- d0.ver %>%
  select(USUBJID, TIME, QTcF.CFB, CONC, QTcF.mB, QTcF.B, TREAT, ACTIVE) %>%
  ## Time as factor 
  mutate(TIME=as.factor(TIME)) %>% 
  filter(!is.na(CONC)) %>% 
  filter(!is.na(QTcF.CFB)) %>%
  ## Ensure a dichotomous variable
  mutate(ACTIVE=as.factor(ACTIVE)) %>%
  # re-scale verapamil concentration (to avoid convergence issues)
  mutate(conc2=CONC/1000)

## Fitting prespecified model

Pre.spec.model.ver <- lmerTest::lmer(QTcF.CFB ~ TIME+
                                       ACTIVE+
                                       I(QTcF.B-QTcF.mB)+
                                       CONC+
                                       (CONC||USUBJID),
                                     data = exp.resp.ver)

# the mnodel using CONC as concentration, doesn't converge
# to fix the issue with non-convergence, concentration was re-scaled (by converting from ng/mL to ug/mL)

# model using conc2 (rescaled concentration):

Pre.spec.model.ver <- lmerTest::lmer(QTcF.CFB ~ TIME+
                                       ACTIVE+
                                       I(QTcF.B-QTcF.mB)+
                                       conc2+
                                       (conc2||USUBJID),
                                     data = exp.resp.ver)

summary(Pre.spec.model.ver)

## Extraction of parameter estimates 
Par.tab.temp.ver<-as.data.frame(coef(summary(Pre.spec.model.ver, 
                                             ddf="Kenward-Roger"))) 

## Formating the table and calculating CI
Par.tab.ver <-Par.tab.temp.ver %>%
  rownames_to_column(var="Fixed effect parameter") %>%
  rename(p  = `Pr(>|t|)`, 
         se = `Std. Error`) %>%
  mutate(`Relative standard error (%)`=se/Estimate*100, 
         `Lower 95% CI` = Estimate + qt(0.025, df=df) * se,
         `Upper 95% CI` = Estimate + qt(0.975, df=df) * se) %>%
  mutate(p=ifelse(p<0.001, "<0.001", signif(p,3)),
         Estimate = round(Estimate,3),
         `Relative standard error (%)`= round(`Relative standard error (%)`,3), 
         `Lower 95% CI` = round(`Lower 95% CI`,3),
         `Upper 95% CI` = round(`Upper 95% CI`,3)
         ) %>%
  select(`Fixed effect parameter`, Estimate, `Lower 95% CI`,
         `Upper 95% CI`, `Relative standard error (%)`, `p`)

Par.tab.ver

## Create table using package flextable
tab2 <- Par.tab.ver %>% flextable() %>%  
  ## Remove borders
  border_remove() %>% 
  ## Set borders for header and bottom
  border(i = NULL,
         j = NULL,
         border = NULL,
         border.top = fp_border(color="black", width=2),
         border.bottom = fp_border(color="black", width=2),
         border.left = NULL,
         border.right = NULL,
         part = "header") %>%
  hline_bottom(j = NULL, border = fp_border(color="black", width=2), part = "body") %>%
  ## Set header style
  bold(bold = TRUE, part = "header") %>%
  align(align="left", part = "all") %>% 
  align(align="right", j=2:6, part="body") %>%
  ## Set font and size
  font(fontname = "Times new roman", part = "all") %>% 
  fontsize(size = 10, part = "all") %>%
  fontsize(size = 8, part="footer") %>%
  autofit()
tab2

## Goodnes of fit plots------------------------------------------------------------------------------------------

GOF.dat<-augment(Pre.spec.model.ver) %>%
  mutate(Sd.rediduals=residuals(Pre.spec.model.ver, scaled = T))

IPRED.DV<-ggplot(data=GOF.dat, aes(x=.fitted, y=QTcF.CFB))+
  geom_point(alpha=0.2)+
  geom_abline(intercept=0, slope=1)+
  geom_smooth(aes(x=.fitted, y=QTcF.CFB),
              method="loess",se=T)+
  labs(x="Model-predicted value (msec)",
       y="Observed Value (msec)")+
  #coord_cartesian(ylim = c(, 140), xlim = c(-40, 140))+
  theme(legend.position="none")

QQ.plot<-ggplot(GOF.dat, aes(sample = Sd.rediduals)) +
  stat_qq(alpha=0.2) +
  #coord_cartesian(ylim = c(-5, 5), xlim = c(-5, 5))+
  geom_abline(slope=1)+
  labs(x="Theoretical Quantiles",
       y="Standardized Residuals") +
  scale_color_gdocs()+
  theme(legend.position="none")

RES.CONC<-ggplot(GOF.dat, aes(x=conc2, y= Sd.rediduals)) +
  geom_point(alpha=0.2) +
  geom_abline(slope=0) +
  labs(x="Drug concentration (ng/mL)",
       y="Standardized Residuals") +
  geom_smooth()+
  geom_abline(intercept=1.96,slope=0,linetype="dashed") +
  geom_abline(intercept=-1.96,slope=0,linetype="dashed") +
  scale_color_gdocs() +
  theme(legend.position="none")

RES.Base.QTcF<-ggplot(GOF.dat, aes(x=`I(QTcF.B - QTcF.mB)`, y= Sd.rediduals)) + 
  geom_point(alpha=0.2) +
  geom_abline(slope=0) +
  labs(x="Centered Baseline QTcF (msec)",
       y="Standardized Residuals") +
  geom_smooth() +
  geom_abline(intercept=1.96,slope=0,linetype="dashed") +
  geom_abline(intercept=-1.96,slope=0,linetype="dashed") +
  scale_color_gdocs() +
  theme(legend.position="none")

RES.TIME<-ggplot(GOF.dat, aes(x=TIME, y= Sd.rediduals)) +
  geom_boxplot(notch = TRUE, outlier.color="white") +
  geom_abline(slope=0, col="#f44336") +
  labs(x="Nominal Time (hours)",
       y="Standardized Residuals") +
  geom_abline(intercept=1.96,slope=0,linetype="dashed") +
  geom_abline(intercept=-1.96,slope=0,linetype="dashed") +
  scale_color_gdocs() +
  theme(legend.position="none")

# add ACTIVE label to indicate Placebo and Drug on the plot
GOF.dat$ACTIVE_label <- " "
GOF.dat$ACTIVE_label[GOF.dat$ACTIVE==0] <- "Placebo\n(ACTIVE=0)"
GOF.dat$ACTIVE_label[GOF.dat$ACTIVE==1] <- "Drug\n(ACTIVE=1)"

RES.ACTIVE<-ggplot(GOF.dat, aes(x=ACTIVE_label, y= Sd.rediduals)) + 
  geom_boxplot(notch = TRUE, outlier.color="white") +
  geom_abline(slope=0, col="#f44336") +
  labs(x="Active treatment",
       y="Standardized Residuals") +
  geom_abline(intercept=1.96,slope=0,linetype="dashed") +
  geom_abline(intercept=-1.96,slope=0,linetype="dashed") +
  scale_color_gdocs() +
  theme(legend.position="none")

allplot.ver <- plot_grid(IPRED.DV, QQ.plot, RES.CONC, RES.Base.QTcF, RES.TIME, RES.ACTIVE, 
                     ncol=3, align = "hv", labels="auto", hjust = 0.1)

allplot.ver

## Specify the concentrations of interest - this can be geo mean Cmax in this example

Cmax.ver <- d0.ver %>% 
  # calculate conc2, to match re-scaled concentration that was used in modeling
  mutate(conc2=CONC/1000) %>%
  # Remove missing concentrations 
  filter(!is.na(conc2)) %>%
  # exclude placebo
  filter(TREAT!="Placebo") %>%                    
  group_by(USUBJID) %>% 
  ## Get Cmax for each ID
  mutate(Cmax.i = max(conc2)) %>%   
  ungroup() %>%   
  # Only select the first observation per ID 
  distinct(USUBJID, .keep_all=T) %>%  
  # Summary statistics 
  group_by(TREAT) %>%                 
  summarise(GMEAN    = exp(mean(log(Cmax.i))),
            GCV      = 100 * sqrt(exp(var(log(Cmax.i))) - 1),
            n.id     = n_distinct(USUBJID))

## Specify the concentrations of interest (here, it is geomean of active arm)
Cmax.ref.grid.ver <- ref.grid(Pre.spec.model.ver, 
                              at=list(conc2=c(0,Cmax.ver$GMEAN)), 
                              ACTIVE=c(0,1)) 

## Generate LS-Means estimated of baseline adjusted QTcF for active and placebo
## by conc2, predictions are averaged over the levels of TIME. 
Cmax.lsm.ver <-lsmeans::lsmeans(Cmax.ref.grid.ver, c("conc2", "ACTIVE"))

## Get DQTCF estimates
Cmax.dQTcF.ver<-summary(Cmax.lsm.ver, level=0.9) %>% filter(ACTIVE==1, conc2!=0)

## Get  DDQTCF estimates
Cmax.ddQTcF.tmp.ver <-summary(contrast(Cmax.lsm.ver, method = "trt.vs.ctrl1"),
                              infer = c(TRUE, FALSE), 
                              level = .90, adjust = "none")
## Select only relevant contrasts
Cmax.ddQTcF.tmp2.ver <-Cmax.ddQTcF.tmp.ver[(nrow(Cmax.ddQTcF.tmp.ver)/2):nrow(Cmax.ddQTcF.tmp.ver)+1,]

## Formating of ddQTcF output table
Cmax.ddQTcF.ver<- Cmax.ddQTcF.tmp2.ver %>% 
  separate(contrast, into = "conc2", sep = "\\s",extra="drop") %>%
  mutate(conc2 = sub('.....', '', conc2)) %>%
  mutate(conc2 = as.numeric(as.character(conc2))) %>%
  filter(conc2!=0) %>%  #remove zero conc
  dplyr::select(-SE, -df) %>%
  round(.,2)

## exposure-response plot with concentration bins

Bin.qtpk.ver <- exp.resp.ver %>%
  filter(ACTIVE==1) %>%
  group_by(Decile=ntile(conc2, 10), ACTIVE) %>% 
  mutate(conc2_min=min(conc2,na.rm=T),
         CONC_median=median(conc2,na.rm=T),
         DQTCF_UCL=mean(QTcF.CFB,na.rm=T)+qnorm(0.95)*(sd(QTcF.CFB,na.rm=T)/sqrt(sum(!is.na(QTcF.CFB)))),
         DQTCF_LCL= mean(QTcF.CFB,na.rm=T)+qnorm(0.05)*(sd(QTcF.CFB,na.rm=T)/sqrt(sum(!is.na(QTcF.CFB)))),
         DQTCF_mean=mean(QTcF.CFB,na.rm=T)
  ) %>%
  mutate(ACTIVE=as.factor(ACTIVE))

## Create data for plotting bin width
plot.bins.ver<-expand.grid(cutpoints=c(Bin.qtpk.ver$conc2_min, max(exp.resp.ver$conc2, na.rm=T)),
                           y.plot=min(exp.resp.ver$QTcF.CFB, na.rm=T)*1.2)

# extract model estimates for whole concentration range

tmp1.ver <- ref.grid(Pre.spec.model.ver,                                                   # Model 
                     at=list(conc2=seq(from=0,to=max(exp.resp.ver$conc2),length.out = 200)),     # Concentration to make prediction
                     ACTIVE=c(0,1))                                                             # Placebo and active

tmp2.ver<-lsmeans::lsmeans(tmp1.ver, c("conc2", "ACTIVE"))
tmp3.ver<-summary(tmp2.ver, level=0.9) %>% filter(ACTIVE==1, conc2!=0)

exp.resp.ver$ACTIVE <- factor(exp.resp.ver$ACTIVE, levels = c(1,0))

Pred.plot.dQTcF.bin.ver <-ggplot()+
  # scatter of observed data
  geom_point(data=exp.resp.ver, aes(x=conc2*1000, y=QTcF.CFB, col=as.factor(ACTIVE)), alpha=0.2) +
  geom_smooth(data=exp.resp.ver, aes(x=conc2*1000, y=QTcF.CFB), method = "loess", color = "red", se = FALSE, linetype = "dashed") +
  # model-predictions
  geom_ribbon(data=tmp3.ver, 
              aes(x=conc2*1000, ymin= lower.CL, ymax=upper.CL), fill="black", alpha=0.1)+
  geom_line(data=tmp3.ver,  aes(x=conc2*1000, y=lsmean), size=1, col="black")+
  labs(title="Verapamil",x="Drug concentration (ng/mL)", y=delta.QTcF.label)+
  # add bins of observed data
  geom_pointrange(data=Bin.qtpk.ver, aes(x=CONC_median*1000,
                                         ymax=DQTCF_UCL,
                                         ymin= DQTCF_LCL,
                                         y=DQTCF_mean,
                                         col=ACTIVE)) +
  geom_point(data=plot.bins.ver, aes(y=y.plot, x=cutpoints*1000), shape="|", size=2 ) +
  geom_line(data=plot.bins.ver, aes(y=y.plot, x=cutpoints*1000), size=0.25) +
  scale_color_gdocs()+
  theme(legend.position="none")
Pred.plot.dQTcF.bin.ver

## summary table with dofetilide and verapamil ddQTcF estimates and 90% CI
# rename verapamin variables to match dofetilides
Cmax.ddQTcF.ver <- Cmax.ddQTcF.ver %>% rename(CONC=conc2)
# change the units for verapamil back to ng/mL
Cmax.ddQTcF.ver$CONC <- Cmax.ddQTcF.ver$CONC*1000 

summ.ddqtcf <- rbind(Cmax.ddQTcF.dof,Cmax.ddQTcF.ver)
summ.ddqtcf[, 'Treatment'] <- c("Dofetilide (500 ug)","Verapamil (120 mg)")

summ.ddqtcf <- summ.ddqtcf %>%
  relocate(Treatment,CONC,estimate,lower.CL,upper.CL)
# output as flextable
tab3 <- summ.ddqtcf %>%
  flextable() %>%  
  set_header_labels(CONC     = "Concentration", 
                    estimate = "Estimate",
                    lower.CL = "Lower 90% CI",
                    upper.CL = "Upper 90% CI") %>%
  ## Remove borders
  border_remove() %>% 
  ## Set borders for header and bottom
  border(i = NULL,
         j = NULL,
         border = NULL,
         border.top = fp_border(color="black", width=2),
         border.bottom = fp_border(color="black", width=2),
         border.left = NULL,
         border.right = NULL,
         part = "header") %>%
  hline_bottom(j = NULL, border = fp_border(color="black", width=2), part = "body") %>%
  ## Set header style
  bold(bold = TRUE, part = "header") %>%
  align(align="left", part = "all") %>% 
  ## Set font and size
  font(fontname = "Times new roman", part = "all") %>% 
  fontsize(size = 10, part = "all") %>%
  fontsize(size = 8, part="footer") %>%
  autofit()
tab3

# plot with dofetilide and verapamil exposure-response
allplot.er <- plot_grid(Pred.plot.dQTcF.bin.dof, Pred.plot.dQTcF.bin.ver, 
                         ncol=2, align = "hv", labels="auto", hjust = 0.1)

allplot.er
# extract model estimates for whole concentration range

tmp1.ver <- ref.grid(Pre.spec.model.ver,                                                         # Model 
                     at=list(conc2=seq(from=0,to=max(exp.resp.ver$conc2),length.out = 200)),     # Concentration to make prediction
                     ACTIVE=c(0,1))                                                              # Placebo and active

tmp2.ver<-lsmeans::lsmeans(tmp1.ver, c("conc2", "ACTIVE"))
tmp3.ver<-summary(tmp2.ver, level=0.9) %>% filter(ACTIVE==1, conc2!=0)
tmp4.ver<-summary(contrast(tmp2.ver, method = "trt.vs.ctrl1"),
                  infer = c(TRUE, FALSE), 
                  level = .90, adjust = "none") # Do not adjust for multiplicity

## Select only relevant contrasts
tmp5.ver<-tmp4.ver[(nrow(tmp4.ver)/2):nrow(tmp4.ver)+1,]
tmp6.ver<- tmp5.ver %>% 
  separate(contrast, into = "conc2", sep = "\\s",extra="drop") %>%
  mutate(conc2 = sub('.....', '', conc2)) %>%
  mutate(conc2 = as.numeric(as.character(conc2))) %>%
  filter(conc2!=0) %>%  #remove zero conc
  dplyr::select(-SE, -df) %>%
  round(.,3)

# exp-resp plot, ddQTc

Pred.plot.ddQTcF.bin.ver <-ggplot()+
  # model-prediction
  geom_ribbon(data=tmp6.ver, 
              aes(x=conc2*1000, ymin= lower.CL, ymax=upper.CL), fill="black", alpha=0.1)+
  geom_line(data=tmp6.ver,  aes(x=conc2*1000, y=estimate), size=1, col="black")+
  geom_hline(aes(yintercept=10), linetype=2)+
  labs(title="Verapamil",x="Drug concentration (ng/mL)", y=deltadelta.QTcF.label) +
  geom_segment(data=Cmax.ddQTcF.ver[1,], 
               aes(x=CONC, xend=CONC,
                   y=-1, yend=upper.CL), col= "#f44336", size=1)+
  geom_segment(data=Cmax.ddQTcF.ver[1,], 
               aes(x=CONC, xend=0,
                   y=upper.CL, yend=upper.CL), col= "#f44336", size=1,
               arrow = arrow(length = unit(0.25,"cm")))+
  geom_label(data=Cmax.ddQTcF.ver[1,], 
             aes(x=0, y=upper.CL), 
             label=paste0("Upper 90% CI at geomean Cmax: ",Cmax.ddQTcF.ver$upper.CL, " (ms)"),
             hjust=0, vjust=-1, fill="#f44336", 
             col="white")
Pred.plot.ddQTcF.bin.ver

# plot ddQTcF vs conc together

# plot with dofetilide and verapamil exposure-response
allplot.ddqt.er <- plot_grid(Pred.plot.ddQTcF.bin.dof, Pred.plot.ddQTcF.bin.ver, 
                        ncol=2, align = "hv", labels="auto", hjust = 0.1)

allplot.ddqt.er
