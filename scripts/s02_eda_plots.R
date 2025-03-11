## Created by: pharmacometrician
## Date: dd/mm/yyyy
## generate exploratory data analysis plots

#Clean environment
rm(list=ls())

library(dplyr)
library(cowplot)
library(ggplot2)
library(ggrepel)
library(readr)
library(tidyr)
library(ggthemes)

## Import datasetS
d0.dof <- read.csv("data/derived/qtpk_dofetilide.csv") 
d0.ver <- read.csv("data/derived/qtpk_verapamil.csv") 

## Settings for plots and labels
conc.label <- "Drug concentration (ng/mL)" 

# QT labels
delta.QTcF.label <- expression(paste(Delta,"QTcF (msec)"))
deltadelta.QTcF.label <- expression(paste(Delta,Delta,"QTcF (msec)"))

# time labels
time.label <- "Time (hours)"

# ggplot themes
theme_set(theme_bw(base_size = 12))

plot_settings <- list(
theme_classic(),
theme(axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    strip.text.x = element_text(size = 15)))

## summary function (calculates arithmetic mean, median, min, max, sd (standard deviation), n (number of samples), 
# se (standard error), LCL (lower limit of 90% CI) and UCL (upper limit of 90% CI) )
my.sum.fun<-funs(mean   = mean(. ,na.rm=T),
                 median = median(. , na.rm=T),
                 min    = min(. ,na.rm=T),
                 max    = max(. ,na.rm=T),
                 sd     = sd(. ,na.rm=T), 
                 n      = sum(!is.na(.)),
                 se     = sd(. ,na.rm=T)/sqrt(sum(!is.na(.))),
                 LCL    = mean(. ,na.rm=T)+qnorm(0.05)*(sd(. ,na.rm=T)/sqrt(sum(!is.na(.)))),
                 UCL    = mean(. ,na.rm=T)+qnorm(0.95)*(sd(. ,na.rm=T)/sqrt(sum(!is.na(.)))))

## Time course of dHR and ddHR: format the data and create exploratory plots for HR

# dofetilide

hr.dt.dof <-d0.dof %>%
  # Make long dataset
  dplyr::select(USUBJID, TREAT, TIME, QTcF.CFB, HR.CFB, ddHR, CONC) %>%
  gather(key = KEY, value = VALUE, -USUBJID, -TREAT, -TIME) %>%
  ## Group by these factors
  group_by(KEY, TREAT, TIME) %>%
  ## Summarize HR, PK, and QT by factors using my summary function 
  summarise_at(vars(VALUE), my.sum.fun) %>%
  ## Show mean and Sd for PK and CI for HR and QT
  mutate(P.UCL = ifelse(KEY  == "CONC", mean+sd, UCL),
         P.LCL = ifelse(KEY  == "CONC", mean-sd, LCL)) %>%
  ungroup() %>%
  ## Labels for facets in plot
  mutate(key2=case_when(KEY=="CONC" ~ conc.label,
                        KEY=="HR.CFB" ~ "\u0394HR (bpm)",
                        KEY=="ddHR" ~ "\u0394\u0394HR (bpm)",
                        KEY=="QTcF.CFB" ~ "\u0394QTcF (msec)"
  ))

## Create cut-offs for time course plots
hline.data <- data.frame(cutoff = c(NA, 10, 10, 10),
                         KEY = c("CONC", "HR.CFB", "ddHR","QTcF.CFB"))

## To arrange group legend in plot
hr.dt.dof$TREAT <- factor(hr.dt.dof$TREAT, levels = c("Dofetilide","Placebo"))

## Make HR by time figure

By.Time.HR.dof <-ggplot(data=hr.dt.dof %>% filter(KEY %in% c("HR.CFB","ddHR")), aes(x=TIME,y=mean, col=TREAT))+
  geom_linerange(aes(ymin=P.LCL, ymax=P.UCL))+
  geom_point()+
  geom_line(aes(x=TIME,y=mean, linetype=TREAT, group=TREAT))+
  labs(title="Dofetilide",x="Time (h)", y="")+
  facet_grid(key2~.)+ 
  geom_hline(data = hline.data %>% filter(KEY %in% c("HR.CFB","ddHR")), aes(yintercept = cutoff), linetype="dashed") +
  theme(strip.text = element_text(
    size = 10)) +
  plot_settings +
  theme(legend.title=element_blank(),
        legend.position="bottom",
        plot.title = element_text(hjust = 0.5,size=15),
        strip.text = element_text(size = 11))+
  scale_color_gdocs()
By.Time.HR.dof

# verapamil
hr.dt.ver <-d0.ver %>%
  # Make long dataset
  dplyr::select(USUBJID, TREAT, TIME, QTcF.CFB, HR.CFB, ddHR, CONC) %>%
  gather(key = KEY, value = VALUE, -USUBJID, -TREAT, -TIME) %>%
  ## Group by these factors
  group_by(KEY, TREAT, TIME) %>%
  ## Summarize HR, PK, and QT by factors using my summary function 
  summarise_at(vars(VALUE), my.sum.fun) %>%
  ## Show mean and Sd for PK and CI for HR and QT
  mutate(P.UCL = ifelse(KEY  == "CONC", mean+sd, UCL),
         P.LCL = ifelse(KEY  == "CONC", mean-sd, LCL)) %>%
  ungroup() %>%
  ## Labels for facets in plot
  mutate(key2=case_when(KEY=="CONC" ~ conc.label,
                        KEY=="HR.CFB" ~ "\u0394HR (bpm)",
                        KEY=="ddHR" ~ "\u0394\u0394HR (bpm)",
                        KEY=="QTcF.CFB" ~ "\u0394QTcF (msec)"
  ))

## To arrange group legend in plot
hr.dt.ver$TREAT <- factor(hr.dt.ver$TREAT, levels = c("Verapamil HCL","Placebo"))

## Make HR by time figure

By.Time.HR.ver <-ggplot(data=hr.dt.ver %>% filter(KEY %in% c("HR.CFB","ddHR")), aes(x=TIME,y=mean, col=TREAT))+
  geom_linerange(aes(ymin=P.LCL, ymax=P.UCL))+
  geom_point()+
  geom_line(aes(x=TIME,y=mean, linetype=TREAT, group=TREAT))+
  labs(title="Verapamil",x="Time (h)", y="")+
  facet_grid(key2~.)+ 
  geom_hline(data = hline.data %>% filter(KEY %in% c("HR.CFB","ddHR")), aes(yintercept = cutoff), linetype="dashed") +
  theme(strip.text = element_text(
    size = 10)) +
  plot_settings +
  theme(legend.title=element_blank(),
        legend.position="bottom",
        plot.title = element_text(hjust = 0.5,size=15),
        strip.text = element_text(size = 11))+
  scale_color_gdocs()
By.Time.HR.ver

By.Time.HR.both <- plot_grid(By.Time.HR.dof, By.Time.HR.ver, align="hv")
By.Time.HR.both

### Generate QTc versus RR plot

## dofetilide
# format the data
mplot.data.dof <- d0.dof %>% 
  dplyr::select(USUBJID, TREAT, ACTIVE, QT=QTm, RRm, QTcF, QTcB) %>%
  gather(value = value, key = type, -RRm, -USUBJID, -TREAT, -ACTIVE) 

## To arrange group legend in plot
mplot.data.dof$TREAT <- factor(mplot.data.dof$TREAT, levels = c("Dofetilide","Placebo"))

# generate the plot
m.plot.dof <- ggplot(data=mplot.data.dof, aes(RRm, value, col=TREAT))+
  geom_point(alpha=0.2)+
  facet_wrap(~type, nrow=2)+          
  geom_smooth(aes(RRm, value, fill=TREAT),
              linewidth=1, method="lm", se=T)+
  labs(title="Dofetilide",x="RR interval (msec)", 
       y="") +
  plot_settings +
  theme(legend.title=element_blank(),
        legend.position="bottom",
        plot.title = element_text(hjust = 0.5,size=15))+
  scale_color_gdocs()+
  scale_fill_gdocs() 
m.plot.dof

## verapamil
# format the data
mplot.data.ver <- d0.ver %>% 
  dplyr::select(USUBJID, TREAT, ACTIVE, QT=QTm, RRm, QTcF, QTcB) %>%
  gather(value = value, key = type, -RRm, -USUBJID, -TREAT, -ACTIVE) 

## To arrange group legend in plot
mplot.data.ver$TREAT <- factor(mplot.data.ver$TREAT, levels = c("Verapamil HCL","Placebo"))

# generate the plot
m.plot.ver <- ggplot(data=mplot.data.ver, aes(RRm, value, col=TREAT))+
  geom_point(alpha=0.2)+
  facet_wrap(~type, nrow=2)+          
  geom_smooth(aes(RRm, value, fill=TREAT),
              linewidth=1, method="lm", se=T)+
  labs(title="Verapamil",x="RR interval (msec)", 
       y="") +
  plot_settings +
  theme(legend.title=element_blank(),
        legend.position="bottom",
        plot.title = element_text(hjust = 0.5,size=15))+
  scale_color_gdocs()+
  scale_fill_gdocs()
m.plot.ver

m.plot.both <- plot_grid(m.plot.dof,m.plot.ver, align="hv")

# extract lm equation and summary, if needed; example code for dofetilide, QTcF vs RR:
fit <- lm(value~RRm,data=mplot.data.dof %>% filter(TREAT=="Dofetilide" & type=="QTcF"))
summary(fit) 

### Generate time course of drug concentration and dQTcF

## dofetilide

# identify time of max concentration and max dQTcF to include in the plot
max.effect.dof <- hr.dt.dof %>% 
  filter(TREAT=="Dofetilide") %>%
  filter(KEY %in% c("CONC","QTcF.CFB")) %>%
  group_by(KEY) %>%
  summarize(max.eff = max(mean))
# extract which timepoint this corresponds to
max.time.dof <- hr.dt.dof %>% 
  filter(KEY=="CONC" & mean==max.effect.dof$max.eff[max.effect.dof$KEY=="CONC"] | 
        KEY=="QTcF.CFB" & mean==max.effect.dof$max.eff[max.effect.dof$KEY=="QTcF.CFB"])

# generate the plot
By.Time.QT.dof <-ggplot(data=hr.dt.dof %>% filter(KEY %in% c("CONC","QTcF.CFB")), aes(x=TIME,y=mean, col=TREAT))+
  geom_linerange(aes(ymin=P.LCL, ymax=P.UCL))+
  geom_point()+
  geom_line(aes(x=TIME,y=mean, linetype=TREAT, group=TREAT))+
  labs(title="Dofetilide",x="Time (h)", y="")+
  facet_grid(key2~.,scales = "free_y")+ 
  theme(strip.text = element_text(
    size = 10)) +
  plot_settings +
  theme(legend.title=element_blank(),
        legend.position="bottom",
        plot.title = element_text(hjust = 0.5,size=15),
        strip.text = element_text(size = 11))+
  scale_color_gdocs() +
  # add vertical lines for dofetilide max effect
  geom_vline(data=max.time.dof,aes(xintercept=TIME), linetype="dashed")
By.Time.QT.dof

## verapamil

# identify time of max concentration and max dQTcF to include in the plot
max.effect.ver <- hr.dt.ver %>% 
  filter(TREAT=="Verapamil HCL") %>%
  filter(KEY %in% c("CONC","QTcF.CFB")) %>%
  group_by(KEY) %>%
  summarize(max.eff = max(mean))
# extract which timepoint this corresponds to
max.time.ver <- hr.dt.ver %>% 
  filter(KEY=="CONC" & mean==max.effect.ver$max.eff[max.effect.ver$KEY=="CONC"] | 
           KEY=="QTcF.CFB" & mean==max.effect.ver$max.eff[max.effect.ver$KEY=="QTcF.CFB"])

# generate the plot
By.Time.QT.ver <-ggplot(data=hr.dt.ver %>% filter(KEY %in% c("CONC","QTcF.CFB")), aes(x=TIME,y=mean, col=TREAT))+
  geom_linerange(aes(ymin=P.LCL, ymax=P.UCL))+
  geom_point()+
  geom_line(aes(x=TIME,y=mean, linetype=TREAT, group=TREAT))+
  labs(title="Verapamil",x="Time (h)", y="")+
  facet_grid(key2~.,scales = "free_y")+ 
  theme(strip.text = element_text(
    size = 10)) +
  plot_settings +
  theme(legend.title=element_blank(),
        legend.position="bottom",
        plot.title = element_text(hjust = 0.5,size=15),
        strip.text = element_text(size = 11))+
  scale_color_gdocs() +
  # add vertical lines for verapamil max effect
  geom_vline(data=max.time.ver,aes(xintercept=TIME), linetype="dashed")
By.Time.QT.ver

By.Time.QT.both <- plot_grid(By.Time.QT.dof,By.Time.QT.ver, align="hv")

### Generate hysteresis plot 

##dofetilide

# format the data
hyst.qtpk.dof <- d0.dof %>%
  ## grouping factors 
  group_by(TREAT, TIME) %>% 
  ## summarize PK and ddQT by factors using my summary function 
  summarise_at(vars(CONC, ddQTcF), my.sum.fun) %>%
  mutate(TIME_LABEL = parse_factor(paste0(as.character(TIME)," h")))

# generate the plot
Hysteresis.plot.dof <- ggplot(data=hyst.qtpk.dof %>%filter(TREAT!="Placebo"), aes(x=CONC_mean, y=ddQTcF_mean))+
  geom_path(alpha=0.7)+
  geom_pointrange(aes(ymin=ddQTcF_LCL, ymax=ddQTcF_UCL), alpha=0.5)+
  geom_text_repel(size =3, aes(label=TIME_LABEL))+
  labs(title="Dofetilide",
       x=conc.label, y=deltadelta.QTcF.label) +
  theme(plot.title = element_text(hjust = 0.5,size=15))
Hysteresis.plot.dof

## verapamil

# format the data
hyst.qtpk.ver <- d0.ver %>%
  ## grouping factors 
  group_by(TREAT, TIME) %>% 
  ## summarize PK and ddQT by factors using my summary function 
  summarise_at(vars(CONC, ddQTcF), my.sum.fun) %>%
  mutate(TIME_LABEL = parse_factor(paste0(as.character(TIME)," h")))

# generate the plot
Hysteresis.plot.ver <- ggplot(data=hyst.qtpk.ver %>%filter(TREAT!="Placebo"), aes(x=CONC_mean, y=ddQTcF_mean))+
  geom_path(alpha=0.7)+
  geom_pointrange(aes(ymin=ddQTcF_LCL, ymax=ddQTcF_UCL), alpha=0.5)+
  geom_text_repel(size =3, aes(label=TIME_LABEL))+
  labs(title="Verapamil",
       x=conc.label, y=deltadelta.QTcF.label) +
  theme(plot.title = element_text(hjust = 0.5,size=15))
Hysteresis.plot.ver

### generate concentration versus ΔΔQTcF bin-plot

## dofetilide

# format the data
Bin.qtpk.dof <- d0.dof %>%
  ## remove missing CONC, and ddQTcF and Placebo
  filter(!is.na(CONC) & !is.na(ddQTcF) & TREAT!="Placebo") %>%
  ## Create Deciles and use that as a factor for summary function 
  group_by(Decile=ntile(CONC, 10), TREAT) %>% 
  summarise_at(vars(CONC, ddQTcF), my.sum.fun)

## Create data for plotting bin width
## Shows the cutpoints of the bins (min and max)
plot.bins.dof<-expand.grid(cutpoints=c(Bin.qtpk.dof$CONC_min, max(d0.dof$CONC, na.rm=T)),
                       y.plot=min(Bin.qtpk.dof$ddQTcF_min, na.rm=T)*1.2)

d0.dof$TREAT <- factor(d0.dof$TREAT, levels = c("Dofetilide","Placebo"))

# generate the plot
Explor.plot.dof <-ggplot()+
  geom_point(data=d0.dof, aes(x=CONC, y=ddQTcF, col=TREAT), alpha=0.2)+
  geom_smooth(data=d0.dof, method="lm",se=F, col="black", size=0.5, linetype=2, aes(x=CONC, y=ddQTcF))+
  geom_smooth(data=d0.dof, method="loess",se=F, col="#f44336", aes(x=CONC, y=ddQTcF))+
  geom_pointrange(data=Bin.qtpk.dof, aes(x=CONC_median, 
                                     ymax=ddQTcF_mean+ddQTcF_sd, 
                                     ymin=ddQTcF_mean-ddQTcF_sd,
                                     y=ddQTcF_mean, 
                                     col=TREAT)) +
  geom_point(data=plot.bins.dof, aes(y=y.plot, x=cutpoints), shape="|", size=2 )+
  geom_line(data=plot.bins.dof, aes(y=y.plot, x=cutpoints), size=0.25)+
  labs(title="Dofetilide",x=conc.label, y=deltadelta.QTcF.label) +
  theme(plot.title = element_text(hjust = 0.5,size=15)) +
  theme(legend.position="none") +
  scale_color_gdocs()
Explor.plot.dof

## verapamil

# format the data
Bin.qtpk.ver <- d0.ver %>%
  ## remove missing CONC, and ddQTcF and Placebo
  filter(!is.na(CONC) & !is.na(ddQTcF) & TREAT!="Placebo") %>%
  ## Create Deciles and use that as a factor for summary function 
  group_by(Decile=ntile(CONC, 10), TREAT) %>% 
  summarise_at(vars(CONC, ddQTcF), my.sum.fun)

## Create data for plotting bin width
## Shows the cutpoints of the bins (min and max)
plot.bins.ver<-expand.grid(cutpoints=c(Bin.qtpk.ver$CONC_min, max(d0.ver$CONC, na.rm=T)),
                           y.plot=min(Bin.qtpk.ver$ddQTcF_min, na.rm=T)*1.2)

## To arrange group legend in plot
d0.ver$TREAT <- factor(d0.ver$TREAT, levels = c("Verapamil HCL","Placebo"))

# generate the plot
Explor.plot.ver <-ggplot()+
  geom_point(data=d0.ver, aes(x=CONC, y=ddQTcF, col=TREAT), alpha=0.2)+
  geom_smooth(data=d0.ver, method="lm",se=F, col="black", size=0.5, linetype=2, aes(x=CONC, y=ddQTcF))+
  geom_smooth(data=d0.ver, method="loess",se=F, col="#f44336", aes(x=CONC, y=ddQTcF))+
  geom_pointrange(data=Bin.qtpk.ver, aes(x=CONC_median, 
                                         ymax=ddQTcF_mean+ddQTcF_sd, 
                                         ymin=ddQTcF_mean-ddQTcF_sd,
                                         y=ddQTcF_mean, 
                                         col=TREAT)) +
  geom_point(data=plot.bins.ver, aes(y=y.plot, x=cutpoints), shape="|", size=2 )+
  geom_line(data=plot.bins.ver, aes(y=y.plot, x=cutpoints), size=0.25)+
  labs(title="Verapamil",x=conc.label, y=deltadelta.QTcF.label) +
  theme(plot.title = element_text(hjust = 0.5,size=15)) +
  theme(legend.position="none") +
  scale_color_gdocs()
Explor.plot.ver