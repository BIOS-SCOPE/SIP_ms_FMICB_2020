#####################
#SIP1614_1712_1726 Experiment
#####################

#This script was written to wrangle, process, and analyze the experimental data (BC using biovolume)

#Shuting Liu, Oct 24, 2019, adapted from Nicholas Huynh, November 6, 2018

######PREP DATA########

install.packages("growthcurver_0.3.0.tar", repos=NULL, type="source")
#Part 1. load necessary packages and data ----
library(tidyverse)
library(data.table)
library(growthcurver)
library(oce)
library(zoo)
library(RColorBrewer)
library(gridExtra)
library(dplyr)
library(ggplot2)

#Part 2. Reading in the dataset ----

#this dataset contains the processed bacterial counts from lab 5 as well as the TOC data
master.df <- read.csv("SIP1614_1712_1726_data_for_growthcurver.csv")

#Part 3. Add derived variables ----

#add a column of the natural logarithm of the cell counts, a column with a bact CCF of 148 fgC/um^3(Gundersen et al,2002_L&O), and bact C based on the CCF of 148

master.df <- 
  master.df %>% 
  mutate(lncellsml = round(log(Cells_mL), 2),
         CCF=148,
         bactC = round(CCF * biovolume *Cells_mL *(1/(12*10^6)),2)) %>%       #bactC unit umole C/L
  group_by(Experiment, Bottle) %>%
  ungroup()

#Part 3a. Interpolate the organic carbon values ----

master.df <- 
  master.df %>%
  group_by(Experiment, Bottle) %>%
  mutate(INTERP_DOC = na.approx(DOC,Days, na.rm=F)) %>%  #use Days to interpolate according to x distance #no need to interp bc because every time points have data
  ungroup()
  


#######GROWTH CURVER FOR BACT AbUNDANCE#####

#Part 4. Fitting the growth curves to a logistic equation ----

#Here, we are using the growthcurver package to fit growth curve data to the standard form of the logistic equation and then return a data table with population-level information. Specifically, the carrying capacity k, specific growth rate r, the initial population size N0, the goodness of fit to the logistic equation sigma as well as the degrees of freedom df, the time at which 1/2 carrying capcicity is reached t_mid, the doubling time t_gen, the area under the logistic curve auc_l, and the emprical area under the experimental curve auc_e. 

#4a. Prepare data for growthcurver ----

gcurve.df <- 
 master.df %>% 
  group_by(Experiment, Bottle) %>%
  #initial cells_ml value for each experiment will be subtracted from the cells_ml values 
  mutate(delta_cells = Cells_mL - Cells_mL[which.min(Days)],
         delta_doc = INTERP_DOC[which.min(Days)] - INTERP_DOC,
         delta_bc = bactC - bactC[which.min(Days)],
         Experiment_Bottle = paste(Experiment,Bottle)) %>% #add this column as a comoon variable for join with gc_out_select data
  ungroup() 

#no dead phase exclusion for this SIP dataset


#split the dataframe above by experiment and bottle 
gcurve.list <- split(gcurve.df,paste(gcurve.df$Experiment,gcurve.df$Bottle))
#store the names of each list object (i.e. ID2)
headers = names(gcurve.list)

#4b. Create a function/for loop that plots each of the curves ----

#first apply the summarize growth function to each of the experiments
gcurveplot.func <- function(elf){
  gc_fit <- SummarizeGrowth(elf$Days, elf$delta_cells) 
}
gcplot.list <- lapply(gcurve.list, gcurveplot.func)

#4c. save the plots as a pdf ----
pdf("growthcurves_BA.pdf")
for (i in 1:length(gcplot.list)) {
  plot(gcplot.list[[i]], main = names(gcplot.list[i]) ) 
}
dev.off()

#4d. Create a function that returns growthcurver output fromt the test curves as a dataframe ----

gcurve.func <- function(narwhal){
  gc.fit <- SummarizeGrowth(narwhal$Days, narwhal$delta_cells)
  gc_output.df  <- as.data.frame(unlist(gc.fit$vals), stringsAsFactors = FALSE)
}

#4e. Apply the function to all the growth curves ----
gc_out.list <- lapply(gcurve.list, gcurve.func)
#save the list as a data frame
gc_out.df <- data.frame(gc_out.list, stringsAsFactors = FALSE)
#transpose the data frame
gc_out.df <- as.data.frame(t(gc_out.df), stringsAsFactors = FALSE)
#replace the character strings
gc_out.df$note <- gsub("cannot fit data", "1", gc_out.df$note)  #SIP1726 control cannot fit or too big stationary
#gc_out.df$note <- gsub("questionable fit", "2", gc_out.df$note)
#coerce the data frame to be one of numerics, not characters
gc_out.df <- as.data.frame(sapply(gc_out.df, as.numeric)) 
#reapply the sample names as the row names
rownames(gc_out.df) <- headers
#make the row names the first column of the data frame and change the column name to "ID2"
gc_out.df <-setDT(gc_out.df, keep.rownames = TRUE)[] #make the rownames into the first column
colnames(gc_out.df)[1] <- "Experiment_Bottle"

#4f. Tidy the growth curve data and find the transition points of the growth phases----

gc_out.df <- 
  gc_out.df %>%
  mutate_at(vars(k, k_se, n0, n0_se, sigma, auc_l, auc_e),funs(round(., 0))) %>% #auc_l area under curve from fitting, auc_e from measurement, these are auc from 0 to end of time, k carrying capacity,se standard error,p value,n0 initial population size, r growth rate,tmid time to reach half carrying capacity, sigma residual standard error from fit,t_gen doubling time
  mutate_at(vars(t_mid, t_gen),funs(round(., 1))) %>%
  mutate_at(vars(k_p, n0_p, r, r_p),funs(round(., 2))) %>% 
  mutate(stationary = t_mid*2)  #stationary phase defined as 2 times half the time it takes reach the carrying capacity (delta).
  
write.csv(gc_out.df,"gc_out_BA.csv")

gc_out_select.df <- gc_out.df %>%
                    select(Experiment_Bottle, stationary,t_gen)


#4g. merge the growthcurve output with the input dataframe

gcurve.df <- 
  gcurve.df %>% 
  full_join(., gc_out_select.df)

write.csv(gcurve.df,"gcurve.csv")


#there we will modify the data frame in excel, adding the stationary time points (if not measured) as new rows and then re-import the file. 
#we will then interpolate data for those missing points
#this is necessary to be able to calculate accurate areas under the curves and thus, bge
########IMPORT DATA##########

bge.df <- read.csv("gcurve_interp.csv")

bge.df <- bge.df %>% 
  select(Experiment, Bottle, Experiment_Bottle,Treatment, Treatment_Btl,Days,stationary,t_gen, Cells_mL, DOC, bactC,lncellsml,INTERP_DOC, delta_cells, delta_bc, delta_doc) %>% 
  group_by(Experiment, Bottle) %>%
  mutate(INTERP_cells = na.approx(Cells_mL, Days, na.rm=F,rule=2),
         INTERP_bC = na.approx(bactC, Days, na.rm=F,rule=2),
         INTERP2_DOC = na.approx(INTERP_DOC, Days, na.rm=F,rule=2),
         INTERP_delta_cells=na.approx(delta_cells,Days, na.rm=F,rule=2),
         INTERP_delta_bc = na.approx(delta_bc, Days,na.rm=F,rule=2),
         INTERP_delta_doc = na.approx(delta_doc, Days,na.rm = F,rule=2)
  )

###########CALCULATE BGEs##########

#Part 5. calc bge using area under the curve ---- using true measured data
auc_bge.df <- 
  bge.df %>% 
  group_by(Experiment, Bottle) %>%
  mutate(
    AUC_cells = integrateTrapezoid(Days, INTERP_delta_cells, type="cA"),
    AUC_bactC = integrateTrapezoid(Days, INTERP_delta_bc, type="cA"), #"cA" give cumulative results
    AUC_doc = integrateTrapezoid(Days, INTERP_delta_doc, type="cA"), #delta_doc already include interpolated all time points data
    AUC_bge = AUC_bactC/AUC_doc
  ) %>% 
  distinct()     #unique rows  
#this gives BGE for every time points

#add column for AUC normalized to days 
auc_bge.df <- 
  auc_bge.df %>% 
  group_by(Experiment, Bottle) %>%
  mutate(
    AUC_cells_norm = AUC_cells/Days,
    AUC_bactC_norm = AUC_bactC/Days,
    AUC_doc_norm = AUC_doc/Days
  ) 

####PLOTTING####

#order the data
order <- unique(gcurve.df$Experiment_Bottle)
gcurve.df$Experiment_Bottle <- as.factor(gcurve.df$Experiment_Bottle)
gcurve.df$Experiment_Bottle <- factor(gcurve.df$Experiment_Bottle, levels = order)


bactC.graph <- ggplot(gcurve.df, aes(x=Days, y=bactC)) +
  geom_point(aes(colour=Experiment_Bottle, size =1)) +
  geom_line(aes(colour=Experiment_Bottle, size=0.8)) +
  guides(colour = FALSE, size=FALSE) +
  labs(x="Days", y="Bacterial Carbon (µM C)") +
  ggtitle("Bacterial Carbon") +
  theme_classic() +
 
  scale_x_continuous(breaks = seq(min(gcurve.df$Days), max(gcurve.df$Days), by = 10)) + 
  #scale_y_continuous(breaks = seq(min(gcurve.df$bactC), max(gcurve.df$bactC)),by=0.2) + 
  theme(title = element_text(size =17),
        axis.title = element_text(size=17),
        axis.text = element_text(size=16),
        legend.title=element_text(size=16),  
        legend.text = element_text(size=16))

doc.graph <- ggplot(gcurve.df, aes(x=Days, y=DOC)) +
  geom_point(aes(colour=Experiment_Bottle, size =1)) +
  geom_line(data=gcurve.df[!is.na(gcurve.df$DOC),],aes(colour=Experiment_Bottle, size=0.8)) +
  guides(size=FALSE) +
  labs(x="Hours", y="DOC (µM C)", colour = "Experiment_Bottle") +
  ggtitle("DOC") +
  theme_classic() +

  scale_x_continuous(breaks = seq(min(gcurve.df$Days), max(gcurve.df$Days), by = 10)) + 
  #scale_y_continuous(breaks = seq(min(gcurve.df$DOC), max(gcurve.df$DOC), by = 1)) + 
  theme(title = element_text(size =17),
        axis.title = element_text(size=17),
        axis.text = element_text(size=16),
        legend.title=element_text(size=16),  
        legend.text = element_text(size=16))



raw <- arrangeGrob(bactC.graph, doc.graph, ncol=2)
ggsave("bactc_doc.pdf", raw, width = 18, height = 8, device = "pdf")

lncells.graph <- ggplot(gcurve.df, aes(x=Days, y=lncellsml)) +
  geom_point(aes(colour=Experiment_Bottle, size =1)) +
  geom_line(aes(colour=Experiment_Bottle, size=0.8)) +
  guides(size=FALSE) +
  labs(x="Days", y="ln(Cells/ml)", colour = "Experiment_Bottle") +
  ggtitle("ln(Cells/ml)") +
  theme_classic() +
  
  scale_x_continuous(breaks = seq(min(gcurve.df$Days), max(gcurve.df$Days), by = 10)) + 
  scale_y_continuous(breaks = seq(min(gcurve.df$lncellsml), max(gcurve.df$lncellsml), by = 0.2)) + 
  theme(title = element_text(size =17),
        axis.title = element_text(size=17),
        axis.text = element_text(size=16),
        legend.title=element_text(size=16),  
        legend.text = element_text(size=16))

ggsave("lncells.pdf", width = 10, height = 8, device = "pdf")



###SPECIFIC GROWTH RATES#####

#view the lncells plot to determine the exponential phases of growth, also check growth rate excel file

#break up the dataframe into the bottles. we will work with each bottle separately
a.df <- filter(auc_bge.df, Experiment_Bottle == "SIP1614 A") 
b.df <- filter(auc_bge.df, Experiment_Bottle == "SIP1614 B") 
c.df <- filter(auc_bge.df, Experiment_Bottle == "SIP1614 E")
d.df <- filter(auc_bge.df, Experiment_Bottle == "SIP1614 F")
e.df <- filter(auc_bge.df, Experiment_Bottle == "SIP1614 G")
f.df <- filter(auc_bge.df, Experiment_Bottle == "SIP1614 H")
g.df <- filter(auc_bge.df, Experiment_Bottle == "SIP1614 K")
h.df <- filter(auc_bge.df, Experiment_Bottle == "SIP1614 L")
i.df <- filter(auc_bge.df, Experiment_Bottle == "SIP1712 A")
j.df <- filter(auc_bge.df, Experiment_Bottle == "SIP1712 B")
k.df <- filter(auc_bge.df, Experiment_Bottle == "SIP1712 E")
l.df <- filter(auc_bge.df, Experiment_Bottle == "SIP1712 F")
m.df <- filter(auc_bge.df, Experiment_Bottle == "SIP1726 A")
n.df <- filter(auc_bge.df, Experiment_Bottle == "SIP1726 B")
o.df <- filter(auc_bge.df, Experiment_Bottle == "SIP1726 C")
p.df <- filter(auc_bge.df, Experiment_Bottle == "SIP1726 D")
q.df <- filter(auc_bge.df, Experiment_Bottle == "SIP1726 E")
r.df <- filter(auc_bge.df, Experiment_Bottle == "SIP1726 F")
s.df <- filter(auc_bge.df, Experiment_Bottle == "SIP1726 G")
t.df <- filter(auc_bge.df, Experiment_Bottle == "SIP1726 H")

#add a new column to each data frame with the specific growth rate
#growth rates are calculated as the slope of the exponential phase for each bottle
a.df <- mutate(a.df, µ = round(lm(lncellsml~Days,a.df[c(1:2),])$coeff[[2]],2))
b.df <- mutate(b.df, µ = round(lm(lncellsml~Days,b.df[c(1:2),])$coeff[[2]],2))
c.df <- mutate(c.df, µ = round(lm(lncellsml~Days,c.df[c(1:4),])$coeff[[2]],2))
d.df <- mutate(d.df, µ = round(lm(lncellsml~Days,d.df[c(1:4),])$coeff[[2]],2))
e.df <- mutate(e.df, µ = round(lm(lncellsml~Days,e.df[c(1:3),])$coeff[[2]],2))
f.df <- mutate(f.df, µ = round(lm(lncellsml~Days,f.df[c(1:3),])$coeff[[2]],2))
g.df <- mutate(g.df, µ = round(lm(lncellsml~Days,g.df[c(1:3),])$coeff[[2]],2))
h.df <- mutate(h.df, µ = round(lm(lncellsml~Days,h.df[c(1:3),])$coeff[[2]],2))
i.df <- mutate(i.df, µ = round(lm(lncellsml~Days,i.df[c(1:7),])$coeff[[2]],2))
j.df <- mutate(j.df, µ = round(lm(lncellsml~Days,j.df[c(1:4),])$coeff[[2]],2))
k.df <- mutate(k.df, µ = round(lm(lncellsml~Days,k.df[c(3:5),])$coeff[[2]],2))
l.df <- mutate(l.df, µ = round(lm(lncellsml~Days,l.df[c(3:6),])$coeff[[2]],2))
m.df <- mutate(m.df, µ = round(lm(lncellsml~Days,m.df[c(1:2),])$coeff[[2]],2))
n.df <- mutate(n.df, µ = round(lm(lncellsml~Days,n.df[c(1:3),])$coeff[[2]],2))
o.df <- mutate(o.df, µ = round(lm(lncellsml~Days,o.df[c(2:5),])$coeff[[2]],2))
p.df <- mutate(p.df, µ = round(lm(lncellsml~Days,p.df[c(2:5),])$coeff[[2]],2))
q.df <- mutate(q.df, µ = round(lm(lncellsml~Days,q.df[c(2:6),])$coeff[[2]],2))
r.df <- mutate(r.df, µ = round(lm(lncellsml~Days,r.df[c(2:4),])$coeff[[2]],2))
s.df <- mutate(s.df, µ = round(lm(lncellsml~Days,s.df[c(3:6),])$coeff[[2]],2))
t.df <- mutate(t.df, µ = round(lm(lncellsml~Days,t.df[c(3:5),])$coeff[[2]],2))


#put all of the dataframes back together
spgrowth <- rbind(a.df, b.df, c.df, d.df, e.df, f.df, g.df, h.df,i.df,j.df,k.df,l.df,m.df,n.df,o.df,p.df,q.df,r.df,s.df,t.df)


bge_calc <- spgrowth %>%
  select(Treatment_Btl,Experiment_Bottle, Days, stationary, t_gen,bactC, DOC,INTERP_DOC, delta_cells,delta_bc, delta_doc, INTERP_cells, INTERP_bC, INTERP2_DOC, INTERP_delta_cells,INTERP_delta_bc,INTERP_delta_doc,AUC_cells, AUC_bactC, AUC_doc, AUC_bge,AUC_cells_norm,AUC_bactC_norm,AUC_doc_norm, µ) 


#save the data frame as a csv file 
write.csv(bge_calc, "Calculated_Master.csv")




