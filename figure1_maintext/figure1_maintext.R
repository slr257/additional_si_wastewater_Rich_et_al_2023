#remove all variables and set working directory
rm(list=ls(all=TRUE)) 
setwd("~/repos/additional_si_wastewater_Rich_et_al_2023/figure1_maintext")

#libraries needed
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(matrixStats)
library(reshape2)
library(RColorBrewer)
library(wesanderson)
library(ggrepel)
library(gridExtra)


#-----INPUT DATA AND DATA CLEANING-------------------

#Input concentration data
input.inf <- "input//Nov2020_Targets_INF_all_no_yellow_221205_IBU.csv"
input.eff <- "input//Nov2020_Targets_EFF_all_no_yellow_221205_IBU.csv"


data.inf.wide <- read.csv(file = paste(input.inf), sep = ",", header = T)
data.eff.wide <- read.csv(file = paste(input.eff), sep = ",", header = T)

#delete first column
data.inf.wide$X <- NULL
data.eff.wide$X <- NULL

#Input compount names for looping
compound_select <- "input//Compounds.txt"
compound_select <- read.delim(file = paste(compound_select), sep="\t", header=T)

#this needs to be in vector format to subset the Nov.2020 dataset
compound_select <- compound_select$Compound

#clean data in wide format - only includes compounds with good analytics and cal curves with R2 > 0.85
data.inf.wide.cleaned <- data.inf.wide[data.inf.wide$Compound %in% compound_select, ]
data.eff.wide.cleaned <- data.eff.wide[data.eff.wide$Compound %in% compound_select, ]


#only include good data points for further analysis
data.inf.wide.trunc <- data.inf.wide.cleaned[, c(1:3,25:66)]
data.eff.wide.trunc <- data.eff.wide.cleaned[, c(1:3,25:66)]

#-----------------------------------------------------------------------

#-----------------------------------------------------------------------

#define columns to perform a function on later with dplyr
gathercols <- c("Nov18_1","Nov18_2","Nov18_3",
                "Nov19_1","Nov19_2","Nov19_3",
                "Nov20_1","Nov20_2","Nov20_3",
                "Nov21_1","Nov21_2","Nov21_3",
                "Nov22_1","Nov22_2","Nov22_3",
                "Nov23_1","Nov23_2","Nov23_3",
                "Nov24_1","Nov24_2","Nov24_3",
                "Nov25_1","Nov25_2","Nov25_3",
                "Nov26_1","Nov26_2","Nov26_3",
                "Nov27_1","Nov27_2","Nov27_3",
                "Nov28_1","Nov28_2","Nov28_3",
                "Nov29_1","Nov29_2","Nov29_3",
                "Nov30_1","Nov30_2","Nov30_3",
                "Dec1_1","Dec1_2","Dec1_3")


#Input list of dates
Dates <- "input//Dates_trunc.txt"
Dates <- read.delim(file = paste(Dates), sep="\t", header=T)
Dates <- Dates$Dates

#gather function in "tidyr" gathers column names included in vector titled "gathercols" - converts data to long format
data.inf.long <-gather(data.inf.wide.trunc,gathercols, key = "Date", value = "Concentration")
data.eff.long <-gather(data.eff.wide.trunc,gathercols, key = "Date", value = "Concentration")



#--------HISTOGRAM OF INF DETECTIONS BEFORE LOQ REPLACEMENT---------------------------

#separate date triplicate subsets
INF.data <- data.inf.long %>% separate(Date, c("Date", NA))
#makes dates in correct order
INF.data$Date <- factor(INF.data$Date, levels=Dates)

# replaces non-numeric entries with NAs
INF.data$Concentration <- as.numeric(INF.data$Concentration)

# NEED TO REPLACE NEGATIVE VALUES WITH NAS OR ELSE THEY GET COUNTED IN AVERAGING ADDED 12/04/2022---------------------------------------------------
INF.data$Concentration[INF.data$Concentration <= 0] <- NA

#get mean values out of triplicates
INF.means <- INF.data %>%
  group_by(Compound, Date, LOQ) %>% # group_by groups data by listed values
  summarise(m=mean(Concentration, na.rm=T))

colnames(INF.means) <- c("Compound","Date","LOQ","Concentration")
INF.means$Concentration <- replace(INF.means$Concentration, which(is.nan(INF.means$Concentration)), NA)


#------counting occurances--------------
INF.data.occ <- INF.means
colnames(INF.data.occ) <- c("Compound","Date","LOQ","Detection")

#binary data
INF.data.occ$Detection <- replace(INF.data.occ$Detection, which(INF.data.occ$Detection > 0), 1) # concentration greater than zero (detected)
INF.data.occ$Detection <- replace(INF.data.occ$Detection, which(INF.data.occ$Detection < 0), 0) # concentration less than zero (not-detected)
INF.data.occ$Detection <- replace(INF.data.occ$Detection, which(is.na(INF.data.occ$Detection)), 0) # NA - also not detected


detections.inf <- data.frame(matrix(ncol = 2, nrow = 152))
colnames(detections.inf) <- c("Compound","Detections")
detections.inf$Compound <- as.character(compound_select)


for (i in 1:length(compound_select)){
  data.current <- INF.data.occ[INF.data.occ$Compound == as.character(compound_select[i]),]
  c <- sum(data.current$Detection)
  detections.inf$Detections[i] <- c
  
}

#--------HISTOGRAM OF EFF DETECTIONS BEFORE LOQ REPLACEMENT---------------------------

#separate date triplicate subsets
EFF.data <- data.eff.long %>% separate(Date, c("Date", NA))
#makes dates in correct order
EFF.data$Date <- factor(EFF.data$Date, levels=Dates)

# replaces non-numeric entries with NAs
EFF.data$Concentration <- as.numeric(EFF.data$Concentration)

# NEED TO REPLACE NEGATIVE VALUES WITH NAS OR ELSE THEY GET COUNTED IN AVERAGING ADDED 12/04/2022---------------------------------------------------
EFF.data$Concentration[EFF.data$Concentration <= 0] <- NA

#get mean values out of triplicates
EFF.means <- EFF.data %>%
  group_by(Compound, Date, LOQ) %>% # group_by groups data by listed values
  summarise(m=mean(Concentration, na.rm=T))

colnames(EFF.means) <- c("Compound","Date","LOQ","Concentration")
EFF.means$Concentration <- replace(EFF.means$Concentration, which(is.nan(EFF.means$Concentration)), NA)

#------counting occurances--------------
EFF.data.occ <- EFF.means
colnames(EFF.data.occ) <- c("Compound","Date","LOQ","Detection")

#binary data
EFF.data.occ$Detection <- replace(EFF.data.occ$Detection, which(EFF.data.occ$Detection > 0), 1) # concentration greater than zero (detected)
EFF.data.occ$Detection <- replace(EFF.data.occ$Detection, which(EFF.data.occ$Detection < 0), 0) # concentration less than zero (not-detected)
EFF.data.occ$Detection <- replace(EFF.data.occ$Detection, which(is.na(EFF.data.occ$Detection)), 0) # NA - also not detected


detections.eff <- data.frame(matrix(ncol = 2, nrow = 152))
colnames(detections.eff) <- c("Compound","Detections")
detections.eff$Compound <- as.character(compound_select)


for (i in 1:length(compound_select)){
  data.current <- EFF.data.occ[EFF.data.occ$Compound == as.character(compound_select[i]),]
  c <- sum(data.current$Detection)
  detections.eff$Detections[i] <- c
  
}

# format data for figure 1A
detections.inf$Location <- rep("INF",152)
detections.eff$Location <- rep("EFF",152)

df.detections <- rbind(detections.inf,detections.eff)
colnames(df.detections) <- c("Compound","Detections","Location")

str(df.detections)
df.detections$Location <- factor(df.detections$Location, levels = c("INF","EFF"))


#---------------------------------FIGURE 1A-------------------------------------
p1 <- ggplot(df.detections, aes(x = Detections, fill = Location))+
  geom_histogram(position = "dodge", binwidth = 1, color = "grey50")+
  stat_bin(bins = 15, geom = "text", size=2, aes(label = ..count..), position=position_dodge(width = 1), vjust = -0.2)+
  theme_bw()+
  theme(text = element_text(size = 12),
        plot.title = element_text(vjust = -7, hjust = 0.01,size = 12),
        legend.position = c(0.3, 0.85),
        legend.text = element_text(size=rel(0.5)),
        legend.title = element_text(size=rel(0.5)),
        legend.key.size = unit(0.05, 'in'),
        plot.margin=unit(c(-0.25,0.25,0.25,0.1), "cm")
  )+
  ylab("Count")+
  xlab("Number of Days Observed")+
  scale_x_continuous(breaks=0:14)+
  scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70,80), limits = c(0,80))+
  labs(title = "[A]")


#-----------------------------------------------------------------------------------


#---------------------prepare dataset for other analyses------------------------
data.inf.long$Location <- "Influent"
data.eff.long$Location <- "Effluent"

#combine inf and eff datasets
Nov.2020.data <- rbind(data.inf.long,data.eff.long)

#separate date triplicate subsets
Nov.2020.data <- Nov.2020.data %>%
  separate(Date, c("Date", NA))

#makes dates in correct order
Nov.2020.data$Date <- factor(Nov.2020.data$Date, levels=Dates)

#change data type to numeric!
Nov.2020.data$Concentration <- as.numeric(Nov.2020.data$Concentration)


#-------------DATASET WITHOUT LOQ REPLACEMENT----------------------------------------------------------------------------
#replace negative values with NA (this also replaces <LOQ entries with NA)
Nov.2020.data$Concentration <- replace(Nov.2020.data$Concentration, which(Nov.2020.data$Concentration <= 0), NA)

#replace values that are less than reported LOQ with NA
Nov.2020.data$Concentration <- replace(Nov.2020.data$Concentration, which(Nov.2020.data$Concentration < Nov.2020.data$LOQ), NA)

#-----------------------------------------------------------------------------------------------------------------------------
#get mean values out of triplicates
df.means <- Nov.2020.data %>% 
  group_by(Compound, Date, Location,LOQ) %>% # group_by groups data by listed values
  summarise(M=mean(Concentration, na.rm=T)) # summary statistic on concentration (can also be st.dev. or min/max)
#spread(Location, m)

#somehow NA values turn into NaNs, which might be a data type thing
df.means$M <- replace(df.means$M, which(is.nan(df.means$M)), NA)
df.means$Location <- factor(df.means$Location, levels = c("Influent","Effluent"))

#------REPLACE nondetects WITH LOQ--------------------------------------------------------------------------------------------------------------------------------------------------------------

df.means$M[is.na(df.means$M) == T] <- df.means$LOQ[is.na(df.means$M) == T]
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# remove values set to LOQ from the matrix used for the histogram
df.means.histogram <- df.means[df.means$M != df.means$LOQ,]


#-------------------------DATA FOR CALCULATING RATE CONSTANTS------------------------------
colnames(df.means) <- c("Compound","Date","Location","LOQ","Concentration")

mv.means.calcs <- df.means %>%
  spread(Location, Concentration)

#Input list of flowrates----------------
Flow <- "input//spdes_flowrates.csv"
Flow <- read.delim(file = paste(Flow), sep=",", header=T)
Flow <- Flow[,-2]

mgd_avg <-rep(Flow,152)

#add column to total dataset
mv.means.flow <- cbind(mv.means.calcs,mgd = mgd_avg) #this dataframe will be good for making loading plots!!


#calcualte rate constants------------------
vol.tank <- 0.491*3 # million gallons - this is one tank multiplied by 3 because there are three tanks online and the HRT changes with this

mv.means.flow$HRT <- vol.tank/mv.means.flow$mgd # calculates HRT in days - HRT = V/Q

mv.means.flow$k <-log(mv.means.flow$Effluent/mv.means.flow$Influent)*(-1/mv.means.flow$HRT) # rate constants in d-1 using -1/HRT*ln(cout/cin)

mv.means.flow$Rem <- (mv.means.flow$Influent- mv.means.flow$Effluent)*100/mv.means.flow$Influent

mv.means.flow.pos <- mv.means.flow[mv.means.flow$Rem > 0,]

sum_stats_pos <- mv.means.flow.pos%>%
  group_by(Compound) %>%
  summarise(I.m = mean(Influent), I.sd = sd(Influent),R.m = mean(Rem), R.sd = sd(Rem), k.m = mean(k), k.sd = sd(k), n = length(Compound))

sum_stats_pos$R.CoV <- sum_stats_pos$R.sd/sum_stats_pos$R.m
sum_stats_pos$k.CoV <- sum_stats_pos$k.sd/sum_stats_pos$k.m

total.k.mean <- mean(sum_stats_pos$k.m)
total.k.sd <- sd(sum_stats_pos$k.m)

sum_stats_pos$k.z <-(sum_stats_pos$k.m - total.k.mean)/total.k.sd


#-------------GET %R data------------------
rem.dat <- (mv.means.calcs$Influent- mv.means.calcs$Effluent)*100/mv.means.calcs$Influent
removal.data <- cbind(mv.means.calcs,Removal = rem.dat)
removal.data.plot <- removal.data[removal.data$Removal > 0,]


#change labels for plot legend
df.means.histogram.plot <- df.means.histogram

# way to change this from google
levels(df.means.histogram.plot$Location)[match("Influent",levels(df.means.histogram.plot$Location))] <- "INF"
levels(df.means.histogram.plot$Location)[match("Effluent",levels(df.means.histogram.plot$Location))] <- "EFF"

#------------------------FIGURE 1B------------------------------------------------

p2 <- ggplot(df.means.histogram.plot, aes(x = log10(M), fill = Location)) +
  geom_histogram(position = "identity", alpha = 0.3)+
  labs(
    title = "[B]",
    x = expression("Concentration ("~ mu*"g/L)"),
    y = "Count"
  )+
  theme_bw()+
  # values on x-axis for concentration are log(ng/L), but I just changed the labels to reflect ug/L by shifting by 10^-3
  scale_x_continuous(breaks=seq(0,6,1),limits = c(0,6), labels = c(expression(10^-3),expression(10^-2),expression(10^-1),expression(10^0),expression(10^1),expression(10^2),expression(10^3)))+
  #scale_x_continuous(breaks=seq(0,7,1))+
  scale_y_continuous(breaks=seq(0,200,20))+
  theme(
    plot.title = element_text(vjust = -7, hjust = 0.01,size = 12),
    #axis.text = element_text(size = 12),
    #axis.title = element_text(size = 12),
    legend.text = element_text(size=rel(0.5)),
    legend.title = element_text(size=rel(0.5)),
    legend.position = c(0.8, 0.8),
    legend.key.size = unit(0.05, 'in'),
    plot.margin=unit(c(-0.5,0.25,0.25,0.1), "cm")
  )


# remove points with no SD before plotting----------------
sum_stats_pos_atleast3 <- sum_stats_pos[sum_stats_pos$n > 2,]
names.MP78.atleast3 <- sum_stats_pos_atleast3$Compound


# test out dataset before replacing NAs with LOQ for the 34 pos MPs
MP_plotlabels <- read.delim(file = paste("input//MP_plot_list_Fig1C.txt"), sep="\t", header=F)

#------------------FIGURE 1C---------------------------------------------------
p3 <- ggplot(sum_stats_pos_atleast3, aes(x = log10(I.m), y = R.m))+
  geom_point(shape=21, stroke=1, aes(fill=R.sd),color="grey50",size=2)+ 
  #scale_fill_gradient(low = "white", high = "darkgreen")+
  scale_fill_gradientn(colors = rev(heat.colors(6)))+
  theme_bw()+
  scale_x_continuous(breaks=seq(0,6,1), limits = c(0,6),labels = c(expression(10^-3),expression(10^-2),expression(10^-1),expression(10^0),expression(10^1),expression(10^2),expression(10^3)))+
  scale_y_continuous(limits = c(0,105))+
  #scale_x_continuous(breaks=seq(0,6,1))+
  theme(
    plot.title = element_text(vjust = -7, hjust = 0.01,size = 12),
    # axis.text = element_text(size = 12),
    # axis.title = element_text(size = 12),
    legend.text = element_text(size=rel(0.5)),
    legend.position = c(0.9, 0.45),
    legend.title = element_text(size=rel(0.5)),
    # panel.grid.major.x = element_line(color = "gray96"),
    # panel.grid.major.y = element_line(color = "gray50"),
    legend.key.size = unit(0.06, 'in'),
    plot.margin=unit(c(-0.5,0.25,0.1,0.1), "cm")
  )+
  labs(
    title = "[C]",
    x =  expression("INF Concentration ("~mu*"g/L)"),
    y = "% Removal",
    fill = "SD",
  )+
  geom_text_repel(
    data = sum_stats_pos_atleast3[sum_stats_pos_atleast3$Compound %in% MP_plotlabels$V1,], aes(label = Compound),
    #label=sum_stats_pos_atleast3$Compound, 
    # nudge_x = 0.5, nudge_y = 0.5, 
    # check_overlap = T,
    size = 2,
    fontface = "bold"
  )


tiff("FIGURE1_ABC.tiff",width=3.25,height=9,units="in",res=300)
grid.arrange(p1,p2,p3, nrow=3)
dev.off()


#### END ####