#remove all variables and set working directory
rm(list=ls(all=TRUE)) 
setwd("~/repos/additional_si_wastewater_Rich_et_al_2023/figure3_maintext")

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
library(grid)


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

#-----------------GENERATE HISTOGRAM BEFORE REPLACING AVERAGE WITH LOQ-------------------------------------
#png("histogram_means_NAs_in_triplicates.png", width = 1000, height = 1000)
df.means$Location <- factor(df.means$Location, levels = c("Influent","Effluent"))


#replace average with LOQ------------------------------------------------------------------------------------------------
#df.means$M <- replace(df.means$M, which(is.na(df.means$M)), df.means$LOQ) # this is in previous scripts and is the wrong way to replace NAs with LOQ

df.means$M[is.na(df.means$M) == T] <- df.means$LOQ[is.na(df.means$M) == T]
#---------------------------------------------------------------------------------------------------------------


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

set_date_levels <- c("Nov18","Nov19","Nov20","Nov21","Nov22","Nov23","Nov24","Nov25","Nov26","Nov27","Nov28","Nov29","Nov30","Dec1")


# End data cleaning--------- mv.means.flow is a good dataframe-------------

# 2. input text files used to subset the data
MPBC <- read.delim(file = paste("input//MP_binBC.txt"), sep="\t", header=F)

# generate boxplot of rate constants for SI
plot.MPBC <- mv.means.flow[mv.means.flow$Compound %in% MPBC[,1],]

#format data for heatmap
get_data <- function(names){
  data_sub1 <- mv.means.flow[mv.means.flow$Compound %in% names[,1],]
  data_sub1 <- data_sub1[,-c(3:7)]
  matrix <-  data_sub1 %>%
    spread(Date,k)
  return(matrix)
}

# 3. change "heatmap_data" variable to get matrix for different lists of MPs-------------------------------------------------------------------
heatmap_data <- get_data(MPBC)

# rename rows
rn <- heatmap_data$Compound
heatmap_data <- heatmap_data[,-c(1)]
rownames(heatmap_data) <- rn


#changing colors of heatmap with diverging scale----------------------------
ncol <- 100

## Make a vector with n colors
cols <- RColorBrewer:::brewer.pal(11,"RdBu")  # OR c("purple","white","orange")
#cols <- c("blue","white","red")
rampcols <- colorRampPalette(colors = rev(cols), space="Lab")(ncol)

## Make a vector with n+1 breaks
rampbreaks <- seq(0, 100, length.out = ncol+1)

#-------- FIGURE 3 - CLUSTERED HEATMAP --------------------
rn <- rownames(heatmap_data)
rn[22] <- "Ritalinic acid"
rn[3] <- "Benzotriazole-methyl-1H"
rn[9] <- "Dimethyl phthalate"
rn[14] <- "Gabapentin-lactam"
rn[24] <- "Tributyl phosphate"

rownames(heatmap_data) <- rn

fig3 <- pheatmap(heatmap_data, cluster_cols = FALSE, scale = "row",
                 clustering_method = "ward.D2",fontsize_row = 12,fontsize_col = 12,cutree_rows = 4, cellwidth = 10, cellheight = 10,
                 color = rampcols, width = 2, height = 5)
dev.off()

colors <- rep(NA, 24)
colors[1:24] <- c("black")
for (i in c(1,2,6,7,8,10,11,12) ){
  colors[i] <- "red"
}

fig3$gtable$grobs[[4]]$gp <- gpar(col=colors)

png("Figure3.tiff",width = 5.5,height=5,units="in",res=300)
fig3
dev.off()
