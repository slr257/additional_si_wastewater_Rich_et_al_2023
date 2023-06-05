#remove all variables and set working directory
rm(list=ls(all=TRUE)) 
# setwd("~/Nov2020_Project2/Target Screening/221205_Add_Ibuprofen/concentration_profiles")

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
library(gridExtra)
library(grid)
library(scales)
library(gtable)


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

# test out dataset before replacing NAs with LOQ for the 34 pos MPs
MP35 <- read.delim(file = paste("input//MP35_all_names.txt"), sep="\t", header=F)

Nov.2020.MP35 <- Nov.2020.data[Nov.2020.data$Compound %in% MP35$V1,]


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
#--------------------PLOT LOADINGS USING DF.MEANS-----------------------------


##------ Plotting mean values for concentration profiles-----------------
compound_select.df <- as.data.frame(compound_select)
colnames(compound_select.df) <- c("Compound")

#-----------------------get %R using mv.means.calcs-------------------
mv.means.calcs$Removal <- (mv.means.calcs$Influent - mv.means.calcs$Effluent)*100/mv.means.calcs$Influent

add.Rem <- as.data.frame(cbind(mv.means.calcs$Compound,mv.means.calcs[,2],rep("%R",2128), mv.means.calcs$LOQ, mv.means.calcs$Removal))

# rename %R values to concentrations to match df.means format for easy plotting
colnames(add.Rem) <- c("Compound","Date","Location","LOQ","Concentration")

df.means.plot <- rbind(df.means,add.Rem)

# rename dataframe for simplicity
df.means <- df.means.plot

#-----------------------formatted grid plot of profiles---------------------------------------
levels(df.means$Location)[match("Influent",levels(df.means$Location))] <- "INF"
levels(df.means$Location)[match("Effluent",levels(df.means$Location))] <- "EFF"

df.means$Location <- factor(df.means$Location, levels = c("Influent","Effluent","%R"))

plots.keep <- c(34,124,70)
plot.vec <- list()


for (i in plots.keep){
  data.current <- df.means[df.means$Compound == as.character(compound_select.df$Compound[i]),]
  
  data.current.inf <- data.current[data.current$Location == "Influent",]
  data.current.eff <- data.current[data.current$Location == "Effluent",]
  data.current.rem <- data.current[data.current$Location == "%R",]
  
  plot.data <- cbind(data.current.inf[,c(1,2,5)], data.current.eff[,5],data.current.rem[,5])
  colnames(plot.data) <- c("Compound","Date","INF","EFF","Removal")
  
  if(i == 34){
    
    coeff1 = max(data.current.inf$Concentration[data.current.inf$Compound == "Caffeine"])
    
    plot.vec[[i]]<- ggplot(data = plot.data, aes(x = Date))+
      #ggtitle(compound_select.df$Compound[i])+
      geom_line(data = plot.data, aes(y = Removal/(100/coeff1), color = "%R", size = 0.4), group = 1, lwd = 2)+
      geom_line(data = plot.data, aes(y = INF, color = "INF", size = 0.4), group = 1, lwd = 2)+
      geom_line(data = plot.data, aes(y = EFF, color = "EFF", size = 0.4), group = 1, lwd = 2)+
      scale_size(guide = "none")+
      scale_fill_discrete(labels = factor(c("INF","EFF","%R"), levels = c("INF","EFF","%R")))+
      scale_color_manual(values = c("grey", "#00BFC4","#F8766D"))+
      # scale_fill_manual(labels=c("INF","EFF","%R"), values=c("dodgerblue4", "firebrick4","grey"))+
      ggtitle(compound_select.df$Compound[i])+
      theme(axis.text = element_text(size = 8),
            plot.title = element_text(vjust = -7, hjust = 0.01,size = 12),
            legend.text = element_text(size = 8),
            #axis.text.x = element_text(size = 9, angle = 45, vjust = 0.5),
            axis.text.x =  element_blank(),
            axis.text.y = element_text(size = 8.5),
            legend.title = element_blank(),
            legend.key.height = unit(0.4, "cm"), 
            legend.key.width = unit(0.2,"cm"),
            legend.position = c(0.2,0.8),
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            panel.background = element_blank(), 
            panel.grid.major = element_line(color = "gray97"), 
            panel.grid.minor = element_line(color = "gray97"),
            panel.border = element_rect(color = "lightgrey", fill = NA),
            text = element_text(size=12)
      )+
      scale_x_discrete(limits = set_date_levels)+
      scale_y_continuous(limits = c(0,coeff1), labels = comma,
                         # Features of the first axis
                         name = "Concentration",
                         # Add a second axis and specify its features
                         sec.axis = sec_axis(~.*(100/coeff1), name="%R"))+
      labs(fill="Location", title = "[A]")
    #ylab("Concentration (ng/L)")
  }else if(i == 124){
    
    coeff2 = max(data.current.inf$Concentration[data.current.inf$Compound == "Propranolol"])
    
    plot.vec[[i]]<- ggplot(data = plot.data, aes(x = Date))+
      #ggtitle(compound_select.df$Compound[i])+
      geom_line(data = plot.data, aes(y = Removal/(100/coeff2), color = "%R", size = 0.4), group = 1, lwd = 2)+
      geom_line(data = plot.data, aes(y = INF, color = "INF", size = 0.4), group = 1, lwd = 2)+
      geom_line(data = plot.data, aes(y = EFF, color = "EFF", size = 0.4), group = 1, lwd = 2)+
      scale_size(guide = "none")+
      scale_fill_discrete(labels = factor(c("INF","EFF","%R"), levels = c("INF","EFF","%R")))+
      scale_color_manual(values = c("grey", "#00BFC4","#F8766D"))+
      theme(axis.text = element_text(size = 8),
            plot.title = element_text(vjust = -7, hjust = 0.01,size = 12),
            # legend.text = element_text(size = 8),
            #axis.text.x = element_text(size = 9, angle = 45, vjust = 0.5),
            axis.text.x =  element_blank(),
            axis.text.y = element_text(size = 9),
            legend.title = element_blank(),
            # legend.key.height = unit(0.05, "cm"), 
            # legend.key.width = unit(0.05,"cm"),
            # legend.position = c(0.85,0.85),
            legend.position = "none",
            # axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            panel.background = element_blank(), 
            panel.grid.major = element_line(color = "gray97"), 
            panel.grid.minor = element_line(color = "gray97"),
            panel.border = element_rect(color = "lightgrey", fill = NA),
            text = element_text(size=12),
      )+
      scale_x_discrete(limits = set_date_levels)+
      labs(fill="Location", title = "[B]")+
      scale_y_continuous(limits = c(0,coeff2),labels = comma,
                         # Features of the first axis
                         name = "Concentration (ng/L)",
                         
                         # Add a second axis and specify its features
                         sec.axis = sec_axis(~.*(100/coeff2), name="%R"))
    #ylab("Concentration (ng/L)")
  }else if(i == 70){
    
    coeff3 = max(data.current.inf$Concentration[data.current.inf$Compound == "Gabapentin"])
    
    plot.vec[[i]]<- ggplot(data = plot.data, aes(x = Date))+
      #ggtitle(compound_select.df$Compound[i])+
      geom_line(data = plot.data, aes(y = Removal/(100/coeff3), color = "%R", size = 0.4), group = 1, lwd = 2)+
      geom_line(data = plot.data, aes(y = INF, color = "INF", size = 0.4), group = 1, lwd = 2)+
      geom_line(data = plot.data, aes(y = EFF, color = "EFF", size = 0.4), group = 1, lwd = 2)+
      scale_size(guide = "none")+
      scale_fill_discrete(labels = factor(c("INF","EFF","%R"), levels = c("INF","EFF","%R")))+
      scale_color_manual(values = c("grey", "#00BFC4","#F8766D"))+
      theme(axis.text = element_text(size = 8),
            plot.title = element_text(vjust = -7, hjust = 0.01,size = 12),
            legend.text = element_text(size = 8),
            axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5),
            axis.text.y = element_text(size = 9),
            legend.title = element_blank(),
            legend.key.height = unit(0.05, "cm"), legend.key.width = unit(0.05,"cm"),
            #legend.position = c(0.85,0.85),
            legend.position = "none",
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            panel.background = element_blank(), 
            panel.grid.major = element_line(color = "gray97"), 
            panel.grid.minor = element_line(color = "gray97"),
            panel.border = element_rect(color = "lightgrey", fill = NA),
            text = element_text(size=12),
      )+
      scale_x_discrete(limits = set_date_levels)+
      scale_y_continuous(limits = c(0,coeff3),labels = comma,
                         # Features of the first axis
                         name = "Concentration",
                         
                         # Add a second axis and specify its features
                         sec.axis = sec_axis(~.*(100/coeff3), name="%R"))+
      labs(fill="Location", title = "[C]")
    
  }
}

p1 <- plot.vec[[34]]
p2 <- plot.vec[[124]]
p3 <- plot.vec[[70]]


grid.arrange(p1,p2,p3, nrow = 1)
g1 <- list(p1,p2,p3)

ga <- ggplotGrob(p1)
gb <- ggplotGrob(p2)
gc <- ggplotGrob(p3)
g <- rbind(ga,gb, gc, size = "first")
g$widths <- unit.pmax(ga$widths, gb$widths, gc$widths)
grid.newpage()
grid.draw(g)

png("Figure2.tiff",width=3.5,height=10,units="in",res=300)
grid.draw(g)
dev.off()
