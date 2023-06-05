#remove all variables and set working directory
rm(list=ls(all=TRUE))
# setwd("~/repos/additional_si_wastewater_Rich_et_al_2023/figure4_maintext")

# input csv files
obs_data <- read.csv(file = paste("input/MPBC_clust4_obsdata.csv"), sep = ",", header = T)
pred_data <- read.csv(file = paste("input/MPBC_clust4_preddata.csv"), sep = ",", header = T)

# libraries
library(ggplot2)
library(dplyr)
library(scales)


obs_data$Type <- rep("Observed",nrow(obs_data))
colnames(obs_data) <- c("MP","Cluster","btrule","Reaction","Type")
pred_data$Type <- rep("Predicted", nrow(pred_data))
colnames(pred_data) <- c("MP","Cluster","btrule","Reaction","Type")
pred_data <- na.omit(pred_data)
MP.TP.combined <- rbind(obs_data, pred_data)
MP.TP.combined$Type <- factor(MP.TP.combined$Type, levels = c("Predicted","Observed"))


color.values.stepped <- c("#B263EC","#653EB3", # purple
                          "#BE0032", # red
                          "#0F8299","#3E9FB3","#7ABECC","#B8DEE6", #blues
                          "#F38400","#FEAF16", # oranges
                          "#999999")


MP.TP.combined$btrule <- factor(MP.TP.combined$btrule, levels = c("bt0001","bt0002",
                                                                  "bt0005",
                                                                  "bt0011","bt0012","bt0013","bt0014",
                                                                  "bt0241","bt0242"))


MP.TP.combined.noNA <- na.omit(MP.TP.combined)


# Way to relabel facet grid names with a function and labeller----------------------------

Types <- list(
  'Predicted'="[A] Predicted",
  'Observed'="[B] Observed"
)

facet_labeller <- function(variable,value){
  return(Types[value])
}
# 


# -------------------FIGURE 4-------------------------------------------
tiff("figure4.tiff", units="in", width=7, height=4, res=300)
ggplot(data.frame(MP.TP.combined.noNA), aes(x=Cluster, fill = btrule))+
  geom_bar(position = "stack")+
  scale_fill_manual(values = color.values.stepped)+
  theme(text = element_text(size=16))+
  facet_wrap(~Type, labeller = facet_labeller)+
  theme(text = element_text(size=18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.key.size = unit(0.3, 'cm'),
        axis.text.x = element_text(vjust = 1),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "gray97"),
        panel.grid.minor = element_line(color = "gray97"),
        panel.border = element_rect(color = "grey", fill = NA))+
  scale_x_discrete(limits=c("1","2","3","4"))+
  ylab("Count")
dev.off()
