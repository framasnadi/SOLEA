#############################
#   SOLEA: Survial toOL on a scEnario bAsis
#############################
# Developed by Enrico Nicola Armelloni & Francesco Masnadi (CNR-IRBIM, Ancona)
##############################
#
#
# Include some brief explanation
#
#
##############################
rm(list=ls())
set.seed(42)


##############
#
#  Setup input parameters and folders
#
############## 

### The working directory must be the folder were the two codes and input data are stored. Please insert here the working directory, or press ctr+shift+h and navigate to the right folder
setwd("~/CNR/SOS/github/SOLEA") 
surv_data<-"Relative" ## Absolute if KM, relative if not
n_scenarios<-as.numeric(4)
censor<-as.numeric(120)
tree_best<-as.numeric(500)
try_best<-as.numeric(2)
filename<-"Input_data.csv"


##############
#
#  Install missing packages and get data
#
############## 
list.of.packages <- c("tidyverse", "caret","rpart","rpart.plot","e1071","randomForest","survival","survminer","data.table","Boruta")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
# Load Packages
library(tidyverse); library(caret);library(rpart);library(rpart.plot);library("e1071");library(randomForest);library(survival);library(survminer);library(data.table);library(Boruta)
dir.create(file.path(".", paste0("Plots & Graphs", surv_data)))
plotdir <- (paste0("./Plots & Graphs", surv_data, "/"))
source("support_functions.R") # source external functions
Data <- import(filename)

##############
#
# Variables exploration
#
##############  
#  Step 1: Explore collinearity between variables with pairs plots (Zuur et al., 2009) and remove (manually) the variables that correlates (>0.7)
Collinearity(Data) ## After run, please check the plot in the folder Plots&Graphics"
Data<-Data %>% dplyr::select(-Towing_speed) # After evaluating the collinearity, we decided to eliminate the variable "Towing speed"

# Step 2: Explore meaningfulness of remaining variables with Boruta model
BS<-Boruta_screen(Data); print(BS)
remove<-dplyr::pull(BS %>%dplyr::filter(value!="Confirmed")%>%dplyr::select(name))
Data<-Data%>%dplyr::select(-remove) # Remove variables not meaningful

##############
#
#  Random Forest (RF) model to identify the stressors affecting immediate survival
#
############## 
rf<-RF_screen(Data)
# Partial dependence plot
imp <- randomForest::importance(rf)
impvar<- rownames(imp)
db <- Data %>% dplyr::select(-Survivability_days, -Vitality_class)
tiff(paste0(plotdir,"PartialDep_RF.tif"),width = 85, height = 85, units = "mm", res = 1200, pointsize = 5)
op<-par(mfrow=c(2, 3))
for (j in seq_along(rownames(imp))) {
  partialPlot(rf, as.data.frame(db), rownames(imp)[j], xlab=rownames(imp)[j],ylab ="Centered Log Odds of Dead",main=paste("Partial Dependence on", rownames(imp)[j]))
}
par(op)
dev.off()

##############
#
#  Classification Tree (CT) to split the dataset into "fishing scenarios"
#
############## 
tree.rpart<-CT_create(Data)

##############
#
#  Fishing scenarios
#
############## 
Scenario_0<-Data %>% dplyr::mutate(scen_set = "Aggregate")
Scenario<-Data %>%dplyr::mutate(scen_set = as.character(tree.rpart[["where"]])) %>%bind_rows(., Scenario_0)
Scenario_list<-split(Scenario,Scenario$scen_set)

#   SI by scenario
SI_list <- lapply(Scenario_list, SI_calculation)
if(surv_data=="Absolute")
  {
  Surv_list<-lapply(Scenario_list, KM_model)
  }else{
  Surv_list<-lapply(Scenario_list, Cox_model)
  
}
#   SO by Scenario


##############
#
#  Final result
#
############## 
Final_result<-mapply(SR,SI_list,Surv_list)


##############
#
#  More plots
#
############## 

# plot SI histogram
plot_scenarios<-bind_rows(SI_list) #%>%dplyr::mutate(Scenar= factor(Scenar, levels=titles))
Vitplot<-ggplot(plot_scenarios, aes(fill=Vitality_class, y=Percentage, x=Scenar ,color=Vitality_class)) +  geom_bar( stat="identity",position="fill")+fill_palette(palette = c("green4", "blue", "red", "grey"))+color_palette(palette = c("black", "black", "black", "black")) +   
  theme(  text=element_text(family="Tahoma") ,legend.text = element_text(size = 8),legend.title = element_blank(),axis.title.x = element_blank(),axis.title.y = element_text( size = 7),
          axis.text.x=element_text(colour="black", size = 8, angle = 90),
          axis.text.y=element_text(colour="black", size = 8)) + ggtitle("Onboard Vitality Assessment")
ggsave("ABCDproportion.tiff", Vitplot, path = plotdir)

# plot SO histogram
plot_SR<-as.data.frame(t(as.data.frame(Final_result)))%>%dplyr::rename("Scenar"="Scenario")%>%dplyr::mutate(Scen_x=as.numeric(seq(1:nrow(.))))%>%dplyr::mutate(SR=as.numeric(SR)) %>%dplyr::mutate(clr=ifelse(SR==max(SR), "green", ifelse(SR==min(SR), "red", "#999999")))

colrs<-c("red" = "red", "green" = "green", "#999999" = "#999999")

SO_PLOT<-ggplot(data=plot_SR, aes(x=as.character(Scenar), y=as.numeric(SR)))+geom_col(aes(fill=clr), show.legend = FALSE)+ geom_linerange(aes(ymin = as.numeric(low_ci), ymax = as.numeric(upper_ci)))+ylab("SO")+xlab("Scenario")+ggtitle(paste0(surv_data, " SO among Scenarios"))+theme_bw()+ scale_fill_manual(values=colrs)

ggsave("SO_PLOT.tiff", SO_PLOT,path = plotdir )
