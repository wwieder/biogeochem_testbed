############################################################
##   R code for plotting soil bgc testbed experimental    ##
##    results of N perturbations                          ##
##                                                        ##
##   Created by: Brooke Eastman, 2023                     ##
#----------------------------------------------------------#

# first run import_POINT_FEF_final.R and **RESTART R** #


print(list.filenames.CASA)
print(list.filenames.MIMICS)

# plot vegetation results        #
##  annual fluxes ##
r <- c(1:31) # select years (all years in 31-year experiment) #
# create color palette #
cl <- c("darkorange4","wheat3","yellow","orange2","darkolivegreen3","navy","darkorchid4","lightseagreen")
par(mfrow=c(3,2), mar=c(2,4.4,0.5,0.5))
for (i in c(2,3,5:17)){    # select variables to plot. default 2,3,5:17 #
  # option to set y axis range:
  ymin=min(Veg_fluxes_2plot[[1]][r,i],Veg_fluxes_2plot[[2]][r,i],Veg_fluxes_2plot[[3]][r,i],Veg_fluxes_2plot[[4]][r,i],Veg_fluxes_2plot[[5]][r,i],Veg_fluxes_2plot[[6]][r,i],Veg_fluxes_2plot[[7]][r,i],Veg_fluxes_2plot[[8]][r,i])
  ymax=max(Veg_fluxes_2plot[[1]][r,i],Veg_fluxes_2plot[[2]][r,i],Veg_fluxes_2plot[[3]][r,i],Veg_fluxes_2plot[[4]][r,i],Veg_fluxes_2plot[[5]][r,i],Veg_fluxes_2plot[[6]][r,i],Veg_fluxes_2plot[[7]][r,i],Veg_fluxes_2plot[[8]][r,i])
  plot(Veg_fluxes_2plot[[1]][r,i],type="l",lwd=3,col=cl[1],ylim=c(ymin,ymax),ylab=colnames(Veg_fluxes_2plot[[1]][i]),cex.lab=2,cex.axis=1.8)
  for (j in 2:length(Veg_fluxes_2plot)) {
    lines(Veg_fluxes_2plot[[j]][r,i],col=cl[j],lwd=3,cex.lab=2,cex.axis=1.8)
    legend("topleft", legend=c(list.filenames.CASA[1:j]),
          col=c(cl[1:j]),
          lty=1, cex=0.8,  bg=NULL, box.lty=0)
  }
}

## mean annual pools ##
r <- c(1:31) # select years (all years in 31-year experiment) #
# create color palette #
cl <- c("darkorange4","wheat3","yellow","orange2","darkolivegreen3","navy","darkorchid4","lightseagreen")
par(mfrow=c(4,2), mar=c(1.5,4.4,0.5,0.5))
for (i in c(2:20)){ # select variables to plot. default 2:20
  # option to set y axis range:
  ymin=min(Veg_pools_2plot[[1]][r,i],Veg_pools_2plot[[2]][r,i],Veg_pools_2plot[[3]][r,i],Veg_pools_2plot[[4]][r,i],Veg_pools_2plot[[5]][r,i],Veg_pools_2plot[[6]][r,i],Veg_pools_2plot[[7]][r,i],Veg_pools_2plot[[8]][r,i])
  ymax=max(Veg_pools_2plot[[1]][r,i],Veg_pools_2plot[[2]][r,i],Veg_pools_2plot[[3]][r,i],Veg_pools_2plot[[4]][r,i],Veg_pools_2plot[[5]][r,i],Veg_pools_2plot[[6]][r,i],Veg_pools_2plot[[7]][r,i],Veg_pools_2plot[[8]][r,i])
  plot(Veg_pools_2plot[[1]][r,i],type="l",col=cl[1],lwd=3,ylim=c(ymin,ymax),ylab=colnames(Veg_pools_2plot[[1]][i]),cex.lab=2,cex.axis=1.8)
  for (j in 2:length(Veg_pools_2plot)) {
    lines(Veg_pools_2plot[[j]][r,i],col=cl[j],lwd=3,cex.lab=1.2,cex.axis=1.8)
    legend("topleft", legend=c(list.filenames.CASA[1:j]),
           col=c(cl[1:j]),
           lty=1, cex=.8,  bg=NULL, box.lty=0)
  }
}



### plotting soil so as to line up casa and mimics results (from two different list items) ###
## mean annual soil pools ##
r <- c(1:31)
lc <- length(list.filenames.CASA) - length(Msoil_fluxes_2plot)   # length of casa list that is actually casa output, not mimics #
# create color palette #
cl <- c("darkorange4","wheat3","yellow","orange2","darkolivegreen3","navy","darkorchid4","lightseagreen")
par(mfrow=c(3,2), mar=c(1.5,4,0.5,0.5))
for (i in c(2:23)){       # select variables to plot.  default 2:23 #
  # option to set y axis range:
  ymin=min(CSoil_pools_2plot[[1]][r,i],CSoil_pools_2plot[[2]][r,i],CSoil_pools_2plot[[3]][r,i],CSoil_pools_2plot[[4]][r,i],Msoil_pools_2plot[[1]][r,i],Msoil_pools_2plot[[2]][r,i],Msoil_pools_2plot[[3]][r,i],Msoil_pools_2plot[[4]][r,i])
  ymax=max(CSoil_pools_2plot[[1]][r,i],CSoil_pools_2plot[[2]][r,i],CSoil_pools_2plot[[3]][r,i],CSoil_pools_2plot[[4]][r,i],Msoil_pools_2plot[[1]][r,i],Msoil_pools_2plot[[2]][r,i],Msoil_pools_2plot[[3]][r,i],Msoil_pools_2plot[[4]][r,i])
  plot(CSoil_pools_2plot[[1]]$iYrCnt,CSoil_pools_2plot[[1]][r,i],type="l",lwd=3,col=cl[1],ylim=c(ymin,ymax),ylab=colnames(CSoil_pools_2plot[[1]][i]),cex.lab=1.2,cex.axis=1.8)
  for (j in 2:lc) {
    lines(CSoil_pools_2plot[[j]][r,i],col=cl[j],lwd=3,cex.lab=1.2,cex.axis=1.8) 
    }
  for (j in 1:length(Msoil_pools_2plot)) {
    lines(Msoil_pools_2plot[[j]][r,i],col=cl[j+lc],lwd=3,cex.lab=1.2,cex.axis=1.8)
    }
  legend("topleft", legend=c(list.filenames.CASA[1:length(list.filenames.CASA)]),
        col=c(cl[1:length(list.filenames.CASA)]),
        lty=1, cex=.8,  bg=NULL, box.lty=0)
}

## annual soil fluxes ##
r <- c(1:31)
lc <- length(list.filenames.CASA) - length(Msoil_fluxes_2plot)   # length of casa list that is actually casa output, not mimics #
# create color palette #
cl <- c("darkorange4","wheat3","yellow","orange2","darkolivegreen3","navy","darkorchid4","lightseagreen")
par(mfrow=c(4,2), mar=c(1.5,4,0.5,0.5))
for (i in c(2:9)){           # select variables you want to plot. default 2:9 #
  # option to set y axis range:
  ymin=min(Csoil_fluxes_2plot[[1]][r,i],Csoil_fluxes_2plot[[2]][r,i],Csoil_fluxes_2plot[[3]][r,i],Csoil_fluxes_2plot[[4]][r,i],Msoil_fluxes_2plot[[1]][r,i],Msoil_fluxes_2plot[[2]][r,i],Msoil_fluxes_2plot[[3]][r,i],Msoil_fluxes_2plot[[4]][r,i])
  ymax=max(Csoil_fluxes_2plot[[1]][r,i],Csoil_fluxes_2plot[[2]][r,i],Csoil_fluxes_2plot[[3]][r,i],Csoil_fluxes_2plot[[4]][r,i],Msoil_fluxes_2plot[[1]][r,i],Msoil_fluxes_2plot[[2]][r,i],Msoil_fluxes_2plot[[3]][r,i],Msoil_fluxes_2plot[[4]][r,i])
  plot(Csoil_fluxes_2plot[[1]]$iYrCnt,Csoil_fluxes_2plot[[1]][r,i],type="l",lwd=3,col=cl[1],ylim=c(ymin,ymax),ylab=colnames(Csoil_fluxes_2plot[[1]][i]),cex.lab=1.2,cex.axis=1.8)
  for (j in 2:lc) {
    lines(Csoil_fluxes_2plot[[j]][r,i],col=cl[j],lwd=3) 
  }
  for (j in 1:length(Msoil_fluxes_2plot)) {
    lines(Msoil_fluxes_2plot[[j]][r,i],col=cl[j+lc],lwd=3)
  }
  j <- length(list.filenames.CASA)
  legend("topleft", legend=c(list.filenames.CASA[1:length(list.filenames.CASA)]),
         col=c(cl[1:length(list.filenames.CASA)]),
         lty=1, cex=0.8, box.lty=0, bg=NULL) 
}

# plot MIMICS soil microbial results #
par(mfrow=c(3,2), mar=c(1.5,4,0.5,0.5))
cl <- c("darkolivegreen3","navy","darkorchid4","lightseagreen") # create mimics-only color palette #
for (i in c(24:28)){
  ymin=min(Msoil_pools_2plot[[1]][r,i],Msoil_pools_2plot[[2]][r,i],Msoil_pools_2plot[[3]][r,i],Msoil_pools_2plot[[4]][r,i])
  ymax=max(Msoil_pools_2plot[[1]][r,i],Msoil_pools_2plot[[2]][r,i],Msoil_pools_2plot[[3]][r,i],Msoil_pools_2plot[[4]][r,i])
  plot(Msoil_pools_2plot[[1]]$iYrCnt,Msoil_pools_2plot[[1]][r,i],type="l",lwd=3,col=cl[1],ylim=c(ymin,ymax),ylab=colnames(Msoil_pools_2plot[[1]][i]))
  for(j in 2:length(Msoil_pools_2plot)) {
    lines(Msoil_pools_2plot[[j]][r,i],col=cl[j],lwd=3)
  }
  legend("topleft", legend=c(list.filenames.MIMICS[1:length(list.filenames.MIMICS)]),
         col=c(cl[1:length(list.filenames.MIMICS)]),
         lty=1, cex=0.8, box.lty=0, bg=NULL) 
}



###############################################
####      DAILY OUTPUT                  #######

# plot 10-year mean daily variability #
## NOTE: best data visualization when 2 simulations are selected   ##
##       here, we select the default calibrated experimental (+N)  ##
##       for CASA (CASAtenmean[[1]]) and MIMICS (CASAtenmean[[5]]) ##

par(mfrow=c(4,1), mar=c(1.5,4,0.5,0.5))
y <- c(7665:11315)
cl <- c("darkorange4","darkolivegreen3")
axis <- c("Min Soil N","Nmin Uptake","CrSoil","NPP")
for (j in c(46,47,25,34)){  # these refer to Nsoilpool, Nuptake, soil respiration, and NPP
  ymin=min(CASAtenmean[[1]][[j]],CASAtenmean[[5]][[j]])
  ymax=max(CASAtenmean[[1]][[j]],CASAtenmean[[5]][[j]])
  plot(CASAtenmean[[1]][[j]],type="l",lwd=2,col=cl[1],ylim=c(ymin,ymax),ylab=colnames(CASAtenmean[[1]][j]))
  for (i in 2:2) {
    lines(CASAtenmean[[i]][[j]],type="l",lwd=2,col=cl[i])
    abline(h=0, col= "gray60", lty=3)
    legend("topleft", legend=c(list.filenames.CASA[1:length(list.filenames.CASA)]),
          col=c(cl[1:length(list.filenames.CASA)]),
           lty=1, cex=0.8, box.lty=0) 
  }
}

# comnpare 10 year mean (above) to final year
CASA2019 <- lapply(list.CASA2, '[', c(10950:11315),)            # subset last year

par(mfrow=c(4,1), mar=c(2.1,4,0.5,0.5))
y <- c(7665:11315)
cl <- c("darkorange4","darkolivegreen3")
axis <- c("Min Soil N","Nmin Uptake","Nleaching","NPP")
for (j in c(46,47,51,34)){  # these refer to Nsoilpool, Nuptake, xkNlimiting, and NPP
  ymin=min(CASA2019[[1]][[j]],CASA2019[[5]][[j]])
  ymax=max(CASA2019[[1]][[j]],CASA2019[[5]][[j]])
  plot(CASA2019[[1]][[j]],type="l",lwd=2,col=cl[1],ylim=c(ymin,ymax),ylab=colnames(CASA2019[[1]][j]),
       cex.axis=2)
  for (i in 2:length(CASA2019)) {
    lines(CASA2019[[i]][[j]],type="l",lwd=2,col=cl[i])
    abline(h=0, col= "gray60", lty=3)
    #legend("topleft", legend=c(list.filenames.CASA[1:length(list.filenames.CASA)]),
    #col=c(cl[1:length(list.filenames.CASA)]),
    #lty=1, cex=0.8, box.lty=0) 
  }
}



##############################################################################
### estimate mean annual pools and fluxes for last ten years               ### 
#####  for experimental N addition response ratio estimates and figures  #####

library(tidyverse)
library(data.table)
library(stringr)

veg_pools_summary <- map(Veg_pools_2plot, ~.x %>% filter(iYrCnt>=22) %>%
                           summarise_if(is.numeric,list(mean,sd)))
veg_pools_summary <- do.call(rbind.data.frame,veg_pools_summary)
veg_pools_summary <- veg_pools_summary %>% rownames_to_column('simulation')
#pull out columns of interest:
veg_pools_summary1 <- veg_pools_summary[,c(1,3:8,14:18)]

veg_flux_summary <- map(Veg_fluxes_2plot, ~.x %>% filter(iYrCnt>=22) %>%
                           summarise_if(is.numeric,list(mean,sd)))
veg_flux_summary <- do.call(rbind.data.frame,veg_flux_summary)
#pull out columns of interest:
veg_flux_summary1 <- veg_flux_summary[,c(3,8,9,13)]

write.csv(veg_flux_summary,"veg_flux_summary.csv")

#set number of casa runs
cruns <- 3
casa_spools_summary <- map(CSoil_pools_2plot, ~.x %>% filter(iYrCnt>=22) %>%
                           summarise_if(is.numeric,list(mean,sd)))
casa_spools_summary <- do.call(rbind.data.frame,casa_spools_summary)
casa_spools_summary <- casa_spools_summary %>% rownames_to_column('simulation')
#pull out columns of interest:
casa_spools_summary1 <- casa_spools_summary[1:cruns,c(1,21:24)]

casa_sfluxes_summary <- map(Csoil_fluxes_2plot, ~.x %>% filter(iYrCnt>=22) %>%
                           summarise_if(is.numeric,list(mean,sd)))
casa_sfluxes_summary <- do.call(rbind.data.frame,casa_sfluxes_summary)
#pull out columns of interest:
casa_sfluxes_summary1 <- casa_sfluxes_summary[1:cruns,c(2)]

# extract soil mineral N
casa_NsoilMin <- map(list.CASA2, ~.x %>% filter(iYrCnt>=22) %>%
                             summarise_if(is.numeric,mean))
casa_NsoilMin <- do.call(rbind.data.frame,casa_NsoilMin)
casa_NsoilMin <- casa_NsoilMin %>% rownames_to_column('simulation')
casa_NsoilMin <- casa_NsoilMin[1:cruns,c(46)]

#set number of mimics runs
mruns <- 3
mimics_spools_summary <- map(Msoil_pools_2plot, ~.x %>% filter(iYrCnt>=22) %>%
                             summarise_if(is.numeric,list(mean,sd)))
mimics_spools_summary <- do.call(rbind.data.frame,mimics_spools_summary)
mimics_spools_summary <- mimics_spools_summary %>% rownames_to_column('simulation')
#pull out columns of interest:
# include mic r:k ratio (Column 29)
mimics_spools_summary1_r2k <- mimics_spools_summary[1:mruns,c(1,21:24,29)]
#exclude r:k ratio (to be able to match and bind to CASA soil values)
mimics_spools_summary1 <- mimics_spools_summary[1:mruns,c(1,21:24)]


mimics_sfluxes_summary <- map(Msoil_fluxes_2plot, ~.x %>% filter(iYrCnt>=22) %>%
                              summarise_if(is.numeric,list(mean,sd)))
mimics_sfluxes_summary <- do.call(rbind.data.frame,mimics_sfluxes_summary)
#pull out columns of interest:
mimics_sfluxes_summary1 <- mimics_sfluxes_summary[1:mruns,c(2)]

#add root flux (to compare to TBCF)
rootflux_summary <- veg_flux_summary[,c(17)]


#combine data frames to one
veg_summary <- cbind(veg_pools_summary1,veg_flux_summary1)
casa_ssummary <- cbind(casa_spools_summary1,casa_sfluxes_summary1)
mimics_ssummary <- cbind(mimics_spools_summary1,mimics_sfluxes_summary1)
#rename casaflux Crsoil to match the column name in mimics df for binding
colnames(casa_ssummary)[colnames(casa_ssummary) == 'casaflux%Crsoil_fn1'] <- 'hresp_fn1'
soil_summary <- rbind(casa_ssummary,mimics_ssummary)
soil_summary <- cbind(soil_summary,rootflux_summary)
  
all_summary <- cbind(veg_summary,soil_summary)

#  NOTE:  I also calculated aboveground NPP from the sum of leaf C flux (Cleaf2met + Cleaf2str)
## and wood C increment from previous year  ### using the allpoolsCASA wood C pool data
write.csv(allpoolsCASA,"allpoolsCASA.csv",row.names = F)
write.csv(allannualCASA,"allfluxCASA.csv",row.names=F)
# calculated wood C increment, and added that to leaf flux
CASA_ANPP <- read.csv("allpoolsCASA.csv")
CASA_ANPP <- CASA_ANPP %>% group_by(column_label)  %>% filter(iYrCnt>20) %>% summarise_at('ANPP',mean)

##############################################################
# write summary files to csvs.  #

write.csv(all_summary,"all_summary_forRR_012423.csv",row.names = FALSE)

#### !!! NOTES: I edited the above in excel, where I calculate response ratios,  ####
### !!! add in the observational values, and put into a format ###
## !!! usable for the lattice forest plots below ##


#######################################################################
#---------------------------------------------------------------------#
## Create latticed forest plots of response ratio for obs & sims     ##
#######################################################################
# 
library(tidyverse)
library(scales)

#  read in data in proper format. see: https://datascienceplus.com/lattice-like-forest-plot-using-ggplot2-in-r/ 


###################################################################################
## !!!!!!!!!!!!!! FINAL VERSION !!!!!!!!!!!! ##
## selected key veg responses and soil responses to plot all models in two figures
#
# import vegetation response ratio data: "RR_veg_FINAL_data_ANPP_111822.csv" #
RR_veg_data <- read.csv("C:/Users/Brooke/Desktop/SOIL_BGC_TESTBED/FEF_Naddexp2/RR_veg_FINAL_data_ANPP_111822.csv")
# lock in factor level order #
RR_veg_data$Ecosystem.Property <- factor(RR_veg_data$Ecosystem.Property, levels = unique(RR_veg_data$Ecosystem.Property))
override.shape <- c(16,14,17,8,14,17)

# set label names for legend and axes #
RR_veg_data <- RR_veg_data %>% mutate(Model = factor(Model, levels=c("CASA enzyme inhibition","CASA allocation shift","CASA","MIMICS enzyme inhibition","MIMICS allocation shift","MIMICS","Observed")))
RR_veg_data <- RR_veg_data %>% mutate(Ecosystem.Property = factor(Ecosystem.Property, levels=c("ANPP","Leaf C","Wood C","Fine root C")))


######## PLOT VEGETATION RR #################

RR_FINAL_veg = ggplot(data=RR_veg_data, aes(x = Ecosystem.Property,y = Response.Ratio, ymin = Lower.Limit, ymax = Upper.Limit ))+
  geom_pointrange(aes(col=Model,shape=Model), size=1, position=position_dodge(width=0.6))+ 
  scale_color_manual(name = "Model",
                     labels= c("CASA enzyme inhibition","CASA allocation shift","CASA","MIMICS enzyme inhibition","MIMICS allocation shift","MIMICS","Observed"),
                     values=c("darkorange4","darkorange4","darkorange4","blue","blue","blue","black")) +
  scale_shape_manual(name = "Model",
                     labels= c("CASA enzyme inhibition","CASA allocation shift","CASA","MIMICS enzyme inhibition","MIMICS allocation shift","MIMICS","Observed"),
                     values=c(8,7,17,8,7,17,16)) +
  geom_hline(aes(fill=Model),yintercept =1, linetype=2)+
  geom_vline(xintercept = seq(from=1.5, to=5.5, by = 1)) +
  scale_x_discrete(limits=rev) + ylim(0.5,1.75) +
  labs(color = "Model", shape ="Model", linetype = "Model") +
  xlab('Ecosystem Property\n')+ ylab("N Addition Response Ratio")+
  geom_errorbar(aes(ymin=Lower.Limit, ymax=Upper.Limit,col=Model),width=0.8,cex=1, 
                show.legend = FALSE, position=position_dodge(width=0.6))+ 
  #facet_wrap(~Ecosystem.Property,strip.position="left",nrow=3,scales = "free_y") +
  theme_bw() +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.ticks.y=element_blank(),
        axis.text=element_text(size=20,color="black"),
        axis.title=element_text(size=22,face="bold"),
        legend.text=element_text(size=16),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
  #guides(col = guide_legend(reverse = TRUE)) +
  coord_flip()
RR_FINAL_veg


## repeat for soil response ratios ##

# read in soil response ratio data, "RR_soil_FINAL_data_ALL_110422.csv" #
RR_soil_data <- read.csv("C:/Users/Brooke/Desktop/SOIL_BGC_TESTBED/FEF_Naddexp2/RR_soil_FINAL_data_ALL_110422.csv")
#lock in factor level order
RR_soil_data$Ecosystem.Property <- factor(RR_soil_data$Ecosystem.Property, levels = unique(RR_soil_data$Ecosystem.Property))
override.shape <- c(16,14,17,8,14,17)
# define label names for legend and axes #
RR_soil_data <- RR_soil_data %>% mutate(Model = factor(Model, levels=c("CASA enzyme inhibition","CASA allocation shift","CASA","MIMICS enzyme inhibition","MIMICS allocation shift","MIMICS","Observed")))


####################### PLOT SOIL RR ###########################

RR_FINAL_soil = ggplot(data=RR_soil_data, aes(x = Ecosystem.Property,y = Response.Ratio, ymin = Lower.Limit, ymax = Upper.Limit ))+
  geom_pointrange(aes(col=Model,shape=Model), size=1, position=position_dodge(width=0.6))+ 
  scale_color_manual(name = "Model",
                     labels= c("CASA enzyme inhibition","CASA allocation shift","CASA","MIMICS enzyme inhibition","MIMICS allocation shift","MIMICS","Observed"),
                     values=c("darkorange4","darkorange4","darkorange4","blue","blue","blue","black")) +
  scale_shape_manual(name = "Model",
                     labels= c("CASA enzyme inhibition","CASA allocation shift","CASA","MIMICS enzyme inhibition","MIMICS allocation shift","MIMICS","Observed"),
                     values=c(8,7,17,8,7,17,16)) +
  geom_hline(aes(fill=Model),yintercept =1, linetype=2)+
  geom_vline(xintercept = seq(from=1.5, to=5.5, by = 1)) +
  #scale_y_continuous(labels=label_number(accuracy=0.01)) +
  scale_x_discrete(limits=rev) + ylim(0.625,1.5) +
  labs(color = "Model", shape ="Model", linetype = "Model") +
  xlab('Ecosystem Property')+ ylab("N Addition Response Ratio")+
  geom_errorbar(aes(ymin=Lower.Limit, ymax=Upper.Limit,col=Model),width=0.8,cex=1, 
                show.legend = FALSE, position=position_dodge(width=0.6))+ 
  #facet_wrap(~Ecosystem.Property,strip.position="left",nrow=3,scales = "free_y") +
  theme_bw() +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.ticks.y=element_blank(),
        axis.text=element_text(size=20,color="black"),
        axis.title=element_text(size=22,face="bold"),
        legend.text=element_text(size=16),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
  #guides(col = guide_legend(reverse = TRUE)) +
  annotate("text",x=0.1,y=13,label="(b)") +
  coord_flip()
RR_FINAL_soil

# save plots #
ggsave("RR_veg_ALL_111822.png", plot=RR_FINAL_veg, device = "png",
       height=7, width=10, units="in", dpi=500)
ggsave("RR_soil_ALL_111622.png", plot=RR_FINAL_soil, device = "png",
       height=7, width=10, units="in", dpi=500)








######################################################################
## Plot comparing observed relationship between fraction C in light ##
##     fraction (POM) vs bulk mineral soil C:N ratio                ##
######################################################################

  ##########################################################
  
  #import csv with data from observed AND modeled, combined: "FractionsBarPlotCNvLF_110522.csv"
  CNvLF_data <- read.csv("FractionsBarPlotCNvLF_110522.csv")
  
  my_theme2 <- theme(axis.line.x = element_line(size = 0.5, colour = "black"), axis.line.y = element_line(size = 0.5, colour = "black"),
                     axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), panel.border = element_blank(),  panel.background = element_blank(), text=element_text(size = 14),
                     axis.text.x=element_text(colour="black", size = 20),
                     axis.text.y=element_text(colour="black", size = 20), 
                     axis.title=element_text(size=20,face="bold"),
                     legend.title = element_blank(),
                     legend.text=element_text(size=14),
                     legend.position = c(0.85,0.2),
                     legend.key=element_blank())
  
  LF_v_CN <- ggplot(CNvLF_data, aes(FractionTotalC,SoilCN)) + 
    geom_smooth(data=CNvLF_data[10:29,],method = "lm", col = "black") +
    geom_smooth(data=CNvLF_data[1:4,],method = "lm", col = "darkorange4", linetype = "dashed", se = FALSE) +
    geom_smooth(data=CNvLF_data[5:8,],method = "lm", col = "blue", linetype = "dashed", se = FALSE) +
    geom_point(aes(shape=Model,col=Model),size=5) +
    scale_shape_manual(values=c(2,17,7,8,2,17,7,8,19,1)) +
    scale_color_manual(values=c("darkorange4","darkorange4","darkorange4","darkorange4",
                                "blue","blue","blue","blue",
                                "black","black")) +
    my_theme2 +  theme(legend.position="bottom",
                       axis.text=element_text(size=20)) + 
    scale_y_continuous(name="C:N ratio of bulk soil\n") +
    xlab("\nLight fraction proportion of total soil C") + 
    labs(fill="Model")
  LF_v_CN
  
  # save figure #
  ggsave("LFvCN_110522.png", LF_v_CN, device = "png", dpi=500, height=8, width = 9, units="in")

  

###########################################################
  #################### END ###########################
###########################################################



