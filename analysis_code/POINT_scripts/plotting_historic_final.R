############################################################
##   R code for plotting soil bgc testbed historic        ##
##    results of site-calibrated models                   ##
##                                                        ##
##   Created by: Brooke Eastman, 2023                     ##
#----------------------------------------------------------#

# first run import_POINT_FEF_final.R and **RESTART R** #

print(list.filenames.CASA)
print(list.filenames.MIMICS)

# plot vegetation results        #
##  annual fluxes ##
r <- c(1:89) # all years #
cl <- c("orange2","lightseagreen")          # create color palette #
par(mfrow=c(3,2), mar=c(2,4.4,0.5,0.5))
for (i in c(2,3,5:15)){    # variables of interest included 2,3,5:15, though more can be selected
  # option to set y axis range:
  ymin=min(Veg_fluxes_2plot[[1]][r,i],Veg_fluxes_2plot[[2]][r,i])
  ymax=max(Veg_fluxes_2plot[[1]][r,i],Veg_fluxes_2plot[[2]][r,i])
  plot(Veg_fluxes_2plot[[1]][r,i],type="l",lwd=3,col=cl[1],ylim=c(ymin,ymax),ylab=colnames(Veg_fluxes_2plot[[1]][i]),cex.lab=2,cex.axis=1.8)
  for (j in 2:length(Veg_fluxes_2plot)) {
    lines(Veg_fluxes_2plot[[j]][r,i],col=cl[j],lwd=3,cex.lab=2,cex.axis=1.8)
    legend("topleft", legend=c(list.filenames.CASA[1:j]),
          col=c(cl[1:j]),
          lty=1, cex=0.8, box.lty=0)
  }
}

## mean annual vegetation pools ##
r <- c(1:89) # all years #
cl <- c("orange2","lightseagreen")          # create color palette #
par(mfrow=c(4,2), mar=c(1.5,4.4,0.5,0.5))
for (i in c(2:19)){ # default 2:19, select variables of interest #
  # option to set y axis range #
  ymin=min(Veg_pools_2plot[[1]][r,i],Veg_pools_2plot[[2]][r,i])
  ymax=max(Veg_pools_2plot[[1]][r,i],Veg_pools_2plot[[2]][r,i])
  plot(Veg_pools_2plot[[1]][r,i],type="l",col=cl[1],lwd=3,ylim=c(ymin,ymax),ylab=colnames(Veg_pools_2plot[[1]][i]),cex.lab=2,cex.axis=1.8)
  for (j in 2:length(Veg_pools_2plot)) {
    lines(Veg_pools_2plot[[j]][r,i],col=cl[j],lwd=3,cex.lab=1.2,cex.axis=1.8)
    legend("topleft", legend=c(list.filenames.CASA[1:j]),
           col=c(cl[1:j]),
           lty=1, cex=.8, box.lty=0)
  }
}



# plotting soil output so as to line up casa and mimics results (from different lists) #
## mean annual pools ##
r <- c(1:89) # all years #
lc <- length(list.filenames.CASA) - length(Msoil_fluxes_2plot)   # length of casa list that is actually casa output, not mimics #
cl <- c("orange2","lightseagreen")          # create color palette #
par(mfrow=c(3,2), mar=c(1.5,4,0.5,0.5))
for (i in c(2:23)){       #2:23 default variables of interest / select variables you want to plot #
  # option to set y axis range including all list items for CASA and MIMICS:
  ymin=min(CSoil_pools_2plot[[1]][r,i],Msoil_pools_2plot[[1]][r,i])
  ymax=max(CSoil_pools_2plot[[1]][r,i],Msoil_pools_2plot[[1]][r,i])
  plot(CSoil_pools_2plot[[1]]$iYrCnt,CSoil_pools_2plot[[1]][r,i],type="l",lwd=3,col=cl[1],ylim=c(ymin,ymax),ylab=colnames(CSoil_pools_2plot[[1]][i]),cex.lab=1.2,cex.axis=1.8)
  # plot CASA soil output #
  for (j in 1:lc) {
    lines(CSoil_pools_2plot[[j]][r,i],col=cl[j],lwd=3,cex.lab=1.2,cex.axis=1.8) 
  }
  # add MIMICS soil output #
  for (j in 1:length(Msoil_pools_2plot)) {
    lines(Msoil_pools_2plot[[j]][r,i],col=cl[j+lc],lwd=3)
  }
  j <- length(list.filenames.CASA)
  legend("topleft", legend=c(list.filenames.CASA[1:length(list.filenames.CASA)]),
         col=c(cl[1:length(list.filenames.CASA)]),
         lty=1, cex=0.8, box.lty=0, bg=NULL) 
}


## annual fluxes ##
r <- c(1:89) # all years in historical simulation #
lc <- length(list.filenames.CASA) - length(Msoil_fluxes_2plot) # length of casa list that is actually a CASA output (not MIMICS) #
cl <- c("orange2","lightseagreen")          # create color palette #
par(mfrow=c(4,2), mar=c(1.5,4,0.5,0.5))
for (i in c(2:9)){                            # select variables you want to plot
  # option to set y axis range:
  ymin=min(Csoil_fluxes_2plot[[1]][r,i],Msoil_fluxes_2plot[[1]][r,i])
  ymax=max(Csoil_fluxes_2plot[[1]][r,i],Msoil_fluxes_2plot[[1]][r,i])
  # plot CASA soil fluxes #
  plot(Csoil_fluxes_2plot[[1]]$iYrCnt,Csoil_fluxes_2plot[[1]][r,i],type="l",lwd=3,col=cl[1],ylim=c(ymin,ymax),ylab=colnames(Csoil_fluxes_2plot[[1]][i]),cex.lab=1.2,cex.axis=1.8)
  # add MIMICS soil fluxes #
  for (j in 1:length(Msoil_fluxes_2plot)) {
    lines(Msoil_fluxes_2plot[[j]][r,i],col=cl[j+lc],lwd=3)
  }
  j <- length(list.filenames.CASA)
  legend("topleft", legend=c(list.filenames.CASA[1:length(list.filenames.CASA)]),
         col=c(cl[1:length(list.filenames.CASA)]),
         lty=1, cex=0.8, box.lty=0, bg=NULL) 
}


# plot MIMICS microbial output data #
par(mfrow=c(3,2), mar=c(1.5,4,0.5,0.5))
cl <- c("navy","lightseagreen","purple1","hotpink","limegreen")  # select new color palette - MIMICS only #
for (i in c(22:28)){
  # define y axis range #
  ymin=min(Msoil_pools_2plot[[1]][r,i])
  ymax=max(Msoil_pools_2plot[[1]][r,i])
  plot(Msoil_pools_2plot[[1]]$iYrCnt,Msoil_pools_2plot[[1]][r,i],type="l",lwd=3,col=cl[1],ylim=c(ymin,ymax),ylab=colnames(Msoil_pools_2plot[[1]][i]))
  for(j in 1:length(Msoil_pools_2plot)) {
    lines(Msoil_pools_2plot[[j]][r,i],col=cl[j],lwd=3)
  }
  legend("topleft", legend=c(list.filenames.MIMICS[1:length(list.filenames.MIMICS)]),
         col=c(cl[1:length(list.filenames.MIMICS)]),
         lty=1, cex=0.8, box.lty=0, bg=NULL) 
}




###############################################
####             DAILY OUTPUT           #######
# 10-year mean intrannual variability #
par(mfrow=c(4,1), mar=c(1.5,4,0.5,0.5))
y <- c(28836:32485)
cl <- c("orange2","lightseagreen") # create color palette #
axis <- c("Min Soil N","Nmin Uptake","CrSoil","NPP")
for (j in c(46,47,25,34)){  # select variables. default (46,47,25,34) refer to Nsoilpool, Nuptake, Crsoil, and NPP #
  ymin=min(CASAtenmean[[1]][[j]],CASAtenmean[[2]][[j]])
  ymax=max(CASAtenmean[[1]][[j]],CASAtenmean[[2]][[j]])
  plot(CASAtenmean[[1]][[j]],type="l",lwd=2,col=cl[1],ylim=c(ymin,ymax),ylab=colnames(CASAtenmean[[1]][j]))
  for (i in 2:length(CASAtenmean)) {
    lines(CASAtenmean[[i]][[j]],type="l",lwd=2,col=cl[i])
    #legend("topleft", legend=c(list.filenames.CASA[1:length(list.filenames.CASA)]),
     #      col=c(cl[1:length(list.filenames.CASA)]),
      #     lty=1, cex=0.8, box.lty=0) 
  }
}







######################################################################
### estimate mean annual pools and fluxes for last ten years to make a 
##### table comparing model and observations at steady state (ish)

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
cruns <- 1
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
mruns <- 1
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



###################################################################################
# compare observed soil T and moisture with GSWP3 forcing data #

# import observed met data: soil T and Soil Moisture: #
library(readxl)
# import observed soil temperature data (see "soilThourly1yr.xslx" in analysis_code/ directory) #
soilThourly1yr <- read_excel("~/Desktop/NCAR-CLM/DATA/soilThourly1yr.xlsx")
# take daily average over 24 hourly measurements #
soilTdaily <- tapply(soilThourly1yr$Mean,rep(seq_along(soilThourly1yr$Mean),each=24,length.out=length(soilThourly1yr$Mean)),mean)
DOY <- c(1:365)
ObsSoilT <- as.data.frame(cbind(DOY,soilTdaily))

#import observed soil moisture data (see "ObsSoilMoisture.xslx" in analysis_code/ directory) #
SoilMoistureObs <- read_excel("~/Desktop/NCAR-CLM/DATA/ObsSoilMoisture.xlsx")
moistmean <- tapply(SoilMoistureObs$Moisture,SoilMoistureObs$DOY,mean,na.rm=TRUE)
DOYmois <- tapply(SoilMoistureObs$DOY,SoilMoistureObs$DOY,mean)



################  PLOT  ##########################################
# select/plot ten year average DOY annual cycles from GSW3P #
## met data: soil temperature, moisture ##
dev.off()
r <- c(5,6,26)
cl <- c("maroon","dodgerblue3")            # create a color palette #
par(mar = c(5, 4, 2, 4) + 0.3)  # Leave space for z axis #
plot(CASAtenmean[[1]][[1]],CASAtenmean[[1]][[5]]-273.15,type="l",lwd=3,col=cl[1],xlab="DOY",ylab="Soil Temperature (C)",cex.lab=2,cex.axis=1.8) 
points(ObsSoilT$DOY, ObsSoilT$soilTdaily,col="maroon")
par(new = TRUE)
plot(DOYmois,moistmean/100, col=cl[2], pch=16, axes = FALSE, bty = "n", xlab = "", ylab = "")
lines(CASAtenmean[[1]][[1]], CASAtenmean[[1]][[6]],lwd=3,col=cl[2])
axis(side=4, at = c(0,0.1,0.2,0.3,0.4),tick=TRUE,cex.lab=2,cex.axis=1.8)
mtext("Soil Moisture", side=4, line=3,cex = 2)





#########  DATA EXPLORATION  ########################################################################
###############################################################
## plotting C and N balances from balance .csv output files ###
###############################################################


################## C BALANCE #######################################################
# plot  10-year mean intrannual variability #
par(mfrow=c(4,1), mar=c(1.5,4,0.5,0.5))
cl <- c("orange2","lightseagreen")            # create a color palette #
for (j in 4:12){  # these refer to the output variables in the C balance
  ymin=min(CBALANCEtenmean[[1]][[j]],CBALANCEtenmean[[2]][[j]])#,CBALANCEtenmean[[3]][[j]],CBALANCEtenmean[[5]][[j]],CBALANCEtenmean[[6]][[j]])
  ymax=max(CBALANCEtenmean[[1]][[j]],CBALANCEtenmean[[2]][[j]])#,CBALANCEtenmean[[3]][[j]],CBALANCEtenmean[[5]][[j]],CBALANCEtenmean[[6]][[j]])
  plot(CBALANCEtenmean[[1]][[j]],type="l",lwd=2,col=cl[1],ylim=c(ymin,ymax),ylab=colnames(CBALANCEtenmean[[1]][j]))
  for (i in 2:length(CBALANCEtenmean)) {
    lines(CBALANCEtenmean[[i]][[j]],type="l",lwd=2,col=cl[i])
    legend("topleft", legend=c(list.filenames.CBALANCE[1:length(list.filenames.CBALANCE)]),
           col=c(cl[1:length(list.filenames.CBALANCE)]),
           lty=1, cex=0.8, box.lty=0, bg=NULL) 
  }
}

# plot daily output from last ten years of historic run 
par(mfrow=c(4,1), mar=c(1.5,4,0.5,0.5))
y <- c(28836:32485)
for (j in c(4:12)){  # plot all output variables (columns)
  ymin=min(list.CBALANCE[[1]][[j]][y],list.CBALANCE[[2]][[j]][y])
  ymax=max(list.CBALANCE[[1]][[j]][y],list.CBALANCE[[2]][[j]][y])
  plot(list.CBALANCE[[1]][[j]][y],type="l",lwd=2,col=cl[1],ylim=c(ymin,ymax),ylab=colnames(list.CBALANCE[[1]][j]))
  for (i in 2:length(list.CBALANCE)) {
    lines(list.CBALANCE[[i]][[j]][y],type="l",lwd=2,col=cl[i])
    legend("topleft", legend=c(list.filenames.CBALANCE[1:length(list.filenames.CBALANCE)]),
           col=c(cl[1:length(list.filenames.CBALANCE)]),
           lty=1, cex=0.8, box.lty=0, bg=NULL) 
  }
}




################## N BALANCE #######################################################
# plot  10-year mean intrannual variability #
par(mfrow=c(4,1), mar=c(1.5,4,0.5,0.5))
cl <- c("orange2","lightseagreen")            # create a color palette #
for (j in 4:12){  # these refer to the output variables in the C balance
  ymin=min(NBALANCEtenmean[[1]][[j]],NBALANCEtenmean[[2]][[j]])
  ymax=max(NBALANCEtenmean[[1]][[j]],NBALANCEtenmean[[2]][[j]])
  plot(NBALANCEtenmean[[1]][[j]],type="l",lwd=2,col=cl[1],ylim=c(ymin,ymax),ylab=colnames(NBALANCEtenmean[[1]][j]))
  for (i in 2:length(NBALANCEtenmean)) {
    lines(NBALANCEtenmean[[i]][[j]],type="l",lwd=2,col=cl[i])
    legend("topleft", legend=c(list.filenames.NBALANCE[1:length(list.filenames.NBALANCE)]),
           col=c(cl[1:length(list.filenames.NBALANCE)]),
           lty=1, cex=0.8, box.lty=0, bg=NULL) 
  }
}

# plot daily output from last ten years of historic run 
par(mfrow=c(4,1), mar=c(1.5,4,0.5,0.5))
y <- c(28836:32485)
for (j in c(4:12)){  # plot all output variables (columns)
  ymin=min(list.NBALANCE[[1]][[j]][y],list.NBALANCE[[2]][[j]][y])
  ymax=max(list.NBALANCE[[1]][[j]][y],list.NBALANCE[[2]][[j]][y])
  plot(list.NBALANCE[[1]][[j]][y],type="l",lwd=2,col=cl[1],ylim=c(ymin,ymax),ylab=colnames(list.NBALANCE[[1]][j]))
  for (i in 2:length(list.NBALANCE)) {
    lines(list.NBALANCE[[i]][[j]][y],type="l",lwd=2,col=cl[i])
    legend("topleft", legend=c(list.filenames.NBALANCE[1:length(list.filenames.NBALANCE)]),
           col=c(cl[1:length(list.filenames.NBALANCE)]),
           lty=1, cex=0.8, box.lty=0, bg=NULL) 
  }
}




