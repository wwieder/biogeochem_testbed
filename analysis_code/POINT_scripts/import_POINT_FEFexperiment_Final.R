############################################################
## R code for importing soil bgc testbed result           ##
##    POINT*.csv files as lists from EXPERIMENTAL         ##
##    N perturbation simulations.                         ##
##                      AND                               ##
##    Calculate                                           ##
##    additional variables of interest for analysis       ##
##    and figure-making                                   ##
##                                                        ##
##   Created by: Brooke Eastman, 2023                     ##
#----------------------------------------------------------#

remove(list=ls()) # clear environment #

# set working directory #
## NOTE: I used different directories for the historic simulations ("FEF_cal2") and the experimental simulations ("FEF_Naddexp2") ##
setwd("C:/Users/Brooke/Desktop/SOIL_BGC_TESTBED/FEF_Naddexp2")
library(tidyverse)
library(data.table)
library(stringr)

# read in all .csv files into two lists: list.CASA (CASA and CASA-MIMICS), and list.MIMICS  #
# note: "POINT" files are output from testbed simulations, and were renamed in terminal after each simulation #
##       CASA-MIMICS point files were renamed from POINT_casa*.csv to POINT_casamimics*.csv to distinguish from CASA-only simulations #
list.filenames.CASA <- list.files(pattern="POINT_casa",full.names = F)
list.CASA <- list()

## for some reason, R is detecting an extra column ##
## drop and fill is best workaround to kep column names and all data, it adds a column (V55) of NAs at the end  ##
y <- for (i in 1:length(list.filenames.CASA)){
  list.CASA[[i]] <- fread(list.filenames.CASA[i],drop=1, fill=TRUE,colClasses = list(numeric = 1:66))}    
names(list.CASA) <- list.filenames.CASA
##  so delete the last column from the items with an extra column                                        ##
for (i in 1:length(list.CASA)){
  list.CASA[[i]] <- list.CASA[[i]][,1:65]}

# print name of all files imported #
print(list.filenames.CASA)
#   

# repeat for MIMICS output #
list.filenames.MIMICS <- list.files(pattern=c("POINT_mimics"),full.names = F)
list.MIMICS <- list()
for(i in 1:length(list.filenames.MIMICS)){
  list.MIMICS[[i]] <- fread(list.filenames.MIMICS[i],drop=1,fill=TRUE,colClasses = list(numeric=1:141))}
names(list.MIMICS) <- list.filenames.MIMICS
##  delete the last column from these items                                        ##
for (i in 1:length(list.filenames.MIMICS)){
  list.MIMICS[[i]] <- list.MIMICS[[i]][,1:141]}
print(list.filenames.MIMICS)
#    


# mutate items in lists to calculate additional variables of interest:  #
list.CASA2 <- map(list.CASA, ~.x %>% mutate(AutoResp=`casaflux%Crmplant(LEAF)`+`casaflux%Crmplant(WOOD)`+`casaflux%Crmplant(FROOT)`+`casaflux%Crgplant`,
                                            LEAFCN=`casapool%Cplant(LEAF)`/`casapool%Nplant(LEAF)`, # leaf C:N ratio #
                                            WOODCN=`casapool%Cplant(WOOD)`/`casapool%Nplant(WOOD)`, # wood C:N ratio #
                                            FROOTCN=`casapool%Cplant(FROOT)`/`casapool%Nplant(FROOT)`, # fine root C:N ratio #
                                            LitMET.CN=`casapool%Clitter(MET)`/`casapool%Nlitter(MET)`, # metabolic litter C:N ratio #
                                            LitSTR.CN=`casapool%Clitter(STR)`/`casapool%Nlitter(STR)`, # structural litter C:N ratio #
                                            SoilMIC.CN=`casapool%Csoil(MIC)`/`casapool%Nsoil(MIC)`, # soil microbial pool C:N ratio #
                                            SoilSLOW.CN=`casapool%Csoil(SLOW)`/`casapool%Nsoil(SLOW)`, # slow soil pool C:N ratio #
                                            SoilPASS.CN=`casapool%Csoil(PASS)`/`casapool%Nsoil(PASS)`, # passive soil pool C:N ratio # 
                                            NetNMin=`casaflux%Nsmin`+`casaflux%Nsimm`+`casaflux%Nlittermin`, # net N mineralization flux #
                                            TotalSoilC= `casapool%Csoil(MIC)`+`casapool%Csoil(SLOW)`+`casapool%Csoil(PASS)`, # total mineral soil C pool #
                                            TotalSoilN= `casapool%Nsoil(MIC)`+`casapool%Nsoil(SLOW)`+`casapool%Nsoil(PASS)`, # total mineral soil N pool #
                                            TotalVegC= `casapool%Cplant(LEAF)`+`casapool%Cplant(WOOD)`+`casapool%Cplant(FROOT)`, # total vegetation C pool #
                                            TotalVegN= `casapool%Nplant(LEAF)`+`casapool%Nplant(WOOD)`+`casapool%Nplant(FROOT)`, # total vegetation N pool #
                                            NlitInptTot= `casaflux%NlitInptMet`+`casaflux%NlitInptStruc`, # total litter N flux to soil #
                                            Shoot2Root= (`casapool%Cplant(LEAF)`+`casapool%Cplant(WOOD)`)/`casapool%Cplant(FROOT)`, # shoot 2 root ratio #
                                            SoilC_minOH=`casapool%Csoil(MIC)`+`casapool%Csoil(SLOW)`+`casapool%Csoil(PASS)`+`casapool%Clitter(MET)`+`casapool%Clitter(STR)`, # total soil C pool (mineral and organic) #
                                            SoilN_minOH=`casapool%Nsoil(MIC)`+`casapool%Nsoil(SLOW)`+`casapool%Nsoil(PASS)`+`casapool%Nlitter(MET)`+`casapool%Nlitter(STR)`, # total soil N pool (mineral and organic) #
                                            ))

list.CASA2 <- map(list.CASA2, ~.x %>% mutate(BulkSoilCN=TotalSoilC/TotalSoilN, TotalVegCN=TotalVegC/TotalVegN, Soil_minOH_CN=SoilC_minOH/SoilN_minOH))
## where BulkSoilCN is mineral soil C:N ratio; TotalVegCN is total vegetation C:N ratio; Soil_minOH_CN is the total soil C:N (mineral and organic) ##

list.MIMICS2 <- map(list.MIMICS, ~.x %>% mutate(inptMet.CN= inptMetC/inptMetN, # metabolic litter inputs to soil C:N ratio #
                                                inptStr.CN= inputStrC/inputStrN, # structural litter inputs to soil C:N ratio #
                                                LitMET.CN=LITm/LITmN, # metabolic litter pool C:N ratio #
                                                LitSTR.CN=LITs/LITsN, # structural litter pool C:N ratio #
                                                SOMa.CN=SOMa/SOMaN, # SOMa pool C:N ratio #
                                                SOMc.CN=SOMc/SOMcN, # SOMc pool C:N ratio #
                                                SOMp.CN=SOMp/SOMpN, # SOMp pool C:N ratio #
                                                dLITs.CN=dLITs/dLITs, # C:N ratio of loss/gain in structural litter #
                                                dLITm.CN=dLITm/dLITmN, # C:N ratio of loss/gain in metabolic litter #
                                                MICr2MICk=MICr/MICk, # ratio of MICr biomass C to MICk biomass C #
                                                TotalSoilC=SOMa+SOMc+SOMp+MICr+MICk, # total mineral soil C pool (including microbial biomass C) #
                                                TotalSoilN=SOMaN+SOMcN+SOMpN+MICrN+MICkN, # total mineral soil N pool (including microbial biomass N) #
                                                SoilC_minOH=SOMa+SOMc+SOMp+LITm+LITs+MICr+MICk, # total soil C pool (mineral, organic and microbial biomass) #
                                                SoilN_minOH=SOMaN+SOMcN+SOMpN+LITmN+LITsN+MICrN+MICkN)) # total soil N pool (mineral, organic, and microbial biomass) #
list.MIMICS2 <- map(list.MIMICS2, ~.x %>% mutate(BulkSoilCN=TotalSoilC/TotalSoilN, Soil_minOH_CN=SoilC_minOH/SoilN_minOH))
## where BulkSoilCN is the mineral soil C:N ratio; Soil_minOH_CN is the total soil C:N ratio (mineral and organic) ##

#################### last ten years mean annual trends #######################
CASAten <- lapply(list.CASA2, '[', c(7665:11315),)            # subset last ten years
MIMICSten <- lapply(list.MIMICS2, '[',c(7665:11315),)         # subset last ten years

# average each day of year (DOY) over the last ten years for a mean annual cycle #
CASAtenmean <- map(CASAten, ~.x %>% group_by(idoy) %>% summarise_all(mean))
MIMISCtenmean <- map(MIMICSten, ~.x %>% group_by(doy) %>% summarise_all(mean))

######################## Annual Summaries ####################################
# calculate annual estimates of fluxes (sums) and ten year annual means for pools
CASAannflux <- map(list.CASA2, ~.x %>% select(contains("flux"),NlitInptTot,AutoResp,NetNMin,iYrCnt) %>% group_by(iYrCnt) %>% summarise_all(sum))
CASAmeanpool <- map(list.CASA2, ~.x %>% select(contains("pool"),contains("CN"),TotalSoilC,TotalSoilN,SoilC_minOH,SoilN_minOH,iYrCnt,xkNlimiting,Shoot2Root) %>% group_by(iYrCnt) %>% summarise_all(mean))

# create lists of variables from each model in an order that will line up with contrasting model (CASA vs MIMICS) #
## note: the order here is important! when selecting columns, they are set up to match columns in contrasting model! #
Veg_pools_2plot <- map(list.CASA2, ~.x %>% select(c(2,7:12,40:42,46,53,67:69,78,79,81,85)) %>% group_by(iYrCnt) %>% summarise_all(mean))
Veg_fluxes_2plot <- map(list.CASA2, ~.x %>% select(c(2,33:38,47,51,64,65,66,75,60,61,80)) %>% group_by(iYrCnt) %>% summarise_all(sum))
CSoil_pools_2plot <- map(list.CASA2, ~.x %>% select(c(2,19:21,22:24,13:17,70:74,76,77,82:84,86,46)) %>% group_by(iYrCnt) %>% summarise_all(mean))
Msoil_pools_2plot <- map(list.MIMICS2, ~.x %>% select(c(2,13:15,98:101,9:10,70,94,95,144:148,152:157,11:12,96:97,151,101)) %>% group_by(iYrCnt) %>% summarise_all(mean))
Csoil_fluxes_2plot <- map(list.CASA2, ~.x %>% select(c(2,25,26:32)) %>% group_by(iYrCnt) %>% summarise_all(sum))
Msoil_fluxes_2plot <- map(list.MIMICS2, ~.x %>% select(c(2,6,72:78)) %>% group_by(iYrCnt) %>% summarise_all(sum))

# convert MIMICS units (mg/cm3) to match CASA units (g/m2), in this case, with a soil depth of 45 cm, conversion factor = * 450 #
## NOTE: 45 cm is not default model soil depth, but is what we used here because we have observations to 45 cm depth to compare with ##
M2Cunits <- function(x){ (x*450)}
Msoil_pools_2plot <- map(Msoil_pools_2plot, ~.x %>% mutate_at(c(2:9,11:12,18:21,24:27,29),M2Cunits))
Msoil_fluxes_2plot <- map(Msoil_fluxes_2plot, ~.x %>% mutate_at(2,M2Cunits))


# create one dataframe with all listed dataframes:     #
allannualCASA <- bind_rows(CASAannflux, .id= "column_label")     # turn list into one big df #
allpoolsCASA <- bind_rows(CASAmeanpool, .id="column_label")      # turn list into one big df #

# repeat process of summarizing for MIMICS output files  #
MIMICSannflux <- map(list.MIMICS2, ~.x %>% select(hresp,iYrCnt,Nspill_r,Nspill_k,Overflow_r,Overflow_k,dLITm,
                                                  dLITs,dLITmN,dLITsN,dSOMa,dSOMc,dSOMp,dSOMaN,dSOMcN,dSOMpN,inptMetC,
                                                  inputStrC,inptMetN,inputStrN) %>% group_by(iYrCnt) %>% summarise_all(sum))
MIMICSmeanpool <- map(list.MIMICS2, ~.x %>% select(NPPan,LITm,LITs,CWD,LITmN,LITsN,SOMa,SOMp,SOMc,SOMaN,SOMpN,SOMcN,SOMa.CN,SOMc.CN,SOMp.CN,
                                                   DIN,MICr,MICk,MICrN,MICkN,MICr2MICk,ratLigN,iYrCnt,TotalSoilC,TotalSoilN,BulkSoilCN) %>% group_by(iYrCnt) %>% summarise_all(mean))
allpoolsMIMICS <- bind_rows(MIMICSmeanpool, .id="column_label")
allannualMIMICS <- bind_rows(MIMICSannflux, .id = "column_label")

####################################################################################################################
##    Import csv C and N balance outputs .csvs
list.filenames.CBALANCE <- list.files(pattern="Cbalance",full.names = F)
# read in all .csv balance files into one list  #
list.CBALANCE <- list()
y <- for (i in 1:length(list.filenames.CBALANCE)){
  list.CBALANCE[[i]] <- fread(list.filenames.CBALANCE[i],drop=1, fill=TRUE,colClasses = list(numeric = 1:10))}    ## for some reason, R is detecting an extra column ##
## drop and fill is best workaround to kep column names and all data, it adds a column (V55) of NAs at the end  ##
names(list.CBALANCE) <- list.filenames.CBALANCE
##  so delete the last column from the items with an extra column                                        ##
for (i in 1:length(list.CBALANCE)){
  list.CBALANCE[[i]] <- list.CBALANCE[[i]][,1:10]}
print(list.filenames.CBALANCE)
list.CBALANCE <- map(list.CBALANCE, ~.x %>% mutate(TotalC=Cplant+Clitter+Cmic+Csom+Clabile, Cbalance=Cin - Cout))

# summarize by year
#CBALANCE_annual <- map(list.CASA2, ~.x %>% select(contains("flux"),NlitInptTot,AutoResp,NetNMin,iYrCnt) %>% group_by(iYrCnt) %>% summarise_all(sum))
#CASAmeanpool <- map(list.CASA2, ~.x %>% select(contains("pool"),contains("CN"),TotalSoilC,TotalSoilN,SoilC_minOH,SoilN_minOH,iYrCnt,xkNlimiting,Shoot2Root) %>% group_by(iYrCnt) %>% summarise_all(mean))


# 

list.filenames.NBALANCE <- list.files(pattern="Nbalance",full.names = F)
# read in all .csv balance files into one list  #
list.NBALANCE <- list()
y <- for (i in 1:length(list.filenames.NBALANCE)){
  list.NBALANCE[[i]] <- fread(list.filenames.NBALANCE[i],drop=1, fill=TRUE,colClasses = list(numeric = 1:10))}    ## for some reason, R is detecting an extra column ##
## drop and fill is best workaround to kep column names and all data, it adds a column (V55) of NAs at the end  ##
names(list.NBALANCE) <- list.filenames.CBALANCE
##  so delete the last column from the items with an extra column                                        ##
for (i in 1:length(list.NBALANCE)){
  list.NBALANCE[[i]] <- list.NBALANCE[[i]][,1:10]}
print(list.filenames.NBALANCE)
list.NBALANCE <- map(list.NBALANCE, ~.x %>% mutate(TotalN=Nplant+Nlitter+Nmic+Nsom+Nsoilminrl, Nbalance=Nin-Nout))
# 


#################### last ten years mean annual trends (DOY means) #######################
CBALANCEten <- lapply(list.CBALANCE, '[', c(7665:11315),)            # subset last ten years
NBALANCEten <- lapply(list.NBALANCE, '[',c(7665:11315),)         # subset last ten years

CBALANCEtenmean <- map(CBALANCEten, ~.x %>% group_by(doy) %>% summarise_all(mean))
NBALANCEtenmean <- map(NBALANCEten, ~.x %>% group_by(doy) %>% summarise_all(mean))

#################### annual means (n=89) #################################################
TOTALCannual <- map(list.CBALANCE, ~.x  %>% group_by(iYrCnt) %>% summarise_all(mean))
CBALANCEannual <- map(list.CBALANCE, ~.x  %>% group_by(iYrCnt) %>% summarise_all(sum))

NBALANCEannual <- map(list.NBALANCE, ~.x  %>% group_by(iYrCnt) %>% summarise_all(mean))


###########################################################################
## see file: "plotting_Naddexperiment.R"                                 ##
##    for data visualization scripts                                     ##
###########################################################################
## NOTE: R must be restarted after importing files and before plotting   ##
#------------------------------------END----------------------------------#
###########################################################################






