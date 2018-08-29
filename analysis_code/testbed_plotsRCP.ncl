;-------------------------------------------------------------------
;
;   File:       testbed_maps.ncl
;
;   Author:     Melannie Hartman & Will Wieder
;
;   Date:       March 14, 2013
;               October 20, 2014
;               March 23, 2015
;               November 9, 2016
;
;   Description:    
;       Reads a netCDF file that contains gridded output
;       variables (pools and fluxes) from the end of a casaclm simulation.
;
;       THE USER MUST MODIFY filedir and file1 below!
;-------------------------------------------------------------------

load "$NCARG_LIB/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_LIB/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_LIB/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_LIB/ncarg/nclscripts/csm/shea_util.ncl"

;----------------------------------------------------------------------
; Function to attach a labelbar outside of an XY plot
;----------------------------------------------------------------------
function attach_labelbar(wks,plot,labels,colors)
local lbres, vph, vpw, nboxes
begin
  nboxes = dimsizes(labels) - 1
  getvalues plot                      ; Get plot size for use in
    "vpHeightF" : vph                 ; creating labelbar.
    "vpWidthF"  : vpw
  end getvalues

  lbres                    = True          ; labelbar only resources
  lbres@lbAutoManage       = False          ; Necessary to control sizes
  lbres@lbPerimOn          = False
  lbres@vpWidthF           = 0.35* vpw     ; labelbar width
  lbres@vpHeightF          = 0.8 * vph     ; labelbar height
  lbres@lbFillColors       = colors(::-1)  ; reverse color bar, again
  lbres@lbMonoFillPattern  = True          ; Solid fill pattern
  lbres@lbLabelFontHeightF = 0.022         ; font height. default is small
  lbres@lbLabelJust        = "CenterLeft"  ; left justify labels
  lbres@lbBoxLinesOn       = False
  lbres@lbLabelAlignment   ="ExternalEdges"
  lbres@lbTitleString      = "% Saturation"
  lbres@lbTitleFontHeightF = 0.025         ; default 0.025 ?
  lbres@lbLeftMarginF      = 0.10         ; default 0.05 
  lbres@lbTitleJust        = "BottomCenter"
  lbid = gsn_create_labelbar(wks,nboxes,labels(::-1),lbres) ; reverse labels
;
; Now, create some annotation resources indicating how we want to
; attach the labelbar to the plot.
;                                                                 
  amres                = True
  amres@amParallelPosF = 0.70     ; Move away from plot
  annoid = gsn_add_annotation(plot,lbid,amres)

  return(annoid)
end

;----------------------------------------------------------------------
; Main code
;----------------------------------------------------------------------
begin
 
;-------------------------------------------------------------------
;  Read a netCDF file using addfile.
;  USER MUST MODIFY filedir, file1!
;-------------------------------------------------------------------
 
 dir     = "/project/tss/wwieder/CASACLM/GRID/OUTPUT_GRID_CRU_NCEP/"
 exper   = "RCP45"
 
 mim_case = exper+"/"

 case0   = "CASACNP_"+exper+"/"
 case1   = "MIMICS_"+mim_case
 case2   = "CORPSE_"+exper+"_minFW/"

 if (exper .eq. "RCP45") then
   file0      = addfile(dir + case0 + "ann_casaclm_rcp45.nc","R")
   file1      = addfile(dir + case1 + "ann_mimics_rcp45.nc" ,"R")
   corpsedata = addfile(dir + case2 + "ann_corpse_rcp45.nc" ,"R")
 end if

 if (exper .eq. "RCP85") then
   file0      = addfile(dir + case0 + "ann_casaclm_rcp85.nc","R")
   file1      = addfile(dir + case1 + "ann_mimics_rcp85.nc" ,"R")
   corpsedata = addfile(dir + case2 + "ann_corpse_rcp85.nc" ,"R")
 end if

 out_dir    = dir + case1 
 print("out dir "+out_dir)
    
 obs_dir    = "/project/tss/wwieder/soilCN/global_run/input_data/"
 obs_data   = addfile(obs_dir + "Soil_obs.nc", "R")
 area_data  = addfile("/project/tss/wwieder/soilCN/clm_vars/clm_LITFALL_data.nc", "R")
 arc_data   = addfile(obs_dir + "NCSCD_Circumarctic_100cmSOCC_025deg.nc","R")
 print("files added")

 dims    = dimsizes(file1->cMICr)
 ntim    = dims(0)
 nlat    = dims(1)
 nlon    = dims(2)
;print(out_dir)
;print(dims)

 ;print(ntim)
 ;print(nlat)
 ;print(nlon)
 fyear = 1901
 lyear = 2100 
 time  = ispan(fyear,lyear,1)

;-------------------------------------------------------------------
;  Read in data      
;-------------------------------------------------------------------

; ------ Read in OBS variables -----------
HWSD_SOC_all = obs_data->HWSD_SOC
HWSD_SOC     = dim_sum_n_Wrap(HWSD_SOC_all, 0)
;printVarSummary(HWSD_SOC)

NCSCD         = arc_data->NCSCD_SOC
NCSCD         = NCSCD / 10    ; kg/m2 (orig units unclear)
NCSCD@units   = "kgC/m2"
;printVarSummary(NCSCD)

area_1deg    = area_data->areaupsc * area_data->landfrac ;km2

; ------ Read in MIMICS variables -----------
landarea = file1->landarea
cLITm    = file1->cLITm
cLITs    = file1->cLITs
cMICr    = file1->cMICr
cMICk    = file1->cMICk
cSOMa    = file1->cSOMa
cSOMc    = file1->cSOMc
cSOMp    = file1->cSOMp
mim_HR   = file1->cHresp
mim_fW   = file1->fW
mim_HR   = mim_HR * 365
mim_HR@units = "gC m-2 y-1"

mim_IN   = file1->cLitInput_metbc + file1->cLitInput_struc
mim_INmet= file1->cLitInput_metbc 
mim_INstr= file1->cLitInput_struc
mim_IN   = mim_IN * 365
mim_INmet= mim_INmet * 365
mim_INstr= mim_INstr * 365
lat      = mim_HR&lat

; -- sum pools --
; units = gC/m2
mim_LIT = cLITm + cLITs
mim_MIC = cMICr + cMICk
mim_SOM = cSOMa + cSOMp + cSOMc
mim_TOT = mim_LIT + mim_MIC + mim_SOM
mim_SOM_pro = cSOMp 

copy_VarMeta(cLITm , mim_LIT)
copy_VarMeta(cMICr , mim_MIC)
copy_VarMeta(cSOMa , mim_SOM)
copy_VarMeta(cLITm , mim_TOT)
copy_VarMeta(mim_HR, mim_IN)

mim_LIT@long_name = "Total Litter C, MIMICS"
mim_MIC@long_name = "Total Microbial C, MIMICS"
mim_SOM@long_name = "Total Soil C, MIMICS"
mim_TOT@long_name = "Total C, MIMICS (litter + microbe + SOM)"
mim_IN@long_name  = "MIMICS litter inputs"
mim_INmet@long_name  = "MIMICS metabolic  inputs"
mim_INstr@long_name  = "MIMICS structural inputs"

mim_fracMIC = where(mim_TOT.gt.0, mim_MIC / mim_TOT, mim_MIC@_FillValue)

; ------ Read in CASA variables -----------
; units, pools  = gC/m^2   & 
;        fluxes = gC/m^2/d

cas_VEG  = file0->cleaf    + file0->cwood   + file0->cfroot
cas_CWD  = file0->clitcwd
cas_LIT  = file0->clitmetb + file0->clitstr + file0->clitcwd
cas_SOM  = file0->csoilmic + file0->csoilslow + file0->csoilpass
cas_SOM_pro = file0->csoilpass 
; passive pool only for CASA
;cas_SOM_pro = cas_SOM_pro + file0->csoilslow
cas_TOT  = cas_LIT + cas_SOM
cas_HR   = file0->cresp
cas_GPP  = file0->cgpp
cas_NPP  = file0->cnpp
cas_IN   = file0->litInptStruc + file0->litInptMet
cas_INmet= file0->litInptMet
cas_INstr= file0->litInptStruc 
tair     = file0->tairC 
pft      = file0->IGBP_PFT
cas_HR   = cas_HR  * 365
cas_GPP  = cas_GPP * 365
cas_NPP  = cas_NPP * 365
cas_IN   = cas_IN  * 365
cas_INmet= cas_INmet* 365
cas_INstr= cas_INstr* 365
cas_AR   = cas_GPP - cas_NPP
cas_fW   = file0->fW
cas_fT   = file0->fT

cas_HR@units  = "gC m-2 y-1"
cas_GPP@units = "gC m-2 y-1"
cas_NPP@units = "gC m-2 y-1"
cas_IN@units  = "gC m-2 y-1"
cas_INstr@units  = "gC m-2 y-1" 
cas_INmet@units  = "gC m-2 y-1"

copy_VarMeta(cas_HR, cas_AR)
copy_VarMeta( mim_LIT, cas_VEG)
copy_VarMeta( mim_LIT, cas_CWD)
copy_VarMeta( mim_LIT, cas_LIT)
copy_VarMeta( mim_SOM, cas_SOM)
copy_VarMeta( mim_TOT, cas_TOT)

cas_AR@long_name  = "Total AR, CASA GPP-NPP"
cas_VEG@long_name = "Total VEG C, CASA"
cas_CWD@long_name = "Total CWD C, CASA"
cas_LIT@long_name = "Total Litter C, CASA"
cas_SOM@long_name = "Total Soil C, CASA"
cas_TOT@long_name = "Total C, CASA (litter + SOM)"
cas_HR@long_name  = "Total HR, CASA"

; ------ Read in CORPSE variables -----------
LitterLayer_C1 = corpsedata->LitterLayer_C1
LitterLayer_C2 = corpsedata->LitterLayer_C2
LitterLayer_C3 = corpsedata->LitterLayer_C3

litLiveMicrobeC  = corpsedata->LitterLayer_LiveMicrobeC
soilLiveMicrobeC = corpsedata->Soil_LiveMicrobeC        

Soil_C1 = corpsedata->Soil_C1
Soil_C2 = corpsedata->Soil_C2
Soil_C3 = corpsedata->Soil_C3

SoilProtected_C1 = corpsedata->SoilProtected_C1
SoilProtected_C2 = corpsedata->SoilProtected_C2
SoilProtected_C3 = corpsedata->SoilProtected_C3

LitterLayer_CO2 = corpsedata->LitterLayer_CO2
Soil_CO2        = corpsedata->Soil_CO2
TSOI            = corpsedata->Ts
TSOI            = TSOI - 273.15
; read in mimics thetaLiq calculation, for now
Theta           = file1->thetaLiq
Theta           = Theta * 100 ; convert from fraction to % saturation
;printVarSummary(Soil_CO2)

; -- sum CORPSE pools --
cor_LIT  = LitterLayer_C1  + LitterLayer_C2 + LitterLayer_C3 
cor_MIC  = litLiveMicrobeC + soilLiveMicrobeC
cor_SOM  = Soil_C1 + Soil_C2 + Soil_C3 + SoilProtected_C1 + SoilProtected_C2 + SoilProtected_C3
cor_TOT  = cor_LIT + cor_MIC + cor_SOM
cor_HR   = (LitterLayer_CO2 + Soil_CO2) *365
cor_fW   = corpsedata->fW
cor_SOM_pro  = SoilProtected_C1 + SoilProtected_C2 + SoilProtected_C3 
cor_TOT2 = where(cor_TOT.gt.0, cor_TOT, cor_TOT@_FillValue)
;printMinMax(cor_TOT ,1)
;printMinMax(cor_TOT2,1)

cor_fracMIC = cor_MIC/cor_TOT2

min_val  = 20    ; mask out low productivity cells
cas_TOT = where(cas_NPP.le.min_val, cas_TOT@_FillValue, cas_TOT)
mim_TOT = where(cas_NPP.le.min_val, mim_TOT@_FillValue, mim_TOT)
cor_TOT = where(cas_NPP.le.min_val, cor_TOT@_FillValue, cor_TOT)
landarea = where(cas_NPP(0,:,:).le.min_val, landarea@_FillValue, landarea)

cas_hr = where(cas_HR.le.min_val, cas_HR@_FillValue, cas_HR)
mim_hr = where(mim_HR.le.min_val, mim_HR@_FillValue, mim_HR)
cor_hr = where(cor_HR.le.min_val, cor_HR@_FillValue, cor_HR)
cas_TAU = where(cas_hr.gt.min_val, cas_TOT/cas_hr, cas_hr@_FillValue)
mim_TAU = where(mim_hr.gt.min_val, mim_TOT/mim_hr, mim_hr@_FillValue)
cor_TAU = where(cor_hr.gt.min_val, cor_TOT/cor_hr, cor_hr@_FillValue)
cor_TAU = where(cor_TAU.le.1e-3   , cor_TAU@_FillValue, cor_TAU)
;printMinMax(cor_TAU,1)
;printVarSummary(cor_HR)
;printVarSummary(cas_HR)

copy_VarMeta(LitterLayer_C1  , cor_LIT )
copy_VarMeta(litLiveMicrobeC , cor_MIC )
copy_VarMeta(SoilProtected_C1, cor_SOM )
copy_VarMeta(SoilProtected_C1, cor_TOT )
copy_VarMeta(SoilProtected_C1, cor_SOM_pro )
copy_VarMeta(LitterLayer_CO2 , cor_HR  )


cor_LIT@long_name = "Total Litter C, CORPSE"
cor_MIC@long_name = "Total Microbial C, CORPSE"
cor_SOM@long_name = "Total Soil C, CORPSE"
cor_TOT@long_name = "Total C, CORPSE (litter + microbe + SOM)"
cor_SOM_pro@long_name = "Total protected SOM, CORPSE"
cor_HR@long_name  = "Total HR, CORPSE"
cor_HR@units      = "gC m-2 y-1"


; ------------------------------------------------------------
;          calculate global stockstestbed_plots3.ncl
; ------------------------------------------------------------
landareaC = conform_dims(dims, landarea, (/1,2/))
landareaWGT = landareaC / sum(landarea)
copy_VarCoords(cor_SOM, landareaC)
copy_VarCoords(cor_SOM, landareaWGT)

; -- HWSD_OBS --
gHWSD_SOC = sum(HWSD_SOC * area_1deg * 1e6 ) / 1e12 ; convert kgC/m2 to PgC globally

; -- CASA --   
Gcas_VEG  = dim_sum_n(cas_VEG * landareaC * 1e6, (/1,2/) ) / 1e15  ;convert gC/m2 to Pg globally
Gcas_CWD  = dim_sum_n(cas_CWD * landareaC * 1e6, (/1,2/) ) / 1e15  
Gcas_LIT  = dim_sum_n(cas_LIT * landareaC * 1e6, (/1,2/) ) / 1e15  
Gcas_SOM  = dim_sum_n(cas_SOM * landareaC * 1e6, (/1,2/) ) / 1e15
Gcas_TOT  = dim_sum_n(cas_TOT * landareaC * 1e6, (/1,2/) ) / 1e15
Gcas_HR   = dim_sum_n(cas_HR  * landareaC * 1e6, (/1,2/) ) / 1e15
Gcas_AR   = dim_sum_n(cas_AR  * landareaC * 1e6, (/1,2/) ) / 1e15
Gcas_IN   = dim_sum_n(cas_IN  * landareaC * 1e6, (/1,2/) ) / 1e15
Gcas_INstr= dim_sum_n(cas_INstr*landareaC * 1e6, (/1,2/) ) / 1e15
Gcas_INmet= dim_sum_n(cas_INmet*landareaC * 1e6, (/1,2/) ) / 1e15
Gcas_GPP  = dim_sum_n(cas_GPP * landareaC * 1e6, (/1,2/) ) / 1e15
Gcas_NPP  = dim_sum_n(cas_NPP * landareaC * 1e6, (/1,2/) ) / 1e15
Gcas_TAU  = dim_sum_n(cas_TAU * landareaWGT, (/1,2/))
Gtair     = dim_sum_n(tair    * landareaC, (/1,2/) ) / sum(landarea)  ;weighted annual air tem
GTSOI     = dim_sum_n(TSOI    * landareaC, (/1,2/) ) / sum(landarea)
GTheta    = dim_sum_n(Theta   * landareaC, (/1,2/) ) / sum(landarea)
Gcas_fT   = dim_sum_n(cas_fT  * landareaC, (/1,2/) ) / sum(landarea)
Gcas_fW   = dim_sum_n(cas_fW  * landareaC, (/1,2/) ) / sum(landarea)
Gmim_fW   = dim_sum_n(mim_fW  * landareaC, (/1,2/) ) / sum(landarea)
Gcor_fW   = dim_sum_n(cor_fW  * landareaC, (/1,2/) ) / sum(landarea)

Acas_NPP  = dim_sum_n(cas_NPP(:,{50:90},:) * landareaC(:,{50:90},:) * 1e6, (/1,2/) ) / 1e15
Acas_GPP  = dim_sum_n(cas_GPP(:,{50:90},:) * landareaC(:,{50:90},:) * 1e6, (/1,2/) ) / 1e15
Acas_VEG  = dim_sum_n(cas_VEG(:,{50:90},:) * landareaC(:,{50:90},:) * 1e6, (/1,2/) ) / 1e15  
Acas_TOT  = dim_sum_n(cas_TOT(:,{50:90},:) * landareaC(:,{50:90},:) * 1e6, (/1,2/) ) / 1e15
Amim_TOT  = dim_sum_n(mim_TOT(:,{50:90},:) * landareaC(:,{50:90},:) * 1e6, (/1,2/) ) / 1e15
Acor_TOT  = dim_sum_n(cor_TOT(:,{50:90},:) * landareaC(:,{50:90},:) * 1e6, (/1,2/) ) / 1e15
ATSOI     = dim_sum_n(TSOI(:,{50:90},:)    * landareaC(:,{50:90},:), (/1,2/) ) / sum(landarea({50:90},:))
ATheta    = dim_sum_n(Theta(:,{50:90},:)   * landareaC(:,{50:90},:), (/1,2/) ) / sum(landarea({50:90},:))
Acas_fT    = dim_sum_n(cas_fT(:,{50:90},:)   * landareaC(:,{50:90},:), (/1,2/) ) / sum(landarea({50:90},:))
Acas_fW    = dim_sum_n(cas_fW(:,{50:90},:)   * landareaC(:,{50:90},:), (/1,2/) ) / sum(landarea({50:90},:))
Amim_fW    = dim_sum_n(mim_fW(:,{50:90},:)   * landareaC(:,{50:90},:), (/1,2/) ) / sum(landarea({50:90},:))
Acor_fW    = dim_sum_n(cor_fW(:,{50:90},:)   * landareaC(:,{50:90},:), (/1,2/) ) / sum(landarea({50:90},:))


Tcas_NPP  = dim_sum_n(cas_NPP(:,{-90:50},:) * landareaC(:,{-90:50},:) * 1e6, (/1,2/) ) / 1e15
Tcas_GPP  = dim_sum_n(cas_GPP(:,{-90:50},:) * landareaC(:,{-90:50},:) * 1e6, (/1,2/) ) / 1e15
Tcas_VEG  = dim_sum_n(cas_VEG(:,{-90:50},:) * landareaC(:,{-90:50},:) * 1e6, (/1,2/) ) / 1e15
Tcas_TOT  = dim_sum_n(cas_TOT(:,{-90:50},:) * landareaC(:,{-90:50},:) * 1e6, (/1,2/) ) / 1e15
Tmim_TOT  = dim_sum_n(mim_TOT(:,{-90:50},:) * landareaC(:,{-90:50},:) * 1e6, (/1,2/) ) / 1e15
Tcor_TOT  = dim_sum_n(cor_TOT(:,{-90:50},:) * landareaC(:,{-90:50},:) * 1e6, (/1,2/) ) / 1e15
TTSOI     = dim_sum_n(TSOI(:,{-90:50},:)    * landareaC(:,{-90:50},:), (/1,2/) ) / sum(landarea({-90:50},:))
TTheta    = dim_sum_n(Theta(:,{-90:50},:)   * landareaC(:,{-90:50},:), (/1,2/) ) / sum(landarea({-90:50},:))


; -- MIMICS --
Gmim_LIT  = dim_sum_n(mim_LIT * landareaC * 1e6, (/1,2/) ) / 1e15  ;convert gC/m2 to Pg globally
Gmim_MIC  = dim_sum_n(mim_MIC * landareaC * 1e6, (/1,2/) ) / 1e15
Gmim_SOM  = dim_sum_n(mim_SOM * landareaC * 1e6, (/1,2/) ) / 1e15
Gmim_TOT  = dim_sum_n(mim_TOT * landareaC * 1e6, (/1,2/) ) / 1e15
Gmim_HR   = dim_sum_n(mim_HR  * landareaC * 1e6, (/1,2/) ) / 1e15
Gmim_IN   = dim_sum_n(mim_IN  * landareaC * 1e6, (/1,2/) ) / 1e15
Gmim_INmet= dim_sum_n(mim_INmet*landareaC * 1e6, (/1,2/) ) / 1e15
Gmim_INstr= dim_sum_n(mim_INstr*landareaC * 1e6, (/1,2/) ) / 1e15
Gmim_TAU = dim_sum_n(mim_TAU * landareaWGT, (/1,2/))
Gmim_fracMIC = dim_sum_n(mim_fracMIC * landareaC / sum(landarea), (/1,2/))
Gmim_fracMIC = dim_avg_n(mim_fracMIC, (/1,2/))
Gmim_fracMICsd = dim_stddev_n(mim_fracMIC, (/1,2/))

 print("MIM fracMIC " +Gmim_fracMIC(1)+", "+ Gmim_fracMICsd(1))

; -- CORPSE -- 
Gcor_LIT  = dim_sum_n(cor_LIT * landareaC * 1e6, (/1,2/) ) / 1e15  ;convert gC/m2 to Pg globally
Gcor_MIC  = dim_sum_n(cor_MIC * landareaC * 1e6, (/1,2/) ) / 1e15
Gcor_SOM  = dim_sum_n(cor_SOM * landareaC * 1e6, (/1,2/) ) / 1e15
Gcor_TOT  = dim_sum_n(cor_TOT * landareaC * 1e6, (/1,2/) ) / 1e15

printMinMax(cor_HR, 0)
cor_HR = where(cor_HR .gt. 1e5, cor_HR@_FillValue, cor_HR)
printMinMax(cor_HR, 0)

Gcor_HR   = dim_sum_n(cor_HR  * landareaC * 1e6, (/1,2/) ) / 1e15
Gcor_TAU  = dim_sum_n(cor_TAU * landareaWGT , (/1,2/))
Gcor_fracMIC = dim_sum_n(cor_fracMIC * landareaC / sum(landarea), (/1,2/))
Gcor_fracMIC = dim_avg_n(cor_fracMIC, (/1,2/))
Gcor_fracMICsd = dim_stddev_n(cor_fracMIC, (/1,2/))

print("COR fracMIC " + Gcor_fracMIC(1)+", "+ Gcor_fracMICsd(1))

Gcas_SOM_pro  = dim_sum_n(cas_SOM_pro * landareaC * 1e6, (/1,2/) ) / 1e15
Gmim_SOM_pro  = dim_sum_n(mim_SOM_pro * landareaC * 1e6, (/1,2/) ) / 1e15
Gcor_SOM_pro  = dim_sum_n(cor_SOM_pro * landareaC * 1e6, (/1,2/) ) / 1e15
print("Gcas_SOM_pro init = "+Gcas_SOM_pro(0)+" final = "+Gcas_SOM_pro(199))
print("Gmim_SOM_pro init = "+Gmim_SOM_pro(0)+" final = "+Gmim_SOM_pro(199))
print("Gcor_SOM_pro init = "+Gcor_SOM_pro(0)+" final = "+Gcor_SOM_pro(199))

print("")
print("mimics Inputs = "+avg(Gmim_IN( 0:9))   + \
               ", MET= "+avg(Gmim_INmet( 0:9))+ \
               ", STR= "+avg(Gmim_INstr( 0:9)))
print("casa   Inputs = "+avg(Gcas_IN( 0:9))   + \
               ", MET= "+avg(Gcas_INmet( 0:9))+ \
               ", STR= "+avg(Gcas_INstr( 0:9)))

print("")
print("mean initial NPP = "+avg(Gcas_NPP(0:9))+ \
                 ", GPP = "+avg(Gcas_GPP(0:9))+ \
                 ", VEG = "+avg(Gcas_VEG(0:9))+ \
                 ", CWD = "+avg(Gcas_CWD(0:9))+ \
               ", mimIN = "+avg(Gmim_IN( 0:9))+ \ 
               ", casIN = "+avg(Gcas_IN( 0:9)))
print("stddev       NPP = "+stddev(Gcas_NPP(0:9))+ \
                 ", GPP = "+stddev(Gcas_GPP(0:9))+ \
                 ", VEG = "+stddev(Gcas_VEG(0:9))+ \
                 ", CWD = "+stddev(Gcas_CWD(0:9))+ \
               ", mimIN = "+stddev(Gmim_IN( 0:9))+ \
               ", casIN = "+stddev(Gcas_IN( 0:9)))
print("")
print("mean 2010  NPP = "+avg(Gcas_NPP(100:109))+ \
               ", GPP = "+avg(Gcas_GPP(100:109))+ \
               ", VEG = "+avg(Gcas_VEG(100:109))+ \
               ", CWD = "+avg(Gcas_CWD(100:109))+ \
             ", mimIN = "+avg(Gmim_IN( 100:109))+ \
             ", casIN = "+avg(Gcas_IN( 100:109)))
print("")
print("mean final NPP = "+avg(Gcas_NPP(ntim-10:ntim-1))+ \
               ", GPP = "+avg(Gcas_GPP(ntim-10:ntim-1))+ \
               ", VEG = "+avg(Gcas_VEG(ntim-10:ntim-1))+ \
               ", CWD = "+avg(Gcas_CWD(ntim-10:ntim-1))+ \
             ", mimIN = "+avg(Gmim_IN( ntim-10:ntim-1))+ \
             ", casIN = "+avg(Gcas_IN( ntim-10:ntim-1)))
print("stddev     NPP = "+stddev(Gcas_NPP(ntim-10:ntim-1))+ \
               ", GPP = "+stddev(Gcas_GPP(ntim-10:ntim-1))+ \
               ", VEG = "+stddev(Gcas_VEG(ntim-10:ntim-1))+ \
               ", CWD = "+stddev(Gcas_CWD(ntim-10:ntim-1))+ \
             ", mimIN = "+stddev(Gmim_IN( ntim-10:ntim-1))+ \
             ", casIN = "+stddev(Gcas_IN( ntim-10:ntim-1)))
print("")
print("mean initial TSOI   = "+avg(GTSOI(0:9))  +  \
                 ", Theta  = "+avg(GTheta(0:9)) +  \
                 ", cas_fT = "+avg(Gcas_fT(0:9))+  \
                 ", cas_fW = "+avg(Gcas_fW(0:9))+  \
                 ", mim_fW = "+avg(Gmim_fW(0:9))+  \
                 ", cas_fW = "+avg(Gcas_fW(0:9)) )
print("mean 2010    TSOI   = "+avg(GTSOI(100:109))  +  \
			 ", Theta  = "+avg(GTheta(100:109)) +  \
			 ", cas_fT = "+avg(Gcas_fT(100:109))+  \
                 ", cas_fW = "+avg(Gcas_fW(100:109))+  \
                 ", mim_fW = "+avg(Gmim_fW(100:109))+  \
                 ", cas_fW = "+avg(Gcas_fW(100:109)) )
print("mean final   TSOI   = "+avg(GTSOI(ntim-10:ntim-1))  +  \
                 ", Theta  = "+avg(GTheta(ntim-10:ntim-1)) +  \
                 ", cas_fT = "+avg(Gcas_fT(ntim-10:ntim-1))+  \
                 ", cas_fW = "+avg(Gcas_fW(ntim-10:ntim-1))+  \
                 ", mim_fW = "+avg(Gmim_fW(ntim-10:ntim-1))+  \
                 ", cas_fW = "+avg(Gcas_fW(ntim-10:ntim-1)) )
print("")

print("                 CASA   MIMICS  CORPSE     ")
print("initial tot  C = "+decimalPlaces(avg(Gcas_TOT(  0:9   )),1,True) + \
                     ", "+decimalPlaces(avg(Gmim_TOT(  0:9   )),1,True) + \
                     ", "+decimalPlaces(avg(Gcor_TOT(  0:9   )),1,True) )
print("2010   tot  C = "+decimalPlaces(avg(Gcas_TOT(100:109)),1,True) + \
                     ", "+decimalPlaces(avg(Gmim_TOT(100:109)),1,True) + \
                     ", "+decimalPlaces(avg(Gcor_TOT(100:109)),1,True) )
print("final   tot  C = "+decimalPlaces(avg(Gcas_TOT(ntim-10:ntim-1)),1,True) + \
                     ", "+decimalPlaces(avg(Gmim_TOT(ntim-10:ntim-1)),1,True) + \
                     ", "+decimalPlaces(avg(Gcor_TOT(ntim-10:ntim-1)),1,True) )

print("")
print("initial pro  C = "+decimalPlaces(avg(Gcas_SOM_pro(  0:9   )),1,True) + \
                     ", "+decimalPlaces(avg(Gmim_SOM_pro(  0:9   )),1,True) + \
                     ", "+decimalPlaces(avg(Gcor_SOM_pro(  0:9   )),1,True) )
print("final ProSoilC = "+decimalPlaces(avg(Gcas_SOM_pro(ntim-10:ntim-1)),1,True) + \
                     ", "+decimalPlaces(avg(Gmim_SOM_pro(ntim-10:ntim-1)),1,True) + \
                     ", "+decimalPlaces(avg(Gcor_SOM_pro(ntim-10:ntim-1)),1,True) )

print("")
print("initial RH       "+decimalPlaces(avg(Gcas_HR(  0:9   )),1,True) + \ 
                     ", "+decimalPlaces(avg(Gmim_HR(  0:9   )),1,True) + \
                     ", "+decimalPlaces(avg(Gcor_HR(  0:9   )),1,True) )
print("final   RH       "+decimalPlaces(avg(Gcas_HR(ntim-10:ntim-1)),1,True) + \ 
                     ", "+decimalPlaces(avg(Gmim_HR(ntim-10:ntim-1)),1,True) + \
                     ", "+decimalPlaces(avg(Gcor_HR(ntim-10:ntim-1)),1,True) )

print("")
print("initial turnover "+decimalPlaces(avg(Gcas_TAU(  0:9   )),1,True) + \
                     ", "+decimalPlaces(avg(Gmim_TAU(  0:9   )),1,True) + \
                     ", "+decimalPlaces(avg(Gcor_TAU(  0:9   )),1,True) )
print("final   turnover "+decimalPlaces(avg(Gcas_TAU(ntim-10:ntim-1)),1,True) + \ 
                     ", "+decimalPlaces(avg(Gmim_TAU(ntim-10:ntim-1)),1,True) + \
                     ", "+decimalPlaces(avg(Gcor_TAU(ntim-10:ntim-1)),1,True) )

print("initial MIC C = 0.00"  + \
                     ", "+decimalPlaces(avg(Gmim_MIC(  0:9   )),2,True) + \
                     ", "+decimalPlaces(avg(Gcor_MIC(  0:9   )),2,True) )
print("final   MIC C = 0.00 " + \
                     ", "+decimalPlaces(avg(Gmim_MIC(ntim-10:ntim-1)),2,True) + \
                     ", "+decimalPlaces(avg(Gcor_MIC(ntim-10:ntim-1)),2,True) )

print("")


;---------------------------
; Arctic only results
;---------------------------
print("")
print("--------Arctic only (> 50 N) -----------")
print("mean initial NPP = "+avg(Acas_NPP(0:9))+ \
                 ", GPP = "+avg(Acas_GPP(0:9))+ \
                 ", VEG = "+avg(Acas_VEG(0:9))+ \
                 ", TSOI = "+avg(ATSOI(0:9))  + \
                 ", Theta = "+avg(ATheta(0:9)) )
print("mean initial, cas_fT = "+avg(Acas_fT(0:9))+  \
                  ", cas_fW = "+avg(Acas_fW(0:9))+  \
                  ", mim_fW = "+avg(Amim_fW(0:9))+  \
                  ", cas_fW = "+avg(Acas_fW(0:9)) )

print("")
print("mean 2010    NPP = "+avg(Acas_NPP(100:109))+ \
                 ", GPP = "+avg(Acas_GPP(100:109))+ \
                 ", VEG = "+avg(Acas_VEG(100:109))+ \
                 ", TSOI = "+avg(ATSOI(100:109))  + \
                 ", Theta = "+avg(ATheta(100:109)) )
print("mean 2010 , cas_fT = "+avg(Acas_fT(100:109))+  \
                ", cas_fW = "+avg(Acas_fW(100:109))+  \
                ", mim_fW = "+avg(Amim_fW(100:109))+  \
                ", cas_fW = "+avg(Acas_fW(100:109)) )
print("")
print("mean final   NPP = "+avg(Acas_NPP(ntim-10:ntim-1))+ \
                 ", GPP = "+avg(Acas_GPP(ntim-10:ntim-1))+ \
                 ", VEG = "+avg(Acas_VEG(ntim-10:ntim-1))+ \ 
                 ", TSOI = "+avg(ATSOI(ntim-10:ntim-1))  + \
                 ", Theta = "+avg(ATheta(ntim-10:ntim-1)) )
print("mean final, cas_fT = "+avg(Acas_fT(ntim-10:ntim-1))+  \
                ", cas_fW = "+avg(Acas_fW(ntim-10:ntim-1))+  \
                ", mim_fW = "+avg(Amim_fW(ntim-10:ntim-1))+  \
                ", cas_fW = "+avg(Acas_fW(ntim-10:ntim-1)) )

print("")
print("                 CASA   MIMICS  CORPSE     ")
print("initial Soil C = "+decimalPlaces(avg(Acas_TOT(  0:9   )),1,True) + \
                     ", "+decimalPlaces(avg(Amim_TOT(  0:9   )),1,True) + \
                     ", "+decimalPlaces(avg(Acor_TOT(  0:9   )),1,True) )
print("2010    Soil C = "+decimalPlaces(avg(Acas_TOT(100:109)),1,True) + \
                     ", "+decimalPlaces(avg(Amim_TOT(100:109)),1,True) + \
                     ", "+decimalPlaces(avg(Acor_TOT(100:109)),1,True) )
print("final   Soil C = "+decimalPlaces(avg(Acas_TOT(ntim-10:ntim-1)),1,True) + \
                     ", "+decimalPlaces(avg(Amim_TOT(ntim-10:ntim-1)),1,True) + \
                     ", "+decimalPlaces(avg(Acor_TOT(ntim-10:ntim-1)),1,True) )
;---------------------------
; Temperate & tropical only results
;---------------------------
print("")
print("--------Temperate & tropical (< 50 N) -----------")
print("mean initial NPP = "+avg(Tcas_NPP(0:9))+ \
                 ", GPP = "+avg(Tcas_GPP(0:9))+ \
                 ", VEG = "+avg(Tcas_VEG(0:9))+ \
                 ", TSOI = "+avg(TTSOI(0:9))  + \
                 ", Theta = "+avg(TTheta(0:9)) )
print("")
print("mean final   NPP = "+avg(Tcas_NPP(ntim-10:ntim-1))+ \
                 ", GPP = "+avg(Tcas_GPP(ntim-10:ntim-1))+ \
                 ", VEG = "+avg(Tcas_VEG(ntim-10:ntim-1))+ \
                 ", TSOI = "+avg(TTSOI(ntim-10:ntim-1))  + \
                 ", Theta = "+avg(TTheta(ntim-10:ntim-1)) )

print("")
print("                 CASA   MIMICS  CORPSE     ")
print("initial Soil C = "+decimalPlaces(avg(Tcas_TOT(  0:9   )),1,True) + \
                     ", "+decimalPlaces(avg(Tmim_TOT(  0:9   )),1,True) + \
                     ", "+decimalPlaces(avg(Tcor_TOT(  0:9   )),1,True) )
print("final   Soil C = "+decimalPlaces(avg(Tcas_TOT(ntim-10:ntim-1)),1,True) + \
                     ", "+decimalPlaces(avg(Tmim_TOT(ntim-10:ntim-1)),1,True) + \
                     ", "+decimalPlaces(avg(Tcor_TOT(ntim-10:ntim-1)),1,True) )


; --- calculate change in stocks ---------
; global stocks - average of first 10 year
; ----------------------------------------

d_cas_GPP = new(dims,"float")  ; define arrays
copy_VarCoords(cas_GPP,d_cas_GPP)
d_cas_NPP = d_cas_GPP         
d_cas_VEG = d_cas_GPP         
d_cas_CWD = d_cas_GPP         
d_cas_LIT = d_cas_GPP          
d_cas_SOM = d_cas_GPP         
d_cas_SOM_pro = d_cas_GPP         
d_cas_TOT = d_cas_GPP         
d_cas_HR  = d_cas_GPP         
d_cas_AR  = d_cas_GPP         
d_TAIR    = d_cas_GPP         
d_TSOI    = new(dims,"double")
copy_VarCoords(TSOI, d_TSOI)
d_Theta   = d_TSOI            

d_mim_LIT = d_cas_GPP          ; define arrays
d_mim_MIC = d_cas_GPP         
d_mim_SOM = d_cas_GPP         
d_mim_SOM_pro = d_cas_GPP         
d_mim_TOT = d_cas_GPP         
d_mim_HR  = d_cas_GPP         
d_mim_IN  = d_cas_GPP         

d_cor_LIT = d_TSOI             ; define arrays
d_cor_MIC = d_TSOI              
d_cor_SOM = d_TSOI             
d_cor_SOM_pro = d_TSOI             
d_cor_TOT = d_TSOI             
d_cor_HR  = d_TSOI             

do i = 0,ntim - 1               ; 

   d_cas_GPP(i,:,:) = cas_GPP(i,:,:) - dim_avg_n(cas_GPP(0:9,:,:),0)
   d_cas_NPP(i,:,:) = cas_NPP(i,:,:) - dim_avg_n(cas_NPP(0:9,:,:),0)
   d_cas_VEG(i,:,:) = cas_VEG(i,:,:) - dim_avg_n(cas_VEG(0:9,:,:),0)
   d_cas_CWD(i,:,:) = cas_CWD(i,:,:) - dim_avg_n(cas_CWD(0:9,:,:),0)
   d_cas_LIT(i,:,:) = cas_LIT(i,:,:) - dim_avg_n(cas_LIT(0:9,:,:),0)
   d_cas_SOM(i,:,:) = cas_SOM(i,:,:) - dim_avg_n(cas_SOM(0:9,:,:),0)
   d_cas_TOT(i,:,:) = cas_TOT(i,:,:) - dim_avg_n(cas_TOT(0:9,:,:),0)
   d_cas_HR(i,:,:)  = cas_HR (i,:,:) - dim_avg_n(cas_HR (0:9,:,:),0)
   d_cas_AR(i,:,:)  = cas_AR (i,:,:) - dim_avg_n(cas_AR (0:9,:,:),0)
   d_TAIR(i,:,:)    = tair(i,:,:)    - dim_avg_n(tair   (0:9,:,:),0)
   d_TSOI(i,:,:)    = TSOI(i,:,:)    - dim_avg_n(TSOI   (0:9,:,:),0)
   d_Theta(i,:,:)   = Theta(i,:,:)   - dim_avg_n(Theta  (0:9,:,:),0)

   d_mim_LIT(i,:,:) = mim_LIT(i,:,:) - dim_avg_n(mim_LIT(0:9,:,:),0) 
   d_mim_MIC(i,:,:) = mim_MIC(i,:,:) - dim_avg_n(mim_MIC(0:9,:,:),0)     
   d_mim_SOM(i,:,:) = mim_SOM(i,:,:) - dim_avg_n(mim_SOM(0:9,:,:),0)     
   d_mim_TOT(i,:,:) = mim_TOT(i,:,:) - dim_avg_n(mim_TOT(0:9,:,:),0)     
   d_mim_HR(i,:,:)  = mim_HR (i,:,:) - dim_avg_n(mim_HR (0:9,:,:),0)     
   d_mim_IN(i,:,:)  = mim_IN (i,:,:) - dim_avg_n(mim_IN (0:9,:,:),0)     

   d_cor_LIT(i,:,:) = cor_LIT(i,:,:) - dim_avg_n(cor_LIT(0:9,:,:),0)
   d_cor_MIC(i,:,:) = cor_MIC(i,:,:) - dim_avg_n(cor_MIC(0:9,:,:),0)
   d_cor_SOM(i,:,:) = cor_SOM(i,:,:) - dim_avg_n(cor_SOM(0:9,:,:),0)
   d_cor_TOT(i,:,:) = cor_TOT(i,:,:) - dim_avg_n(cor_TOT(0:9,:,:),0)
   d_cor_HR(i,:,:)  = cor_HR (i,:,:) - dim_avg_n(cor_HR (0:9,:,:),0)

   d_cas_SOM_pro(i,:,:) = cas_SOM_pro(i,:,:) - dim_avg_n(cas_SOM_pro(0:9,:,:),0)
   d_mim_SOM_pro(i,:,:) = mim_SOM_pro(i,:,:) - dim_avg_n(mim_SOM_pro(0:9,:,:),0) 
   d_cor_SOM_pro(i,:,:) = cor_SOM_pro(i,:,:) - dim_avg_n(cor_SOM_pro(0:9,:,:),0)

end do

;-------------------------------------------------------------------
; -- global sum of changes --
;-------------------------------------------------------------------
d_Gcas_GPP = dim_sum_n(d_cas_GPP* landareaC * 1e6, (/1,2/))/1e15
d_Gcas_NPP = dim_sum_n(d_cas_NPP* landareaC * 1e6, (/1,2/))/1e15
d_Gcas_VEG = dim_sum_n(d_cas_VEG* landareaC * 1e6, (/1,2/))/1e15
d_Gcas_CWD = dim_sum_n(d_cas_CWD* landareaC * 1e6, (/1,2/))/1e15
d_Gcas_LIT = dim_sum_n(d_cas_LIT* landareaC * 1e6, (/1,2/))/1e15
d_Gcas_SOM = dim_sum_n(d_cas_SOM* landareaC * 1e6, (/1,2/))/1e15
d_Gcas_TOT = dim_sum_n(d_cas_TOT* landareaC * 1e6, (/1,2/))/1e15
d_Gcas_HR  = dim_sum_n(d_cas_HR * landareaC * 1e6, (/1,2/))/1e15
d_Gcas_AR  = dim_sum_n(d_cas_AR * landareaC * 1e6, (/1,2/))/1e15
d_Gtair    = dim_sum_n(d_TAIR   * landareaC, (/1,2/))/ sum(landarea)
d_GTSOI    = dim_sum_n(d_TSOI   * landareaC, (/1,2/))/ sum(landarea)
d_GTheta   = dim_sum_n(d_Theta  * landareaC, (/1,2/))/ sum(landarea)

d_Gmim_LIT = dim_sum_n(d_mim_LIT* landareaC * 1e6, (/1,2/))/1e15
d_Gmim_MIC = dim_sum_n(d_mim_MIC* landareaC * 1e6, (/1,2/))/1e15
d_Gmim_SOM = dim_sum_n(d_mim_SOM* landareaC * 1e6, (/1,2/))/1e15
d_Gmim_TOT = dim_sum_n(d_mim_TOT* landareaC * 1e6, (/1,2/))/1e15
d_Gmim_HR  = dim_sum_n(d_mim_HR * landareaC * 1e6, (/1,2/))/1e15
d_Gmim_IN  = dim_sum_n(d_mim_IN * landareaC * 1e6, (/1,2/))/1e15

d_Gcor_LIT = dim_sum_n(d_cor_LIT* landareaC * 1e6, (/1,2/))/1e15
d_Gcor_MIC = dim_sum_n(d_cor_MIC* landareaC * 1e6, (/1,2/))/1e15
d_Gcor_SOM = dim_sum_n(d_cor_SOM* landareaC * 1e6, (/1,2/))/1e15
d_Gcor_TOT = dim_sum_n(d_cor_TOT* landareaC * 1e6, (/1,2/))/1e15
d_Gcor_HR  = dim_sum_n(d_cor_HR * landareaC * 1e6, (/1,2/))/1e15

d_Gcas_SOM_pro = dim_sum_n(d_cas_SOM_pro* landareaC * 1e6, (/1,2/))/1e15
d_Gmim_SOM_pro = dim_sum_n(d_mim_SOM_pro* landareaC * 1e6, (/1,2/))/1e15
d_Gcor_SOM_pro = dim_sum_n(d_cor_SOM_pro* landareaC * 1e6, (/1,2/))/1e15

cas_NPP2 = where(cas_NPP .eq. 0, cas_NPP@_FillValue, cas_NPP)
mim_MIC2 = where(mim_MIC .eq. 0, mim_MIC@_FillValue, mim_MIC)
cor_MIC2 = where(cor_MIC .eq. 0, cor_MIC@_FillValue, cor_MIC)

pd_cas_NPP = dim_avg_n(cas_NPP2((ntim-10):(ntim-1),:,:),0) \
           / dim_avg_n(cas_NPP2(0:9,:,:),0) 
pd_mim_MIC = dim_avg_n(mim_MIC2((ntim-10):(ntim-1),:,:),0) \
           / dim_avg_n(mim_MIC2(0:9,:,:),0) 
pd_cor_MIC = dim_avg_n(cor_MIC2((ntim-10):(ntim-1),:,:),0) \
          /  dim_avg_n(cor_MIC2(0:9,:,:),0) 

pd_cas_NPP = (pd_cas_NPP -1 ) * 100
pd_mim_MIC = (pd_mim_MIC -1 ) * 100
pd_cor_MIC = (pd_cor_MIC -1 ) * 100

copy_VarCoords(landarea, pd_cas_NPP)
copy_VarCoords(landarea, pd_mim_MIC)
copy_VarCoords(landarea, pd_cor_MIC)

print("")
print("finished calculating chages in global pools")
print("")

;-------------------------------------------------------------------
; -- zonal mean (lat x day) --
;-------------------------------------------------------------------

 d_cas_NPP_zone = transpose(dim_avg_n_Wrap(d_cas_NPP,2))    ; lat x day
 d_TSOI_zone    = transpose(dim_avg_n_Wrap(d_TSOI,2))                   
 d_Theta_zone   = transpose(dim_avg_n_Wrap(d_Theta,2))                  

 d_cas_TOT_zone = transpose(dim_avg_n_Wrap(d_cas_TOT,2))    
 d_mim_TOT_zone = transpose(dim_avg_n_Wrap(d_mim_TOT,2))    
 d_cor_TOT_zone = transpose(dim_avg_n_Wrap(d_cor_TOT,2))    


;-------------------------------------------------------------------
; plot changes in stocks & fluxes over time
;-------------------------------------------------------------------
years         = ispan(1901,2100,1)
HR_all        = new((/4,ntim/),float)
HR_all!1      = "years"
HR_all&years  =  years

HR_runAvg      = HR_all
GPP_all        = HR_all
GPP_runAvg     = HR_all
NEE_all        = HR_all
LIT_all        = HR_all
SOM_all        = LIT_all
TOT_all        = LIT_all
ENV_all        = LIT_all

HR_all(0,:)   = 0.
HR_all(1,:)   = (/Gcas_HR/)
HR_all(2,:)   = (/Gmim_HR/)
HR_all(3,:)   = (/tofloat(Gcor_HR)/)

GPP_all(0,:)   = 0.
GPP_all(1,:)   = (/Gcas_AR /)
GPP_all(2,:)   = (/Gcas_GPP/)
GPP_all(3,:)   = (/Gcas_NPP/)

; cacluate 5 year running average for C fluxes
HR_runAvg(0,:) = 0.
HR_runAvg(1,:) = (/runave(Gcas_HR,  5, 0)/)
HR_runAvg(2,:) = (/runave(Gmim_HR,  5, 0)/)
HR_runAvg(3,:) = (/tofloat(runave(Gcor_HR,  5, 0))/)


GPP_runAvg(0,:) = 0.
GPP_runAvg(1,:) = (/runave(Gcas_AR , 5, 0)/)
GPP_runAvg(2,:) = (/runave(Gcas_GPP, 5, 0)/)
GPP_runAvg(3,:) = (/runave(Gcas_NPP, 5, 0)/)

; anomaly from 5 year running mean 
HR_anom       = HR_all  - HR_runAvg
GPP_anom      = GPP_all - GPP_runAvg
TAIR_anom     = Gtair    - runave(Gtair, 5,0)

NEE_all(0,:)  = 0.
NEE_all(1,:)  = Gcas_HR - Gcas_NPP       
NEE_all(2,:)  = Gmim_HR - Gcas_NPP      
NEE_all(3,:)  = tofloat(Gcor_HR) - Gcas_NPP       

LIT_all(0,:) = 0.
LIT_all(1,:) = (/d_Gcas_LIT/)
LIT_all(2,:) = (/d_Gmim_LIT/)
LIT_all(3,:) = (/tofloat(d_Gcor_LIT)/)

SOM_all(0,:) = 0.
SOM_all(1,:) = (/d_Gcas_SOM/)
SOM_all(2,:) = (/d_Gmim_SOM/)
SOM_all(3,:) = (/tofloat(d_Gcor_SOM)/)

TOT_all(0,:) = 0.
TOT_all(1,:) = (/d_Gcas_TOT/)
TOT_all(2,:) = (/d_Gmim_TOT/)
TOT_all(3,:) = (/tofloat(d_Gcor_TOT)/)

ENV_all(0,:) = 0.
ENV_all(1,:) = (/tofloat(d_GTSOI)/)
ENV_all(2,:) = (/tofloat(d_GTheta)/)  
ENV_all(3,:) = (/d_Gcas_NPP/)           
; -- write out csv files of annual values -- 
        years = toint((/HR_all&years/))
; anomaly fluxes
        fil_out = out_dir+"anom_fluxes.txt"
        system("/bin/rm -f "+fil_out)   ; remove any pre-existing file
        write_table(fil_out, "a", [/"YEAR","TAIR","GPP","NPP","AR","CAS_HR","MIM_HR","COR_HR"/], \
                                     "%s %s %s %s %s %s %s %s ")
        write_table(fil_out, "a", [/years        , TAIR_anom,                       \
                                     GPP_anom(2,:), GPP_anom(3,:), GPP_anom(1,:),    \
                                      HR_anom(1,:),  HR_anom(2,:),  HR_anom(3,:) /], \
                                   "%i %f %f %f %f %f %f %f")
       print("wrote "+ fil_out)    

; annual fluxes
        fil_out2 = out_dir+"ann_fluxes.txt"
        system("/bin/rm -f "+fil_out2)   ; remove any pre-existing file
        write_table(fil_out2, "a", [/"YEAR","TAIR","GPP","NPP","AR","CAS_HR","MIM_HR","COR_HR"/], \
                                     "%s %s %s %s %s %s %s %s ")
        write_table(fil_out2, "a", [/years        , (/tofloat(Gtair)/),                       \
                                     (/tofloat(Gcas_GPP)/), (/tofloat(Gcas_NPP)/), (/tofloat(Gcas_AR)/),    \
                                     (/tofloat(Gcas_HR)/),  (/tofloat(Gmim_HR )/), (/tofloat(Gcor_HR)/) /], \
                                     "%i %f %f %f %f %f %f %f")
       print("wrote "+ fil_out2)    

;--------------------------------------------
; -- make time series plots --
;--------------------------------------------
  fout  = out_dir+"Testbed_TS_plots_"+exper
  wks  = gsn_open_wks("ps",fout)

  res                        = True
;  res@gsnMaximize            = True
  res@gsnFrame               = False
  res@gsnDraw                = True
  res@tiXAxisString          = "Year"
  res@trXMinF                = 1901                   ; use this and next line to change timescale that is graphed
  res@trXMaxF                = 2100  
  res@pmLegendDisplayMode    = "Always"
  res@pmLegendSide           = "Top"
  res@pmLegendParallelPosF   = .29                 ; move units right
  res@pmLegendOrthogonalPosF = -0.325              ; move units down
  res@pmLegendWidthF         = 0.15                ; Change width and
  res@pmLegendHeightF        = 0.18                ; height of legend.
  res@lgLabelFontHeightF     = .015                 ; label font height
  res@lgPerimOn              = False               ; turn off box around

  res@xyLineThicknesses      = (/1,7,7,7/)
 ;set dashed lines for rcp45, solid for rcp58
 if (exper .eq. "RCP45") then
   res@xyDashPatterns        = (/1,16,16,16/)
 end if
 if (exper .eq. "RCP85") then
   res@xyDashPatterns        = (/1,0,0,0/)
 end if

; color blind friendly colors, these dddon't have enough contrast
;  colors = (/"black", \ 
;             "(/  0. , 0.286, 0.286 /)", \
;             "(/  0. , 0.427, 0.859 /)", \
;             "(/0.572, 0.286,  0.   /)"/)

  colors = (/"black", \
             "(/0.008, 0.498, 0.498 /)", \
;            "(/  0. , 0.286, 0.286 /)", \
             "(/0.286, 0.   , 0.572 /)", \
             "(/0.572, 0.286,  0.   /)"/)
;             "(/0.859, 0.819,  0.   /)"/)

  res@xyLineColors           = colors
  res@xyExplicitLegendLabels = (/""," CASA", " MIMICS"," CORPSE"/)
  res@lgItemOrder            = (/3,2,1,0/)           ; reverse legend order
  res@pmLegendSide           = "Bottom"
  res@pmLegendOrthogonalPosF = -0.450              ; move units down
  res@trYMinF          = -80.
  res@trYMaxF          =  20.
  res@tiYAxisString    = "~F33~D ~F21~Litter C (Pg C)"
  plot                 = gsn_csm_xy(wks,LIT_all&years,LIT_all,res)
  frame(wks)

  res@pmLegendSide     = "Top"
  res@pmLegendOrthogonalPosF = -0.325              ; move units down
  res@trYMinF          = -80.
  res@trYMaxF          =  80.
  res@tiYAxisString    = "~F33~D ~F21~SOC (Pg C)"
  plot                 = gsn_csm_xy(wks,SOM_all&years,SOM_all,res)
  frame(wks)

  delete(plot)

; -- make Fig_5, change in drivers and total soil C -- 
  plot = new(2, "graphic")
  res@gsnDraw          = False
  res@gsnFrame         = False
  res@pmLegendDisplayMode    = "Never"
  res@tiYAxisString    = "~F33~D ~F21~Total C (Pg C)"
  res@trYMinF       = -80.
  res@trYMaxF       =  80.
  plot(1)              = gsn_csm_xy(wks,TOT_all&years,TOT_all,res)
print("hist. TOT change, "+dim_avg_n(TOT_all(:,100:109),1))
print("RCP   TOT change, "+dim_avg_n(TOT_all(:,190:199),1))

  res1             = True
  res2             = True
  res3             = True

  res1@tmXBLabelsOn  = False                  ; turn off x labels
  res1@trXMinF       = 1901                   ; use this and next line to change timescale that is graphed
  res1@trXMaxF       = 2100
  res1@trYMinF       = -0.5
  res1@trYMaxF       =  4.0
  res1@tiXAxisString = ""
  res1@gsnFrame      = False
  res1@gsnDraw       = False
  res1@xyLineColors   = (/"black","red","blue"/)
  res1@xyLineThicknesses      = (/1,4,4/)
 ;set dashed lines for rcp45, solid for rcp58
 if (exper .eq. "RCP45") then
   res1@xyDashPatterns        = (/1,16,16/)
 end if
 if (exper .eq. "RCP85") then
   res1@xyDashPatterns        = (/1,0,0/)
 end if


  res2@tmXBLabelsOn  = False                  ; turn off x labels
  res2@trYMinF       = -5. 
  res2@trYMaxF       =  40.
  res2@xyLineColor   = "black"
  res2@xyLineThicknesses  = 5.                  ; thicker line
  if (exper .eq. "RCP45") then
    res2@xyDashPattern        = 16
  end if
  if (exper .eq. "RCP85") then
    res2@xyDashPattern        = 0
  end if

  res1@tiYAxisString    = "~F33~D ~F21~Soil Temperature & Moisture "      
  res2@tiYAxisString    = "~F33~D ~F21~Litter Inputs (Pg C)"

  plot(0) = gsn_csm_xy2(wks,ENV_all&years,ENV_all(0:2,:),ENV_all(3,:),\ 
                                          res1        ,res2          ) 
  resP                     = True               ; modify the panel plot
  resP@gsnFrame            = False              ; don't advance panel plot
; resP@gsnMaximize         = True               ; use full page
  resP@gsnPanelFigureStringsJust = "TopLeft"
  resP@gsnPanelFigureStringsFontHeightF = 0.02
  resP@gsnPanelFigureStringsBackgroundFillColor = -1 ;transparent
  resP@gsnPanelFigureStrings= (/"(a)","(b)"/)
  resP@gsnPanelFigureStringsPerimOn = False      ; turns off box around fig labels
  resP@gsnStringFont  = "helvetica-bold"
  resP@amJust         = "TopLeft"
  gsn_panel(wks,plot,(/2,1/),resP)               ; now draw as one plot
; frame(wks)

; --- end Fig_5 ---

  
  delete([/plot,res, res1, res2/])

  print("------------------------")
  print("wrote "+ fout+".ps")
  print("")

;-------------------------------------------------------------------
; Calculate initial pools / fluxes
;-------------------------------------------------------------------

  ; average over first decade & convert pools to kgC/m2
  TAIR_i    = dim_avg_n_Wrap(tair(0:9,:,:), 0) 
  TSOI_i    = dim_avg_n_Wrap(TSOI(0:9,:,:), 0)
  Theta_i   = dim_avg_n_Wrap(Theta(0:9,:,:), 0)

  cas_GPP_i = dim_avg_n_Wrap(cas_GPP(0:9,:,:), 0)
  cas_GPP_f = dim_avg_n_Wrap(cas_GPP(ntim-10:ntim-1,:,:), 0)
  cas_NPP_i = dim_avg_n_Wrap(cas_NPP(0:9,:,:), 0)
  cas_NPP_f = dim_avg_n_Wrap(cas_NPP(ntim-10:ntim-1,:,:), 0)
  cas_HR_i  = dim_avg_n_Wrap(cas_HR(0:9,:,:), 0)
  mim_HR_i  = dim_avg_n_Wrap(mim_HR(0:9,:,:), 0)
  cor_HR_i  = dim_avg_n_Wrap(cor_HR(0:9,:,:), 0)
  cas_HR_f  = dim_avg_n_Wrap(cas_HR(ntim-10:ntim-1,:,:), 0)
  mim_HR_f  = dim_avg_n_Wrap(mim_HR(ntim-10:ntim-1,:,:), 0)
  cor_HR_f  = dim_avg_n_Wrap(cor_HR(ntim-10:ntim-1,:,:), 0)

  cas_LIT_i = dim_avg_n_Wrap(cas_LIT(0:9,:,:), 0) *1e-3
  mim_LIT_i = dim_avg_n_Wrap(mim_LIT(0:9,:,:), 0) *1e-3
  cor_LIT_i = dim_avg_n_Wrap(cor_LIT(0:9,:,:), 0) *1e-3

  mim_MIC_i = dim_avg_n_Wrap(mim_MIC(0:9,:,:), 0) *1e-3 ;kgC/m2
  cor_MIC_i = dim_avg_n_Wrap(cor_MIC(0:9,:,:), 0) *1e-3
  mim_MIC_f = dim_avg_n_Wrap(mim_MIC(ntim-10:ntim-1,:,:),0)*1e-3
  cor_MIC_f = dim_avg_n_Wrap(cor_MIC(ntim-10:ntim-1,:,:),0)*1e-3

  cas_SOM_i = dim_avg_n_Wrap(cas_SOM(0:9,:,:), 0) *1e-3
  mim_SOM_i = dim_avg_n_Wrap(mim_SOM(0:9,:,:), 0) *1e-3
  cor_SOM_i = dim_avg_n_Wrap(cor_SOM(0:9,:,:), 0) *1e-3

  cas_SOM_pro_i = dim_avg_n_Wrap(cas_SOM_pro(0:9,:,:), 0) *1e-3
  mim_SOM_pro_i = dim_avg_n_Wrap(mim_SOM_pro(0:9,:,:), 0) *1e-3
  cor_SOM_pro_i = dim_avg_n_Wrap(cor_SOM_pro(0:9,:,:), 0) *1e-3
  cas_SOM_pro_f = dim_avg_n_Wrap(cas_SOM_pro(ntim-10:ntim-1,:,:), 0) *1e-3
  mim_SOM_pro_f = dim_avg_n_Wrap(mim_SOM_pro(ntim-10:ntim-1,:,:), 0) *1e-3
  cor_SOM_pro_f = dim_avg_n_Wrap(cor_SOM_pro(ntim-10:ntim-1,:,:), 0) *1e-3

  cas_TOT_i = dim_avg_n_Wrap(cas_TOT(0:9,:,:), 0) *1e-3
  mim_TOT_i = dim_avg_n_Wrap(mim_TOT(0:9,:,:), 0) *1e-3
  cor_TOT_i = dim_avg_n_Wrap(cor_TOT(0:9,:,:), 0) *1e-3
  cas_TOT_f = dim_avg_n_Wrap(cas_TOT(ntim-10:ntim-1,:,:), 0) *1e-3
  mim_TOT_f = dim_avg_n_Wrap(mim_TOT(ntim-10:ntim-1,:,:), 0) *1e-3
  cor_TOT_f = dim_avg_n_Wrap(cor_TOT(ntim-10:ntim-1,:,:), 0) *1e-3

  cas_TAU_i = dim_avg_n_Wrap(cas_TAU(0:9,:,:), 0)
  mim_TAU_i = dim_avg_n_Wrap(mim_TAU(0:9,:,:), 0)
  cor_TAU_i = dim_avg_n_Wrap(cor_TAU(0:9,:,:), 0)
  cas_TAU_f = dim_avg_n_Wrap(cas_TAU(ntim-10:ntim-1,:,:), 0)
  mim_TAU_f = dim_avg_n_Wrap(mim_TAU(ntim-10:ntim-1,:,:), 0)
  cor_TAU_f = dim_avg_n_Wrap(cor_TAU(ntim-10:ntim-1,:,:), 0)

  ; calculate global mic fraction of total SOM
  mim_fracMIC_i = dim_avg_n_Wrap(mim_fracMIC(0:9,:,:), 0) * 100
  cor_fracMIC_i = dim_avg_n_Wrap(cor_fracMIC(0:9,:,:), 0) * 100

  copy_VarCoords(landarea, cas_GPP_i)
  copy_VarCoords(landarea, cas_GPP_f)
  copy_VarCoords(landarea, cas_NPP_i)
  copy_VarCoords(landarea, cas_NPP_f)
  copy_VarCoords(landarea, cas_LIT_i)
  copy_VarCoords(landarea, mim_LIT_i)
  copy_VarCoords(landarea, cor_LIT_i)
  copy_VarCoords(landarea, mim_MIC_i)
  copy_VarCoords(landarea, cor_MIC_i)
  copy_VarCoords(landarea, mim_MIC_f)
  copy_VarCoords(landarea, cor_MIC_f)

  copy_VarCoords(landarea, cas_SOM_i)
  copy_VarCoords(landarea, mim_SOM_i)
  copy_VarCoords(landarea, cor_SOM_i)
  copy_VarCoords(landarea, cas_SOM_pro_i)
  copy_VarCoords(landarea, mim_SOM_pro_i)
  copy_VarCoords(landarea, cor_SOM_pro_i)
  copy_VarCoords(landarea, cas_SOM_pro_f)
  copy_VarCoords(landarea, mim_SOM_pro_f)
  copy_VarCoords(landarea, cor_SOM_pro_f)

  copy_VarCoords(landarea, cas_TOT_i)
  copy_VarCoords(landarea, mim_TOT_i)
  copy_VarCoords(landarea, cor_TOT_i)
  copy_VarCoords(landarea, cas_TOT_f)
  copy_VarCoords(landarea, mim_TOT_f)
  copy_VarCoords(landarea, cor_TOT_f)

  copy_VarCoords(landarea, cas_TAU_i)
  copy_VarCoords(landarea, mim_TAU_i)
  copy_VarCoords(landarea, cor_TAU_i)
  copy_VarCoords(landarea, cas_TAU_f)
  copy_VarCoords(landarea, mim_TAU_f)
  copy_VarCoords(landarea, cor_TAU_f)

  copy_VarCoords(landarea, mim_fracMIC_i)
  copy_VarCoords(landarea, cor_fracMIC_i)

  cas_LIT_i@units = "kg C m~S~-2~N~"
  mim_LIT_i@units = "kg C m~S~-2~N~"
  cor_LIT_i@units = "kg C m~S~-2~N~"

  cas_SOM_i@units = "kg C m~S~-2~N~"
  mim_SOM_i@units = "kg C m~S~-2~N~"
  cor_SOM_i@units = "kg C m~S~-2~N~"
  cas_SOM_pro_i@units = "kg C m~S~-2~N~"
  mim_SOM_pro_i@units = "kg C m~S~-2~N~"
  cor_SOM_pro_i@units = "kg C m~S~-2~N~"
  cas_SOM_pro_f@units = "kg C m~S~-2~N~"
  mim_SOM_pro_f@units = "kg C m~S~-2~N~"
  cor_SOM_pro_f@units = "kg C m~S~-2~N~"

  cas_TOT_i@units = "kg C m~S~-2~N~"
  mim_TOT_i@units = "kg C m~S~-2~N~"
  cor_TOT_i@units = "kg C m~S~-2~N~"
  cas_TOT_f@units = "kg C m~S~-2~N~"
  mim_TOT_f@units = "kg C m~S~-2~N~"
  cor_TOT_f@units = "kg C m~S~-2~N~"

  cas_TAU_i@units = "y"
  mim_TAU_i@units = "y"
  cor_TAU_i@units = "y"
  cas_TAU_f@units = "y"
  mim_TAU_f@units = "y"
  cor_TAU_f@units = "y"

  cas_GPP_i@long_name = cas_GPP@long_name
  cas_NPP_i@long_name = cas_NPP@long_name
  cas_NPP_f@long_name = cas_NPP@long_name
  cas_HR_i@long_name = cas_HR@long_name
  mim_HR_i@long_name = mim_HR@long_name
  cor_HR_i@long_name = cor_HR@long_name

  cas_LIT_i@long_name = cas_LIT@long_name
  mim_LIT_i@long_name = mim_LIT@long_name
  cor_LIT_i@long_name = cor_LIT@long_name

  cas_SOM_i@long_name = cas_SOM@long_name
  mim_SOM_i@long_name = mim_SOM@long_name
  cor_SOM_i@long_name = cor_SOM@long_name
  cas_SOM_pro_i@long_name = cas_SOM_pro@long_name
  mim_SOM_pro_i@long_name = mim_SOM_pro@long_name
  cor_SOM_pro_i@long_name = cor_SOM_pro@long_name

  cas_TOT_i@long_name = cas_TOT@long_name
  mim_TOT_i@long_name = mim_TOT@long_name
  cor_TOT_i@long_name = cor_TOT@long_name

  ; changes in environmental scalars
  cas_fT_i  = dim_avg_n_Wrap(cas_fT(0:9,:,:), 0)
  cas_fW_i  = dim_avg_n_Wrap(cas_fW(0:9,:,:), 0)
  mim_fW_i  = dim_avg_n_Wrap(mim_fW(0:9,:,:), 0)
  cor_fW_i  = dim_avg_n_Wrap(cor_fW(0:9,:,:), 0)
  d_cas_fT  = dim_avg_n_Wrap(cas_fT(ntim-10:ntim-1,:,:), 0)
  d_cas_fW  = dim_avg_n_Wrap(cas_fW(ntim-10:ntim-1,:,:), 0)
  d_mim_fW  = dim_avg_n_Wrap(mim_fW(ntim-10:ntim-1,:,:), 0)
  d_cor_fW  = dim_avg_n_Wrap(cor_fW(ntim-10:ntim-1,:,:), 0)

  d_cas_fT  = d_cas_fT - cas_fT_i 
  d_cas_fW  = d_cas_fW - cas_fW_i 
  d_mim_fW  = d_mim_fW - mim_fW_i 
  d_cor_fW  = d_cor_fW - cor_fW_i 

  r_cas_fT  = dim_avg_n_Wrap(cas_fT(ntim-10:ntim-1,:,:), 0)
  r_cas_fW  = dim_avg_n_Wrap(cas_fW(ntim-10:ntim-1,:,:), 0)
  r_mim_fW  = dim_avg_n_Wrap(mim_fW(ntim-10:ntim-1,:,:), 0)
  r_cor_fW  = dim_avg_n_Wrap(cor_fW(ntim-10:ntim-1,:,:), 0)

  cas_fT_d  = where(cas_fT_i.le.0, cas_fT_i@_FillValue, cas_fT_i)
  cas_fW_d  = where(cas_fW_i.le.0, cas_fW_i@_FillValue, cas_fW_i)
  mim_fW_d  = where(mim_fW_i.le.0, mim_fW_i@_FillValue, mim_fW_i)
  cor_fW_d  = where(cor_fW_i.le.0, cor_fW_i@_FillValue, cor_fW_i)

  r_cas_fT  = r_cas_fT / cas_fT_d
  r_cas_fW  = r_cas_fW / cas_fW_d
  r_mim_fW  = r_mim_fW / mim_fW_d
  r_cor_fW  = r_cor_fW / cor_fW_d

; ------------------------------------------------------
;   Write out .nc file & maps   
; ------------------------------------------------------
; ------------------------------------------------------

   cdffile = out_dir+ "init_TAU.nc"
   system("rm -f " +  cdffile)
   a = addfile(cdffile,"c")
   a@source = systemfunc("pwd")+"/"+get_script_name()

   a->cas_TAU   = cas_TAU_i
   a->mim_TAU   = mim_TAU_i
   a->cor_TAU   = cor_TAU_i
   a->r_cas_fT  = r_cas_fT
   a->r_cas_fW  = r_cas_fW
   a->r_mim_fW  = r_mim_fW
   a->r_cor_fW  = tofloat(r_cor_fW)
print(" Wrote " + cdffile)
print("")

;-------------------------------------------------------------------
;-------------------------------------------------------------------
; plotting turnover times of initial stocks
; start Fig_4 
;-------------------------------------------------------------------
;-------------------------------------------------------------------
 fout      = out_dir+"Testbed_initial_turnover_"+exper 
 wks       = gsn_open_wks("pdf",fout)
 ny        = 3
 nx        = nlat*nlon
 ; array for initial turnover times
 data2d      = new ((/ny,nx/), "float")
 data2d(0,:) = (/ndtooned(cas_TAU_i)/)
 data2d(1,:) = (/ndtooned(mim_TAU_i)/)
 data2d(2,:) = (/ndtooned(tofloat(cor_TAU_i))/)

 ; array for change turnover times
 delta2d      = new ((/ny,nx/), "float")
 delta2d(0,:) = (/ndtooned(cas_TAU_f - cas_TAU_i)/)
 delta2d(1,:) = (/ndtooned(mim_TAU_f - mim_TAU_i)/)
 delta2d(2,:) = (/ndtooned(tofloat(cor_TAU_f - cor_TAU_i))/)

 theta1d     = (/ndtooned(Theta_i)/)
 levels    = fspan(0,75, 6)
 nlevels   = dimsizes(levels)
  labels   = new((nlevels),string)
  labels   = ""
  labels   = "" + levels
; print(labels)

;---Create new 3D array to hold groupings of values
  data3d = new((/nlevels-1,ny,nx/),typeof(data2d))
  delta3d = new((/nlevels-1,ny,nx/),typeof(delta2d))

;---Group the values by level and put in 3D array.
  do n=0,nlevels-2
    do j=0,ny-1
      data3d(n,j,:) = where(theta1d .ge.levels(n).and.theta1d.lt.levels(n+1), \
                            data2d(j,:), data2d@_FillValue)
      delta3d(n,j,:) = where(theta1d .ge.levels(n).and.theta1d.lt.levels(n+1), \
                            delta2d(j,:), delta2d@_FillValue)
      ;printMinMax(data3d(n,j,:), 0)
    end do
  end do

; -- select color map --
  gsn_define_colormap(wks, "MPL_gist_earth")
  gsn_reverse_colormap(wks)
  nstep   = 25; 99 /(nlevels-1) 
  colors2 = ispan(28,128,nstep)

 res                     = True                   ; plot mods desired
 res@gsnFrame            = False
 res@gsnDraw             = False
 res@gsnMaximize         = False                  ; maximize plot in frame
 res@xyMarkLineMode      = "Markers"              ; choose which have markers
 res@xyMarker            = 16                      ; choose type of marker
 res@xyMarkerSizeF       = 0.008                   ; Marker size (default 0.01)
 res@xyMarkerOpacities   = 0.80                    ; make markers semi-transparent
 res@trYMinF             =  1                      ; axis max
 res@trYMaxF             = 5000                    ; axis min
 res@xyYStyle            = "Log"
 res@trXMinF             = -22                     ; axis max
 res@trXMaxF             = 30                      ; axis min
 res@tiYAxisString       = "Turnover time (y)" 
 res@tiXAxisString       = "Mean Air Temp (~S~o~N~C)"
 res@tmEqualizeXYSizes   = "True"                  ; X & Y labels equal
 res@tmXBLabelFontHeightF= 0.03                    ; default 0.02 for a viewport width of 0.6
;res@vpXF             = 0.12                ; default is 0.2 (aspect ratio)
;res@vpYF             = 0.8                 ; default is 0.8
;res@vpHeightF        = 0.4                 ; default is 0.6
;res@vpWidthF         = 0.4                 ; default is 0.6

 TAIR1d                  = ndtooned(TAIR_i)

 res@tiYAxisString       = ""
 res@xyMarkerColor       = colors2(0) 
 overlays = new(nlevels-1,graphic)
 plot0 = gsn_csm_xy (wks,TAIR1d,data3d(0,0,:),res)
 do n=1,nlevels-2
     res@xyMarkerColor       = colors2(n)
     overlays(n-1) = gsn_csm_xy (wks,TAIR1d,data3d(n,0,:),res)
     overlay(plot0,overlays(n-1))
 end do
 delete(overlays) 

 res@tiYAxisString       = "Inferred Turnover Time (y)"
 overlays = new(nlevels-1,graphic)
 res@xyMarkerColor       = colors2(0)
 plot1 = gsn_csm_xy (wks,TAIR1d,data3d(0,1,:),res)
 do n=1,nlevels-2
     res@xyMarkerColor       = colors2(n)
     overlays(n-1) = gsn_csm_xy (wks,TAIR1d,data3d(n,1,:),res)
     overlay(plot1,overlays(n-1))
 end do
 delete(overlays)

 overlays = new(nlevels-1,graphic)
 res@tiYAxisString       = ""
 res@xyMarkerColor       = colors2(0)
 plot2 = gsn_csm_xy (wks,TAIR1d,data3d(0,2,:),res)
 do n=1,nlevels-2
     res@xyMarkerColor       = colors2(n)
     overlays(n-1) = gsn_csm_xy (wks,TAIR1d,data3d(n,2,:),res)
     overlay(plot2,overlays(n-1))
 end do
 delete(overlays)

;---Attach a labelbar
  labelbar = attach_labelbar(wks,plot1,labels,colors2)

; Add pannel ID to each plot
   txres               = True
   txres@txFontHeightF = 0.033
   txres@txFont        = 22 ; helvetica bold
;  txres@txFontThicknessF = 2.0
   txres@txFontColor   = "black"
;  txres@txBackgroundFillColor = "White"
   text0 = gsn_add_text(wks,plot0,"a",-18.0,3150,txres)
   text1 = gsn_add_text(wks,plot1,"b",-18.0,3150,txres)
   text2 = gsn_add_text(wks,plot2,"c",-18.0,3150,txres)


; calculate error, based on equation for Charlie's line
; from excel... y = 0.0008x2 - 0.0401x + 1.5636

pred  = (0.0008 * TAIR1d^2 - 0.0401*TAIR1d + 1.5636)
predC = conform_dims(dimsizes(data2d), pred, 1)
error = log10(data2d) - predC

 delete([/res@trYMinF, res@trYMaxF,labels/])

 Koven = "/project/tss/wwieder/CASACLM/GRID/"
 xdata = readAsciiTable(Koven+"x_confint.txt", 1, "float", 0)
 ydata = asciiread(Koven+"y_confint.txt", (/3,51/), "float")

 delete([/res@xyMarker,res@xyMarkLineMode /])
 res@trYMinF             =  1                      ; axis max
 res@trYMaxF             = 5000                    ; axis min
 res@xyYStyle            = "Log"
 res@pmLegendDisplayMode    = "Never"            ; legend resources
 res@xyMarkLineMode    = "Lines"        ; choose which have markers
 res@xyDashPatterns    = (/0,1,1/)      ; choose type of lines
 res@xyLineColor       = "black"        ; Line color
 res@xyLineThicknesses = (/4,3,3/)
 plotA0 = gsn_csm_xy (wks,xdata(:,0),ydata(:,:),res)        ; create plot
 plotA1 = gsn_csm_xy (wks,xdata(:,0),ydata(:,:),res)        ; create plot
 plotA2 = gsn_csm_xy (wks,xdata(:,0),ydata(:,:),res)        ; create plot

 overlay(plot0, plotA0)
 overlay(plot1, plotA1)
 overlay(plot2, plotA2)

;now attach plots all together
  attachres1 = True

  attachres1@gsnAttachPlotsXAxis = True ; attaches along x-axis 

  attachres2 = True

  attachid1  = gsn_attach_plots(plot0,(/plot1,plot2/),attachres1,attachres2)

  maximize_output(wks,False)

  print("------------------------")
  print("wrote "+ fout+".ps")
  print("")

  delete([/ res@xyDashPatterns, res@xyLineThicknesses, \
            plot0, plot1, plot2, plotA0, plotA1, plotA2 /])

; -- end   Fig_4 --

; make Fig_S4, average line of tau vs. temp for each model

  CAS_coef = lspoly(TAIR1d, log10(data2d(0,:)), 1, 3)
  MIM_coef = lspoly(TAIR1d, log10(data2d(1,:)), 1, 3)
  COR_coef = lspoly(TAIR1d, log10(data2d(2,:)), 1, 3)
; print(MIM_coef)

  CAS_mean = CAS_coef(0) + CAS_coef(1) * xdata(:,0) + CAS_coef(2) * xdata(:,0)^2
  MIM_mean = MIM_coef(0) + MIM_coef(1) * xdata(:,0) + MIM_coef(2) * xdata(:,0)^2
  COR_mean = COR_coef(0) + COR_coef(1) * xdata(:,0) + COR_coef(2) * xdata(:,0)^2
  
  res@xyDashPattern     = 0
  res@xyLineThicknesses = 6

  if (exper .eq. "_loam/") then
      res@xyDashPattern     = 16  ; wider spaces between dashes
      res@xyLineThicknesses = 4
  end if

  res@xyLineColor       = colors(1)
  plot0    = gsn_csm_xy (wks,xdata(:,0),10^CAS_mean,res)
  res@tiYAxisString       = "Inferred Turnover Time (y)"
  res@xyLineColor  = colors(2)
  plot1    = gsn_csm_xy (wks,xdata(:,0),10^MIM_mean,res)
  res@tiYAxisString       = ""
  res@xyLineColor  = colors(3)
  plot2    = gsn_csm_xy (wks,xdata(:,0),10^COR_mean,res)

   text0 = gsn_add_text(wks,plot0,"(b) CASA-CNP",-11.0,3300,txres)
   text1 = gsn_add_text(wks,plot1,"(d) MIMICS",-13.0,3300,txres)
   text2 = gsn_add_text(wks,plot2,"(f) CORPSE",-12.5,3300,txres) 
  delete(res@xyDashPattern)
  delete(res@xyLineThicknesses)

  res@xyDashPatterns    = (/0,1,1/)      ; choose type of lines
  res@xyLineColor       = "black"        ; Line color
  res@xyLineThicknesses = (/4,3,3/)

  plotA0 = gsn_csm_xy (wks,xdata(:,0),ydata(:,:),res)        ; create plot
  plotA1 = gsn_csm_xy (wks,xdata(:,0),ydata(:,:),res)        ; create plot
  plotA2 = gsn_csm_xy (wks,xdata(:,0),ydata(:,:),res)        ; create plot

  overlay(plot0, plotA0)
  overlay(plot1, plotA1)
  overlay(plot2, plotA2)

  attachid1  = gsn_attach_plots(plot0,(/plot1,plot2/),attachres1,attachres2)

  maximize_output(wks,False)

  delete([/res,plot0, plot1, plot2/])


;-------------------------------------------------------------------
;-------------------------------------------------------------------
; plot the map of Total soil carbon
;-------------------------------------------------------------------
;-------------------------------------------------------------------

 res                     = True            ; plot mods desired
 res@gsnMaximize         = True  
 res@gsnPaperOrientation = "portrait"
 res@gsnDraw             = False
 res@gsnFrame            = False
 res@cnSmoothingOn       = False
 res@mpProjection        = "Robinson"
 res@mpLimitMode         = "LatLon"
 res@mpMinLatF           = -60.             ; get rid of antarctica
;res@mpOutlineOn         = True
;res@mpPerimOn           = False
;res@mpGridAndLimbOn     = True         ;draws Limb line around plots, w/o lat-lon lines
;res@mpGridLatSpacingF   = 180
;res@mpGridLonSpacingF   = 180
;res@mpGridLineThicknessF= 0.
;res@mpGridLineColor     = "transparent"
;res@mpGridMaskMode      = "MaskLand"
 res@cnInfoLabelOn       = False           ; turn off cn info label
 res@cnFillOn            = True            ; turn on color
 res@cnLinesOn           = False           ; turn off contour lines
 res@cnLineLabelsOn      = False           ; turn off contour labels
 res@gsnSpreadColors     = True            ; spread out color table
 res@lbOrientation       = "Horizontal"
 res@lbLabelBarOn        = False           ; turn off individual cb's

 res@vpXF             = 0.12                ; default is 0.2 (aspect ratio)
 res@vpYF             = 0.8                 ; default is 0.8
 res@vpHeightF        = 0.4                 ; default is 0.6
 res@vpWidthF         = 0.75                ; default is 0.6

  res@lbLabelFontHeightF   =.022                   ; make labels larger

  resP                  = True               ; modify the panel plot
  resP@gsnFrame         = False              ; don't advance panel plot
  resP@gsnPanelLabelBar = True               ; add common colorbar
  resP@lbOrientation    =   "horizontal"     ; vertical label bar


;-------------------------------------------------------------------
;-------------------------------------------------------------------
; plot the map of Total litter + microbe + soil carbon
;-------------------------------------------------------------------
;-------------------------------------------------------------------
  fout  = out_dir+"TestbedInitial_SOMplots_"+exper 
  wks3  = gsn_open_wks("ps",fout)
  gsn_define_colormap(wks3, "WhViBlGrYeOrRe")
  plot = new(3, "graphic")
  print(" starting initial map")
  ; Fig_S2 (a)
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels             = fspan(0,1600,17)
  res@gsnLeftString        = ""
  res@gsnRightString       = ""
  res@lbLabelBarOn         = True            ; turn off individual cb's

  plot(0) = gsn_csm_contour_map(wks3,cas_NPP_i,res)
  plot(1) = gsn_csm_contour_map(wks3,cas_NPP_f,res)

  resP@gsnPanelLabelBar = False              ; add common colorbar
  resP@gsnPanelFigureStringsJust = "TopLeft"
  resP@gsnPanelFigureStringsFontHeightF = 0.022 
  resP@gsnPanelFigureStrings= (/"(a)","(b)"/)
  resP@gsnPanelFigureStringsBackgroundFillColor = -1 ;transparent
; resP@gsnPanelFigureStrings= (/decimalPlaces(avg(Gcas_NPP(0:9           )),1,True)+ " Pg C y~S~-1~N~", \
;                               decimalPlaces(avg(Gcas_NPP(ntim-10:ntim-1)),1,True)+ " Pg C y~S~-1~N~"/)
  resP@gsnPanelFigureStringsPerimOn = False      ; turns off box around fig labels
  resP@gsnStringFont  = "helvetica-bold"
  resP@amJust         = "TopLeft" 
  resP@gsnPanelMainString = "inital, final NPP"             
  gsn_panel(wks3,plot,(/2,2/),resP)               ; now draw as one plot
  frame(wks3)
  delete([/resP@gsnPanelFigureStrings,res@cnLevels/])

; -- Microbial biomass (% SOC)--
  res@lbLabelBarOn         = False           ; turn off individual cb's
  res@cnLevels             = fspan(0,5.,11)
  plot(0) = gsn_csm_contour_map(wks3,mim_fracMIC_i,res)
  plot(1) = gsn_csm_contour_map(wks3,cor_fracMIC_i,res)

  resP@gsnPanelLabelBar = True               ; add common colorbar
  resP@gsnPanelMainString = "microbial biomass, % Total C"  
  resP@gsnPanelFigureStringsFontHeightF = 0.022
  resP@gsnPanelFigureStrings= (/"(a)","(b)"/)
; resP@gsnPanelFigureStrings= (/decimalPlaces(avg(Gmim_fracMIC(0:9)*100),3,True)+ "% ", \
;                               decimalPlaces(avg(Gcor_fracMIC(0:9)*100),3,True)+ "% "  /)
  gsn_panel(wks3,plot,(/2,2/),resP)               ; now draw as one plot
  frame(wks3)
  delete([/res@cnLevels /])

; -- Microbial biomass (stock)--
  res@cnLevels             = fspan(0,0.3,13)
  plot(0) = gsn_csm_contour_map(wks3,mim_MIC_i,res)
  plot(1) = gsn_csm_contour_map(wks3,cor_MIC_i,res)

  resP@gsnPanelLabelBar = True               ; add common colorbar
  resP@gsnPanelMainString = "microbial biomass MIM, COR (kgC/m2)"
  resP@gsnPanelFigureStringsFontHeightF = 0.022
  gsn_panel(wks3,plot,(/2,2/),resP)               ; now draw as one plot
  frame(wks3)

  delete([/plot,resP@gsnPanelFigureStrings, res@cnLevels/]);, res@cnLevelSelectionMode /])

; -- initial specific respiration rates
; this should likely be calculated instantenously, 
; current units gCO2 / gMIC / year
  plot = new(4, "graphic")

  xMIMi = where(mim_MIC_i.le.0., mim_MIC_i@_FillValue, mim_MIC_i*1e3) ;convert to gC/m2 from kgC/m2     
  xCORi = where(cor_MIC_i.le.0., cor_MIC_i@_FillValue, cor_MIC_i*1e3)
  xMIMi = mim_HR_i / xMIMi     ;gCO2/gMIC/y
  xCORi = cor_HR_i / xCORi
  copy_VarCoords(mim_HR_i, xMIMi)
  copy_VarCoords(cor_HR_i, xCORi)
; printVarSummary(mim_HR_i)
; printVarSummary(mim_MIC_i)

  xMIMf = where(mim_MIC_f.le.0., mim_MIC_f@_FillValue, mim_MIC_f*1e3) ;convert to gC/m2 (from kgC/m2)
  xCORf = where(cor_MIC_f.le.0., cor_MIC_f@_FillValue, cor_MIC_f*1e3)
  xMIMf = mim_HR_f / xMIMf
  xCORf = cor_HR_f / xCORf
  copy_VarCoords(mim_HR_f, xMIMf)
  copy_VarCoords(cor_HR_f, xCORf)

  dxMIM  = xMIMf-xMIMi
  dxCOR  = xCORf-xCORi
  copy_VarCoords(xMIMi,dxMIM)
  copy_VarCoords(xCORi,dxCOR)

  res@cnLevels             = fspan(0,12,13)
  res@lbLabelBarOn         = True            ; turn off individual cb's
  plot(0) = gsn_csm_contour_map(wks3,xMIMi, res)
  plot(1) = gsn_csm_contour_map(wks3,xCORi, res)


  res@cnLevels             = (fspan(-2,2,13))
  plot(2) = gsn_csm_contour_map(wks3,dxMIM, res)
  plot(3) = gsn_csm_contour_map(wks3,dxCOR, res)

  resP@gsnPanelLabelBar = False              ; add common colorbar
  resP@gsnPanelMainString = "specific respiration (gC/gMIC/y)"
  resP@gsnPanelFigureStringsFontHeightF = 0.022
  resP@gsnPanelFigureStrings= (/"(a)","(b)","(c)","(d)"/)
  gsn_panel(wks3,plot,(/2,2/),resP)               ; now draw as one plot
  frame(wks3)
  delete([/plot,resP@gsnPanelFigureStrings,res@cnLevels,res@cnLevelSelectionMode /])


; -- initial temp & moisture plots
; Fig_S1 a, b
  plot = new(3, "graphic")
  res@lbLabelBarOn         = True            ; turn off individual cb's
  plot(0) = gsn_csm_contour_map(wks3,TSOI_i       ,res)
  plot(1) = gsn_csm_contour_map(wks3,Theta_i      ,res)

  resP@gsnPanelLabelBar = False              ; add common colorbar
  resP@gsnPanelMainString = "initial tsoi, moisture" 
  resP@gsnPanelFigureStringsFontHeightF = 0.022
  resP@gsnPanelFigureStrings= (/"(a)","(b)"/)
  gsn_panel(wks3,plot,(/2,2/),resP)               ; now draw as one plot
  frame(wks3)
; delete([/resP@gsnPanelFigureStrings /])
; print("-- init temp, theta")

;repeat above, but with new color par for Fig_S1b
  gsn_define_colormap(wks3, "MPL_gist_earth")
  gsn_reverse_colormap(wks3)
; nstep   = 25; 99 /(nlevels-1)
; print(nstep)
; colors2 = ispan(28,128,nstep)

  plot(0) = gsn_csm_contour_map(wks3,TSOI_i       ,res)
  plot(1) = gsn_csm_contour_map(wks3,Theta_i      ,res)

  resP@gsnPanelLabelBar = False              ; add common colorbar
  resP@gsnPanelMainString = "initial tsoi, moisture"
  resP@gsnPanelFigureStringsFontHeightF = 0.022
  resP@gsnPanelFigureStrings= (/"(a)","(b)"/)
  gsn_panel(wks3,plot,(/2,2/),resP)               ; now draw as one plot
  frame(wks3)
  delete([/resP@gsnPanelFigureStrings /])
; print("-- init temp, theta")



;-- initial HR plots -- 
; not plotted for expediency
  gsn_define_colormap(wks3, "WhViBlGrYeOrRe")

  res@lbLabelBarOn         = False           ; turn off individual cb's
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels             = fspan(0,1600,17)
; plot(0) = gsn_csm_contour_map(wks3,cas_HR_i,res)
; plot(1) = gsn_csm_contour_map(wks3,mim_HR_i,res)
; plot(2) = gsn_csm_contour_map(wks3,cor_HR_i,res)

  resP@gsnPanelLabelBar = True               ; add common colorbar
  resP@gsnPanelMainString = "Initail HR, gC/m2/y"
  resP@gsnPanelFigureStrings= (/"(a)","(b)","(c)"/)
; resP@gsnPanelFigureStrings= (/decimalPlaces(avg(Gcas_HR(0:9)),1,True)+ " Pg C y~S~-1~N~", \
;                               decimalPlaces(avg(Gmim_HR(0:9)),1,True)+ " Pg C y~S~-1~N~", \
;                               decimalPlaces(avg(Gcor_HR(0:9)),1,True)+ " Pg C y~S~-1~N~"  /)
  resP@gsnPanelFigureStringsPerimOn = False      ; turns off box around fig labels
  resP@gsnStringFont  = "helvetica-bold"
  resP@amJust         = "TopLeft"
; gsn_panel(wks3,plot,(/3,1/),resP)               ; now draw as one plot
; frame(wks3)
  delete([/resP@gsnPanelFigureStrings,res@cnLevels/])
; print("-- init HR")

; -- initial TURNOVER TIME plots --
  res@cnMinLevelValF       =  0.01
  res@cnMaxLevelValF       =  1000.  
  res@cnLevels             = ((/0.1,1,2,5,10,20,50,100,200,300,500,1000/))
  plot(0) = gsn_csm_contour_map(wks3,cas_TAU_i,res)
  plot(1) = gsn_csm_contour_map(wks3,mim_TAU_i,res)
  plot(2) = gsn_csm_contour_map(wks3,cor_TAU_i,res)

  resP@gsnPanelMainString = "Initail Turnover time (y)"
  resP@gsnPanelFigureStrings= (/"(a)","(b)","(c)"/)
; resP@gsnPanelFigureStrings= (/decimalPlaces(avg(Gcas_TAU(0:9)),1,True)+ " y ", \
;                               decimalPlaces(avg(Gmim_TAU(0:9)),1,True)+ " y", \
;                               decimalPlaces(avg(Gcor_TAU(0:9)),1,True)+ " y"  /)
  resP@gsnPanelFigureStringsPerimOn = False      ; turns off box around fig labels
  resP@gsnStringFont  = "helvetica-bold"
  resP@amJust         = "TopLeft"
  gsn_panel(wks3,plot,(/3,1/),resP)               ; now draw as one plot
  frame(wks3)
  
; -- final   TURNOVER TIME plots --
  plot(0) = gsn_csm_contour_map(wks3,cas_TAU_f,res)
  plot(1) = gsn_csm_contour_map(wks3,mim_TAU_f,res)
  plot(2) = gsn_csm_contour_map(wks3,cor_TAU_f,res)

  resP@gsnPanelMainString = "Final Turnover time (y)"
  gsn_panel(wks3,plot,(/3,1/),resP)               ; now draw as one plot
  frame(wks3)


  delete([/resP@gsnPanelFigureStrings,res@cnLevels/])
; print("-- turnover")  
; -- Initial litter -- 
  res@cnMinLevelValF       =  0.001
  res@cnMaxLevelValF       =  200.
; res@cnLevels             = (/0.1,5,10,15,20,25,30,35,40,45,50,55,60,65/)
  res@cnLevels             = ((/0.1,1,2,3,5,10,20,30,50,100,200/))
  plot(0) = gsn_csm_contour_map(wks3,cas_LIT_i,res)
  plot(1) = gsn_csm_contour_map(wks3,mim_LIT_i,res)
; plot(2) = gsn_csm_contour_map(wks3,cor_LIT_i,res)

  resP@gsnPanelMainString = "Initail LIT, kgC/m2"
  resP@gsnPanelFigureStrings= (/"(a)","(b)","(c)"/)
; resP@gsnPanelFigureStrings= (/decimalPlaces(avg(Gcas_LIT(0:9)),1,True)+ " Pg C", \
;                               decimalPlaces(avg(Gmim_LIT(0:9)),1,True)+ " Pg C", \
;                               decimalPlaces(avg(Gcor_LIT(0:9)),1,True)+ " Pg C"  /)
  resP@gsnPanelFigureStringsPerimOn = False      ; turns off box around fig labels
  resP@gsnStringFont  = "helvetica-bold"
  resP@amJust         = "TopLeft" 
  gsn_panel(wks3,plot,(/3,1/),resP)               ; now draw as one plot
  frame(wks3)

; -- Initial soil C --
  plot(0) = gsn_csm_contour_map(wks3,cas_SOM_i,res)
  plot(1) = gsn_csm_contour_map(wks3,mim_SOM_i,res)
  plot(2) = gsn_csm_contour_map(wks3,cor_SOM_i,res)
  resP@gsnPanelMainString = "Initail Soil C, kgC/m2"
  gsn_panel(wks3,plot,(/3,1/),resP)               ; now draw as one plot
  frame(wks3)

; -- Initial protected soil C
  plot(0) = gsn_csm_contour_map(wks3,cas_SOM_pro_i,res)
  plot(1) = gsn_csm_contour_map(wks3,mim_SOM_pro_i,res)
  plot(2) = gsn_csm_contour_map(wks3,cor_SOM_pro_i,res)
  resP@gsnPanelMainString = "Initail protected Soil C, kgC/m2"
  gsn_panel(wks3,plot,(/3,1/),resP)               ; now draw as one plot
  frame(wks3)

; -- fraction protected soi C
cas_SOM_pro_frac = where(cas_SOM_i.gt. 0.1, cas_SOM_i, cas_SOM_i@_FillValue)
mim_SOM_pro_frac = where(mim_SOM_i.gt. 0.1, mim_SOM_i, mim_SOM_i@_FillValue)
cor_SOM_pro_frac = where(cor_SOM_i.gt. 0.1, cor_SOM_i, cor_SOM_i@_FillValue)
cas_SOM_pro_frac = cas_SOM_pro_i/cas_SOM_pro_frac*100
mim_SOM_pro_frac = mim_SOM_pro_i/mim_SOM_pro_frac*100
cor_SOM_pro_frac = cor_SOM_pro_i/cor_SOM_pro_frac*100
copy_VarCoords(cas_SOM_pro_i, cas_SOM_pro_frac)
copy_VarCoords(mim_SOM_pro_i, mim_SOM_pro_frac)
copy_VarCoords(cor_SOM_pro_i, cor_SOM_pro_frac)

  res@cnMinLevelValF       =  0.001
  res@cnMaxLevelValF       =  100.
  res@cnLevels             = ((/1.,5,10,20,30,40,50,60,70,80,90/))
  plot(0) = gsn_csm_contour_map(wks3,cas_SOM_pro_frac,res)
  plot(1) = gsn_csm_contour_map(wks3,mim_SOM_pro_frac,res)
  plot(2) = gsn_csm_contour_map(wks3,cor_SOM_pro_frac,res)
  resP@gsnPanelMainString = "Initail protected Soil C, %"
  gsn_panel(wks3,plot,(/3,1/),resP)               ; now draw as one plot
  frame(wks3)
  delete([/resP@gsnPanelFigureStrings, plot/])


; -- Total Soil C stocks --
; --        Fig_2        -- 
  plot = new(4, "graphic")
  res@cnMinLevelValF       =  0.001
  res@cnMaxLevelValF       =  200.
  res@cnLevels             = ((/0.1,1,2,3,5,10,20,30,50,100,200/))
  plot(0) = gsn_csm_contour_map(wks3,cas_TOT_i,res)  
  plot(1) = gsn_csm_contour_map(wks3,mim_TOT_i,res)
  plot(2) = gsn_csm_contour_map(wks3,cor_TOT_i,res)
  plot(3) = gsn_csm_contour_map(wks3,HWSD_SOC,res)
  resP@gsnPanelMainString = "Total Soil C, "+cas_TOT_i@units
  resP@gsnPanelFigureStrings= (/"(a)","(b)","(c)","(d)"/)
; resP@gsnPanelFigureStrings= (/decimalPlaces(avg(Gcas_TOT(0:9)),1,True)+ " Pg C", \
;                               decimalPlaces(avg(Gmim_TOT(0:9)),1,True)+ " Pg C", \
;                               decimalPlaces(avg(Gcor_TOT(0:9)),1,True)+ " Pg C", \
;                               decimalPlaces(gHWSD_SOC,1,True)         + " Pg C"  /)
  gsn_panel(wks3,plot,(/2,2/),resP)               ; now draw as one plot
  frame(wks3)

  res@cnFillMode = "RasterFill"
  plot(0) = gsn_csm_contour_map(wks3,cas_TOT_i,res)
  plot(1) = gsn_csm_contour_map(wks3,mim_TOT_i,res)
  plot(2) = gsn_csm_contour_map(wks3,cor_TOT_i,res)
  plot(3) = gsn_csm_contour_map(wks3,HWSD_SOC,res)
  gsn_panel(wks3,plot,(/2,2/),resP)               ; now draw as one plot
  frame(wks3)

  delete([/plot,resP@gsnPanelFigureStrings,res@cnFillMode/])

; -- make arctic only map, showing NCSCD too -- 
; ------ SI Fig. 3 ------
  Plot = new(6, "graphic")
  RES  = True
  RES@mpMinLatF  = 50                           ; minimum lat to plot
  RES@mpGridLatSpacingF    = 15                 ; change latitude  line spacing
  RES@mpGridLonSpacingF    = 90                 ; change longitude line spacing
  RES@mpFillOn   = False
  RES@gsnDraw    = False
  RES@gsnFrame   = False

  RES@cnFillOn             = True             ; color fill
  RES@cnFillMode           = "RasterFill"
  RES@cnLinesOn            = False            ; no contour lines
  RES@lbLabelBarOn         = False            ; turn off individual cb's
  RES@cnLevelSelectionMode = "ExplicitLevels"
  RES@cnMinLevelValF       =  0.001
  RES@cnMaxLevelValF       =  200.
  RES@cnLevels             = ((/0.1,1,2,5,10,20,50,100,200/))
  RES@gsnPolar             = "NH"            ; select the N. Hem.
  RES@gsnLeftString   = ""               ; turn off units and long name from plots
  RES@gsnRightString  = ""

  RES@gsnLeftString   = "CASA-CNP arctic total C"
  i_arc_cas_tot       = avg(Acas_TOT(0:9))
  RES@gsnRightString  = round(i_arc_cas_tot,1) + " Pg C"
  Plot(0) = gsn_csm_contour_map_polar(wks3,cas_TOT_i,RES)

; decimalPlaces(avg(Acas_TOT(ntim-10:ntim-1)) - avg(Acas_TOT(0:9)),1,True)


  RES@gsnLeftString   = "MIMICS arctic total C"
  i_arc_mim_tot       = avg(Amim_TOT(0:9))
  RES@gsnRightString  = round(i_arc_mim_tot,1) + " Pg C"
  Plot(2) = gsn_csm_contour_map_polar(wks3,mim_TOT_i,RES)

  RES@gsnLeftString   = "CORPSE   arctic total C"
  i_arc_cor_tot       = avg(Acor_TOT(0:9))
  RES@gsnRightString  = round(i_arc_cor_tot,1) + " Pg C"
  Plot(4) = gsn_csm_contour_map_polar(wks3,cor_TOT_i,RES)

  RES@gsnLeftString   = ""               ; turn off units and long name from plots
  RES@gsnRightString  = ""
  Plot(1) = gsn_csm_contour_map_polar(wks3,HWSD_SOC ,RES)
  Plot(3) = gsn_csm_contour_map_polar(wks3,NCSCD    ,RES)

  RESP     = True
  RESP@gsnMaximize                  = True        ; fill the page
  RESP@gsnPanelFigureStringsPerimOn = False       ; turns off box around fig labels
  RESP@gsnStringFont                = "helvetica-bold"
  RESP@amJust                       = "TopLeft"
  RESP@gsnPanelFigureStrings = (/"(a)","(d)","(b)","(e)","(c)","(f)"/)
  RESP@gsnPanelLabelBar      = True               ; add common colorbar
; RESP@gsnPanelMainString    = "Total Soil C, "+cas_TOT_i@units

  gsn_panel(wks3,Plot,(/3,2/),RESP)               ; now draw as one plot
  delete([/RES, RESP, Plot/])

  print("------------------------")
  print("wrote "+ fout+".ps")
  print("")

;-------------------------------------------------------------------
;-------------------------------------------------------------------
; plot the map of change in total soil carbon (last - first)
;-------------------------------------------------------------------
;-------------------------------------------------------------------
  ; ----------- average chagne in stocks --------------
  D_cas_GPP = dim_avg_n_Wrap(d_cas_GPP(ntim-10:ntim-1,:,:), 0)
  D_cas_NPP = dim_avg_n_Wrap(d_cas_NPP(ntim-10:ntim-1,:,:), 0)
  D_cas_HR = dim_avg_n_Wrap(d_cas_HR(ntim-10:ntim-1,:,:), 0)
  D_mim_HR = dim_avg_n_Wrap(d_mim_HR(ntim-10:ntim-1,:,:), 0)
  D_cor_HR = dim_avg_n_Wrap(d_cor_HR(ntim-10:ntim-1,:,:), 0)

  D_cas_LIT = dim_avg_n_Wrap(d_cas_LIT(ntim-10:ntim-1,:,:), 0)
  D_mim_LIT = dim_avg_n_Wrap(d_mim_LIT(ntim-10:ntim-1,:,:), 0)
  D_cor_LIT = dim_avg_n_Wrap(d_cor_LIT(ntim-10:ntim-1,:,:), 0)

  D_mim_MIC = dim_avg_n_Wrap(d_mim_MIC(ntim-10:ntim-1,:,:), 0)
  D_cor_MIC = dim_avg_n_Wrap(d_cor_MIC(ntim-10:ntim-1,:,:), 0)

  D_cas_SOM = dim_avg_n_Wrap(d_cas_SOM(ntim-10:ntim-1,:,:), 0)
  D_mim_SOM = dim_avg_n_Wrap(d_mim_SOM(ntim-10:ntim-1,:,:), 0)
  D_cor_SOM = dim_avg_n_Wrap(d_cor_SOM(ntim-10:ntim-1,:,:), 0)

  D_cas_SOM_pro = dim_avg_n_Wrap(d_cas_SOM_pro(ntim-10:ntim-1,:,:), 0)
  D_mim_SOM_pro = dim_avg_n_Wrap(d_mim_SOM_pro(ntim-10:ntim-1,:,:), 0)
  D_cor_SOM_pro = dim_avg_n_Wrap(d_cor_SOM_pro(ntim-10:ntim-1,:,:), 0)
  printMinMax(D_cas_SOM_pro, 1)

  D_cas_TOT = dim_avg_n_Wrap(d_cas_TOT(ntim-10:ntim-1,:,:), 0) 
  D_mim_TOT = dim_avg_n_Wrap(d_mim_TOT(ntim-10:ntim-1,:,:), 0)
  D_cor_TOT = dim_avg_n_Wrap(d_cor_TOT(ntim-10:ntim-1,:,:), 0)

  D_cas_TOT_2000 = dim_avg_n_Wrap(d_cas_TOT(ntim-100:ntim-91,:,:), 0)
  D_mim_TOT_2000 = dim_avg_n_Wrap(d_mim_TOT(ntim-100:ntim-91,:,:), 0)
  D_cor_TOT_2000 = dim_avg_n_Wrap(d_cor_TOT(ntim-100:ntim-91,:,:), 0)

  D_TSOI    = dim_avg_n_Wrap(d_TSOI(ntim-10:ntim-1,:,:), 0)
  D_Theta   = dim_avg_n_Wrap(d_Theta(ntim-10:ntim-1,:,:), 0)


  ; -------- add attribues ------------
  D_cas_GPP@units = "g C m~S~-2~N~ y~S~-1~N~"
  D_cas_NPP@units = "g C m~S~-2~N~ y~S~-1~N~"
  D_cas_HR@units  = "g C m~S~-2~N~ y~S~-1~N~"
  D_mim_HR@units  = "g C m~S~-2~N~ y~S~-1~N~"
  D_cor_HR@units  = "g C m~S~-2~N~ y~S~-1~N~"

  D_cas_LIT@units = "g C m~S~-2~N~"
  D_mim_LIT@units = "g C m~S~-2~N~"
  D_cor_LIT@units = "g C m~S~-2~N~"

  D_mim_MIC@units = "g C m~S~-2~N~"
  D_cor_MIC@units = "g C m~S~-2~N~"

  D_cas_SOM@units = "g C m~S~-2~N~"
  D_mim_SOM@units = "g C m~S~-2~N~"
  D_cor_SOM@units = "g C m~S~-2~N~"

  D_cas_SOM_pro@units = "g C m~S~-2~N~"
  D_mim_SOM_pro@units = "g C m~S~-2~N~"
  D_cor_SOM_pro@units = "g C m~S~-2~N~"

  D_cas_TOT@units = "g C m~S~-2~N~"
  D_mim_TOT@units = "g C m~S~-2~N~"
  D_cor_TOT@units = "g C m~S~-2~N~"

  D_cas_TOT_2000@units = "g C m~S~-2~N~"
  D_mim_TOT_2000@units = "g C m~S~-2~N~"
  D_cor_TOT_2000@units = "g C m~S~-2~N~"

  D_TSOI@units  = "~S~o~N~C"
  D_Theta@units = Theta@units

  D_cas_GPP@long_name = "change "+ cas_GPP@long_name
  D_cas_NPP@long_name = "change "+ cas_NPP@long_name
  D_cas_HR@long_name  = "change "+ cas_HR@long_name
  D_mim_HR@long_name  = "change "+ mim_HR@long_name
  D_cor_HR@long_name  = "change "+ cor_HR@long_name
  
  D_cas_LIT@long_name = "change "+ cas_LIT@long_name
  D_mim_LIT@long_name = "change "+ mim_LIT@long_name
  D_cor_LIT@long_name = "change "+ cor_LIT@long_name

  D_mim_MIC@long_name = "change "+ mim_MIC@long_name
  D_cor_MIC@long_name = "change "+ cor_MIC@long_name

  D_cas_SOM@long_name = "change "+ cas_SOM@long_name
  D_mim_SOM@long_name = "change "+ mim_SOM@long_name
  D_cor_SOM@long_name = "change "+ cor_SOM@long_name

  D_cas_SOM_pro@long_name = "change "+ cas_SOM_pro@long_name
  D_mim_SOM_pro@long_name = "change "+ mim_SOM_pro@long_name
  D_cor_SOM_pro@long_name = "change "+ cor_SOM_pro@long_name

  D_cas_TOT@long_name = "change "+ cas_TOT@long_name
  D_mim_TOT@long_name = "change "+ mim_TOT@long_name
  D_cor_TOT@long_name = "change "+ cor_TOT@long_name

  D_cas_TOT_2000@long_name = "Y2K change "+ cas_TOT@long_name
  D_mim_TOT_2000@long_name = "Y2K change "+ mim_TOT@long_name
  D_cor_TOT_2000@long_name = "Y2K change "+ cor_TOT@long_name

  D_TSOI@long_name = "change "+ TSOI@long_name
  D_Theta@long_name = "change "+ Theta@long_name

; printVarSummary(D_cas_TOT)
;-- Write out .nc file & maps --

   cdffile = out_dir+ "init_final_stocks.nc"
   system("rm -f " +  cdffile)
   a = addfile(cdffile,"c")
   a@source = systemfunc("pwd")+"/"+get_script_name()

   a->cas_TOT   = cas_TOT_i
   a->mim_TOT   = mim_TOT_i
   a->cor_TOT   = doubletofloat(cor_TOT_i)
   a->D_cas_TOT   = D_cas_TOT
   a->D_mim_TOT   = D_mim_TOT
   a->D_cor_TOT   = doubletofloat(D_cor_TOT)
   a->D_cas_TOT_2000   = D_cas_TOT_2000
   a->D_mim_TOT_2000   = D_mim_TOT_2000
   a->D_cor_TOT_2000   = doubletofloat(D_cor_TOT_2000)

  print(" Wrote " + cdffile)
  print("")
  print("------------------------")

  ; ------ plot figures ----------------
  fout  = out_dir+"TestbedChange_SOMplots_"+exper
  wks2  = gsn_open_wks("ps",fout)
  gsn_define_colormap(wks2,"CBR_drywet")
  plot = new(6, "graphic")
;  print(" starting change SOM plots ")
; -- change in GPP,NPP --------  
  res@lbLabelBarOn        = True            ; turn off individual cb's
  res@cnMinLevelValF       =  -400. 
  res@cnMaxLevelValF       =  400.  
  res@cnLevels             = fspan(-400.,400.,11)

  plot(0) = gsn_csm_contour_map(wks2,D_cas_GPP,res)
  plot(1) = gsn_csm_contour_map(wks2,D_cas_NPP,res)

  resP@gsnPanelLabelBar = False              ; add common colorbar
  resP@gsnPanelMainString = "change GPP, NPP "+D_cas_GPP@units
  resP@gsnPanelFigureStrings= (/"(a)","(b)","","","",""/)
  gsn_panel(wks2,plot,(/3,2/),resP)               ; now draw as one plot
  frame(wks2)

  delete([/plot,res@cnMinLevelValF,resP@gsnPanelFigureStrings, \
                res@cnMaxLevelValF,res@cnLevels/])
; print(" d GPP, NPP")
; -- change in TSOI, Theta ----
  plot = new(6, "graphic")
  gsn_define_colormap(wks2,"NCV_blu_red")

  res@lbLabelBarOn        = True            ; turn off individual cb's
  res@cnMinLevelValF       =  -5.
  res@cnMaxLevelValF       =  5.
  res@cnLevels             = fspan(-6.,6.,13)
  plot(2) = gsn_csm_contour_map(wks2,D_TSOI,res)

  resP@gsnPanelMainString = "change TSOI, Theta"
  resP@gsnPanelLabelBar = False              ; add common colorbar
  resP@gsnPanelFigureStrings= (/"","","(c)","","(e)",""/)
  gsn_panel(wks2,plot,(/3,2/),resP)               ; now draw as one plot
  frame(wks2)

  gsn_define_colormap(wks2,"CBR_drywet")
  res@cnMinLevelValF       =  -13.
  res@cnMaxLevelValF       =  13.
  delete(res@cnLevels)
  res@cnLevels             = fspan(-10. ,10. ,21)
  plot(4) = gsn_csm_contour_map(wks2,D_Theta,res)

  gsn_panel(wks2,plot,(/3,2/),resP)               ; now draw as one plot
  frame(wks2)
; print(" d TSOI, THETA ")
  delete([/res@cnLevels,resP@gsnPanelFigureStrings,plot/])


; -- change in HR --------
; commented out to speed up plotting
  plot = new(3, "graphic")
  gsn_define_colormap(wks2,"CBR_drywet")

  res@lbLabelBarOn         = False           ; turn off individual cb's
  res@cnMinLevelValF       =  -200.
  res@cnMaxLevelValF       =  200.
  res@cnLevels             = fspan(-200.,200.,11)
; plot(0) = gsn_csm_contour_map(wks2,D_cas_HR,res)
; plot(1) = gsn_csm_contour_map(wks2,D_mim_HR,res)
; plot(2) = gsn_csm_contour_map(wks2,D_cor_HR,res)
  resP@gsnPanelMainString = "change HR "+D_cas_HR@units
  resP@gsnPanelFigureStrings= (/"(a)","(b)","(c)"/)
  resP@gsnPanelLabelBar = True             ; add common colorbar
; resP@gsnPanelFigureStrings= (/decimalPlaces(avg(d_Gcas_HR(ntim-10:ntim-1)),1,True)+ " Pg C", \
;                               decimalPlaces(avg(d_Gmim_HR(ntim-10:ntim-1)),1,True)+ " Pg C", \
;                               decimalPlaces(avg(d_Gcor_HR(ntim-10:ntim-1)),1,True)+ " Pg C"  /)
; gsn_panel(wks2,plot,(/3,1/),resP)               ; now draw as one plot
; frame(wks2)

; -- change in Litter pools  --------  
  res@cnMinLevelValF       =  -1500.
  res@cnMaxLevelValF       =  1500.
  res@cnLevels             = fspan(-1500,1500,11  )
  plot(0) = gsn_csm_contour_map(wks2,D_cas_LIT,res)
  plot(1) = gsn_csm_contour_map(wks2,D_mim_LIT,res)
; plot(2) = gsn_csm_contour_map(wks2,D_cor_LIT,res)
  resP@gsnPanelMainString = "change Litter "  +D_cas_LIT@units
; resP@gsnPanelFigureStrings= (/decimalPlaces(avg(d_Gcas_LIT(ntim-10:ntim-1)),1,True)+ " Pg C", \
;                               decimalPlaces(avg(d_Gmim_LIT(ntim-10:ntim-1)),1,True)+ " Pg C", \
;                               decimalPlaces(avg(d_Gcor_LIT(ntim-10:ntim-1)),1,True)+ " Pg C"  /)
  gsn_panel(wks2,plot,(/3,1/),resP)               ; now draw as one plot
  frame(wks2)
; print(" d litter"      )

; -- change in SOM pools  --------
  plot(0) = gsn_csm_contour_map(wks2,D_cas_SOM,res)
  plot(1) = gsn_csm_contour_map(wks2,D_mim_SOM,res)
  plot(2) = gsn_csm_contour_map(wks2,D_cor_SOM,res)
  resP@gsnPanelMainString = "change soil C "  +D_cas_SOM@units
; resP@gsnPanelFigureStrings= (/decimalPlaces(avg(d_Gcas_SOM(ntim-10:ntim-1)),1,True)+ " Pg C", \
;                               decimalPlaces(avg(d_Gmim_SOM(ntim-10:ntim-1)),1,True)+ " Pg C", \
;                               decimalPlaces(avg(d_Gcor_SOM(ntim-10:ntim-1)),1,True)+ " Pg C"  /)
  gsn_panel(wks2,plot,(/3,1/),resP)               ; now draw as one plot
  frame(wks2)
; print(" d soil C "      )

; -- change in SOM_pro pools  --------
  res@cnMinLevelValF       =  -1000.
  res@cnMaxLevelValF       =  1000.
  delete(res@cnLevels)
  res@cnMinLevelValF       = -1.5e3
  res@cnMaxLevelValF       =  1.5e3
  res@cnLevels = (/-1.5e3,-1e3,-750,-500,-250,-100,-50,0,50,100,250,500,750,1e3,1.5e3/)

  plot(0) = gsn_csm_contour_map(wks2,D_cas_SOM_pro,res)
  plot(1) = gsn_csm_contour_map(wks2,D_mim_SOM_pro,res)
  plot(2) = gsn_csm_contour_map(wks2,D_cor_SOM,res)
  resP@gsnPanelMainString = "change soil C "  +D_cas_SOM@units
; resP@gsnPanelFigureStrings= (/decimalPlaces(avg(d_Gcas_SOM(ntim-10:ntim-1)),1,True)+ " Pg C", \
;                               decimalPlaces(avg(d_Gmim_SOM(ntim-10:ntim-1)),1,True)+ " Pg C", \
;                               decimalPlaces(avg(d_Gcor_SOM(ntim-10:ntim-1)),1,True)+ " Pg C"  /)
  resP@gsnPanelMainString = "change protected soil C "  +D_cas_SOM_pro@units
  gsn_panel(wks2,plot,(/3,1/),resP)               ; now draw as one plot
  frame(wks2)
  delete(res@cnLevels)
; print(" d protected soil C "      )



; -- change in total soil C pools ---------------
; --   Fig_6    ---
  res@cnMinLevelValF       = (-2.5e3*2)
  res@cnMaxLevelValF       = (2.5e3*2)
  res@cnLevels = ((/-2.5e3,-1.5e3,-1e3,-750,-500,-100,0,100,500,750,1e3,1.5e3,2.5e3/)*2)

  plot(0) = gsn_csm_contour_map(wks2,D_cas_TOT,res)
  plot(1) = gsn_csm_contour_map(wks2,D_mim_TOT,res)
  plot(2) = gsn_csm_contour_map(wks2,D_cor_TOT,res)
  resP@gsnPanelMainString = "change total C "  +D_cas_TOT@units
  gsn_panel(wks2,plot,(/3,1/),resP)               ; now draw as one plot
  frame(wks2)

  res@cnFillMode = "RasterFill"
  plot(0) = gsn_csm_contour_map(wks2,D_cas_TOT,res)
  plot(1) = gsn_csm_contour_map(wks2,D_mim_TOT,res)
  plot(2) = gsn_csm_contour_map(wks2,D_cor_TOT,res)
  gsn_panel(wks2,plot,(/3,1/),resP)               ; now draw as one plot
  frame(wks2)

  delete(res@cnLevels)
  resP@gsnPanelMainString = "Y2K Change total C"  +D_cas_TOT@units
  res@cnFillMode = "RasterFill"
  res@cnMinLevelValF       = (-1.0e3  )
  res@cnMaxLevelValF       = (1.0e3  )
  res@cnLevels = ((/-800. ,-500. ,-300,-150,-50.,0,50, 150, 300, 500, 800/))

  plot(0) = gsn_csm_contour_map(wks2,D_cas_TOT_2000,res)
  plot(1) = gsn_csm_contour_map(wks2,D_mim_TOT_2000,res)
  plot(2) = gsn_csm_contour_map(wks2,D_cor_TOT_2000,res)
  gsn_panel(wks2,plot,(/3,1/),resP)               ; now draw as one plot
  frame(wks2)

  delete(res@cnLevels)
  res@cnMinLevelValF       = (-2.5e3*2)
  res@cnMaxLevelValF       = (2.5e3*2)
  res@cnLevels = ((/-2.5e3,-1.5e3,-1e3,-750,-500,-100,0,100,500,750,1e3,1.5e3,2.5e3/)*2)

; print(" d total C "      )

; -- changes in arctic only (> 50N, roughly corresponding to domain in McGuire 2016, 2018 papers) --
; ------ SI Fig. 3 ------
  Plot = new(6, "graphic")
  RES  = True
  RES@mpMinLatF  = 50                           ; minimum lat to plot
  RES@mpGridLatSpacingF    = 15                 ; change latitude  line spacing
  RES@mpGridLonSpacingF    = 90                 ; change longitude line spacing
  RES@gsnPolar             = "NH"            ; select the N. Hem.
  RES@mpFillOn   = False
  RES@gsnDraw    = False
  RES@gsnFrame   = False

  RES@cnFillOn             = True             ; color fill
  RES@cnFillMode           = "RasterFill"
  RES@cnLinesOn            = False            ; no contour lines
  RES@lbLabelBarOn         = True             ; turn off individual cb's

  gsn_define_colormap(wks2,"NCV_blu_red")
  RES@cnLevelSelectionMode = "ExplicitLevels"
  RES@cnMinLevelValF       =  -5.
  RES@cnMaxLevelValF       =  5.
  RES@cnLevels             = fspan(-6.,6.,13)
  RES@gsnLeftString   = "TSOI change "
  d_arc_TSOI          = avg(ATSOI(ntim-10:ntim-1)) - avg(ATSOI(0:9))
  RES@gsnRightString  = decimalPlaces(d_arc_TSOI, 1,True) + " ~S~o~N~C"
  print("TSOI arctic    change "+RES@gsnRightString)
  Plot(1) = gsn_csm_contour_map(wks2,D_TSOI,RES)
  Plot(0) = gsn_csm_contour_map(wks2,D_TSOI,RES)
  Plot(2) = gsn_csm_contour_map(wks2,D_TSOI,RES)
  Plot(3) = gsn_csm_contour_map(wks2,D_TSOI,RES)
  Plot(4) = gsn_csm_contour_map(wks2,D_TSOI,RES)
  Plot(5) = gsn_csm_contour_map(wks2,D_TSOI,RES)
  RESP     = True
  RESP@gsnMaximize                  = True        ; fill the page
; RESP@gsnPanelFigureStringsPerimOn = False       ; turns off box around fig labels
; RESP@gsnStringFont                = "helvetica-bold"
; RESP@amJust                       = "TopLeft"
; RESP@gsnPanelFigureStrings = (/"(a)","(d)","(b)","(e)","(c)","(f)"/)
  RESP@gsnPanelLabelBar      = False              ; add common colorbar
  RESP@gsnPanelMainString    = "Arctic changes, "+exper

  gsn_panel(wks2,Plot,(/3,2/),RESP)               ; now draw as one plot

  ; add rest of plots w/ right color bar
  delete(RES@cnLevels)
  gsn_define_colormap(wks2,"CBR_drywet")
  RES@cnMinLevelValF       =  res@cnMinLevelValF 
  RES@cnMaxLevelValF       =  res@cnMaxLevelValF
  RES@cnLevels             = res@cnLevels                    

  ;dim_sum_n not working for some reason?     
  RES@gsnLeftString   = "CASA-CNP change total C"  
  d_arc_cas_tot       = avg(Acas_TOT(ntim-10:ntim-1)) - avg(Acas_TOT(0:9)) 
  RES@gsnRightString  = round(d_arc_cas_tot,3) + " Pg C" 
  print("CASA arctic TOT change " + RES@gsnRightString)
  Plot(0) = gsn_csm_contour_map_polar(wks2,D_cas_TOT,RES)

  RES@gsnLeftString   = "MIMICS change total C"
  d_arc_mim_tot       = avg(Amim_TOT(ntim-10:ntim-1)) - avg(Amim_TOT(0:9)) 
  RES@gsnRightString  = round(d_arc_mim_tot,3) + " Pg C" 
  print("MIM  arctic TOT change " + RES@gsnRightString)
  Plot(2) = gsn_csm_contour_map_polar(wks2,D_mim_TOT,RES)

  RES@gsnLeftString   = "CORPSE change total C"
  d_arc_cor_tot       = avg(Acor_TOT(ntim-10:ntim-1)) - avg(Acor_TOT(0:9)) 
  RES@gsnRightString  = round(d_arc_cor_tot,3) + " Pg C" 
  print("COR  arctic TOT change " + RES@gsnRightString)
  Plot(4) = gsn_csm_contour_map_polar(wks2,D_cor_TOT,RES)

  delete(RES@cnLevels)
  RES@cnMinLevelValF       =  -13.
  RES@cnMaxLevelValF       =  13.
  RES@cnLevels             = fspan(-10. ,10. ,21)
  RES@gsnLeftString   = "theta change "
  d_arc_Theta          = avg(ATheta(ntim-10:ntim-1)) - avg(ATheta(0:9))
  RES@gsnRightString  = decimalPlaces(d_arc_Theta, 1,True) + " "
  print("Theta arctic   change "+RES@gsnRightString)
  Plot(3) = gsn_csm_contour_map(wks2,D_Theta,RES)

  ; plot changes in NPP
  delete(RES@cnLevels)
  RES@cnMinLevelValF       =  -400. 
  RES@cnMaxLevelValF       =  400.  
  RES@cnLevels             = fspan(-400.,400.,11)
  RES@gsnLeftString   = "NPP change" 
  d_arc_npp_tot       = avg(Acas_NPP(ntim-10:ntim-1)) - avg(Acas_NPP(0:9))
  RES@gsnRightString  = round(d_arc_npp_tot,3) + " Pg C y^-1"
  print("NPP  arctic TOT change " + RES@gsnRightString)
  Plot(5) = gsn_csm_contour_map(wks2,D_cas_NPP,RES)

  gsn_panel(wks2,Plot,(/3,2/),RESP)               ; now draw as one plot

  delete([/RES, RESP, Plot/])
  delete([/plot, resP@gsnPanelFigureStrings,res@cnFillMode,res@cnLevels/])


; -- change in NPP & MIC     --------
  plot = new(4, "graphic")
  res@lbLabelBarOn        = True            ; turn off individual cb's
  res@cnMinLevelValF      =  -400.
  res@cnMaxLevelValF      =  400.
  res@cnLevels            = fspan(-400.,400.,11)

; plot(0) = gsn_csm_contour_map(wks2,D_cas_GPP,res)
  plot(1) = gsn_csm_contour_map(wks2,D_cas_NPP,res)

  resP@gsnPanelMainString = "abs change NPP, mimMic, corMic"
  resP@gsnPanelFigureStrings= (/"(a)","(b)","(c)","(d)"/)
  resP@gsnPanelLabelBar = False              ; add common colorbar

  res@cnMinLevelValF       =  -50.
  res@cnMaxLevelValF       =  50.
  res@cnLevels             = fspan(-50, 50, 11  )
  plot(2) = gsn_csm_contour_map(wks2,D_mim_MIC,res)
  plot(3) = gsn_csm_contour_map(wks2,D_cor_MIC,res)
  gsn_panel(wks2,plot,(/2,2/),resP)               ; now draw as one plot
  frame(wks2)

  plot(1) = gsn_csm_contour_map(wks2,pd_cas_NPP,res)
  plot(2) = gsn_csm_contour_map(wks2,pd_mim_MIC,res)
  plot(3) = gsn_csm_contour_map(wks2,pd_cor_MIC,res)
  resP@gsnPanelMainString = "% change NPP, mimMic, corMic"      
  gsn_panel(wks2,plot,(/2,2/),resP)               ; now draw as one plot
  frame(wks2)


  delete([/res,resP,wks2,plot/])
  print("")
  print("wrote "+ fout+".ps")
  print("")


;------------------------------------------------------
; generate plot of change in SOM, averaged by latitude
;    Make Fig_XX
;-----------------------------------------------------
  fout  = out_dir+"dStocks_zonal"+exper
  wks   = gsn_open_wks("ps" ,fout)                      ; open ps file
; gsn_define_colormap(wks,"CBR_drywet")

; plot  = new(2, "graphic")
  res                      = True            ; plot mods desired
  res@gsnMaximize          = False
  res@gsnDraw              = False
  res@gsnFrame             = False
  res@cnSmoothingOn        = False
  res@cnFillOn             = True            ; turn on color
  res@cnInfoLabelOn        = False           ; turn on cn  info label
  res@cnLinesOn            = False           ; no contour lines
  res@cnLineLabelsOn       = False           ; turn off contour labels
  res@vpXF                 = 0.12            ; default is 0.2 (aspect ratio)
  res@vpYF                 = 0.9             ; default is 0.8
  res@gsnLeftString        = ""
  res@gsnRightString       = ""
  res@gsnCenterString      = ""
  res@lbLabelFontHeightF   =.03                    ; make labels larger
  res@tmEqualizeXYSizes    = "True"                ; X & Y labels equal
  res@tmXBLabelFontHeightF = 0.03                  ; default 0.02 for a viewport width of 0.6
  res@tiXAxisString        = ""

  res@trYMinF              = -50                   ; set min and max lat
  res@trYMaxF              = 80 
; res@cnLevelSelectionMode = "ManualLevels"        ; manual levels
; res@cnMinLevelValF       = (-1000.*3)
; res@cnMaxLevelValF       = (1000.*3)

  res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
; irregular color bar
  res@cnLevels = (/-5e3 ,-3e3  ,-2e3,-1e3,-200,0,200,500,750,1.e3,1.5e3/)
; more reglular color bar
; res@cnLevels = (/-3e3 ,-1.5e3,-800,-500,-200,0,200,500,800,1.5e3,3e3/)
  cmap       = read_colormap_file("CBR_drywet")
  cmap(5,:)  = 0.

  res@cnFillPalette = cmap
  res@lbOrientation        = "vertical"      ; vertical label bar
  res@lbLabelBarOn         = False           ; turn off individual cb's
  plot1b = gsn_csm_lat_time(wks, d_cas_TOT_zone,res )
  res@lbLabelBarOn         = True            ; turn off individual cb's
  plot3b = gsn_csm_lat_time(wks, d_mim_TOT_zone,res )
  res@lbLabelBarOn         = False           ; turn off individual cb's
  plot5b = gsn_csm_lat_time(wks, d_cor_TOT_zone,res )

  text1b= gsn_add_text(wks,plot1b,"(a)",1915,70,txres)
  text3b= gsn_add_text(wks,plot3b,"(b)",1915,70,txres)
  text5b= gsn_add_text(wks,plot5b,"(c)",1915,70,txres)

  gres                  = True                ; polyline mods desired
  gres@gsLineThicknessF = 1.0                 ; line thickness
  gres@gsLineDashPattern= 1 
  line1b = gsn_add_polyline(wks, plot1b, (/1900,2100/),(/0,0/), gres)   
  line3b = gsn_add_polyline(wks, plot3b, (/1900,2100/),(/0,0/), gres)         
  line5b = gsn_add_polyline(wks, plot5b, (/1900,2100/),(/0,0/), gres)         

              
  attachid1  = gsn_attach_plots(plot1b,(/plot3b,plot5b/),attachres1,attachres2)
  plot       = plot1b
  resP                  = True               ; modify the panel plot
; resP@gsnPanelLabelBar = True               ; add common colorbar
; resP@lbOrientation    = "vertical"         ; vertical label bar
; resP@pmLabelBarParallelPosF   = 0.17       ; move right?default 0.5
  gsn_panel(wks,plot,(/1,1/),resP)           ; now draw as one plot
  frame(wks)

  delete([/res,resP,wks/])

  print("------------------------")
  print("wrote "+ fout+".ps")
  print("")
  ; ---------------- plot zonal means --------------------------

  zTOT         = new((/4,nlat/),float)
  zTOT(0,:)    = dim_avg(cas_TOT_i)
  zTOT(1,:)    = dim_avg(mim_TOT_i)
  zTOT(2,:)    = tofloat(dim_avg(cor_TOT_i))
  zTOT(3,:)    = 0.

  zTAU         = new((/4,nlat/),float)
  zTAU(0,:)    = dim_avg(cas_TAU_i)
  zTAU(1,:)    = dim_avg(mim_TAU_i)
  zTAU(2,:)    = tofloat(dim_avg(cor_TAU_i))
  zTAU(3,:)    = 0.

  zTAU_f         = new((/4,nlat/),float)
  zTAU_f(0,:)    = dim_avg(cas_TAU_f)
  zTAU_f(1,:)    = dim_avg(mim_TAU_f)
  zTAU_f(2,:)    = tofloat(dim_avg(cor_TAU_f))
  zTAU_f(3,:)    = 0.


  ; look at zonal mean and variance of HR over time
  zHR_mean        = new((/3,ntim,nlat/), float)
  zHR_mean!0      = "model"
  zHR_mean!1      = "time"
  zHR_mean!2      = "lat"
  zHR_mean&time   = years
  zHR_mean&lat    = lat
  zHR_runStdev    = zHR_mean
  zHR_runMean     = zHR_mean
  zHR_runStdev@long_name = "HR, 20 y running stdev"
  zHR_runMean@long_name  = "HR, 20 y running mean"
 
  zNPP_mean       = new((/ntim,nlat/), float)
  zNPP_mean!0     = "time"
  zNPP_mean!1     = "lat"
  zNPP_mean&time  = years
  zNPP_mean&lat   = lat

  zGPP_mean       = zNPP_mean
  zTAIR_mean      = zNPP_mean
  zNPP_runStdev   = zNPP_mean
  zGPP_runStdev   = zNPP_mean
  zTAIR_runStdev  = zNPP_mean
  zNPP_runMean    = zNPP_mean
  zGPP_runMean    = zNPP_mean
  zTAIR_runMean   = zNPP_mean

  zNPP_runStdev@lonn_name  = "NPP, 20 y running stdev"
  zGPP_runStdev@lonn_name  = "GPP, 20 y running stdev"
  zTAIR_runStdev@long_name = "TAIR, 20 y running stdev"
  zNPP_runMean@lonn_name   = "NPP, 20 y running mean"
  zGPP_runMean@lonn_name   = "GPP, 20 y running mean"    
  zTAIR_runMean@long_name  = "TAIR, 20 y running mean"    

  zHR_mean(0,:,:) = dim_avg_n(cas_HR, 2)
  zHR_mean(1,:,:) = dim_avg_n(mim_HR, 2)
  zHR_mean(2,:,:) = tofloat(dim_avg_n(cor_HR, 2))
  zNPP_mean(:,:)  = dim_avg_n(cas_NPP, 2)
  zGPP_mean(:,:)  = dim_avg_n(cas_GPP, 2)
  zTAIR_mean(:,:) = dim_avg_n(tair+273.15, 2) ; convert to K


  ; cacluate CV (NOT Variance) to standardize against mean
  do i= 10,ntim-10
     zHR_runStdev(:,i,:) = dim_stddev_n(zHR_mean(:,i-10:i+9,:),1)  
     zHR_runMean(:,i,:)  = dim_avg_n(zHR_mean(:,i-10:i+9,:),1)
     zNPP_runStdev(i,:)  = dim_stddev_n(zNPP_mean(i-10:i+9,:), 0)
     zNPP_runMean(i,:)   = dim_avg_n(zNPP_mean(i-10:i+9,:), 0)
     zGPP_runStdev(i,:)  = dim_stddev_n(zGPP_mean(i-10:i+9,:), 0)
     zGPP_runMean(i,:)   = dim_avg_n(zGPP_mean(i-10:i+9,:), 0)
     zTAIR_runStdev(i,:) = dim_stddev_n(zTAIR_mean(i-10:i+9,:),0) 
     zTAIR_runMean(i,:)  = dim_avg_n(zTAIR_mean(i-10:i+9,:),0)
  end do

  zHR_runMean   = where(zHR_runMean  .ge.1., zHR_runMean,   zHR_runMean@_FillValue) ;avoid high CV with low productivity
  zNPP_runMean  = where(zNPP_runMean .ge.1., zNPP_runMean,  zNPP_runMean@_FillValue)
  zGPP_runMean  = where(zGPP_runMean .ge.1., zGPP_runMean,  zGPP_runMean@_FillValue)
  zTAIR_runMean = where(zTAIR_runMean.ne.0., zTAIR_runMean, zTAIR_runMean@_FillValue)

  zHR_runcv   = zHR_runStdev  / zHR_runMean
  zNPP_runcv  = zNPP_runStdev / zNPP_runMean
  zGPP_runcv  = zGPP_runStdev / zGPP_runMean
  zTAIR_runcv = zTAIR_runStdev/ zTAIR_runMean

  copy_VarCoords(zHR_runMean,   zHR_runcv  )
  copy_VarCoords(zNPP_runMean,  zNPP_runcv )
  copy_VarCoords(zGPP_runMean,  zGPP_runcv )
  copy_VarCoords(zTAIR_runMean, zTAIR_runcv)

  zHR_runcv@long_name   = "HR, 20 y running CV  "
  zTAIR_runcv@long_name = "TAIR, 20 y running CV  "
  zNPP_runcv@long_name  = "NPP, 20 y running CV  "
  zGPP_runcv@long_name  = "GPP, 20 y running CV  "


;  printVarSummary(zHR_runcv)

  zOBS         = dim_avg(HWSD_SOC)
  nlat2        = dimsizes(HWSD_SOC&lat)          
  zNCSCD       = dim_avg(NCSCD)
  nlatNCSCD    = dimsizes(NCSCD&lat)

;----------------------------------------------------------
; define a polygon centered with width of 1 sigma from obs
;-----------------------------------------------------------
 xp    = new( (/2*nlat2/), double)
 yp    = new( (/2*nlat2/), double)
 dx    = dim_stddev(HWSD_SOC)

 do k=0,nlat2-1
   xp(k)           = zOBS(k) + dx(k)
   yp(k)           = HWSD_SOC&lat(k)
   yp(2*nlat2-1-k) = HWSD_SOC&lat(k)
   xp(2*nlat2-1-k) = zOBS(k) - dx(k)
 end do

  zdTOT         = new((/4,nlat/),float)
  zdTOT(0,:)    = 0.
  zdTOT(1,:)    = dim_avg(D_cas_TOT)
  zdTOT(2,:)    = dim_avg(D_mim_TOT)
  zdTOT(3,:)    = tofloat(dim_avg(D_cor_TOT))

  ; mask out zeros for normalized chagne in C stocks, convert to gC/m2 (from kgC)
  cas_i     = where(cas_TOT_i .le. 0.1, cas_TOT_i@_FillValue, cas_TOT_i*1e3)
  mim_i     = where(mim_TOT_i .le. 0.1, mim_TOT_i@_FillValue, mim_TOT_i*1e3)
  cor_i     = where(cor_TOT_i .le. 0.1, cor_TOT_i@_FillValue, cor_TOT_i*1e3)
  gpp_i     = where(cas_GPP_i .le. 0.1, cas_GPP_i@_FillValue, cas_GPP_i)
  npp_i     = where(cas_NPP_i .le. 0.1, cas_NPP_i@_FillValue, cas_NPP_i)
  gpp_d     = where(D_cas_GPP .le. 10., D_cas_GPP@_FillValue, D_cas_GPP)
  npp_d     = where(D_cas_NPP .le. 10., D_cas_NPP@_FillValue, D_cas_NPP)
  
  zdTOT_norm         = new((/4,nlat/),float)
  zdTOT_norm(0,:)    = 0.
; normalized to initial stocks
  zdTOT_norm(1,:)    = dim_avg(D_cas_TOT/cas_i) * 100 
  zdTOT_norm(2,:)    = dim_avg(D_mim_TOT/mim_i) * 100
  zdTOT_norm(3,:)    = tofloat(dim_avg(D_cor_TOT/cor_i)) * 100
  zdGPP_norm         = dim_avg(D_cas_GPP / gpp_i) * 100
  zdNPP_norm         = dim_avg(D_cas_NPP / npp_i) * 100

; normalized to change in npp
; zdTOT_norm(1,:)    = dim_avg(D_cas_TOT/npp_d) 
; zdTOT_norm(2,:)    = dim_avg(D_mim_TOT/npp_d) 
; zdTOT_norm(3,:)    = tofloat(dim_avg(D_cor_TOT/npp_d)) 


  print("normalized D_cas_TOT = "+sum((D_cas_TOT/cas_i)* landarea)/ sum(landarea) )
  print("normalized D_mim_TOT = "+sum((D_mim_TOT/mim_i)* landarea)/ sum(landarea) )
  print("normalized D_cor_TOT = "+sum((D_cor_TOT/cor_i)* landarea)/ sum(landarea) )
  print("--")
  print("normalized D_cor_GPP = "+sum((D_cas_GPP / gpp_i)* landarea)/ sum(landarea) )
  print("normalized D_cor_NPP = "+sum((D_cas_NPP / npp_i)* landarea)/ sum(landarea) )

  fout  = out_dir+"Testbed_zonal_plots_"+exper
  wks  = gsn_open_wks("ps" ,fout)

; --- plot zonal mean of initial stocks ---
  res                        = True
  res@gsnFrame               = False
  res@gsnDraw                = False
  res@tiYAxisString          = "Latitude"
  ymin                       = -50
  ymax                       = 83
  xmin                       = 0
  xmax                       = 100
  res@tiXAxisString          = "Initial stocks ("+cas_TOT_i@units+")"
;---Make sure NCL doesn't pick different min/max for axes 
  res@trXMinF = xmin 
  res@trXMaxF = xmax 
  res@trYMinF = ymin 
  res@trYMaxF = ymax 

  res@xyLineThicknessF = 5.0
  res@xyLineColor      = colors(1)
  plot1                = gsn_csm_xy(wks,zTOT(0,:),D_cas_TOT&lat, res)
  res@xyLineColor      = colors(2)
  plot2                = gsn_csm_xy(wks,zTOT(1,:),D_cas_TOT&lat, res)
  res@xyLineColor      = colors(3)
  plot3                = gsn_csm_xy(wks,zTOT(2,:),D_cas_TOT&lat, res)
  res@xyLineColor      = colors(0)
  plot0                = gsn_csm_xy(wks,zTOT(3,:),D_cas_TOT&lat, res)
  res@xyLineThicknessF = 5.0
  res@xyLineColor      = colors(0)
  plot4                = gsn_csm_xy(wks,zOBS,HWSD_SOC&lat,  res)
  res@pmLegendOrthogonalPosF = -0.685              ; move units down
  
  gsres                   = True                        ; poly res
  gsres@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres@gsFillColor       = "Grey"                 ; color chosen
  dummy = gsn_add_polygon (wks,plot4,xp,yp,gsres)
 
  res@xyLineThicknessF = 3.0
  res@xyLineColor      = colors(0)
  res@xyDashPattern    = 1
  plot5                = gsn_csm_xy(wks,zNCSCD,NCSCD&lat,  res)

  overlay(plot0, plot1)
  overlay(plot0, plot2)
  overlay(plot0, plot3)
  overlay(plot0, plot4)
;try including NCSCD too?
  overlay(plot0, plot5)

; delete(x)
  x         = (/0,5,10,15,20,25,30,60,85,110/)
 irreg_axis = create "irregular" irregularPlotClass wks 
    "trXMinF" : xmin 
    "trYMinF" : ymin 
    "trXMaxF" : xmax 
    "trYMaxF" : ymax 
    "trXCoordPoints" : x 
  end create 
; 
; By overlaying linear plot on irregular X axis, (Y axis 
; is still linear), we put it through transformation, 
; forcing x axis to be irregular. 
; 
  overlay(irreg_axis,plot0)

; Attach a legend
  labels = (/" HWSD"," NCSCD"," CASA", " MIMICS"," CORPSE"/)
  lgres                    = True
  lgres@lgLineColors       = (/"black","black","(/  0. , 0.286, 0.286 /)", \
             "(/0.286, 0.   , 0.572 /)", "(/0.572, 0.286,  0.   /)"/)
  lgres@lgItemType         = "Lines"        ; show lines only (default)
  lgres@lgLabelFontHeightF = .025           ; legend label font thickness
  lgres@lgDashIndexes      = (/0,1,0,0,0/)
  lgres@pmLegendDisplayMode    = "Never"             ; turn on legend
  lgres@pmLegendOrthogonalPosF = -0.625              ; move units down
  lgres@pmLegendSide           = "Top"
  lgres@pmLegendParallelPosF   = 0.83                ; move units right
  lgres@pmLegendWidthF         = 0.15                ; Change width and
  lgres@pmLegendHeightF        = 0.18                ; height of legend.
  lgres@lgLabelFontHeightF     = .025                ; label font height
  lgres@lgPerimOn              = False               ; turn off box around

  legend = gsn_create_legend (wks, 5, labels,lgres)

  draw(irreg_axis) 
  frame(wks)                 

  delete([/labels, lgres /])


; -- now try w/ NCSCD --
  res@trYMinF                = 45                    ; use this and next line to change timescale that is graphed
  res@trYMaxF                = 90
  res@tiXAxisString          = "Initial stocks ("+cas_TOT_i@units+")"
; res@trXMinF                = 0
; res@trXMaxF                = 60
    res@xyDashPattern  = 0
  res@xyLineThicknessF = 5.0
  res@xyLineColor      = colors(1)
  plot1                = gsn_csm_xy(wks,zTOT(0,:),D_cas_TOT&lat, res)

  res@xyLineColor      = colors(2)
  plot2                = gsn_csm_xy(wks,zTOT(1,:),D_cas_TOT&lat, res)

  res@xyLineColor      = colors(3)
  plot3                = gsn_csm_xy(wks,zTOT(2,:),D_cas_TOT&lat, res)

  res@xyLineColor      = colors(0)
  plot0                = gsn_csm_xy(wks,zTOT(3,:),D_cas_TOT&lat, res)

  res@xyLineThicknessF = 5.0
  res@xyLineColor      = colors(0)
  plot4                = gsn_csm_xy(wks,zOBS,HWSD_SOC&lat,  res)

  res@xyLineThicknessF = 2.0
  res@xyLineColor      = colors(0)
  plot5                = gsn_csm_xy(wks,zNCSCD,NCSCD&lat,  res)

  gsres                   = True                        ; poly res
  gsres@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres@gsFillColor       = "Grey"                 ; color chosen
  dummy = gsn_add_polygon (wks,plot4,xp,yp,gsres)

  overlay(plot0, plot1)
  overlay(plot0, plot2)
  overlay(plot0, plot3)
  overlay(plot0, plot4)
  overlay(plot0, plot5)

; Attach a legend
  labels = (/" HWSD"," NCSCD"," CASA", " MIMICS"," CORPSE"/)
  lgres                    = True
  lgres@lgLineColors       = (/"black","black","(/  0. , 0.286, 0.286 /)", \
             "(/0.286, 0.   , 0.572 /)", "(/0.572, 0.286,  0.   /)"/)
  lgres@lgItemType         = "Lines"        ; show lines only (default)
  lgres@lgLabelFontHeightF = .025           ; legend label font thickness
  lgres@lgDashIndexes      = (/0,0,0,0,0/)           
  lgres@pmLegendDisplayMode    = "Never"             ; turn on legend
  lgres@pmLegendOrthogonalPosF = -0.625              ; move units down
  lgres@pmLegendSide           = "Top"               
  lgres@pmLegendParallelPosF   = 0.83                ; move units right
  lgres@pmLegendWidthF         = 0.15                ; Change width and
  lgres@pmLegendHeightF        = 0.18                ; height of legend.
  lgres@lgLabelFontHeightF     = .025                ; label font height
  lgres@lgPerimOn              = False               ; turn off box around
  
  legend = gsn_create_legend (wks, 5, labels,lgres)

  draw(plot0)                
  frame(wks)                 
  delete([/labels, lgres@lgLineColors, lgres@lgDashIndexes /])



  ; --- zonal mean of turnover times ---
  res@trYMinF                = -50                    ; use this and next line to change timescale that is graphed
  res@trYMaxF                = 80  
  res@tiXAxisString          = "Initial turnover time (y)"
  res@trXMinF                = 1   
  res@trXMaxF                = 1000
  res@xyXStyle              = "Log"
  res@xyLineThicknessF = 3.0
  res@xyLineColor      = colors(1)
  plot1                = gsn_csm_xy(wks,zTAU(0,:),D_cas_TOT&lat, res)

  res@xyLineColor      = colors(2)
  plot2                = gsn_csm_xy(wks,zTAU(1,:),D_cas_TOT&lat, res)

  res@xyLineColor      = colors(3)
  plot3                = gsn_csm_xy(wks,zTAU(2,:),D_cas_TOT&lat, res)

  labels = (/" "," CASA", " MIMICS"," CORPSE"/)
  lgres@pmLegendDisplayMode    = "Always"            ; turn on legend
  lgres@lgMonoDashIndex    = True
  lgres@lgDashIndex        = 0
  lgres@lgLineColors       = colors

  overlay(plot1, plot2)
  overlay(plot1, plot3)
  legend = gsn_create_legend (wks, 3, labels,lgres)
  
  draw(plot1)
  frame(wks)


  res@tiXAxisString          = "Final turnover time (y)"
  res@xyLineColor      = colors(1)
  plot1                = gsn_csm_xy(wks,zTAU_f(0,:),D_cas_TOT&lat, res)

  res@xyLineColor      = colors(2)
  plot2                = gsn_csm_xy(wks,zTAU_f(1,:),D_cas_TOT&lat, res)

  res@xyLineColor      = colors(3)
  plot3                = gsn_csm_xy(wks,zTAU_f(2,:),D_cas_TOT&lat, res)
  overlay(plot1, plot2)
  overlay(plot1, plot3)
  legend = gsn_create_legend (wks, 3, labels,lgres)

  draw(plot1)
  frame(wks)


 
  delete([/lgres, gsres, res@xyLineColor, res@xyLineThicknessF,res@xyXStyle,res@trXMinF/]) 
; delete([/res@trXMinF,res@trXMaxF/])

; --- zonal mean, change in soil C socks --- 
  labels = (/" "," CASA", " MIMICS"," CORPSE"/)
  res@xyLineColors           = colors
  res@xyLineThicknesses      = (/1,5,5,5/)
  res@xyDashPatterns         = (/1,0,0,0/)
 
  res@pmLegendDisplayMode    = "Never"             ; turn on legend
  res@pmLegendOrthogonalPosF = -0.625              ; move units down
  res@pmLegendSide           = "Top"
  res@pmLegendParallelPosF   = 0.73                ; move units right
  res@pmLegendWidthF         = 0.15                ; Change width and
  res@pmLegendHeightF        = 0.18                ; height of legend.
  res@lgLabelFontHeightF     = .025                ; label font height
  res@lgPerimOn              = False               ; turn off box around
  res@xyExplicitLegendLabels = labels              ; create explicit labels
  res@lgItemOrder            = (/3,2,1,0/)         ; reverse legend order
  res@tiXAxisString          = "Change stocks ("+D_cas_TOT@units+")"
  res@trXMinF                = -1000
  res@trXMaxF                =  1000
  plot                       = gsn_csm_xy(wks,zdTOT,D_cas_TOT&lat, res)

  draw(plot)
  frame(wks)

; --- normalized zonal mean, change in soil C socks ---
  delete([/res@trXMinF,res@trXMaxF/])
  res@trXMinF                = -10
  res@trXMaxF                =  25
  res@tiXAxisString          = "Normalized stock change(%)"
  plot                       = gsn_csm_xy(wks,zdTOT_norm,D_cas_TOT&lat, res)
  draw(plot)
  frame(wks)

  print("------------------------")
  print("wrote "+ fout+".ps")
  print("")
  delete([/plot, plot2,res, fout, wks/])
;-------------------------------------------------------------------
exit
;-------------------------------------------------------------------


;**************************************
; Plot variance of HR by lat vs. time
;**************************************
 fout  = out_dir+"HR_lat_time_"+exper
 wks  = gsn_open_wks("ps",fout)

 gsn_define_colormap(wks, "WhViBlGrYeOrRe")

 res                     = True                ; plot mods desired
 res@gsnMaximize         = True
;res@gsnPaperOrientation = "landscape"
 res@gsnDraw             = False
 res@gsnFrame            = False
 res@cnSmoothingOn       = False
 res@lbOrientation       = "Horizontal"
 res@cnInfoLabelOn       = False           ; turn off cn info label
 res@cnFillOn            = True            ; turn on color
 res@cnLinesOn           = False           ; turn off contour lines
 res@cnLineLabelsOn      = False           ; turn off contour labels
 res@gsnSpreadColors     = True            ; spread out color table
 res@cnFillOn            = True            ; color on
 res@cnLinesOn           = False           ; turn off contour lines

 res@vpXF             = 0.12                ; default is 0.2 (aspect ratio)
 res@vpYF             = 0.8                 ; default is 0.8
 res@vpHeightF        = 0.4                 ; default is 0.6
 res@vpWidthF         = 0.75                ; default is 0.6
 res@trYMinF          = -50                 ; use this and next line to change timescale that is graphed
 res@trYMaxF          = 85
 res@trXMinF          = 1940
 res@trXMaxF          = 2100

 res@lbLabelBarOn         = False           ; turn off individual cb's
 res@cnLevelSelectionMode = "ExplicitLevels"
 res@cnLevels             = fspan(0,0.1,11)

 zHR_CV   = zHR_runcv(model|:,lat|:,time|:)
 zTAIR_CV = zTAIR_runcv(lat|:,time|:)
 zNPP_CV  = zNPP_runcv(lat|:,time|:)
 zGPP_CV  = zGPP_runcv(lat|:,time|:)

 zHR_runVAR   = zHR_runStdev(model|:,lat|:,time|:)
 zTAIR_runVAR = zTAIR_runStdev(lat|:,time|:)
 zNPP_runVAR  = zNPP_runStdev(lat|:,time|:)
 zGPP_runVAR  = zGPP_runStdev(lat|:,time|:)

 zHR_runVAR   = zHR_runVAR^2
 zTAIR_runVAR = zTAIR_runVAR^2
 zNPP_runVAR  = zNPP_runVAR^2
 zGPP_runVAR  = zGPP_runVAR^2

 zHR_runVAR@long_name   = "HR, 20 y running variance"
 zTAIR_runVAR@long_name = "TAIR, 20 y running variance"
 zNPP_runVAR@long_name  = "NPP, 20 y running variance"
 zGPP_runVAR@long_name  = "GPP, 20 y running variance"

  resP                  = True               ; modify the panel plot
  resP@gsnFrame         = False              ; don't advance panel plot
  resP@gsnPanelLabelBar = True               ; add common colorbar
  resP@lbOrientation    =   "horizontal"     ; vertical label bar

 res@cnLevels             = fspan(0,0.1,11)
 plot    = new(3, "graphic") 
 plot(0) = gsn_csm_lat_time(wks, zHR_CV(0,:,:), res)
 plot(1) = gsn_csm_lat_time(wks, zHR_CV(1,:,:), res)
 plot(2) = gsn_csm_lat_time(wks, zHR_CV(2,:,:), res)

  gsn_panel(wks,plot,(/1,3/),resP)               ; now draw as one plot
  frame(wks)

 delete(res@cnLevels)
 res@cnLevels             = fspan(50,400,15)
 plot    = new(3, "graphic")
 plot(0) = gsn_csm_lat_time(wks, zHR_runVAR(0,:,:), res)
 plot(1) = gsn_csm_lat_time(wks, zHR_runVAR(1,:,:), res)
 plot(2) = gsn_csm_lat_time(wks, zHR_runVAR(2,:,:), res)

  gsn_panel(wks,plot,(/1,3/),resP)               ; now draw as one plot
  frame(wks)

 delete([/res@cnLevelSelectionMode, res@cnLevels,resP@gsnPanelLabelBar,resP@lbOrientation/])

 res@lbLabelBarOn         = True            ; turn off individual cb's
 res@lbOrientation       = "Horizontal"

  plot2    = new(3,"graphic")

  plot2(0) = gsn_csm_lat_time(wks, zTAIR_CV, res)
  plot2(1) = gsn_csm_lat_time(wks, zGPP_CV, res)
  plot2(2) = gsn_csm_lat_time(wks, zNPP_CV, res)
  gsn_panel(wks,plot2,(/1,3/),resP)               ; now draw as one plot
  frame(wks)

  plot2(0) = gsn_csm_lat_time(wks, zTAIR_runVAR, res)
  plot2(1) = gsn_csm_lat_time(wks, zGPP_runVAR, res)
  plot2(2) = gsn_csm_lat_time(wks, zNPP_runVAR, res)
  gsn_panel(wks,plot2,(/1,3/),resP)               ; now draw as one plot
  frame(wks)

  print("------------------------")
  print("wrote "+ fout+".ps")
  print("")


end

