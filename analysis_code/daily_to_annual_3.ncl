; averages and concatinates daily output into annual files

 
begin
; ---- for CASA ----------
;  dir_in  = (/"/project/tss/wwieder/CASACLM/GRID/OUTPUT_GRID_CRU_NCEP/CASACNP_HIST/"/)
;  dir_in  = (/"/project/tss/wwieder/CASACLM/GRID/OUTPUT_GRID_CRU_NCEP/CASACNP_RCP45/"/)
;  model   = "casaclm"

; ---- for MIMICS ---------
;  dir_in  = (/"/project/tss/wwieder/CASACLM/GRID/OUTPUT_GRID_CRU_NCEP/MIMICS_HIST/"/)
;  dir_in  = (/"/project/tss/wwieder/CASACLM/GRID/OUTPUT_GRID_CRU_NCEP/MIMICS_RCP85/"/)

   dir_in  = (/"/project/tss/wwieder/CASACLM/GRID/OUTPUT_GRID_CRU_NCEP/MIMICS_HIST_Vslope2/"/)
;  dir_in  = (/"/project/tss/wwieder/CASACLM/GRID/OUTPUT_GRID_CRU_NCEP/MIMICS_RCP45_Vslope2/"/)
   model   = "mimics"

; ---- for CORPSE ---------
;  dir_in  = (/"/project/tss/wwieder/CASACLM/GRID/OUTPUT_GRID_CRU_NCEP/CORPSE_HIST_minFW/"/)
;  dir_in  = (/"/project/tss/wwieder/CASACLM/GRID/OUTPUT_GRID_CRU_NCEP/CORPSE_RCP45_minFW/"/)
;  model   = "corpse"

  dir_out = dir_in
  yrStrt = 1901                         
  yrLast = 2010
; yrStrt = 2011                         
; yrLast = 2100
  year   = ispan(yrStrt, yrLast, 1)
  nyear  = dimsizes(year)

  print("dir in"+dir_in)
  do y = 0, nyear-1
      YEAR  = year(y)

      ncocmd="ncra -O "+dir_in+model+"_pool_flux_"+YEAR+"_daily.nc "+dir_out+model+"_pool_flux_"+YEAR+"_ann.nc"
      system (ncocmd)
      print("finished "+model+" "+YEAR)

  end do

  ncocmd = "ncrcat -O "+dir_out+model+"_pool_flux_????_ann.nc "+dir_out+"ann_"+model+"_pool_flux_"+yrStrt+"_"+yrLast+".nc"
  system (ncocmd)    

  ncocmd = "rm "+dir_out+"*_pool_flux_????_ann.nc"
  system (ncocmd)
  print("finished annual averaging and cat files")

end

