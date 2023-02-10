; averages and concatinates daily output into annual files

 
begin

  dir    = (/"/project/tss/wwieder/biogeochem_testbed_1.1/GRID_CN/MIMICS_mod5_GSWP3_KO4_exudate0_cwdN/OUTPUT_CN/HIST/"/)
;  dir    = (/"/project/tss/wwieder/biogeochem_testbed_1.1/GRID_CN/MIMICS_mod5_GSWP3_KO4_exudate0_cwdN/OUTPUT_CN/SSP370/"/)
  yrStrt = 1901                         
  yrLast = 2014
  year   = ispan(yrStrt, yrLast, 1)
  nyear  = dimsizes(year)

  print(dir)

  do y = 0, nyear-1
      YEAR  = year(y)
      ncocmd="ncra -O "+dir+"casaclm_pool_flux_"+YEAR+"_daily.nc "+dir+"casaclm_pool_flux_"+YEAR+"_ann.nc"
      system (ncocmd)

      ncocmd="ncra -O "+dir+"mimics_pool_flux_"+YEAR+"_daily.nc "+dir+"mimics_pool_flux_"+YEAR+"_ann.nc"
      system (ncocmd)
      print("finished "+YEAR)
  end do

  ncocmd = "ncrcat -O "+dir+"casaclm_pool_flux_????_ann.nc "+dir+"ann_casaclm_pool_flux_"+yrStrt+"_"+yrLast+".nc"
  system (ncocmd) 

  ncocmd = "ncrcat -O "+dir+"mimics_pool_flux_????_ann.nc "+dir+"ann_mimics_pool_flux_"+yrStrt+"_"+yrLast+".nc"
  system (ncocmd)    

  ncocmd = "rm "+dir+"*_pool_flux_????_ann.nc"
  system (ncocmd)
  print("finished annual averaging and cat files")

end

