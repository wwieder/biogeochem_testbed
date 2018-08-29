#!/usr/local/bin/perl

####################################################################################################
#
# FILE:  run-netcdftools.pl
#
# AUTHOR: Melannie Hartman
# DATE:   April 14, 2014
#         August 10, 2015 (reviewed and updated)
#
# USAGE: perl run-netcdftools.pl
#
# PURPOSE:
#   Run the netcdfTools program from $startYr to $endYr to create met*.nc files 
#   for the CASACNP model. Meterologcial and GPP values written to the met*.nc files 
#   come from for CLM daily climate files.  
#   Daily CLM history files (NetCDF) are available for years 1901-2011.
#
# DESCRIPTION:
#   This script runs netcdfTools repeatedly, one year at a time, from $startYr to $endYr.
#   NetcdfTools generates two files for each calendar year (yyyy):
#     * met_yyyy_yyyy.nc - contains meteorological and GPP drivers for the CASACNP model 
#     * gridinfo_igbpz_yyyy.csv - contains annual N deposition for the CASACNP model 
#   Make sure there is plenty of disk space in $outNcDir. Each "met_yyyy_yyyy.nc" file 
#   is ~ 0.3 GB. NetcdfTools also produces a file named GPP_Daily_Output.yyyy-01-01-00000.nc
#   for each year, where each file is ~19.5 MB.  These files are also moved to $outNcDir.
#   In order to run netcdfTools, this script creates a files.ini file for each year in the 
#   format shown below.
#
# Example format of files.ini (7 lines):
#
# h1         ! file type (h0=monthly, h1=daily)
# 19010101   ! start (yyyymmdd)
# 19011231   ! end (yyyymmdd)
# /project/bgc01/melannie/clm45sp_2deg4506_hist.clm2.h1.1901-01-01-00000.nc
# /project/bgc01/melannie/surfdata_1.9x2.5_simyr1850_c130421.nc
# /project/bgc01/melannie/fndep_clm_hist_simyr1849-2006_1.9x2.5_c100428.nc
# ./CO2_1768-2100.txt
#
# line 1: Set to h1 to ready daily history files
# line 2: Starting date for data to be written to met*.nc
# line 3: Ending date for data to be written to met*.nc 
#   * In this script only one year of data is written to each met*.nc file.
#     However, up to about 7 years can be written to a met*.nc file before it becomes too large.
# line 4: CLM daily history files ($histNcFile)
#  * This is a file name template because the date string in the file name (1901-01-01-00000 here)
#    is replaced with correct date.
# line 5: CLM surface dataset ($surfDataNcFile)
# line 6: CLM N deposition ($ndepNcFile)
# line 7: File containing annual mean atmospheric CO2 concentrations ($co2file)
#
####################################################################################################

# Executable program
$exefile = "./netcdfTools";


# Directory where input NetCDF CLM DAILY history files
#$inNcDir = "/project/bgc01/bonan/casa-cnp/clm/";
#$inNcDir = "/project/tss/wwieder/CASACLM/clm_forcing/CRU_hist/";
#$inNcDir = "/project/tss/wwieder/CASACLM/clm_forcing/GSWP3_hist/";
#$inNcDir = "/project/tss/wwieder/CASACLM/clm_forcing/CRU_RCP45/";
#$inNcDir = "/project/tss/wwieder/CASACLM/clm_forcing/CRU_RCP85/";
$inNcDir = "/project/tss/wwieder/CASACLM/clm_forcing/CLM5sp_GSWP3_hist/";

# Directory where surface data set and N deposition are stored
$inNcDir2 = "/project/tss/wwieder/CASACLM/";
$inNcDir3 = "/project/tss/wwieder/CASACLM/";

# Directory where output NetCDF met*.nc files will be stored
#$outNcDir = "/project/bgc01/melannie/CASACLM/GRID/INPUT_MET_GRID";
#$outNcDir = "/project/bgc01/melannie/CASACLM/GRID/INPUT_MET_GRID_CRU_NCEP_SOILLIQ";
#$outNcDir = "/project/tss/wwieder/CASACLM/GRID/INPUT_MET_GRID_CRU_NCEP_SOILLIQ";
#$outNcDir = "/project/bgc01/melannie/CASACLM/GRID/INPUT_MET_GRID_GSWP3";
#$outNcDir = "/project/tss/bgc01/melannie/CASACLM/GRID/INPUT_MET_GRID_CRU_RCP85";
#$outNcDir = "/project/tss/wwieder/CASACLM/GRID/INPUT_MET_GRID_GSWP3_SOILLIQ";
#$outNcDir = "/project/tss/wwieder/CASACLM/GRID/INPUT_MET_GRID_CRU_RCP45";
$outNcDir = "/project/tss/wwieder/CASACLM/GRID/INPUT_MET_GRID_GSWP3_CLM5_hist";
# Directory where output NetCDF GPP*.nc files will be stored
# $outNcDir2 = "/project/bgc01/melannie/CanopyModel/";

$filesin = "./files.ini";
#$histNcFile = "${inNcDir}clm45sp_2deg4506_hist.clm2.h1.1901-01-01-00000.nc";
#$histNcFile = "${inNcDir}clm4_5_12_r191_CLM45spHIST_CRU.clm2.h1.1901-01-01-00000.nc";
#$histNcFile = "${inNcDir}clm4_5_12_r191_CLM45spHIST_GSWP3.clm2.h1.1901-01-01-00000.nc";
$histNcFile = "${inNcDir}CLM5sp_HIST_GSWP3.clm2.h1.1901-01-01-00000.nc";
#$histNcFile = "${inNcDir}clm4_5_12_r191_CLM45sp_RCP85.clm2.h1.2011-01-01-00000.nc";

$surfDataNcFile = "${inNcDir3}surfdata_1.9x2.5_16pfts_Irrig_CMIP6_simyr1850_c170824.nc";
$ndepNcFile = "${inNcDir2}fndep_clm_rcp4.5_simyr1849-2106_1.9x2.5_c100428.nc";
$co2file = "./CO2_1768-2100.txt";

system("cp clmGrid_IGBP_grid.csv clmGrid_IGBP.csv");

$logfile = "./log.txt";
unlink($logfile);

$startYr = "1901";
$endYr = "2014";

for ($year1 = $startYr; $year1 <= $endYr; $year1++)
{
   $year2 = $year1;
   open(INI, ">$filesin") || die "Can not open file for writing: $filesin\n";
   printf(INI "%s\n", "h1             ! file type (h0=monthly, h1=daily)");
   printf(INI "%s\n", "${year1}0101   ! start (yyyymmdd)");
   printf(INI "%s\n", "${year2}1231   ! end (yyyymmdd)");
   printf(INI "%s\n", $histNcFile);
   printf(INI "%s\n", $surfDataNcFile);
   printf(INI "%s\n", $ndepNcFile);
   printf(INI "%s\n", $co2file);
   close(INI);

   #$cmd = "${exefile} >> ${logfile}";
   $cmd = "${exefile}"; 
   print "$cmd\n";
   system($cmd);

# Moving met*.nc files does not seem to work.  
# Copy them to /project/bgc/ then remove the file from the local directory.
# Melannie 9/29/2014.

   $cmd = "cp met_${year1}_${year2}.nc $outNcDir";
   print "$cmd\n";
   system($cmd);

   $cmd = "rm met_${year1}_${year2}.nc";
   print "$cmd\n";
   system($cmd);

#  $cmd = "mv GPP_Daily_Output.${year1}-01-01-00000.nc $outNcDir2";
#  print "$cmd\n";
#  system($cmd);

#  $cmd = "mv gridinfo_igbpz.csv ${outNcDir}gridinfo_igbpz_${year1}.csv ";
#  print "$cmd\n";
#  system($cmd);

}

