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
#   Daily CLM history files (NetCDF) are available for years 
#      1901-2014 (GSWP3) or 
#      1901-2019 w/ CRUJRA (TRENDY)
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
$inNcDir = "/project/tss/wwieder/biogeochem_testbed_1.1/clm_forcing/clm50_dev110_spHIST/";

# Directory where surface data set and N deposition are stored
$inNcDir2 = "/project/tss/wwieder/biogeochem_testbed_1.1/clm_forcing/";
$inNcDir3 = "/project/tss/wwieder/biogeochem_testbed_1.1/clm_forcing/";

# Directory where output NetCDF met*.nc files will be stored
$outNcDir = "/project/tss/wwieder/biogeochem_testbed_1.1/GRID_CN/INPUT_MET_GRID_GSWP3_CLM5dev110_hist";

$filesin = "./files.ini";
$histNcFile = "${inNcDir}clm50_dev110_spHIST.clm2.h1.1901-01-01-00000.nc";

$surfDataNcFile = "${inNcDir3}surfdata_1.9x2.5_hist_16pfts_Irrig_CMIP6_simyr1850_c190304.nc";
$ndepNcFile = "${inNcDir2}fndep_clm_rcp8.5_simyr1849-2106_1.9x2.5_c100428.nc";
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

