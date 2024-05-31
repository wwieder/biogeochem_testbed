"""utility functions"""

import re
import cftime
import numpy as np
import xarray as xr

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.lines as mlines

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point

#from xr_ds_ex import xr_ds_ex

# sum of all soil C and N pools for MIMICS & CASA
# Slightly different logic when data sets are combined
def sum_pools_combined(ds_in, mod='MIM', CN='True'):
    if CN == 'True': 
        zeroFxn = 1. 
    else:
        zeroFxn = np.nan
    # MIMICS specific pools
    if mod == 'mim':
        ds_in['cTOT'] = ds_in['cLITm']+ds_in['cLITs']+ds_in['cMICr']+ds_in['cMICk']+ \
                        ds_in['cSOMa']+ds_in['cSOMc']+ds_in['cSOMp']
        ds_in['cTOT'].attrs['long_name'] = 'total soil C'
        ds_in['cTOT'].attrs['units'] = ds_in['cLITm'].attrs['units']
        
        ds_in['nTOT'] = (ds_in['nLITm']+ds_in['nLITs']+ds_in['nMICr']+ds_in['nMICk']+ \
                        ds_in['nSOMa']+ds_in['nSOMc']+ds_in['nSOMp']) * zeroFxn
        ds_in['nTOT'].attrs['long_name'] = 'total soil N'
        ds_in['nTOT'].attrs['units'] = ds_in['nLITm'].attrs['units']        

        # CWD fluxes for MIMICS are in cresp, soilC fluxes in cHresp?
        ds_in['cresp'] = ds_in['cHresp'] #+ ds_in['cresp'] 

        ds_in['overRESP'] = ds_in['cOverflow_r'] + ds_in['cOverflow_k']
        ds_in['overRESP'].attrs['long_name'] = 'overflow respiration'
        ds_in['overRESP'].attrs['units'] = ds_in['cOverflow_r'].attrs['units'] 
    
    # CASA specific pools
    if mod == 'cas':
        ds_in['cTOT'] = ds_in['clitmetb']+ds_in['clitstr']+ds_in['csoilmic']+ \
                        ds_in['csoilslow']+ds_in['csoilpass']
        ds_in['cTOT'].attrs['long_name'] = 'total soil C'
        ds_in['cTOT'].attrs['units'] = ds_in['clitmetb'].attrs['units']

        ds_in['nTOT'] = (ds_in['nlitmetb']+ds_in['nlitstr']+ds_in['nsoilmic']+ \
                            ds_in['nsoilslow']+ds_in['nsoilpass'])* zeroFxn
        ds_in['nTOT'].attrs['long_name'] = 'total soil N'
        ds_in['nTOT'].attrs['units'] = ds_in['nlitmetb'].attrs['units']
        
    # output for both models, will be nan for casa or c-only cases
    ds_in['cnLIT']= (ds_in['cLitInput_metb'] + ds_in['cLitInput_struc']) * zeroFxn / \
                    (ds_in['nLitInput_metb'] + ds_in['nLitInput_struc'])
    ds_in['cnLIT'].attrs['long_name'] = 'Litterfall C:N ratio'
    ds_in['cnLIT'].attrs['units'] = ''

    ds_in['cMICtot'] = ds_in['cMICr']+ds_in['cMICk']
    ds_in['cMICtot'].attrs['long_name'] = 'sum of MIC pools'
    ds_in['cMICtot'].attrs['units'] = ds_in['cMICr'].attrs['units']
        
    ds_in['relMIC'] = ds_in['cMICtot']/ds_in['cTOT'] * 100
    ds_in['relMIC'].attrs['long_name'] = 'MIC:SOC'
    ds_in['relMIC'].attrs['units'] = '%'
        
    ds_in['specRESP'] = 1000.* ds_in['cresp'] / ds_in['cMICtot']
    ds_in['specRESP'].attrs['long_name'] = 'specific respiration'
    ds_in['specRESP'].attrs['units'] = 'm'+ds_in['cHresp'].attrs['units'] +'/'+ ds_in['cMICk'].attrs['units']
        
    ds_in['MICrK'] = ds_in['cMICr']/ds_in['cMICk']
    ds_in['MICrK'].attrs['long_name'] = 'MICr:MICK'
    ds_in['MICrK'].attrs['units'] = ''

    ds_in['cnTOT']= ds_in['cTOT'] / ds_in['nTOT']
    ds_in['cnTOT'].attrs['long_name'] = 'Ratio SOM C:N ratio'
    ds_in['cnTOT'].attrs['units'] = ''
            
    ds_in['cnMIC']= (ds_in['cMICr'] + ds_in['cMICk']) * zeroFxn/ \
                    (ds_in['nMICr'] + ds_in['nMICk'])
    ds_in['cnMIC'].attrs['long_name'] = 'MIMICS microbial C:N ratio'
    ds_in['cnMIC'].attrs['units'] = ''

    #ds_in['cnpp'] = ds_in['cnpp']
    #ds_in['cgpp'] = ds_in['cgpp']
    
    ds_in['NEP'] = ds_in['cnpp']-ds_in['cresp']
    ds_in['NEP'].attrs['long_name'] = 'net ecosystem production'
    ds_in['NEP'].attrs['units'] = ds_in['cnpp'].attrs['units']

    ds_in['cVEG'] = ds_in['cleaf']+ds_in['cfroot']+ds_in['cwood']
    ds_in['cVEG'].attrs['long_name'] = 'total veg C'
    ds_in['cVEG'].attrs['units'] = ds_in['cleaf'].attrs['units']

    ds_in['cECO'] = ds_in['cTOT']+ds_in['cVEG'] + ds_in['clitcwd']
    ds_in['cECO'].attrs['long_name'] = 'total ecosystem C'
    ds_in['cECO'].attrs['units'] = ds_in['cleaf'].attrs['units']
        
    ds_in['nVEG'] = (ds_in['nleaf']+ds_in['nfroot']+ds_in['nwood']) * zeroFxn
    ds_in['nVEG'].attrs['long_name'] = 'total veg N'
    ds_in['nVEG'].attrs['units'] = ds_in['nleaf'].attrs['units']

    ds_in['nECO'] = (ds_in['nTOT']+ds_in['nVEG'] + ds_in['nlitcwd']) * zeroFxn
    ds_in['nECO'].attrs['long_name'] = 'total ecosystem N'
    ds_in['nECO'].attrs['units'] = ds_in['cleaf'].attrs['units']
        
    ds_in['cnVEG'] = ds_in['cVEG'] / ds_in['nVEG']  * zeroFxn
    ds_in['cnVEG'].attrs['long_name'] = 'total veg C:N'
    ds_in['cnVEG'].attrs['units'] = ''
            
    ds_in['cnECO'] = ds_in['cECO'] / ds_in['nECO']  * zeroFxn
    ds_in['cnECO'].attrs['long_name'] = 'total ecosystem C:N'
    ds_in['cnECO'].attrs['units'] = ''
        
    return ds_in


# sum of all soil C and N pools for MIMICS & CASA
def sum_pools(ds_in, mod='MIM', CN='True'):
    if mod == 'mim':
        ds_in['cTOT'] = ds_in['cLITm']+ds_in['cLITs']+ds_in['cMICr']+ds_in['cMICk']+ \
                        ds_in['cSOMa']+ds_in['cSOMc']+ds_in['cSOMp']
        ds_in['cTOT'].attrs['long_name'] = 'total soil C'
        ds_in['cTOT'].attrs['units'] = ds_in['cLITm'].attrs['units']

        ds_in['cMICtot'] = ds_in['cMICr']+ds_in['cMICk']
        ds_in['cMICtot'].attrs['long_name'] = 'sum of MIC pools'
        ds_in['cMICtot'].attrs['units'] = ds_in['cMICr'].attrs['units']
        
        ds_in['relMIC'] = ds_in['cMICtot']/ds_in['cTOT'] * 100
        ds_in['relMIC'].attrs['long_name'] = 'MIC:SOC'
        ds_in['relMIC'].attrs['units'] = '%'
        
        # CWD fluxes for MIMICS are in cresp, soilC fluxes in cHresp?
        ds_in['cresp'] = ds_in['cHresp'] #+ ds_in['cresp'] 

        ds_in['specRESP'] = ds_in['cresp'] / ds_in['cMICtot']
        ds_in['specRESP'].attrs['long_name'] = 'specific respiration'
        ds_in['specRESP'].attrs['units'] = ds_in['cHresp'].attrs['units'] +'/'+ ds_in['cMICk'].attrs['units']
        
        ds_in['MICrK'] = ds_in['cMICr']/ds_in['cMICk']
        ds_in['MICrK'].attrs['long_name'] = 'MICr:MICK'
        ds_in['MICrK'].attrs['units'] = ''

        if CN == 'True':
            ds_in['nTOT'] = ds_in['nLITm']+ds_in['nLITs']+ds_in['nMICr']+ds_in['nMICk']+ \
                            ds_in['nSOMa']+ds_in['nSOMc']+ds_in['nSOMp']
            ds_in['nTOT'].attrs['long_name'] = 'total soil N'
            ds_in['nTOT'].attrs['units'] = ds_in['nLITm'].attrs['units']

            ds_in['cnTOT']= ds_in['cTOT'] / ds_in['nTOT']
            ds_in['cnTOT'].attrs['long_name'] = 'MIMICS total SOM C:N ratio'
            ds_in['cnTOT'].attrs['units'] = ''
            
            ds_in['cnMIC']= (ds_in['cMICr'] + ds_in['cMICk'])/ \
                            (ds_in['nMICr'] + ds_in['nMICk'])
            ds_in['cnMIC'].attrs['long_name'] = 'MIMICS microbial C:N ratio'
            ds_in['cnMIC'].attrs['units'] = ''
            
            ds_in['cnLIT']= (ds_in['cLitInput_metb'] + ds_in['cLitInput_struc']) / \
                            (ds_in['nLitInput_metb'] + ds_in['nLitInput_struc'])
            ds_in['cnLIT'].attrs['long_name'] = 'Litterfall C:N ratio'
            ds_in['cnLIT'].attrs['units'] = ''

        
    if mod == 'cas':
        ds_in['cTOT'] = ds_in['clitmetb']+ds_in['clitstr']+ds_in['csoilmic']+ \
                        ds_in['csoilslow']+ds_in['csoilpass']
        ds_in['cTOT'].attrs['long_name'] = 'total soil C'
        ds_in['cTOT'].attrs['units'] = ds_in['clitmetb'].attrs['units']

        if CN == 'True':
            ds_in['nTOT'] = ds_in['nlitmetb']+ds_in['nlitstr']+ds_in['nsoilmic']+ \
                            ds_in['nsoilslow']+ds_in['nsoilpass']
            ds_in['nTOT'].attrs['long_name'] = 'total soil N'
            ds_in['nTOT'].attrs['units'] = ds_in['nlitmetb'].attrs['units']

            ds_in['cnTOT']= ds_in['cTOT'] / ds_in['nTOT']
            ds_in['cnTOT'].attrs['long_name'] = 'CASA total SOM C:N ratio'
            ds_in['cnTOT'].attrs['units'] = ''

    ds_in['NEP'] = ds_in['cnpp']-ds_in['cresp']
    ds_in['NEP'].attrs['long_name'] = 'net ecosystem production'
    ds_in['NEP'].attrs['units'] = ds_in['cnpp'].attrs['units']

    ds_in['cVEG'] = ds_in['cleaf']+ds_in['cfroot']+ds_in['cwood']
    ds_in['cVEG'].attrs['long_name'] = 'total veg C'
    ds_in['cVEG'].attrs['units'] = ds_in['cleaf'].attrs['units']

    ds_in['cECO'] = ds_in['cTOT']+ds_in['cVEG'] + ds_in['clitcwd']
    ds_in['cECO'].attrs['long_name'] = 'total ecosystem C'
    ds_in['cECO'].attrs['units'] = ds_in['cleaf'].attrs['units']
        
    if CN == 'True':
        ds_in['nVEG'] = ds_in['nleaf']+ds_in['nfroot']+ds_in['nwood']
        ds_in['nVEG'].attrs['long_name'] = 'total veg N'
        ds_in['nVEG'].attrs['units'] = ds_in['nleaf'].attrs['units']

        ds_in['nECO'] = ds_in['nTOT']+ds_in['nVEG'] + ds_in['nlitcwd']
        ds_in['nECO'].attrs['long_name'] = 'total ecosystem N'
        ds_in['nECO'].attrs['units'] = ds_in['cleaf'].attrs['units']
        
        ds_in['cnVEG'] = ds_in['cVEG'] / ds_in['nVEG']
        ds_in['cnVEG'].attrs['long_name'] = 'total veg C:N'
        ds_in['cnVEG'].attrs['units'] = ''
            
        ds_in['cnECO'] = ds_in['cECO'] / ds_in['nECO']
        ds_in['cnECO'].attrs['long_name'] = 'total ecosystem C:N'
        ds_in['cnECO'].attrs['units'] = ''
        
    return ds_in


####### Function to calculate global fluxes / stocks ###########################
def globalSum ( dsIn, var, time=-1, conversion=1e-15, units=None, plot=True ):
    area = dsIn.landarea *1e6
    temp = dsIn[var] * area 
    temp = temp.sum(dim=('lat','lon')) * conversion
    if plot==True:
        plt.figure(figsize=[25,6]);
        for i in range(len(var)):
            plt.subplot(1, 4, (1+i))
            plt.ylabel('Global '+var[i]+' (Pg C y^-1)')
            plt.plot(temp['time'], temp[var[i]]);

    else:
        for i in range(len(var)):
            print('global '+var[i]+' '+str(np.round(temp[var[i]].isel(time=time).values,1) )+' Pg C')


######### Generate a function for making panel plots of maps ###############
def map_function(da, cb=0, panel=None, cmap=None, ax=None, 
                 title=None, vmax=None, vmin=None):
    '''a function to make one subplot'''
    wrap_data, wrap_lon = add_cyclic_point(da.values, coord=da.lon)

    if ax is None: ax = plt.gca()
    im = ax.pcolormesh(wrap_lon,da.lat,wrap_data,
                   transform=ccrs.PlateCarree(),
                   vmax=vmax,vmin=vmin,cmap=cmap)
    ax.set_title(title)
    ax.coastlines()
    ax.set_extent([-180,180,-65,80], crs=ccrs.PlateCarree())
    ax.annotate(panel, xy=(0.05, 0.95), xycoords=ax.transAxes,
                ha='center', va='center',fontsize=16)    

    # allows for different colorbars on each plot
    if cb == 1:  # here to right of plots
        fig.colorbar(im, ax=ax,shrink=0.40, pad=0, fraction = 0.1)

    # allows for different colorbars on each plot
    if cb == 2:  # here below plots
        fig.colorbar(im, ax=ax,shrink=0.9, pad=0, fraction = 1, orientation="horizontal")

# strings for labeling panels
panel = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)'] 

################### Function to truncate color map ###################
def truncate_colormap(cmapIn='jet', minval=0.0, maxval=1.0, n=100):
    cmapIn = plt.get_cmap(cmapIn,n)

    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmapIn.name, a=minval, b=maxval),
        cmapIn(np.linspace(minval, maxval, n)))

    arr = np.linspace(0, 50, 100).reshape((10, 10))
    return new_cmap

# add cyclic point
def cyclic_dataarray(da, coord='lon'):
    """ Add a cyclic coordinate point to a DataArray along a specified
    named coordinate dimension.
    >>> from xray import DataArray
    >>> data = DataArray([[1, 2, 3], [4, 5, 6]],
    ...                      coords={'x': [1, 2], 'y': range(3)},
    ...                      dims=['x', 'y'])
    >>> cd = cyclic_dataarray(data, 'y')
    >>> print cd.data
    array([[1, 2, 3, 1],
           [4, 5, 6, 4]])
    """
    assert isinstance(da, xr.DataArray)

    lon_idx = da.dims.index(coord)
    cyclic_data, cyclic_coord = add_cyclic_point(da.values,
                                                 coord=da.coords[coord],
                                                 axis=lon_idx)

    # Copy and add the cyclic coordinate and data
    new_coords = dict(da.coords)
    new_coords[coord] = cyclic_coord
    new_values = cyclic_data

    new_da = xr.DataArray(new_values, dims=da.dims, coords=new_coords)

    # Copy the attributes for the re-constructed data and coords
    for att, val in da.attrs.items():
        new_da.attrs[att] = val
    for c in da.coords:
        for att in da.coords[c].attrs:
            new_da.coords[c].attrs[att] = da.coords[c].attrs[att]

    return new_da

# as above, but for a dataset
# doesn't work because dims are locked in a dataset
'''
def cyclic_dataset(ds, coord='lon'):
    assert isinstance(ds, xr.Dataset)

    lon_idx = ds.dims.index(coord)
    cyclic_data, cyclic_coord = add_cyclic_point(ds.values,
                                                 coord=ds.coords[coord],
                                                 axis=lon_idx)

    # Copy and add the cyclic coordinate and data
    new_coords = dict(ds.coords)
    new_coords[coord] = cyclic_coord
    new_values = cyclic_data

    new_ds = xr.DataSet(new_values, dims=ds.dims, coords=new_coords)

    # Copy the attributes for the re-constructed data and coords
    for att, val in ds.attrs.items():
        new_ds.attrs[att] = val
    for c in ds.coords:
        for att in ds.coords[c].attrs:
            new_ds.coords[c].attrs[att] = ds.coords[c].attrs[att]

    return new_ds
'''
