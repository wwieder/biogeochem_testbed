"""utility functions"""

import re
import cftime
import numpy as np
import xarray as xr

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point

#from xr_ds_ex import xr_ds_ex

# sum of all soil C and N pools for MIMICS & CASA
def sum_pools(ds_in, mod='MIM', CN='True'):
    if mod == 'mim':
        ds_in['cTOT'] = ds_in['cLITm']+ds_in['cLITs']+ds_in['cMICr']+ds_in['cMICk']+ \
                        ds_in['cSOMa']+ds_in['cSOMc']+ds_in['cSOMp']
        ds_in['cTOT'].attrs['long_name'] = 'sum of MIMICS SOC pools'
        ds_in['cTOT'].attrs['units'] = ds_in['cLITm'].attrs['units']

        ds_in['relMIC'] = (ds_in['cMICr']+ds_in['cMICk'])/ds_in['cTOT'] * 100
        ds_in['relMIC'].attrs['long_name'] = 'MIC:SOC'
        ds_in['relMIC'].attrs['units'] = '%'
        
        ds_in['MICrK'] = ds_in['cMICr']/ds_in['cMICk']
        ds_in['MICrK'].attrs['long_name'] = 'MICr:MICK'
        ds_in['MICrK'].attrs['units'] = ''

        if CN == 'True':
            ds_in['nTOT'] = ds_in['nLITm']+ds_in['nLITs']+ds_in['nMICr']+ds_in['nMICk']+ \
                            ds_in['nSOMa']+ds_in['nSOMc']+ds_in['nSOMp']
            ds_in['nTOT'].attrs['long_name'] = 'sum of MIMICS SON pools'
            ds_in['nTOT'].attrs['units'] = ds_in['nLITm'].attrs['units']

            ds_in['cnTOT']= ds_in['cTOT'] / ds_in['nTOT']
            ds_in['cnTOT'].attrs['long_name'] = 'MIMICS total SOM C:N ratio'
            
            ds_in['cnMIC']= (ds_in['cMICr'] + ds_in['cMICk'])/ (ds_in['nMICr'] + ds_in['nMICk'])
            ds_in['cnTOT'].attrs['long_name'] = 'MIMICS microbial C:N ratio'
        
    if mod == 'cas':
        ds_in['cTOT'] = ds_in['clitmetb']+ds_in['clitstr']+ds_in['csoilmic']+ \
                        ds_in['csoilslow']+ds_in['csoilpass']
        ds_in['cTOT'].attrs['long_name'] = 'sum of CASA SOC pools'
        ds_in['cTOT'].attrs['units'] = ds_in['clitmetb'].attrs['units']

        if CN == True:
            ds_in['nTOT'] = ds_in['nlitmetb']+ds_in['nlitstr']+ds_in['nsoilmic']+ \
                            ds_in['nsoilslow']+ds_in['nsoilpass']
            ds_in['nTOT'].attrs['long_name'] = 'sum of CASA SON pools'
            ds_in['nTOT'].attrs['units'] = ds_in['nlitmetb'].attrs['units']

            ds_in['cnTOT']= ds_in['cTOT'] / ds_in['nTOT']
            ds_in['cnTOT'].attrs['long_name'] = 'CASA total SOM C:N ratio'

    return ds_in


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
