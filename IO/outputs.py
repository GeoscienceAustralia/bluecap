"""
Copyright (C) 2020, Monash University, Geoscience Australia
Copyright (C) 2020, Marcus Haynes

Bluecap is released under the Apache License, Version 2.0 (the "License");
you may not use this software except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

The project uses third party components which may have different licenses. 
Please refer to individual components for more details.

"""

import os
from osgeo import gdal,ogr, osr
import numpy as np

def write_point_geojson(filename, layername, pts, dicts, XY=True, in_epsg=4326, out_epsg=4326):
    """
    """
    # create the spatial reference
    target = osr.SpatialReference()
    target.ImportFromEPSG(out_epsg)
    transform = None
    
    # check input spatial reference and reproject, if neccessary
    if (in_epsg != out_epsg):
        source = osr.SpatialReference()
        source.ImportFromEPSG(in_epsg)
        transform = osr.CoordinateTransformation(source, target)
    
    # Create the output Driver
    outDriver = ogr.GetDriverByName('GeoJSON')
    
    # Create the output GeoJSON
    outDataSource = outDriver.CreateDataSource(filename+'.geojson')
    outLayer = outDataSource.CreateLayer(layername, target, geom_type=ogr.wkbMultiPoint, options=["RFC7946=YES"])
    feilds = list(dicts[0].keys())
    for f in feilds:
        fldDef = ogr.FieldDefn(f, ogr.OFTString)
        fldDef.SetWidth(16) #16 char string width
        outLayer.CreateField(fldDef)
    
    # Get the output Layer's Feature Definition
    featureDefn = outLayer.GetLayerDefn()
    
    for i in range(len(dicts)):
        point_geom = ogr.Geometry(ogr.wkbMultiPoint)
        for pt in pts[i]:
            point = ogr.Geometry(ogr.wkbPoint)
            if XY:
                point.AddPoint(pt[0], pt[1])
            else:
                point.AddPoint(pt[1], pt[0])
            if transform:
                point.Transform(transform)
            point_geom.AddGeometry(point)
        
        # create a new feature
        outFeature = ogr.Feature(featureDefn)
        
        # Set new geometry
        outFeature.SetGeometry(point_geom)
        for f in feilds:
            outFeature.SetField(f,dicts[i][f])
        
        # Add new feature to output Layer
        outLayer.CreateFeature(outFeature)
    
        # dereference the feature
        outFeature = None
        point_geom = None
    
    # Save and close DataSources
    outDataSource = None


def array_to_raster(name, array, xpix=0.01, ypix=0.01, xmin=111.995, ymax=-9.995, epsg=4326):
    """Array > Raster
    Save a raster from an array.

    :name array:  str
    :param array: ndarray
    :xpix array:  float
    :ypix array:  float
    :xmin array:  float
    :ymax array:  float
    :epsg array:  int
    """
    
    if len(array.shape) == 2:
      ny,nx = array.shape
      nz = 1
      array = array.reshape((ny,nx,nz))
    else:
      ny,nx,nz = array.shape
    
    proj = osr.SpatialReference()
    proj.ImportFromEPSG(epsg)
    driver = gdal.GetDriverByName('GTiff')
    
    dataset = driver.Create(
        name,
        nx,
        ny,
        nz,
        gdal.GDT_Float32 )

    dataset.SetGeoTransform((
        xmin,
        xpix,
        0,
        ymax,
        0,
        -ypix))

    dataset.SetProjection(proj.ExportToWkt())
    for i in range(nz):
      dataset.GetRasterBand(i+1).WriteArray(array[:,:,i])
      dataset.GetRasterBand(i+1).SetNoDataValue(-99999)
    dataset.FlushCache()  # Write to disk.
    dataset = None

def GetMapOutputStyle(data,type,style="GA"):
  """Sets the output conventions for map data. Returns the colormap, max and min data range and max and min labels."""

  data_min = np.percentile(data[np.isnan(data)==False], 5)
  data_max = np.percentile(data[np.isnan(data)==False], 95)
  if np.isnan(data_min): data_min = 0
  if np.isnan(data_max): data_max = 1
  if data_min == data_max:
    data_min = np.percentile(data[np.isnan(data)==False], 0)
    data_max = np.percentile(data[np.isnan(data)==False], 100)
  label_min = data_min
  label_max = data_max
  cmap = 'viridis'
  
  if(style == "GA"):
    if type == "NPV" or (type == "uncertainNPV"):
      if abs(data_min) > data_max:
        data_max = -1.*data_min
      if data_max > abs(data_min):
        data_min = -1.*data_max
      cmap = 'seismic'
      label_min = data_min * 1e-06
      label_max = data_max * 1e-06
    elif type == "benefit_cost_ratio":
      cmap = 'PiYG'
      #data_min = np.log(data_min)  # log should be dealt with externally - don't want to change data
      #data_max = np.log(data_max)
      if abs(data_min) > data_max:
        data_max = -1.*data_min
      if data_max > abs(data_min):
        data_min = -1.*data_max
      label_min = np.exp(data_min)
      label_max = np.exp(data_max)
      #data = np.log(data)
    elif type == "breakeven_grade":
      cmap = 'magma'
      label_min = data_min
      label_max = data_max
    elif type == "employment":
      data_min = int(data_min)
      data_max = int(data_max)+1
      cmap = 'viridis'
      label_min = data_min
      label_max = data_max
    elif type == "tax":
      cmap = 'inferno'
      label_min = data_min * 1e-06
      label_max = data_max * 1e-06
    elif type == "param1":
      cmap = 'Reds'
    elif type == "param2":
      cmap = 'Greens'
    elif type == "param3":
      cmap = 'Blues'
    elif type == "water":
      cmap = 'Blues'
      data_min = 0
      label_min = data_min * 1e-06
      label_max = data_max * 1e-06
    elif type == "HYDROGEN":
      if abs(data_min) > data_max:
        data_max = -1.*data_min
      if data_max > abs(data_min):
        data_min = -1.*data_max
      cmap = 'seismic'
      label_min = data_min * 1e-06
      label_max = data_max * 1e-06 
    elif type == "HYDROGEN CSP": # Distinct CSP colour map
      values_above_zeros = (x for x in data[np.isnan(data)==False] if x>-350000000000) #this removes values from colour bar below -350k.
      data_min = min(values_above_zeros)
      if abs(data_min) > data_max:
        data_max = -1.*data_min
      if data_max > abs(data_min):
        data_min = -1.*data_max
      cmap = 'seismic'    
      label_min = data_min * 1e-06
      label_max = data_max * 1e-06
    else:
      cmap = 'viridis'
  return cmap,data_max, data_min, label_min, label_max