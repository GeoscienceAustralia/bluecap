"""
Copyright (C) 2019-2021, Monash University, Geoscience Australia
Copyright (C) 2018, Stuart Walsh 

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

import numpy as np

# Check GeoTIFF functionality 
gdal_flag = False
try:
    from osgeo import gdal,osr
    gdal_flag = True
except ImportError:
    pass

def LatLonDegreeToKMCoords(lat,lon):
      avlat = lat*np.pi/180.
      yCoord = 1e-3*( 111132.954*(lat+28) - 559.822 * np.sin( 2 * avlat )*180/(np.pi*2.0) + 1.175 * np.sin( 4 * avlat)*180/(np.pi*4.0) - 0.0023*np.sin( 6 * avlat)*180/(np.pi*6.0) );
      xCoord = 1e-3*( (111412.84 * np.cos( avlat ) - 93.5 * np.cos(3*avlat) + 0.118 * np.cos(5*avlat) ) * (lon - 134 ) )
      return yCoord,xCoord

def UsingGDAL():
  return gdal_flag