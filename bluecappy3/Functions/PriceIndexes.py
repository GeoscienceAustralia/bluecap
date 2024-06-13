"""
Copyright (C) 2019, Monash University, Geoscience Australia
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
from scipy import interpolate  # for interp1d

#Units
from Units.UnitManager import UnitManager



class PriceIndex:

  def __init__(self,years,values, inflator):
    self.indexFunc = interpolate.interp1d(years,values )
    self.minYear = np.min(years)
    self.maxYear = np.max(years)
    self.inflator = inflator

  def IndexedPrice(self,origPrice,origYear,toYear):  

    if toYear < self.minYear:
      rv = self.indexFunc(self.minYear) * (1.0+self.inflator)**(toYear-self.minYear)
    elif  toYear > self.maxYear:
      rv = self.indexFunc(self.maxYear) * (1.0+self.inflator)**(toYear-self.maxYear)
    else:
      rv = self.indexFunc(toYear)
      
    if origYear < self.minYear:
      rv /= self.indexFunc(self.minYear) * (1.0+self.inflator)**(origYear-self.minYear)
    elif  origYear > self.maxYear:
      rv /= self.indexFunc(self.maxYear) * (1.0+self.inflator)**(origYear-self.maxYear)
    else:
      rv /= self.indexFunc(origYear)  
     
    rv = rv*origPrice
    
    return rv
  

MiningEquipmentPriceIndexData = np.array([
	[1990,51.875],
	[1991,53.925],
	[1992,54.3],
	[1993,55.325],
	[1994,55.525],
	[1995,56.65],
	[1996,58.05],
	[1997,58.4],
	[1998,57.9],
	[1999,58.525],
	[2000,63.025],
	[2001,64.775],
	[2002,65.575],
	[2003,66.225],
	[2004,68.75],
	[2005,75.65],
	[2006,81.575],
	[2007,83.925],
	[2008,91.375],
	[2009,89.675],
	[2010,91.425],
	[2011,98.125],
	[2012,102.5],
	[2013,105.25],
	[2014,105.525],
	[2015,99.95],
	[2016,99.225],
	[2017,102.025],
	[2018,105.8]])

MiningEquipmentPriceIndex = PriceIndex(MiningEquipmentPriceIndexData[:,0],MiningEquipmentPriceIndexData[:,1],0.037)
    
