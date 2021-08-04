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



# Common
from Common.Common import BluecapError

# IO
from IO.XML import HasChild,GetChild,AddChild,GetChildren
from IO.XML import HasAttribute
from IO.XML import GetAttributeString,GetAttributeValue,GetAttributeValueOrDefault, SetAttributeString,GetXMLTag


# Units
from Units.UnitManager import UnitManager



# Managers

class OreBodyDataManager():
    def __init__(self,name=""):
      """
      Create an empty orebody manager and default variables. 
      """
      
      self.type = ""
      self.name = name
      self.dip = 0.0
      self.length = 0.0
      self.width = 0.0
      self.height = 0.0
      self.cover = 0.0
      
      self.metalGrades = {}
      self.specificDensity = 2.7 # specific density of the rock
      self.orebodyMass = 0.0
      self.orebodyVolume = 0.0
      self.shapeFactor = 1.0 # recovery based on shape of deposit
      
      self.latLong = np.array([0.0,0.0])

    def ParseXMLNode(self, orebodyDataNode):
      """
      Generate Ore body data from xml tree node. 
      """
      self.type = GetAttributeString(orebodyDataNode,"type")
      
      if(self.type == "K2SO4"):
        self.specificDensity = 1.0  # i.e. brine lake
        
      #self.grade = GetAttributeValue(orebodyDataNode,"grade")
      self.dip = GetAttributeValue(orebodyDataNode,"dip")
      
      self.cover = GetAttributeValue(orebodyDataNode,"cover")
      
      self.latLong[0] = GetAttributeValueOrDefault(orebodyDataNode,"lat",self.latLong[0])
      self.latLong[1] = GetAttributeValueOrDefault(orebodyDataNode,"long",self.latLong[1])
      
      if( HasAttribute(orebodyDataNode,"length") ):
        self.length =GetAttributeValue(orebodyDataNode,"length")
        self.width = GetAttributeValue(orebodyDataNode,"width")
        self.height = GetAttributeValue(orebodyDataNode,"height")
        if( self.width > self.length):
          temp = self.width 
          self.width = self.length
          self.length = temp
      elif( HasAttribute(orebodyDataNode,"mass")  ):
        self.orebodyMass = GetAttributeValue(orebodyDataNode,"mass")
        self.CalculateDepositVolume()
        self.CalculateDepositDimensionsFromVolume()
      else:
        BluecapError("Failed to find orebody mass or dimensions in input.")
      
      for child in GetChildren(orebodyDataNode):
          type = GetXMLTag(child)
          name = GetAttributeString(child,"name")
          grade = GetAttributeValue(child,"grade")
          
          self.metalGrades[name] = grade
      
    def WriteXMLNode(self, node):
      """
      Write ore body to xml node
      """
      SetAttributeString(node,"type",self.type)
      #SetAttributeString(node,"grade",self.grade)
      SetAttributeString(node,"dip",self.dip)
      
      SetAttributeString(node,"length",self.length)
      SetAttributeString(node,"width",self.width)
      SetAttributeString(node,"height",self.height)
      SetAttributeString(node,"cover",self.cover)
      
      # price data
      for name,grade in self.metalGrades.items():
        commNode = AddChild(node,"Commodity")
        SetAttributeString(commNode,"name",name)
        SetAttributeString(commNode,"grade",grade)
      

      return node
      
      
    def CalculateDepositMass(self):
      """
      Calculate mass of ore 
      """
      theUnitManager = UnitManager()
      waterDensity = theUnitManager.ConvertToBaseUnits("1000 kg/m^3") # water density in base units
      self.orebodyVolume = self.shapeFactor*self.length*self.width*self.height
      self.orebodyMass = self.orebodyVolume*self.specificDensity * waterDensity
      
      # print "orebody mass in kg", self.orebodyMass # in kg
      return self.orebodyMass
     
    def CalculateDepositVolume(self):
      """
      Calculate volume of ore given the mass
      """
      theUnitManager = UnitManager()
      waterDensity = theUnitManager.ConvertToBaseUnits("1000 kg/m^3") # water density in base units
      self.orebodyVolume = self.orebodyMass/(self.specificDensity * waterDensity)
      
    def CalculateDepositDimensionsFromVolume(self):
      """
      Calculate dimensions assuming a "cubic" deposit (NB may be slanted)
      """
      if(self.type == "K2SO4"):
        # Assume Potash deposits are tabular
        self.height = 10
        self.length = (self.orebodyVolume/(self.shapeFactor*self.height))**(1.0/2.0)
        self.width = self.length
        assert (self.orebodyVolume/self.shapeFactor) == (self.length * self.width * self.height)
      else:
        self.length = (self.orebodyVolume/self.shapeFactor)**(1.0/3.0)
        self.width = self.length
        self.height = self.length
     
    def ScaleCommodityGrades(self, factor):
      """
      Scale ore grades by a constant factor (used for breakeven analysis)
      """
      for key,value in  self.metalGrades.items():
        self.metalGrades[key] = value*factor