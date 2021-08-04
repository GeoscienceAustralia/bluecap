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

from scipy.optimize import brentq

# IO
from IO.XML import HasChild,GetChild,AddChild
from IO.XML import HasAttribute
from IO.XML import GetAttributeString,GetAttributeValue, SetAttributeString


#Units
from Units.UnitManager import UnitManager

#Functions
from Functions.FunctionManager import FunctionManager


class HydrogenEnergyStorageManager():
    """
        This class holds information representing the hydrogen plant. 
        Note it is a subset of the Hydrogen data manager, which holds information about the hydrogen plant, power plant and infrastructure.
    """
    def __init__(self,plantName=""):
      """
      Create an empty energy storage manager and default variables. 
      """
      
      self.plantName = plantName
      
      self.netOutputCapacityFactor = 1.0
      self.storageOutputCapacityFactor = 0.0
      self.storageOutputCapacity = 0.0 
      self.energyStored = np.array([0.0]) 
      self.energyLoss = np.array([0.0]) 
      
      self.roundTripEfficiency = 0.80   # energy output per unit energy stored
      
      self.capex = np.array([0.0])
      self.opex = np.array([0.0])
      self.startupTime = 1 # in years
      self.actualStartupTime = 1.0
      self.type = "PumpedHydro"
      
      
      #theUnitManager =  UnitManager()
      #self.ugHaulCostPerDepth = theUnitManager.ConvertToBaseUnits("0.0073 AUD/tonne/m")

    def ParseXMLNode(self, storageNode):
      """
      Generate energy storage data from xml tree node.
      """
    
      if(HasAttribute(storageNode,"roundTripEfficiency")):  
        self.roundTripEfficiency = GetAttributeValue(storageNode,"roundTripEfficiency")
      
      self.netOutputCapacity = GetAttributeValue(storageNode,"netOutputCapacity")
      
      if(HasAttribute(storageNode,"type")):  
        self.type = GetAttributeString(storageNode,"type")
        
      return storageNode
    
    
    def WriteXMLNode(self, node):
      """
      Write energy storage data to xml node
      """
      # plant data
      SetAttributeString(node,"roundTripEfficiency",self.roundTripEfficiency)
      SetAttributeString(node,"netOutputCapacity",self.netOutputCapacity)
      SetAttributeString(node,"roundTripEfficiency",self.roundTripEfficiency)
      SetAttributeString(node,"type",self.type)
      
      return node
  

    def CalculateStorageCapacity(self, problemManager, hydrogenManager):
      """Calculate energy storage requirements to meet desired net output capacity."""
      
      powerPlantCF = hydrogenManager.thePowerPlant.capacityFactor
      
      # if energy can be supplied without being stored it will be
      self.storageOutputCapacityFactor = np.maximum(self.netOutputCapacityFactor - powerPlantCF,0.0)
      
      # stored energy (based on amount provided to the hydrogen plant)
      print("self.storageOutputCapacityFactor",self.storageOutputCapacityFactor)
      self.energyStored = self.storageOutputCapacityFactor * hydrogenManager.theHydrogenPlant.energyUse
      
      # amount lost in storage 
      self.energyLoss = self.energyStored*(1.0/self.roundTripEfficiency - 1.0) 
      
      maxEnergyRequired = np.max(self.energyStored)
      
      secondsPerYear = 31536000.
      averagePowerRequired = maxEnergyRequired/secondsPerYear
      
      totalPlantPower = averagePowerRequired/(self.storageOutputCapacityFactor + 1e-64)
      
      self.storageOutputCapacity = totalPlantPower
      
      return self.energyStored
      

    def CalculateCapex(self, hydrogenManager):
      """Calculate capex costs for energy storage facility."""

      self.capex = np.zeros(hydrogenManager.thePowerPlant.projectLife)
      
      theFunctionManager =  FunctionManager()
      theUnitManager =  UnitManager()
      
      capexFunc = theFunctionManager.GetFunction("StorageCapex_" + self.type)
        
      # Capex as a function of output capacity in MW
      self.capex[0] =  capexFunc.f( [self.storageOutputCapacity*theUnitManager.ConvertTo("MW")] )
      
      #print("self.capex", self.capex)
      
      return self.capex
      
      
    def CalculateOpex(self, hydrogenManager):
      """Calculate opex costs for energy storage facility."""
      
      theFunctionManager =  FunctionManager()
      theUnitManager = UnitManager()
      
      opexFunc = theFunctionManager.GetFunction("StorageOpex_" + self.type)
      
      # opex Per MWh stored as a function of output capacity in MW
      opexPerMWh =  opexFunc.f( [self.storageOutputCapacity *theUnitManager.ConvertTo("MW")] )
      self.opex = opexPerMWh*self.energyStored*theUnitManager.ConvertTo("MWh")
      
      return self.opex
