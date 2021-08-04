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


class HydrogenPlantManager():
    """
        This class holds information representing the hydrogen plant. 
        Note it is a subset of the Hydrogen data manager, which holds information about the hydrogen plant, power plant and infrastructure.
    """
    
    def __init__(self,plantName=""):
      """
      Create an empty hydrogen plant data manager and default variables. 
      """
      
      self.plantName = plantName
      self.projectLife = 0
      self.waterCapacity = 0.0
      self.hydrogenProductionCapacity = 0.0   # target output each year
      self.capacityFactor =  1.0
      self.hydrogenPlantCapacity = 0.0 # capacity of hydrogen plant running at 100% capacity
      self.hydrogenProduced = np.array([0.0])  # amount produced each year
      
      self.energyEfficiency = 0.75   # energy input required per unit energy stored
      
      self.waterPerKgH2 = 10.0 # mass water per mass H2 produced
      self.energyPerKgH2 = 0.0
      self.co2PerKgH2 = 0.0 # mass CO2 per mass H2
      self.coalPerKgH2 = 0.0 # mass coal per mass H2
      
      self.waterUse = np.array([0.0])
      self.energyUse = np.array([0.0])
      self.capex = np.array([0.0])
      self.opex = np.array([0.0])
      self.startupTime = 1 # in years
      self.actualStartupTime = 1.0
      self.type = "Electrolysis"
      
      
      #theUnitManager =  UnitManager()
      #self.ugHaulCostPerDepth = theUnitManager.ConvertToBaseUnits("0.0073 AUD/tonne/m")

    def ParseXMLNode(self, hydrogenPlantNode):
      """
      Generate hydrogen plant data from xml tree node.
       - Note this may not be needed in many cases as the mining system is derived from the orebody shape  
      """
    
      if(HasAttribute(hydrogenPlantNode,"energyEfficiency")):  
        self.energyEfficiency = GetAttributeValue(hydrogenPlantNode,"energyEfficiency")
      
      self.hydrogenProductionCapacity = GetAttributeValue(hydrogenPlantNode,"capacity")
      
      if(HasAttribute(hydrogenPlantNode,"waterH2Ratio")):  
        self.waterPerKgH2 = GetAttributeValue(hydrogenPlantNode,"waterH2Ratio")
        
      if(HasAttribute(hydrogenPlantNode,"CO2H2Ratio")):  
        self.co2PerKgH2 = GetAttributeValue(hydrogenPlantNode,"CO2H2Ratio")
        
      if(HasAttribute(hydrogenPlantNode,"type")):  
        self.type = GetAttributeString(hydrogenPlantNode,"type")
        
      
      return hydrogenPlantNode
    
    
      

    def WriteXMLNode(self, node):
      """
      Write hydrogen plant data to xml node
      """
      # plant data
      SetAttributeString(node,"energyEfficiency",self.energyEfficiency)
      SetAttributeString(node,"capacity",self.hydrogenProductionCapacity)
      SetAttributeString(node,"waterH2Ratio",self.waterPerKgH2)
      SetAttributeString(node,"type",self.type)
      
      return node
  
    
    
     
    def DetermineHydrogenSystem(self, problemManager, hydrogenManager):
      """
      Use available data to determine the most likely production method/capacity/opex etc.
      """
    
      theUnitManager =  UnitManager()
      
      self.CalculatePlantCapacity(problemManager, hydrogenManager)
      
      self.DetermineProjectProductionAndCosts(problemManager, hydrogenManager)
      
      return self

      
    def DetermineProjectProductionAndCosts(self, problemManager, hydrogenManager):
      """ Calculate the production type and associated costs for the hydrogen plant."""
      self.CalculateProjectLife(hydrogenManager)
      self.CalculateHydrogenProduced(hydrogenManager)
      
      self.CalculateCapex(hydrogenManager)
      self.CalculateOpex(hydrogenManager)
      
      return self
    
    def SetCapacityFactor(self, cf):
      """ Set the capacity factor for the hydrogen plant."""
      self.capacityFactor =  cf
      return cf
      
    def CalculatePlantCapacity(self,problemManager, hydrogenManager):
      """
      Determine the maximum production rate, water use and energy requirements
      """
      
      if( (self.type == "BrownCoalGasification") or (self.type == "BlackCoalGasification") or (self.type == "SteamMethane") ):
        self.hydrogenPlantCapacity = self.hydrogenProductionCapacity
      else:
        self.hydrogenPlantCapacity = self.hydrogenProductionCapacity/self.capacityFactor
      
      # energy use assumes 120 MJ/kg stored in hydrogen (equivalent to 33.33 kWh/kg)
      # and efficiency is measured relative to the 120 MJ/kg
      theUnitManager = UnitManager()
      hydrogenEnergyDensity = theUnitManager.ConvertToBaseUnits("120 MJ/kg")
      
      self.energyPerKgH2  =  hydrogenEnergyDensity/(self.energyEfficiency+ 1e-64)
      
      self.waterCapacity = self.hydrogenProductionCapacity * self.waterPerKgH2
      
      if(self.co2PerKgH2 == 0.0): # i.e. has not been set manually
        if self.type == "BrownCoalGasification": 
          self.co2PerKgH2 = 22.0  # 22 kg CO2 per kg H - approx based on NREL conversion rates estimates (22-25 kg/kg)
        elif self.type == "BlackCoalGasification":
          self.co2PerKgH2 = 22.0  # 22 kg CO2 per kg H - approx based on NREL conversion rates estimates (22-25 kg/kg)
        elif self.type == "SteamMethane":
          self.co2PerKgH2 = 10.26  # 10-11 kg CO2 per kg H 
          							  # Approx based on CSIRO conversion rates estimates (3.73 kg CH4/ kg H2 -> 3.73 x 44.01 /16 = 10.26 kg CO2 per kg H2)
      
      
      if self.type == "BrownCoalGasification": 
        if(self.coalPerKgH2 == 0.0): # i.e. has not been set manually
          self.coalPerKgH2 = self.co2PerKgH2 * 12.01/44.01     #  mass of carbon required based off amount of CO2 emitted - next we estimate the coal mass
          self.coalPerKgH2 /= 0.65 # assumes typical carbon content for brown coal in Australia (65%) 
      if self.type == "BlackCoalGasification": 
        if(self.coalPerKgH2 == 0.0): # i.e. has not been set manually
          self.coalPerKgH2 = self.co2PerKgH2 * 12.01/44.01     # mass of carbon required based off amount of CO2 emitted - next we estimate the coal mass
          self.coalPerKgH2 /= 0.80 # assumes typical carbon content for black coal in Australia (80%)
      
      return self.hydrogenProductionCapacity
 

    def CalculateProjectLife(self, hydrogenManager):
      """
      Set the life of the project equal to the life of the powerplant + startuptime. 
      """
      
      self.projectLife = hydrogenManager.thePowerPlant.operatingLife + self.startupTime
      self.projectLife = int( np.ceil(self.projectLife ))
      
      return self.projectLife


    def CalculateHydrogenProduced(self, hydrogenManager):
      """ Calculate hydrogen production over plant lifetime and estimate associated water and energy requirements."""
      
      self.hydrogenProduced = np.zeros(self.projectLife)
      
      theUnitManager = UnitManager()
      liters = theUnitManager.ConvertToBaseUnits("L")  # will convert water to kL (m^3)
      
      self.hydrogenProduced[int(self.startupTime):] = self.hydrogenProductionCapacity   #constant in all years
      self.energyUse = self.hydrogenProduced*self.energyPerKgH2
      self.waterUse = self.hydrogenProduced*liters*self.waterPerKgH2
      
      return self.hydrogenProduced


    def CalculateCapex(self, hydrogenManager):
      """ Calculate hydrogen plant startup costs. NB all costs are assumed capitalized."""

      self.capex = np.zeros(self.projectLife)
      
      theFunctionManager =  FunctionManager()
      theUnitManager =  UnitManager()
      
      capexFunc = theFunctionManager.GetFunction("HydrogenCapex_" + self.type)
        
      
      self.capex[0] =  capexFunc.f( [self.hydrogenPlantCapacity*theUnitManager.ConvertTo("1e6 tonne")] )
      
      
      return self.capex
      
    def CalculateOpex(self, hydrogenManager):
      """ Calculate hydrogen plant opex. NB no distinction between opex and sustaining capital."""
      
      theFunctionManager =  FunctionManager()
      theUnitManager = UnitManager()
      
      opexFunc = theFunctionManager.GetFunction("HydrogenOpex_" + self.type)
      
      opexPerTonne =  opexFunc.f( [self.hydrogenPlantCapacity *theUnitManager.ConvertTo("1e6 tonne")] )
      
      if(self.type == "SteamMethane" and hydrogenManager.theEconomicDataManager.calculateGasCosts ):
        AUD2018 = 1.0
        gasPrice = hydrogenManager.theEconomicDataManager.commodityPrices["gas"]*theUnitManager.ConvertTo("1/GJ")*AUD2018
        
        opexPerTonne += (gasPrice-8)*0.1726*1000  
        
        # 0.1726 = change in LCOH/kg per $AUD/GJ price for gas. Based on CSIRO estimate and $8/GJ gas price. 

      if(self.type == "BlackCoalGasification" and hydrogenManager.theEconomicDataManager.calculateBlackCoalCosts ):
        AUD2018 = 1.0
        coalPrice = hydrogenManager.theEconomicDataManager.commodityPrices["blackCoal"]*theUnitManager.ConvertTo("1/GJ")*AUD2018
        
        opexPerTonne += (coalPrice-3)*0.2248*1000 
        # 0.2248 = change in LCOH/kg per $AUD/GJ price for black coal. Based on CSIRO estimate and $3/GJ coal price. 

      if(self.type == "BrownCoalGasification" and hydrogenManager.theEconomicDataManager.calculateBrownCoalCosts ):
        # NB brown coal estimates use CSIRO estimated price for brown coal ($1.50/GJ) but estimate relative change based on black coal price due to a lack of data.  
        AUD2018 = 1.0
        coalPrice = hydrogenManager.theEconomicDataManager.commodityPrices["brownCoal"]*theUnitManager.ConvertTo("1/GJ")*AUD2018
        
        opexPerTonne += (coalPrice-1.5)*0.2248*1000 
        # 0.2248 = change in LCOH/kg per $AUD/GJ price for *black* coal (no equivalent information available for brown coal). Based on CSIRO estimate and $1.5/GJ brown coal price. 


    
      self.opex = opexPerTonne * self.hydrogenProduced  *  theUnitManager.ConvertTo("tonne")   

      return self.opex

      
      
      
