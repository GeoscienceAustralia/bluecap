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


class PowerPlantManager():
    """
        This class holds information representing the power source for hydrogen production. 
        Note it is a subset of the Hydrogen data manager, which holds information about the hydrogen plant, power plant and infrastructure.
    """
    def __init__(self,plantName=""):
      """
      Create an empty power plant data manager and default variables. 
      """
      
      self.plantName = plantName
      self.plantCapacity = 0.0     # capacity of one component of the power plant 
                                   # (i.e. wind or solar *not* whole plant for hybrid plants)
            
      
      self.capacityFactor = 0.35    # capacity factor of whole plant (i.e. not just wind/solar for hybrid plants)
      
      self.waterPerMWh = 0.0 # mass water per MWh
      
      self.energyProduced = np.array([0.0])  # energy produced each year
      
      self.waterUse = np.array([0.0])
      
      self.capex = np.array([0.0])
      self.opex = np.array([0.0])
      
      self.projectLife = 0  # plant life + startup
      self.operatingLife = 25   # life in years before plant must be replaced
      
      self.startupTime = 1  # in years
      
      self.type = "Photovoltaic"
      self.useHybridPlantWithCurtailment = False  # wether to use a hybrid plant with curtailment
      
      self.powerFraction  =  1.0
      self.secondaryPowerSource = None
      

    def ParseXMLNode(self, powerPlantNode):
      """
      Generate power plant data from xml tree node.
       - Note this may not be needed in many cases as the mining system is derived from the orebody shape  
      """
    
      if(HasAttribute(powerPlantNode,"capacityFactor")):  
        self.capacityFactor = GetAttributeValue(powerPlantNode,"capacityFactor")
      
      if(HasAttribute(powerPlantNode,"waterPerMWh")):  
        self.waterPerMWh = GetAttributeValue(powerPlantNode,"waterPerMWh")
        
      if(HasAttribute(powerPlantNode,"operatingLife")):  
        self.operatingLife = GetAttributeValue(powerPlantNode,"operatingLife")
        if (self.operatingLife > 100):
          # assumed life is given in years not integer
          theUnitManager = UnitManager()
          self.operatingLife *= theUnitManager.ConvertTo("year") 
        self.operatingLife = int(np.round( self.operatingLife ) )
        
        
        
      if(HasAttribute(powerPlantNode,"type")):  
        self.type = GetAttributeString(powerPlantNode,"type")
    
      if(HasAttribute(powerPlantNode,"fraction")):  
        self.powerFraction = GetAttributeValue(powerPlantNode,"fraction")
        
      
      if(HasChild(powerPlantNode,"SecondarySystem")):  # and tertiary systems
        secondaryNode = GetChild(powerPlantNode,"SecondarySystem")
        self.secondaryPowerSource = PowerPlantManager()
        self.secondaryPowerSource.ParseXMLNode(secondaryNode) 
        
        if(HasAttribute(powerPlantNode,"useHybridPlantWithCurtailment")):
          self.useHybridPlantWithCurtailment = GetAttributeValue(powerPlantNode,"useHybridPlantWithCurtailment")
          self.secondaryPowerSource.useHybridPlantWithCurtailment = self.useHybridPlantWithCurtailment
      
      return powerPlantNode
    
    
      

    def WriteXMLNode(self, node):
      """
      Write power plant data to xml node
      """
      # plant data
      SetAttributeString(node,"capacityFactor",self.capacityFactor)
      SetAttributeString(node,"capacity",self.hydrogenProductionCapacity)
      SetAttributeString(node,"waterPerMWh",self.waterPerMWh)
      SetAttributeString(node,"operatingLife",self.operatingLife)
      SetAttributeString(node,"type",self.type)
      
      return node
  
    
    
    def SetStartupTimeAndProjectLife(self, startupTime, projectLife): 
      """Set the setup time and project duration."""
      
      self.startupTime = startupTime
      self.projectLife = projectLife
        
      if(self.secondaryPowerSource):
        self.secondaryPowerSource.SetStartupTimeAndProjectLife(startupTime, projectLife)
     
    def DeterminePowerSystem(self, problemManager, hydrogenManager):
      """
      Use available data to determine the most likely production method/capacity/opex etc.
      """
    
      #theUnitManager =  UnitManager()
      
      self.SetStartupTimeAndProjectLife(hydrogenManager.theHydrogenPlant.startupTime,  hydrogenManager.theHydrogenPlant.projectLife)
            
      self.DetermineProjectProductionAndCosts(problemManager, hydrogenManager)
      
      """
      self.CalculateProjectLife(hydrogenManager)
      self.CalculateHydrogenProduction(hydrogenManager)
      self.CalculateCapex(hydrogenManager)
      self.CalculateOpex(hydrogenManager)
      """
      return self

      
    def DetermineProjectProductionAndCosts(self, problemManager, hydrogenManager):
      """Calculate the amount of power required and the power plant startup and ongoing costs."""
      
      self.CalculatePlantCapacity(problemManager, hydrogenManager)
      
      self.CalculateCapex(hydrogenManager)
      self.CalculateOpex(hydrogenManager)
      
      if(hydrogenManager.theEnergyStorage):
        hydrogenManager.theEnergyStorage.CalculateCapex(hydrogenManager)
        hydrogenManager.theEnergyStorage.CalculateOpex(hydrogenManager)      
      
      return self
    
    
    def SetCapacityFactor(self, cf):
      """Set the power plant capacity factor."""
      self.capacityFactor =  cf
      
      if(self.secondaryPowerSource):
        self.secondaryPowerSource.SetCapacityFactor(cf)
    
      return cf
      
      
    def CalculatePlantCapacity(self,problemManager, hydrogenManager):
      """
      Determine the nameplate/installed capacity of the plant
      """
      
      theUnitManager = UnitManager()
      
      
      # energy use assumes 120 MJ/kg stored in hydrogen (equivalent to 33.33 kWh/kg)
      
      
      # without storage the plant must produce as much as max required
      self.energyProduced = np.array(hydrogenManager.theHydrogenPlant.energyUse)
      
      if(hydrogenManager.theEnergyStorage):
        hydrogenManager.theEnergyStorage.CalculateStorageCapacity(problemManager, hydrogenManager)
      
        # with storage we also need to cover the storage losses
        self.energyProduced += hydrogenManager.theEnergyStorage.energyLoss
        print("hydrogenManager.theEnergyStorage.energyLoss",hydrogenManager.theEnergyStorage.energyLoss)
        print("hydrogenManager.theHydrogenPlant.energyUse",hydrogenManager.theHydrogenPlant.energyUse)
      
      
      
      self.waterUse = self.energyProduced*theUnitManager.ConvertTo("MWh")*self.waterPerMWh
      
      liters = theUnitManager.ConvertToBaseUnits("L") 
      self.waterUse *= liters # converts to kL (ie to base unit m^3 from L)
      
      maxEnergyRequired = np.max(self.energyProduced)
      
      secondsPerYear = 31536000.
      averagePowerRequired = maxEnergyRequired/secondsPerYear
      
      totalPlantPower = averagePowerRequired/(self.capacityFactor + 1e-64)
      # nb self.capacityFactor is the weighted capacity factor of the whole plant (not just one component)
      
      if(self.useHybridPlantWithCurtailment):
        if(self.secondaryPowerSource):
           self.secondaryPowerSource.useHybridPlant = True  # so we don't forget to copy all the way down
        
        if( self.powerFraction > 0.5):
          self.plantCapacity = 1.0 *totalPlantPower  # we assume that the hydrogen plant is scaled to the same size as the largest power supply. 
        else:
          self.plantCapacity = totalPlantPower * (self.powerFraction)/(1.0-self.powerFraction+1e-64)
        
        rv= self.plantCapacity
        
        if(self.secondaryPowerSource):
          self.secondaryPowerSource.CalculatePlantCapacity(problemManager,hydrogenManager)
          
      else:
        self.plantCapacity = self.powerFraction*totalPlantPower  # capacity of this component of the plant
      
        rv= self.plantCapacity
      
        if(self.secondaryPowerSource):
          rv += self.secondaryPowerSource.CalculatePlantCapacity(problemManager,hydrogenManager)
      
      return rv




    def CalculateCapex(self, hydrogenManager):
      """Calculate the startup costs for the power plant"""

      self.capex = np.zeros(self.projectLife)
      
      theFunctionManager =  FunctionManager()
      theUnitManager =  UnitManager()
      
      capexFunc = theFunctionManager.GetFunction("EnergyCapex_" + self.type)
        
      self.capex[0] =  capexFunc.f( [self.plantCapacity*theUnitManager.ConvertTo("MW")] )
      
      #print "self.capex", self.capex
      
      rv = np.array(self.capex)
      
      #print "power capex", self.capex
      
      if(self.secondaryPowerSource):
        rv += self.secondaryPowerSource.CalculateCapex(hydrogenManager)
      
      return self.capex
      
    def CalculateOpex(self, hydrogenManager):
      """Calculate the ongoing costs for the power plant"""
    
      #  no distinction between opex and sustaining capital (adjust tax calculation)
      
      theFunctionManager =  FunctionManager()
      theUnitManager = UnitManager()
      
      opexFunc = theFunctionManager.GetFunction("EnergyOpex_" + self.type)
      
      opexPerAnnumPerMW =  opexFunc.f( [self.plantCapacity*theUnitManager.ConvertTo("MW")] )
      
      
      #print "Power plant type: ", self.type
      #print "Power plant capacity (MW)", self.plantCapacity  *  theUnitManager.ConvertTo("MW")
      #print "Capacity factor", self.capacityFactor
      #print "Power plant Opex Per Annum", opexPerAnnum
      
      self.opex = np.zeros(self.projectLife)
      self.opex[self.startupTime:] = opexPerAnnumPerMW  * self.capacityFactor * self.plantCapacity*theUnitManager.ConvertTo("MW")
      
      
      #print "self.plantCapacity",self.plantCapacity
      #print "power opex", self.opex
      
      rv = np.array(self.opex)
      
      if(self.secondaryPowerSource):
        rv += self.secondaryPowerSource.CalculateOpex(hydrogenManager)

      return rv

      
      
      
