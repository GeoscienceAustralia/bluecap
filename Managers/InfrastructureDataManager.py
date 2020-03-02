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


# Common
from Common.Common import Todo

# IO
from IO.XML import HasChild,GetChild,AddChild
from IO.XML import GetAttributeString,GetAttributeValue,GetAttributeValueOrDefault, SetAttributeString, GetAttributeStringOrDefault

# Units
from Units.UnitManager import UnitManager


# Functions
from Functions.FunctionManager import FunctionManager

# Price Indexes
from Functions.PriceIndexes import MiningEquipmentPriceIndex

class InfrastructureDataManager():
    def __init__(self):
      """
      Create an empty processing system data manager and default variables. 
      """
      self.distanceToRoad = 0.0
      self.distanceToRail = 0.0
      self.distanceToPower = 0.0
      self.distanceToWater = 0.0
      self.distanceToGas = 0.0
      
      self.roadTransportationDistance = 0.0
      self.railTransportationDistance = 0.0
      
      self.transportationCapex = np.array([0.0])
      self.transportationOpex = np.array([0.0])
      
      
      self.powerCapex = np.array([0.0])
      self.powerOpex = np.array([0.0])
      
      self.waterCapex = np.array([0.0])
      self.waterOpex = np.array([0.0])
      
      self.infrastructureCapex = np.array([0.0])
      self.infrastructureOpex = np.array([0.0])
      
      #self.powerSupply = "grid"
      self.calculateGas = False
      self.calculateDiesel = False
      

    def ParseXMLNode(self, infrastructureNode):
      """
      Generate processing system data from xml tree node. 
      """
      self.calculateDiesel = GetAttributeStringOrDefault(infrastructureNode,"calculateDiesel", self.calculateDiesel)
      self.calculateGas = GetAttributeValueOrDefault(infrastructureNode,"calculateGas", self.calculateGas)
      
      return infrastructureNode
        

    def WriteXMLNode(self, node):
      """
      Write processing system data to xml node
      """
      # Mine data
      SetAttributeString(node,"roadDist",self.distanceToRoad)
      SetAttributeString(node,"railDist",self.distanceToRail)
      SetAttributeString(node,"powerDist",self.distanceToPower)
      SetAttributeString(node,"waterDist",self.distanceToWater)
      
      if(self.calculateDiesel):
        SetAttributeString(node,"calculateDiesel",self.calculateDiesel)
      if(self.calculateGas):
        SetAttributeString(node,"calculateGas",self.calculateGas)
        SetAttributeString(node,"gasDist",self.distanceToGas)
      
      return node
      
    def DetermineDistanceToInfrastructure(self, problemManager, mineDataManager):
      """
      Use available data to determine the most likely processing method/capacity/opex etc.
      """
      theUnitManager = UnitManager()
      theFunctionManager = FunctionManager()
      
      oneKm = theUnitManager.ConvertToBaseUnits("km")  # assumed that distances are given in km
      self.distanceToRoad = theFunctionManager.GetFunction("DistanceToRoad").f( mineDataManager.mineLatLong  ) *oneKm
      self.distanceToRail = theFunctionManager.GetFunction("DistanceToRail").f( mineDataManager.mineLatLong  )   *oneKm
      self.distanceToWater = theFunctionManager.GetFunction("DistanceToWater").f( mineDataManager.mineLatLong  )  *oneKm
      self.distanceToPower = theFunctionManager.GetFunction("DistanceToPower").f( mineDataManager.mineLatLong  )  *oneKm
           
      self.roadTransportationDistance = theFunctionManager.GetFunction("RoadTransportationDistance").f( mineDataManager.mineLatLong  )  *oneKm
      self.railTransportationDistance = theFunctionManager.GetFunction("RailTransportationDistance").f( mineDataManager.mineLatLong  )  *oneKm 
    
      if(self.calculateGas): 
        self.distanceToGas = theFunctionManager.GetFunction("DistanceToGas").f( mineDataManager.mineLatLong  ) *oneKm
    
    def ZeroDistanceToInfrastructure(self):
      """
      Set distance to infrastructure to 0. 
      Used in regional calculation
      """
      self.distanceToRoad = 0.0
      self.distanceToRail = 0.0
      self.distanceToWater = 0.0
      self.distanceToPower = 0.0
      self.distanceToGas = 0.0
      
      self.roadTransportationDistance = 0.0
      self.railTransportationDistance = 0.0
      
    def CalculateInfrastructureExpenses(self, problemManager, mineDataManager):
      self.CalculateTransportationExpenses(problemManager, mineDataManager)
      self.CalculateWaterExpenses(problemManager, mineDataManager)
      self.CalculatePowerExpenses(problemManager, mineDataManager)
      
      self.infrastructureCapex =  self.transportationCapex +  self.waterCapex + self.powerCapex
      self.infrastructureOpex =  self.transportationOpex +  self.waterOpex + self.powerOpex
      
      
    def CalculateTransportationExpenses(self, problemManager, mineDataManager):
      
      numYears = mineDataManager.theMiningSystem.mineLife
      
      CovertToTodaysPrice = MiningEquipmentPriceIndex.IndexedPrice(1.0,2010,2018)
      
      # System
      #  Road       Rail
      #  Capital costs 1000AUD/km
      #  500-3000 2000-7000
      #  Fleet Costs 1000AUD /(km.Mt/a)
      # 40 100
      # Operating Costs cents/t/km
      #  5-12 1-1.25
      # Table 5: Transportation costs in 2010 AUD [1]
      
      theUnitManager = UnitManager()
      thousandAUDPerkm = 1000. * theUnitManager.ConvertToBaseUnits("AUD/km") # length prices are $1000 per km
      thousandAUDPerMtkmPerYear = 1000. * 1e-6 * theUnitManager.ConvertToBaseUnits("AUD/km/tonne")
      # fleet costs are per 1000 AuD per Mt km/ year
      
      
      fleetCapacity = np.max(mineDataManager.theProcessingSystem.concentrateProduced) 
      
      
      self.roadCapex = self.distanceToRoad * 1750. * CovertToTodaysPrice * thousandAUDPerkm \
                       + 40 * CovertToTodaysPrice * fleetCapacity * self.roadTransportationDistance * thousandAUDPerMtkmPerYear
      
      self.railCapex = self.distanceToRail * 4500. * CovertToTodaysPrice * thousandAUDPerkm \
                       + 100 * CovertToTodaysPrice * fleetCapacity * self.railTransportationDistance * thousandAUDPerMtkmPerYear
      
      self.roadOpex = np.zeros(numYears)  
      self.railOpex = np.zeros(numYears)
      
      centsPerTonnePerKm =  0.01*theUnitManager.ConvertToBaseUnits("AUD/tonne/km")
      self.roadOpex = self.roadTransportationDistance*  \
                          mineDataManager.theProcessingSystem.concentrateProduced * \
                          8.5* CovertToTodaysPrice* centsPerTonnePerKm# 5-12 cents/t/km
      
      self.railOpex = self.railTransportationDistance*  \
                          mineDataManager.theProcessingSystem.concentrateProduced * \
                          1.75* CovertToTodaysPrice*centsPerTonnePerKm # 1-2.5 cents/t/km
      
      # choose between road and rail based on which has lowest NPV (costs are +ve)
      roadNPV = mineDataManager.theEconomicDataManager.CalculateNPV(self.roadOpex) \
                + mineDataManager.theEconomicDataManager.CalculateNPV([self.roadCapex] )
      railNPV = mineDataManager.theEconomicDataManager.CalculateNPV(self.railOpex) \
                + mineDataManager.theEconomicDataManager.CalculateNPV([self.railCapex] )
      
      
      self.transportationCapex = np.zeros(numYears)  
      self.transportationOpex = np.zeros(numYears)
      
      useRoads = roadNPV < railNPV
      if(useRoads):
        self.transportationCapex[0] = self.roadCapex
        self.transportationOpex = self.roadOpex
      else:
        self.transportationCapex[0] = self.railCapex
        self.transportationOpex = self.railOpex
        
    def CalculateRegionalTransportationExpenses(self, distanceToRoad, roadTransportationDistance, \
                                                      distanceToRail, railTransportationDistance, \
                                                      coverMap,concentrateCapacityFunc,discountedTotalConcentrateMassFunc):
        
      AUD2010 = MiningEquipmentPriceIndex.IndexedPrice(1.0,2010,2018)
      
      # System
      #  Road       Rail
      #  Capital costs 1000AUD/km
      #  500-3000 2000-7000
      #  Fleet Costs 1000AUD /(km.Mt/a)
      # 40 100
      # Operating Costs cents/t/km
      #  5-12 1-1.25
      # Table 5: Transportation costs in 2010 AUD [1]
      
      theUnitManager = UnitManager()
      thousandAUDPerkm = 1000. * theUnitManager.ConvertToBaseUnits("AUD") # length prices are $1000 per km  - AUD is AUD/km as regional lengths given in km
      thousandAUDPerMtkmPerYear = 1000. * 1e-6 * theUnitManager.ConvertToBaseUnits("AUD/tonne")  # in AUD/km/tonne - regional lengths given in km
      # fleet costs are per 1000 AuD per Mt km/ year
      
      
      
      # road/rail capex
      roadNPVMap = distanceToRoad * 1750. *AUD2010* thousandAUDPerkm 
      railNPVMap = distanceToRail * 4500. *AUD2010* thousandAUDPerkm 

      # fleet capex
      fleetCapacityMap = concentrateCapacityFunc(coverMap)
      roadNPVMap += fleetCapacityMap * roadTransportationDistance * 40. *AUD2010*  thousandAUDPerMtkmPerYear
      railNPVMap += fleetCapacityMap * railTransportationDistance * 100. *AUD2010 *  thousandAUDPerMtkmPerYear
      
      # discounted continuing costs
      discountedTotalConcentrateMassMap = discountedTotalConcentrateMassFunc(coverMap)
      centsPerTonnePerKm =  0.01*theUnitManager.ConvertToBaseUnits("AUD/tonne")  # actually in AUD/tonne/Km but regional lengths given in km
      roadNPVMap += roadTransportationDistance* discountedTotalConcentrateMassMap * \
                          8.5*AUD2010* centsPerTonnePerKm# 5-12 cents/t/km
      
      railNPVMap += railTransportationDistance*  discountedTotalConcentrateMassMap * \
                          1.75*AUD2010*centsPerTonnePerKm # 1-2.5 cents/t/km
      
      transportationCostNPV = np.minimum(roadNPVMap,railNPVMap)
      
      return transportationCostNPV
        
                
       
    def CalculateWaterExpenses(self, problemManager, mineDataManager):
      """ 
      Calculate water startup infrastructure costs - ongoing costs are assumed included in the mine and processing cost models
      """
      numYears = mineDataManager.theMiningSystem.mineLife
      
      theUnitManager = UnitManager()
      processingCapacity = mineDataManager.theProcessingSystem.processingCapacity*  theUnitManager.ConvertTo("tonne")  # tonnes per year
      
      v = 2.0   # 2 m/s flow rate 
      q = 2.35* processingCapacity # in kL/year assumes 2.35 kL water per tonne 
      sPerYear = 365*24*3600
      diam = np.sqrt( 4* q/ ( np.pi * sPerYear * v) )
      
      distanceScale=  2.6229e6*diam - 125528   # empirical fit to cost of pipeline data per km in 2018 AUD
      
      oneKm = theUnitManager.ConvertToBaseUnits("km")
      
      self.waterCapex = np.zeros(numYears)  

      
      self.waterCapex[0] =  distanceScale*self.distanceToWater/oneKm
      
      #Water Opex not explicitly accounted for as included in mine and processing opex
      self.waterOpex = np.zeros(numYears)
      
    def CalculateRegionalWaterExpenses(self, problemManager,distanceToWater):
          
      theUnitManager = UnitManager()
      processingCapacity = problemManager.theMineDataManager.theProcessingSystem.processingCapacity*  theUnitManager.ConvertTo("tonne")  # tonnes per year
      
      v = 2.0   # 2 m/s flow rate 
      q = 2.35* processingCapacity # in kL/year assumes 2.35 kL water per tonne 
      sPerYear = 365*24*3600
      diam = np.sqrt( 4* q/ ( np.pi * sPerYear * v) )
      
      if diam < 150e-3:  # capped to prevent -ve pipeline costs (approximate diameter of "small" pipes in AUSIMM cost est)
        diam = 150e-3
      
      distanceScale=  2.6229e6*diam - 125528   # empirical fit to pipeline data per km
      
      rv =  distanceScale*distanceToWater
      return rv
  
      
    def CalculatePowerExpenses(self, problemManager, mineDataManager):
      numYears = mineDataManager.theMiningSystem.mineLife
      
      requiredVoltage = 220
      theUnitManager = UnitManager()
      oneKm = theUnitManager.ConvertToBaseUnits("km")
      if(self.distanceToPower < 190*oneKm):  
        requiredVoltage = 132      
        if(self.distanceToPower < 150*oneKm):   # this could also be 100 km
          requiredVoltage = 33      
          if(self.distanceToPower <= oneKm): 
            requiredVoltage = 11 
    
      
      # power capex
      self.powerCapex = np.zeros(numYears)  
      distanceScale = theUnitManager.ConvertToBaseUnits("AUD/km")
      
      AUD2010 = MiningEquipmentPriceIndex.IndexedPrice(1.0,2010,2018)
      AUD2014 = MiningEquipmentPriceIndex.IndexedPrice(1.0,2014,2018)
      
      if( requiredVoltage == 220):
        distanceScale *= 925000*AUD2014
      elif(requiredVoltage == 132):
        distanceScale *= 650000*AUD2014
      elif(requiredVoltage == 33):
        distanceScale *= 300000*AUD2010
      else:
        distanceScale *= 75000*AUD2010
      
      
      self.powerCapex[0] =  distanceScale*self.distanceToPower
      
      
      self.powerOpex = np.zeros(numYears)
      
      #Power opex is ignored as included in mine and processing opex
      

 
    def CalculateRegionalPowerExpenses(self, problemManager, distanceToPower):
      # distanceToPower is assumed given in km
            
      theUnitManager = UnitManager()
      powerCosts = np.zeros(distanceToPower.shape)
      
      AUD2010 = MiningEquipmentPriceIndex.IndexedPrice(1.0,2010,2018) * theUnitManager.ConvertToBaseUnits("AUD")  #in today's dollars
      AUD2014 = MiningEquipmentPriceIndex.IndexedPrice(1.0,2014,2018) * theUnitManager.ConvertToBaseUnits("AUD") 
      AUD2016 = MiningEquipmentPriceIndex.IndexedPrice(1.0,2016,2018) * theUnitManager.ConvertToBaseUnits("AUD") 
     
      
      
      powerCosts[distanceToPower > 190] =  925000*AUD2014
      
      powerCosts[ np.logical_and(distanceToPower <= 190,distanceToPower > 150)  ] = 300000*AUD2014
      
      powerCosts[ np.logical_and(distanceToPower <= 150,distanceToPower > 1)  ] = 300000*AUD2010
      
      powerCosts[distanceToPower <= 1] =  75000*AUD2010
      
      powerCosts *= distanceToPower
      
            
      if(self.calculateDiesel):
        AUD2018 = MiningEquipmentPriceIndex.IndexedPrice(1.0,2018,2018) * theUnitManager.ConvertToBaseUnits("AUD") 
        dieselLCOE = 350*AUD2016   # Renewable energy in Australian mining sector estimate
        retailCOE = 100*AUD2018
        dieselCosts = problemManager.theMineDataManager.theProcessingSystem.processingPower*( dieselLCOE - retailCOE)
        dieselNPC = problemManager.theMineDataManager.theEconomicDataManager.CalculateNPV(dieselCosts)
        powerCosts[ powerCosts > dieselNPC ] = dieselNPC
    
        
      
      return powerCosts
      
      
      #Power opex is ignored as included in mine and processing opex
      
      
       
    def CalculateRegionalGasExpenses(self,problemManager, distanceToGas):
      # distanceToGas is assumed given in km
            
      theUnitManager = UnitManager()
      gasCosts = np.zeros(distanceToGas.shape)
      
      AUD2016 = MiningEquipmentPriceIndex.IndexedPrice(1.0,2016,2018) * theUnitManager.ConvertToBaseUnits("AUD") 
      AUD2014 = MiningEquipmentPriceIndex.IndexedPrice(1.0,2014,2018) * theUnitManager.ConvertToBaseUnits("AUD") 
      
      
      piplelineCostPerKM = 0.315e6*AUD2014  # from core gas production and transmission costs 8 inch class 600 5.6mm wall.
      
      gasCosts = piplelineCostPerKM*distanceToGas
      
      AUD2018 = MiningEquipmentPriceIndex.IndexedPrice(1.0,2018,2018) * theUnitManager.ConvertToBaseUnits("AUD") 
      gasLCOE = 120*AUD2016   # LCOE per Mwh
      retailCOE = 100*AUD2018 # LCOE per Mwh
      
      gasPowerCosts = problemManager.theMineDataManager.theProcessingSystem.processingPower*( gasLCOE - retailCOE)
      gasNPC = problemManager.theMineDataManager.theEconomicDataManager.CalculateNPV(gasPowerCosts)
      gasCosts+= gasNPC
        
      
      return gasCosts
     
                