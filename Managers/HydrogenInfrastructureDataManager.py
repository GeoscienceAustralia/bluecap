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
from scipy import interpolate



# IO
from IO.XML import HasChild,GetChild,AddChild
from IO.XML import GetAttributeString,GetAttributeValue,GetAttributeValueOrDefault, SetAttributeString, GetAttributeStringOrDefault

# Units
from Units.UnitManager import UnitManager


# Functions
from Functions.FunctionManager import FunctionManager
from Functions.BetaDistribution import ScaledBetaDistribution

# Price Indexes
from Functions.PriceIndexes import MiningEquipmentPriceIndex

from .InfrastructureDataManager import InfrastructureDataManager

class HydrogenInfrastructureDataManager(InfrastructureDataManager):
    def __init__(self):
      """
      Create an empty hydrogen infrastructure data manager and default variables. 
      """
      InfrastructureDataManager.__init__(self)
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
      self.totalWaterUse = np.array([0.0]) # total Water Use each year
      
      self.infrastructureCapex = np.array([0.0])
      self.infrastructureOpex = np.array([0.0])
      
      
      # costs to include
      self.includeTransportationCosts = True    # connect plant to the grid
      self.includePowerInfrastructureCosts = True  # connect plant to the grid
      
      

    def ParseXMLNode(self, infrastructureNode):
      """
      Generate infrastructure data from xml tree node. 
      """
      
      self.includeTransportationCosts = GetAttributeValueOrDefault(infrastructureNode,"includeTransportationCosts",self.includeTransportationCosts )
      self.includePowerInfrastructureCosts = GetAttributeValueOrDefault(infrastructureNode,"includePowerInfrastructureCosts",self.includePowerInfrastructureCosts )
      
      return infrastructureNode
        

    def WriteXMLNode(self, node):
      """
      Write infrastructure data to xml node.
      """
      # project data
      SetAttributeString(node,"roadDist",self.distanceToRoad)
      SetAttributeString(node,"railDist",self.distanceToRail)
      SetAttributeString(node,"powerDist",self.distanceToPower)
      SetAttributeString(node,"waterDist",self.distanceToWater)
      
      return node
      
    def DetermineDistanceToInfrastructure(self, problemManager, hydrogenDataManager):
      """
      Use available data to determine distance to infrastructure (for single site calculations).
      """
      theUnitManager = UnitManager()
      theFunctionManager = FunctionManager()
      
      oneKm = theUnitManager.ConvertToBaseUnits("km")  # assumed that distances are given in km
      self.distanceToRoad = theFunctionManager.GetFunction("DistanceToRoad").f( hydrogenDataManager.plantLatLong[::-1]  ) *oneKm
      self.distanceToRail = theFunctionManager.GetFunction("DistanceToRail").f( hydrogenDataManager.plantLatLong[::-1]   )   *oneKm
      self.distanceToWater = theFunctionManager.GetFunction("DistanceToWater").f( hydrogenDataManager.plantLatLong[::-1]   )  *oneKm
      self.distanceToPower = theFunctionManager.GetFunction("DistanceToPower").f( hydrogenDataManager.plantLatLong[::-1]   )  *oneKm
           
      self.roadTransportationDistance = theFunctionManager.GetFunction("RoadTransportationDistance").f( hydrogenDataManager.plantLatLong[::-1]   )  *oneKm
      self.railTransportationDistance = theFunctionManager.GetFunction("RailTransportationDistance").f( hydrogenDataManager.plantLatLong[::-1]   )  *oneKm 
    
      #if(self.calculateGas): 
      #self.distanceToGas = theFunctionManager.GetFunction("DistanceToGas").f( mineDataManager.plantLatLong[::-1]   ) *oneKm
    
    def ZeroDistanceToInfrastructure(self):
      """
      Set distance to infrastructure to 0. Used in regional calculation.
      """
      self.distanceToRoad = 0.0
      self.distanceToRail = 0.0
      self.distanceToWater = 0.0
      self.distanceToPower = 0.0
      self.distanceToGas = 0.0
      
      self.roadTransportationDistance = 0.0
      self.railTransportationDistance = 0.0
      
    def CalculateInfrastructureExpenses(self, problemManager, hydrogenDataManager):
      """Calculate infrastructure costs."""
      self.CalculateTransportationExpenses(problemManager, hydrogenDataManager)
      self.CalculateWaterExpenses(problemManager, hydrogenDataManager)
      self.CalculatePowerExpenses(problemManager, hydrogenDataManager)
      
      self.infrastructureCapex =  self.transportationCapex +  self.waterCapex + self.powerCapex
      self.infrastructureOpex =  self.transportationOpex +  self.waterOpex + self.powerOpex
      
    def ZeroInfrastructureExpenses(self, problemManager, hydrogenDataManager):
      """Set infrastructure costs to 0. Used in regional calculation."""
    
      numYears = hydrogenDataManager.GetDuration()
      
      self.transportationCapex = np.zeros(numYears)  
      self.waterCapex = np.zeros(numYears)   
      self.powerCapex = np.zeros(numYears) 
      
      self.transportationOpex = np.zeros(numYears) 
      self.waterOpex = np.zeros(numYears) 
      self.powerOpex = np.zeros(numYears) 
      
      self.infrastructureCapex = np.zeros(numYears) 
      self.infrastructureOpex = np.zeros(numYears) 
 
      
    def CalculateTransportationExpenses(self, problemManager, hydrogenManager):
      """Calculate transporation costs (Single site calculation)."""
      numYears = hydrogenManager.theHydrogenPlant.projectLife
      
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
      
      
      fleetCapacity = np.max(hydrogenManager.theHydrogenPlant.hydrogenProduced) 
      
      if self.includeTransportationCosts: 
        
        self.roadCapex = self.distanceToRoad * 1750. * CovertToTodaysPrice * thousandAUDPerkm \
                       + 40 * CovertToTodaysPrice * fleetCapacity * self.roadTransportationDistance * thousandAUDPerMtkmPerYear
      
        self.railCapex = self.distanceToRail * 4500. * CovertToTodaysPrice * thousandAUDPerkm \
                       + 100 * CovertToTodaysPrice * fleetCapacity * self.railTransportationDistance * thousandAUDPerMtkmPerYear
      
      else:
        # still need road/rail infrastructure (to build and run the plant)
        self.roadCapex = self.distanceToRoad * 1750. * CovertToTodaysPrice * thousandAUDPerkm 
      
        self.railCapex = self.distanceToRail * 4500. * CovertToTodaysPrice * thousandAUDPerkm 
      
      
      self.roadOpex = np.zeros(numYears)  
      self.railOpex = np.zeros(numYears)
      
      
      if self.includeTransportationCosts: 
          centsPerTonnePerKm =  0.01*theUnitManager.ConvertToBaseUnits("AUD/tonne/km")
          self.roadOpex = self.roadTransportationDistance*  \
                              hydrogenManager.theHydrogenPlant.hydrogenProduced * \
                              8.5* CovertToTodaysPrice* centsPerTonnePerKm # 5-12 cents/t/km
      
          self.railOpex = self.railTransportationDistance*  \
                              hydrogenManager.theHydrogenPlant.hydrogenProduced * \
                              1.75* CovertToTodaysPrice*centsPerTonnePerKm # 1-2.5 cents/t/km
      
      # choose between road and rail based on which has lowest NPV (costs are +ve)
      roadNPV = hydrogenManager.theEconomicDataManager.CalculateNPV(self.roadOpex) \
                + hydrogenManager.theEconomicDataManager.CalculateNPV([self.roadCapex] )
      railNPV = hydrogenManager.theEconomicDataManager.CalculateNPV(self.railOpex) \
                + hydrogenManager.theEconomicDataManager.CalculateNPV([self.railCapex] )
      
      
      self.transportationCapex = np.zeros(numYears)  
      self.transportationOpex = np.zeros(numYears)
      
      useRoads = roadNPV < railNPV
      if(useRoads):
        self.transportationCapex[0] = self.roadCapex
        self.transportationOpex = self.roadOpex
      else:
        self.transportationCapex[0] = self.railCapex
        self.transportationOpex = self.railOpex
    

    def RoadTransportOpex(self, distance, mass):
      CovertTo2018Price = MiningEquipmentPriceIndex.IndexedPrice(1.0,2010,2018)
      
      theUnitManager = UnitManager()
      
      # 8.5 cents per km
      centsPerTonnePerKm =  0.01*theUnitManager.ConvertToBaseUnits("AUD/tonne/km")
      cost = distance* mass * \
                          8.5* CovertTo2018Price* centsPerTonnePerKm
      return cost
    
    def RoadConstructionCapex(self,distance):
      # nb fleet costs are also required
      CovertTo2018Price = MiningEquipmentPriceIndex.IndexedPrice(1.0,2010,2018)
      
      theUnitManager = UnitManager()
      thousandAUDPerkm = 1000. * theUnitManager.ConvertToBaseUnits("AUD/km")
      cost = distance * 1750. * CovertTo2018Price * thousandAUDPerkm
      return cost
        
    def CalculateRegionalTransportationExpenses(self, distanceToRoad, roadTransportationDistance, \
                                                      distanceToRail, railTransportationDistance, \
                                                      coverMap,hydrogenPlantCapacity,discountedTotalConcentrateMassFunc,
                                                      pipelineDistanceToPort =  None):
      """Calculate regional transporation costs (road/rail/pipeline). Note that all projects must have either road or rail access."""
      
      AUD2010 = MiningEquipmentPriceIndex.IndexedPrice(1.0,2010,2018)
      AUD2016 = MiningEquipmentPriceIndex.IndexedPrice(1.0,2016,2018)
      
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
      
      
      
      # road/rail capex - needed just to get road/rail to plant location
      roadNPVMap = distanceToRoad * 1750. *AUD2010* thousandAUDPerkm   # 500 - 3000 
      railNPVMap = distanceToRail * 4500. *AUD2010* thousandAUDPerkm   # 2000 - 7000  # track
      
      pipelineNPVMap =  None
      # transportation costs 
      if self.includeTransportationCosts: 
      
        # pipeline costs must include costs to connect to transportation network - it is assumed the cheapest option (road or rail is used)
        
        if(pipelineDistanceToPort is not None):
          pipelineNPVMap = np.minimum(roadNPVMap,railNPVMap) # still need to build road/rail to the site
          
          #AUD2014 = MiningEquipmentPriceIndex.IndexedPrice(1.0,2014,2018) * theUnitManager.ConvertToBaseUnits("AUD") 
          
          #piplelineCostPerKM = 0.315e6*AUD2014  # from core gas production and transmission costs 8 inch class 600 5.6mm wall.
          #piplelineCostPerKM *= 1.68 # NIST estimates that hydrogen-grade steel pipeline is 68% more costly than standard gas pipeline. 
          
          # From "Evaluation of the Economics of Renewable Hydrogen Supply in the APEC Region", Kan and Shibata, 2018
          AUD2018 = 1.0 #MiningEquipmentPriceIndex.IndexedPrice(1.0,2018,2018)
          netPiplelineCostPerKM = 1.475e6*AUD2018     
          # $0.4 MM USD/km (pipeline) + $50 MM USD /79km (compressors) 
          
          baseCapacity = 194070* 365 * theUnitManager.ConvertToBaseUnits("kg")  
          # estimates based on 194070 kg/day ~ 70 k tonne/year plant
          # assuming linear scaling 
          
          numPipeline = hydrogenPlantCapacity/baseCapacity
          if(numPipeline < 1.0): 
            numPipeline = 1.0
          
          pipelineNPVMap += numPipeline* pipelineDistanceToPort * netPiplelineCostPerKM # pipeline costs
          
          # continuing costs 8% of capex per annum
          opexPerkgkm = 0.08*netPiplelineCostPerKM/baseCapacity
          discountedTotalConcentrateMassMap = discountedTotalConcentrateMassFunc(coverMap)
          pipelineNPVMap += opexPerkgkm* pipelineDistanceToPort * discountedTotalConcentrateMassMap
          
      
        """
        # fleet capex
        roadNPVMap += hydrogenPlantCapacity * roadTransportationDistance * 40. *AUD2010*  thousandAUDPerMtkmPerYear
        railNPVMap += hydrogenPlantCapacity * railTransportationDistance * 100. *AUD2010 *  thousandAUDPerMtkmPerYear
      
        # discounted continuing costs
        discountedTotalConcentrateMassMap = discountedTotalConcentrateMassFunc(coverMap)
        centsPerTonnePerKm =  0.01*theUnitManager.ConvertToBaseUnits("AUD/tonne")  # actually in AUD/tonne/Km but regional lengths given in km
        roadNPVMap += roadTransportationDistance* discountedTotalConcentrateMassMap * \
                          8.5*AUD2010* centsPerTonnePerKm# 5-12 cents/t/km
      
        railNPVMap += railTransportationDistance*  discountedTotalConcentrateMassMap * \
                          1.75*AUD2010*centsPerTonnePerKm # 1-2.5 cents/t/km
        """
        
        # from CSIRO report
        # discounted continuing costs - fleet costs are included in transport costs
        discountedTotalConcentrateMassMap = discountedTotalConcentrateMassFunc(coverMap)
        centsPerTonnePerKm =  0.01*theUnitManager.ConvertToBaseUnits("AUD/tonne")  # actually in AUD/tonne/Km but regional lengths given in km
        roadNPVMap += roadTransportationDistance* discountedTotalConcentrateMassMap * \
                          233*AUD2016* centsPerTonnePerKm # $2.33 AUD/t/km
      
        railNPVMap += railTransportationDistance*  discountedTotalConcentrateMassMap * \
                          55*AUD2016*centsPerTonnePerKm # $0.55 AUD/t/km
        
        # print "here", np.max(railTransportationDistance)
        
      transportationCostNPV = np.minimum(roadNPVMap,railNPVMap) 
      if(pipelineNPVMap is not None):
        transportationCostNPV = np.minimum(transportationCostNPV,pipelineNPVMap)
        
      else:
        transportationCostNPV = np.minimum(roadNPVMap,railNPVMap)
      
      return transportationCostNPV


    def CalculateUncertainRegionalTransportationExpenses(self, distanceToRoad, roadTransportationDistance, \
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
      roadCapexPerKm = ScaledBetaDistribution.FromMeanAndStd(1750.,min=500.,max = 3000.)
      railCapexPerKm = ScaledBetaDistribution.FromMeanAndStd(4500.,min=2000.,max = 7000.)
      roadNPVMap = distanceToRoad * (roadCapexPerKm *AUD2010* thousandAUDPerkm)   # 500 - 3000 
      railNPVMap = distanceToRail * (railCapexPerKm *AUD2010* thousandAUDPerkm)   # 2000 - 7000  # track
      
      
      # fleet capex
      roadFleetCapexPerMtKmYr = ScaledBetaDistribution.FromMeanAndStd(40.,min=28.,max = 52.)  # assuming 10% std dev
      railFleetCapexPerMtKmYr = ScaledBetaDistribution.FromMeanAndStd(100.,min=70.,max = 130.) # assuming 10% std dev
      fleetCapacityMap = concentrateCapacityFunc(coverMap)
      roadNPVMap += fleetCapacityMap * (roadTransportationDistance * roadFleetCapexPerMtKmYr *AUD2010*  thousandAUDPerMtkmPerYear)
      railNPVMap += fleetCapacityMap * (railTransportationDistance * railFleetCapexPerMtKmYr *AUD2010 *  thousandAUDPerMtkmPerYear)
      
      
      # discounted continuing costs
      discountedTotalConcentrateMassMap = discountedTotalConcentrateMassFunc(coverMap)
      centsPerTonnePerKm =  0.01*theUnitManager.ConvertToBaseUnits("AUD/tonne")  # actually in AUD/tonne/Km but regional lengths given in km
      
      
      roadOpexPerTKm = ScaledBetaDistribution.FromMeanAndStd(8.5,min=5.,max = 12.)  # 5-12 cents/t/km
      railOpexPerTKm = ScaledBetaDistribution.FromMeanAndStd(1.75,min=1.,max = 2.5) # 1-2.5 cents/t/km
      
      roadNPVMap += (roadTransportationDistance* discountedTotalConcentrateMassMap) * \
                          (roadOpexPerTKm*AUD2010* centsPerTonnePerKm)# 5-12 cents/t/km
      
      railNPVMap += (railTransportationDistance*  discountedTotalConcentrateMassMap) * \
                          (railOpexPerTKm*AUD2010*centsPerTonnePerKm) # 1-2.5 cents/t/km
      
      
      transportationCostNPV = np.minimum(roadNPVMap,railNPVMap)
      
      
      return transportationCostNPV
      

    def CalculateRegionalCoalTransportationExpenses(self,problemManager,coalTransportationDistance,massCoalPerMassHydrogen, \
                                                     capacityMap,hydrogenPlantCapacity,discountedTotalConcentrateMassFunc):
      """Calculate regional transporation costs for coal transportation."""
        
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
      transportationCostNPV = 0.0  # assumed already paid for (for hydrogen export)
      
      # transportation costs 
      # fleet capex
      transportationCostNPV += massCoalPerMassHydrogen*hydrogenPlantCapacity * coalTransportationDistance * 40. *AUD2010*  thousandAUDPerMtkmPerYear
    
      # discounted continuing costs
      discountedTotalCoalMassMap =  massCoalPerMassHydrogen*discountedTotalConcentrateMassFunc(capacityMap)
      centsPerTonnePerKm =  0.01*theUnitManager.ConvertToBaseUnits("AUD/tonne")  # actually in AUD/tonne/Km but regional lengths given in km
      transportationCostNPV += coalTransportationDistance* discountedTotalCoalMassMap * \
                                  8.5*AUD2010* centsPerTonnePerKm# 5-12 cents/t/km
  
      
      return transportationCostNPV
                
       
    def CalculateWaterExpenses(self, problemManager, hydrogenDataManager):
      """ 
      Calculate water startup infrastructure costs - ongoing costs are assumed included in the mine and processing cost models.
      """
      numYears = hydrogenDataManager.GetProjectDuration()
      
      theUnitManager = UnitManager()
      
      self.totalWaterUse = hydrogenDataManager.theHydrogenPlant.waterUse  +  hydrogenDataManager.thePowerPlant.waterUse
      
      maxWaterRequirement = np.max(self.totalWaterUse)
      
      AUD2018 = MiningEquipmentPriceIndex.IndexedPrice(1.0,2018,2018)
      
      v = 2.0   # 2 m/s flow rate 
      q = maxWaterRequirement  # in kL in kL/year (i.e. m^3/year)
      sPerYear = 365*24*3600
      diam = np.sqrt( 4* q/ ( np.pi * sPerYear * v) )
      
      if diam < 150e-3:  # capped to prevent -ve pipeline costs (approximate diameter of "small" pipes in AUSIMM cost est)
        diam = 150e-3
        
      distanceScale =  (2.6229e6  * diam - 125528) * AUD2018   # empirical fit to cost of pipeline data per km in 2018 AUD
      
      oneKm = theUnitManager.ConvertToBaseUnits("km")
      
      self.waterCapex = np.zeros(numYears)  

      
      self.waterCapex[0] =  distanceScale*self.distanceToWater/oneKm
      
      waterCost = hydrogenDataManager.theEconomicDataManager.commodityPrices["H2O"]  # cost per L/Kg
      waterCost *= 1.0/theUnitManager.ConvertToBaseUnits("L")
      
      #self.waterOpex = np.zeros(numYears)
      self.waterOpex = waterCost*self.totalWaterUse
    
    def CalculateWaterUsePerYear(self, problemManager):
      """Calculate average water use in kL/year."""
      q = np.mean(self.totalWaterUse)
      return q 
      
    def CalculateRegionalWaterExpenses(self, problemManager,distanceToWater):
      """Calculate regional water costs."""
          
      theUnitManager = UnitManager()
      #processingCapacity = problemManager.theMineDataManager.theProcessingSystem.processingCapacity*  theUnitManager.ConvertTo("tonne")  # tonnes per year
      
      
      q = np.max(self.totalWaterUse) # max water use in kL/year
      
      sPerYear = 365*24*3600.
      v = 2.0   # 2 m/s flow rate 
      diam = np.sqrt( 4* q/ ( np.pi * sPerYear * v) )
      
      if diam < 150e-3:  # capped to prevent -ve pipeline costs (approximate diameter of "small" pipes in AUSIMM cost est)
        diam = 150e-3
      
      distanceScale=  2.6229e6*diam - 125528   # empirical fit to pipeline data per km - 2018 dollars
     
      #print "water diam", diam*1e3  # diam in mm
      #print "pipe cost ($/km)", distanceScale  # cost in 2018 dollars
      
      rv =  distanceScale*distanceToWater
      return rv
  
      
    def CalculatePowerExpenses(self, problemManager, hydrogenManager):
      """Calculate power infrastructure expenses (local calculation)."""
      numYears = hydrogenManager.GetProjectDuration()
      self.powerCapex = np.zeros(numYears)  
      
      if(self.includePowerInfrastructureCosts):   
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
      
      #Power opex is ignored as included in plant opex

 
    def CalculateRegionalPowerExpenses(self, problemManager, distanceToPower):
      """ Calculate power transmission infrastructure costs."""
      
      # distanceToPower is assumed given in km
            
      theUnitManager = UnitManager()
      powerCosts = np.zeros(distanceToPower.shape)
      
      if(self.includePowerInfrastructureCosts):  
          AUD2010 = MiningEquipmentPriceIndex.IndexedPrice(1.0,2010,2018) * theUnitManager.ConvertToBaseUnits("AUD")  #in today's dollars
          AUD2014 = MiningEquipmentPriceIndex.IndexedPrice(1.0,2014,2018) * theUnitManager.ConvertToBaseUnits("AUD") 
          AUD2016 = MiningEquipmentPriceIndex.IndexedPrice(1.0,2016,2018) * theUnitManager.ConvertToBaseUnits("AUD") 
      
          powerCosts[distanceToPower > 190] =  925000*AUD2014
      
          powerCosts[ np.logical_and(distanceToPower <= 190,distanceToPower > 150)  ] = 650000*AUD2014
      
          powerCosts[ np.logical_and(distanceToPower <= 150,distanceToPower > 1)  ] = 300000*AUD2010
      
          powerCosts[distanceToPower <= 1] =  75000*AUD2010
      
          powerCosts *= distanceToPower
      
      
      return powerCosts
      
      
      #Power opex is ignored as included in mine and processing opex
      
      
       
    def CalculateRegionalGasExpenses(self,problemManager, distanceToGas):
      """ Calculate gas pipeline transmission costs to hydrogen plant."""
      
      # distanceToGas is assumed given in km
            
      theUnitManager = UnitManager()
      gasCosts = np.zeros(distanceToGas.shape)
      
      #AUD2010 = MiningEquipmentPriceIndex.IndexedPrice(1.0,2010,2018) * theUnitManager.ConvertToBaseUnits("AUD")  #in 2018 dollars
      AUD2016 = MiningEquipmentPriceIndex.IndexedPrice(1.0,2016,2018) * theUnitManager.ConvertToBaseUnits("AUD") 
      AUD2014 = MiningEquipmentPriceIndex.IndexedPrice(1.0,2014,2018) * theUnitManager.ConvertToBaseUnits("AUD") 
      
      
      piplelineCostPerKM = 0.315e6*AUD2014  # from core gas production and transmission costs 8 inch class 600 5.6mm wall.
      
      gasCosts = piplelineCostPerKM*distanceToGas
            
      # The following are excluded as we are using gas lines for transmission - not power (gas costs are included in plant opex). 
      
      #AUD2018 = MiningEquipmentPriceIndex.IndexedPrice(1.0,2018,2018) * theUnitManager.ConvertToBaseUnits("AUD") 
      #gasLCOE = 120*AUD2016   # LCOE per Mwh
      #retailCOE = 100*AUD2018 # LCOE per Mwh

      #gasPowerCosts = problemManager.theMineDataManager.theProcessingSystem.processingPower*( gasLCOE - retailCOE)
      #gasNPC = problemManager.theMineDataManager.theEconomicDataManager.CalculateNPV(gasPowerCosts)
      #gasCosts+= gasNPC
        
      
      return gasCosts
     
     
    def CalculateRegionalCO2Expenses(self,problemManager, distanceToCO2,flowRate):
      """ Calculate CO2 pipeline transmission costs to hydrogen plant."""
      
      # distanceToCO2 is assumed given in km
      
      theUnitManager = UnitManager()
      co2Costs = np.zeros(distanceToCO2.shape)
      
      flowRateInMtYr = flowRate * theUnitManager.ConvertTo("Mt")   # covert to Mt/a
      
      #AUD2010 = MiningEquipmentPriceIndex.IndexedPrice(1.0,2010,2018) * theUnitManager.ConvertToBaseUnits("AUD")  #in 2018 dollars
      #AUD2016 = MiningEquipmentPriceIndex.IndexedPrice(1.0,2016,2018) * theUnitManager.ConvertToBaseUnits("AUD") 
      AUD2015 = MiningEquipmentPriceIndex.IndexedPrice(1.0,2015,2018) * theUnitManager.ConvertToBaseUnits("AUD") 
      
      #powerCosts[distanceToPower > 190] =  925000*AUD2014
      
      # piplelineCostPerKM = 0.315e6*AUD2014  # from core gas production and transmission costs 8 inch class 600 5.6mm wall.
      
      # CCS pipeline cost model from 
      # AUSTRALIAN POWER GENERATION TECHNOLOGY REPORT 
      # Predicts pipeline costs  in the form of a powerlaw a(l)^b where l is in km
      flowRates = np.array([1,3,5,10,15,20,25,30,35,40])   # flow rate Mt/a
      aValues = np.array([0.144,0.243,0.333,0.417,0.466,0.568,0.617,0.693,0.799,0.875] )  # a
      bValues = np.array([1.25,1.25,1.24,1.27,1.29,1.28,1.29,1.29,1.28,1.27]) # b
      
      
      if(flowRateInMtYr < 1.0 ):
        flowRateInMtYr = 1.0

      a_interp = interpolate.interp1d(flowRates, aValues)
      b_interp = interpolate.interp1d(flowRates, bValues)

      aActual = AUD2015*a_interp(flowRateInMtYr)*1e6
      bActual = b_interp(flowRateInMtYr)
      
      co2Costs = aActual*distanceToCO2**bActual
        
      
      return co2Costs
     
                