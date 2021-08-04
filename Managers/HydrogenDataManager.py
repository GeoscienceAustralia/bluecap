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

# external libraries
import numpy as np
import pylab as pl

# Managers
from .HydrogenPlantManager import HydrogenPlantManager
from .PowerPlantManager import PowerPlantManager
from .HydrogenEnergyStorageManager import HydrogenEnergyStorageManager

from .HydrogenEconomicDataManager import HydrogenEconomicDataManager
from .HydrogenInfrastructureDataManager import HydrogenInfrastructureDataManager




# Functions
from Functions.FunctionManager import FunctionManager

# IO
from IO.XML import HasChild,GetChild,GetChildren,AddChild
from IO.XML import GetAttributeValue, GetAttributeString, SetAttributeString
from IO.XML import GetAttributeStringOrDefault

class HydrogenDataManager():
    def __init__(self):
      """
      Create an empty hydrogen data manager and default variables. 
      """
      
      self.theHydrogenPlant = HydrogenPlantManager()
      self.thePowerPlant = PowerPlantManager()
      
      self.theEnergyStorage = None 
      
      self.theEconomicDataManager = HydrogenEconomicDataManager()
      self.theInfrastructureManager = HydrogenInfrastructureDataManager()
      
      
      self.theRehabilitationManager = None  # not yet implemented
      
      
      # location
      self.plantLatLong = np.array([0.0,0.0])
      
      


    def ParseXMLNode(self, hydrogenDataNode):
      """
      Generate Hydrogen Data Manager data from xml tree node. 
      """
      
      # Location
      if(HasChild(hydrogenDataNode,"Location")):
        locNode = GetChild(hydrogenDataNode,"Location")
        self.plantLatLong[0] = GetAttributeValue(locNode,"lat")
        self.plantLatLong[1] = GetAttributeValue(locNode,"long")
      
      
      # Hydrogen plant
      if(HasChild(hydrogenDataNode,"HydrogenPlant")):
        hydrogenPlantNode = GetChild(hydrogenDataNode,"HydrogenPlant")
        self.theHydrogenPlant.ParseXMLNode(hydrogenPlantNode)  
        
      
      # Power plant
      if(HasChild(hydrogenDataNode,"PowerPlant")):
        powerPlantNode = GetChild(hydrogenDataNode,"PowerPlant")
        self.thePowerPlant.ParseXMLNode(powerPlantNode)  
 
       # Energy storage
      if(HasChild(hydrogenDataNode,"EnergyStorage")):
        energyStorageNode = GetChild(hydrogenDataNode,"EnergyStorage")
        self.theEnergyStorage = HydrogenEnergyStorageManager()
        self.theEnergyStorage.ParseXMLNode(energyStorageNode)
        
      # Infrastructure
      if(HasChild(hydrogenDataNode,"Infrastructure")):
        infrastructureNode = GetChild(hydrogenDataNode,"Infrastructure")
        self.theInfrastructureManager.ParseXMLNode(infrastructureNode)      
      
        
      # Economics
      if(HasChild(hydrogenDataNode,"Economics")):
        economicsNode = GetChild(hydrogenDataNode,"Economics")
        self.theEconomicDataManager.ParseXMLNode(economicsNode)

    def WriteXMLNode(self, node):
      """
      Write hydrogen data manager to xml node.
      """
      
      # Location
      locNode = AddChild(node,"Location")
      SetAttributeString(locNode,"lat",self.latLong[0])
      SetAttributeString(locNode,"long",self.latLong[1])
        
      # Economics Node
      economicsNode = AddChild(node,"Economics")
      self.theEconomicDataManager.WriteXMLNode(economicsNode)
      
      # Hydrogen Plant Node
      hydrogenPlantNode = AddChild(node,"HydrogenPlant")
      self.theHydrogenPlant.WriteXMLNode(hydrogenPlantNode)
      
      # Power Plant Node
      powerPlantNode = AddChild(node,"PowerPlant")
      self.thePowerPlant.WriteXMLNode(powerPlantNode)
      
      if(self.theEnergyStorage):
        energyStorageNode = AddChild(node,"EnergyStorage")
        self.theEnergyStorage.WriteXMLNode(energyStorageNode)
      
      # Infrastructure Node
      infrastructureNode = AddChild(node,"Infrastructure")
      self.theInfrastructureManager.WriteXMLNode(infrastructureNode)
      
      #if(self.theRehabilitationManager):
      #  rehabilitationNode = AddChild(node,"Rehabilitation")
      #  self.theRehabilitationManager.WriteXMLNode(rehabilitationNode)
      
      return node
 
 
    def SetHydrogenPlantType(self,plantType):
      self.theHydrogenPlant.plantType = plantType
 
    def SetPowerPlantType(self,plantType):
      self.thePowerPlant.plantType = plantType
    
 
    def CaculateHydrogenProductionAndValue(self,problemManager):
      """
      Determine after tax NPV for the project.
      """
      
      # Hydrogen Model 
      # - nb plant life is determined from the power system, but plant capacity determined based on hydrogen output. 
      # - so power system needs to fix duration 
      # - If so LCOH rather than NPV should be used to compare projects of varying duration. 
      self.DetermineHydrogenSystem(problemManager)
      
      # Power Model
      self.DeterminePowerSystem(problemManager)
      
      # G&A Model
      self.CalculateGandAExpenses(problemManager)
      
      # Infrastructure Model
      self.CalculateInfrastructureCosts(problemManager)
      
      # Rehabilitation Model (if present)
      #self.CalculateRehabilitationCosts(problemManager)
      
      # Cash flow
      self.CalculateBeforeTaxCashFlow(problemManager)
      self.CalculateTaxes(problemManager)
      self.CalculateAfterTaxCashFlow(problemManager)
      
      # EconomicIndicators
      self.CalculateEconomicIndicators(problemManager)
      
      value = self.theEconomicDataManager.atNPV
    
      return value
      
    def GetCapacityFactor(self):
      """ Get the capacity factor for the energy supplied by the power plant(s) to the hydrogen plant (+ storage). """
      
      capacityFactor = self.thePowerPlant.capacityFactor
      
      return capacityFactor      
      
    def SetCapacityFactor(self,capacityFactor):
      """ Set the capacity factor for the energy supplied by the power plant(s) to the hydrogen plant (+ storage). """
      
      self.thePowerPlant.SetCapacityFactor( capacityFactor)
      
      if(self.theEnergyStorage):
        plantCF = self.theEnergyStorage.netOutputCapacity
        self.theHydrogenPlant.SetCapacityFactor( plantCF )
      else:
        self.theHydrogenPlant.SetCapacityFactor( capacityFactor )
      
      return None
   
    def DetermineHydrogenSystem(self,problemManager):
      """Initialize the hydrogen plant variables and production parameters."""
      self.theHydrogenPlant.DetermineHydrogenSystem(problemManager,self)
      
      return self.theHydrogenPlant   
      
    
    def DeterminePowerSystem(self,problemManager):
      """Initialize the power plant variables and production parameters."""
      self.thePowerPlant.DeterminePowerSystem(problemManager,self)
      
      return self.thePowerPlant
      
    def GetProjectDuration(self):
      """Return the total project duration."""
      numYears = self.theHydrogenPlant.startupTime + self.thePowerPlant.operatingLife
      return numYears
      
    # G&A Model
    def CalculateGandAExpenses(self,problemManager):
      """
      General and administrative costs are estimated based on a fixed percentage of the overall project costs.
      """
      
      self.theEconomicDataManager.CalculateGandAExpenses(problemManager,self)
      
      return self.theEconomicDataManager.GandAOpex
    
      
      
    # Infrastructure Model
    def CalculateInfrastructureCosts(self,problemManager):
      """Calculate infrastructure costs associated with a single project."""
      self.theInfrastructureManager.DetermineDistanceToInfrastructure(problemManager,self)
      self.theInfrastructureManager.CalculateInfrastructureExpenses(problemManager,self)
      return None
 
    # Rehabilitation Model
    #def CalculateRehabilitationCosts(self,problemManager):
    #  if(self.theRehabilitationManager):
    #    self.theRehabilitationManager.CalculateRehabilitationExpenses(problemManager,self)
    #  return None
 
 
    # Zero Infrastructure Costs
    def ZeroInfrastructureCosts(self,problemManager):
      """Set all infrastructure costs to zero (used in regional calculation)."""
      self.theInfrastructureManager.ZeroInfrastructureExpenses(problemManager,self)
      return None


    # Zero Processing Capex
    #def ZeroProcessingCapex(self,problemManager):
    #  self.theProcessingSystem.ZeroProcessingCapex(problemManager,self)
    #  return None
      
    # Zero Rehabilitation Costs
    #def ZeroRehabilitationCosts(self,problemManager):
    #  self.theRehabilitationManager.ZeroRehabilitationExpenses(problemManager,self)
    #  return None
      
    # Startup costs
    def CalculateStartupCosts(self):
      """Calculate net startup costs (will be capitalized)."""
      startupCosts = self.theHydrogenPlant.capex \
                      + self.thePowerPlant.capex  \
                      + self.theInfrastructureManager.infrastructureCapex
      
      if(self.theEnergyStorage):
        startupCosts += self.theEnergyStorage.capex
      
      return startupCosts
    
    # Sustaining Costs  
    def CalculateSustainingCosts(self):
      """ Calculate sustaining costs ( nb a proportion will be capitalized - assumed to be sustaining capex)."""
      sustainingCosts = self.theHydrogenPlant.opex \
                      + self.thePowerPlant.opex  \
                      + self.theInfrastructureManager.infrastructureOpex \
                      + self.theEconomicDataManager.GandAOpex
      
      if(self.theEnergyStorage):
        #print("sustainingCosts", sustainingCosts)
        #print("self.theEnergyStorage.opex", self.theEnergyStorage.opex)
        sustainingCosts += self.theEnergyStorage.opex
        
      if(self.theRehabilitationManager):
        sustainingCosts +=  self.theRehabilitationManager.rehabilitationCosts
      
      return sustainingCosts
          
    # Cash flow
    def CalculateBeforeTaxCashFlow(self,problemManager):
      """Calculate the cash flow before state or federal taxes (and any rebates or royalties) are applied."""
      self.theEconomicDataManager.CalculateBeforeTaxCashFlow(problemManager,self)
      return self.theEconomicDataManager.btNCF

    def CalculateTaxes(self,problemManager):
      """Calculate state and federal taxes (and any rebates or royalties) for the project."""
      self.theEconomicDataManager.CalculateTaxes(problemManager,self)
      return self.theEconomicDataManager.taxes

    def CalculateAfterTaxCashFlow(self,problemManager):
      """Calculate cash flow after state and federal taxes (and any rebates or royalties) are accounted for."""
      self.theEconomicDataManager.CalculateAfterTaxCashFlow(problemManager,self)
      return self.theEconomicDataManager.atNCF
      
    # EconomicIndicators
    def CalculateEconomicIndicators(self,problemManager):
      """Calculate economic indicators associated with the project. Returns after tax net present value."""
      self.theEconomicDataManager.CalculateBeforeTaxNPV(problemManager,self)
      self.theEconomicDataManager.CalculateAfterTaxNPV(problemManager,self)
      return self.theEconomicDataManager.atNPV

    def CalculateEconomicIndicator(self,problemManager,type):
      """Record before and after tax NPV and return after tax NPV for the project"""
      rv = 0.0
      if(type == "atNPV"):
        self.theEconomicDataManager.CalculateBeforeTaxNPV(problemManager,self)
        self.theEconomicDataManager.CalculateAfterTaxNPV(problemManager,self)
        rv = self.theEconomicDataManager.atNPV 
      elif(type == "btNPV"):
        self.theEconomicDataManager.CalculateBeforeTaxNPV(problemManager,self)
        rv = self.theEconomicDataManager.btNPV 
      elif(type == "NetRealCost"):
        rv = np.sum(self.theEconomicDataManager.btNCF - self.theEconomicDataManager.revenue)
      return rv
      
      
    def PlotResults(self):
      """Plot the results of the single-site calculation and display on screen."""
      pl.plot(self.theEconomicDataManager.btNCF/1e6,"b-",label="btNCF")
      
      #pl.plot(self.theEconomicDataManager.royalties/1e6,"r--",label="Royalties")
      pl.plot(self.theEconomicDataManager.taxes/1e6,"r-",label="Income tax")
      pl.plot(self.theEconomicDataManager.atNCF/1e6,'g-',label="atNCF")
      pl.legend()
      pl.figure()
      
      inflation =  (1.0+self.theEconomicDataManager.inflation)**np.array( range( self.theHydrogenPlant.projectLife ) )
     
      pl.plot(self.theEconomicDataManager.btNCF/1e6,"b-",label="btNCF(MMAUD)")
      pl.plot(self.theEconomicDataManager.revenue/1e6,"g-",label="Revenue")
      pl.plot(inflation*self.theHydrogenPlant.capex/1e6,"r--",label="Hydrogen Startup")
      pl.plot(inflation*self.theHydrogenPlant.opex/1e6,"r-",label="Hydrogen Sustaining")
      pl.plot(inflation*self.thePowerPlant.capex/1e6,"m--",label="Power Startup")
      pl.plot(inflation*self.thePowerPlant.opex/1e6,"m-",label="Power Sustaining")
      pl.plot(inflation*self.theEconomicDataManager.GandAOpex/1e6,"r:",label="G&A")
      pl.plot(inflation*self.theInfrastructureManager.infrastructureCapex/1e6,"k--",label="Infrastructure Startup")
      pl.plot(inflation*self.theInfrastructureManager.infrastructureOpex/1e6,"k-",label="Infrastructure Sustaining")
      pl.legend(loc="lower right")
      
      pl.show()
      
    def RecordResults(self,problemManager,fmt='txt'):
      """Record the results of the single-site calculation to file."""
      if fmt == 'txt':
        fid = open(problemManager.outputPrefix + ".txt","wt")
        fid.write("#BTNCF:\n")
        np.savetxt(fid,self.theEconomicDataManager.btNCF,newline=" ")
        fid.write("\n")
        
        fid.write("#taxes:\n")
        np.savetxt(fid,self.theEconomicDataManager.taxes,newline=" ")
        fid.write("\n")
        fid.write("#ATNCF:\n")
        np.savetxt(fid,self.theEconomicDataManager.atNCF,newline=" ")
        fid.write("\n")
        
        fid.write("#Revenue:\n")
        np.savetxt(fid,self.theEconomicDataManager.revenue,newline=" ")
        fid.write("\n")
        
        fid.write("#Hydrogen Startup:\n")
        np.savetxt(fid,self.theHydrogenPlant.capex,newline=" ")
        fid.write("\n")
        fid.write("#Hydrogen Sustaining:\n")
        np.savetxt(fid,self.theHydrogenPlant.opex,newline=" ")
        fid.write("\n")
        
        fid.write("#Power Startup:\n")
        np.savetxt(fid,self.thePowerPlant.capex,newline=" ")
        fid.write("\n")
        fid.write("#Power Sustaining:\n")
        np.savetxt(fid,self.thePowerPlant.opex,newline=" ")
        fid.write("\n")
        fid.write("#G&A:\n")
        np.savetxt(fid,self.theEconomicDataManager.GandAOpex,newline=" ")
        fid.write("\n")
        fid.write("#Infrastructure Startup:\n")
        np.savetxt(fid,self.theInfrastructureManager.infrastructureCapex,newline=" ")
        fid.write("\n")
        fid.write("#Infrastructure Sustaining:\n")
        np.savetxt(fid,self.theInfrastructureManager.infrastructureOpex,newline=" ")
        
        fid.close()
      
      elif fmt == 'csv':
        years = np.arange(self.theHydrogenPlant.projectLife, dtype=int)
        with open(problemManager.outputPrefix + ".csv","wt") as f:
          f.write('YEAR,' + ','.join((years + 1).astype(str)) + '\n')
          f.write('BTNCF,' + ','.join(self.theEconomicDataManager.btNCF.astype(str)) + '\n')
          f.write('Royalties,' + ','.join(self.theEconomicDataManager.royalties.astype(str)) + '\n')
          f.write('taxes,' + ','.join(self.theEconomicDataManager.taxes.astype(str)) + '\n')
          f.write('ATNCF,' + ','.join(self.theEconomicDataManager.atNCF.astype(str)) + '\n')
          f.write('Total Mined,' + ','.join(self.theMiningSystem.materialMined.astype(str)) + '\n')
          f.write('Ore Mined,' + ','.join(self.theMiningSystem.oreMined.astype(str)) + '\n')
          f.write('Waste Mined,' + ','.join(self.theMiningSystem.wasteMined.astype(str)) + '\n')
          f.write('Ore processed,' + ','.join(self.theProcessingSystem.oreProcessed.astype(str)) + '\n')
          f.write('Concentrate,' + ','.join(self.theProcessingSystem.concentrateProduced.astype(str)) + '\n')
          f.write('Revenue,' + ','.join(self.theEconomicDataManager.revenue.astype(str)) + '\n')
          f.write('Mining Startup,' + ','.join(self.theMiningSystem.miningCapex.astype(str)) + '\n')
          f.write('Mining Sustaining,' + ','.join(self.theMiningSystem.miningOpex.astype(str)) + '\n')
          f.write('Processing Startup,' + ','.join(self.theProcessingSystem.processingCapex.astype(str)) + '\n')
          f.write('Processing Sustaining,' + ','.join(self.theProcessingSystem.processingOpex.astype(str)) + '\n')
          f.write('G&A,' + ','.join(self.theEconomicDataManager.GandAOpex.astype(str)) + '\n')
          f.write('Infrastructure Startup,' + ','.join(self.theInfrastructureManager.infrastructureCapex.astype(str)) + '\n')
          f.write('Infrastructure Sustaining,' + ','.join(self.theInfrastructureManager.infrastructureOpex.astype(str)) + '\n')
          f.write('Inflation,' + ','.join(((1.0+self.theEconomicDataManager.inflation)**(years)).astype(str)) + '\n')
          f.write('Discount,' + ','.join((1./((1.0+self.theEconomicDataManager.discountRate)**(years+1))).astype(str)) + '\n')