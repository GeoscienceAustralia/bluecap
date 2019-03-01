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

# external libraries
import numpy as np
import pylab as pl

# Managers
from OreBodyDataManager import OreBodyDataManager
from MiningSystemDataManager import MiningSystemDataManager
from ProcessingSystemDataManager import ProcessingSystemDataManager
from EconomicDataManager import EconomicDataManager
from InfrastructureDataManager import InfrastructureDataManager


# Functions
from Functions.FunctionManager import FunctionManager

# IO
from IO.XML import HasChild,GetChild,AddChild
from IO.XML import GetAttributeValue, SetAttributeString

class MineDataManager():
    def __init__(self):
      """
      Create an empty mine data manager and default variables. 
      """
      
      self.mineLatLong = np.array([0.0,0.0])  
      self.theOreBody = OreBodyDataManager()
      self.theMiningSystem = MiningSystemDataManager()
      self.theProcessingSystem = ProcessingSystemDataManager()
      self.theEconomicDataManager = EconomicDataManager()
      self.theInfrastructureManager = InfrastructureDataManager()
      


    def ParseXMLNode(self, mineDataNode):
      """
      Generate Mine Data Manager data from xml tree node. 
      """
      
      # Location
      if(HasChild(mineDataNode,"Location")):
        locNode = GetChild(mineDataNode,"Location")
        self.mineLatLong[0] = GetAttributeValue(locNode,"lat")
        self.mineLatLong[1] = GetAttributeValue(locNode,"long")
      
      # Orebody
      if(HasChild(mineDataNode,"Orebody")):
        orebodyNode = GetChild(mineDataNode,"Orebody")
        self.theOreBody.ParseXMLNode(orebodyNode)
       
      
      theFunctionManager = FunctionManager()
      if( (self.theOreBody.cover < 0.0 ) and (theFunctionManager.HasFunction("DepthOfCover") ) ):
        self.theOreBody.cover = theFunctionManager.GetFunction("DepthOfCover").f( self.mineLatLong  )
        print "Cover set to: ", self.theOreBody.cover 
        
        
        
      # Economics
      if(HasChild(mineDataNode,"Economics")):
        economicsNode = GetChild(mineDataNode,"Economics")
        self.theEconomicDataManager.ParseXMLNode(economicsNode)

    def WriteXMLNode(self, node):
      """
      Write problem to xml node
      """
      
      # Location
      locNode = AddChild(node,"Location")
      SetAttributeString(locNode,"lat",self.mineLatLong[0])
      SetAttributeString(locNode,"long",self.mineLatLong[1])
      
      # Orebody
      orebodyNode = AddChild(node,"Orebody")
      self.theOreBody.WriteXMLNode(orebodyNode)
        
      # Economics Node
      economicsNode = AddChild(node,"Economics")
      self.theEconomicDataManager.WriteXMLNode(economicsNode)
      
      
      # Mining System Node
      miningNode = AddChild(node,"Mining")
      self.theMiningSystem.WriteXMLNode(miningNode)
      
      # Processing System Node
      processingNode = AddChild(node,"Processing")
      self.theProcessingSystem.WriteXMLNode(processingNode)
      
      # Infrastructure Node
      infrastructureNode = AddChild(node,"Infrastructure")
      self.theInfrastructureManager.WriteXMLNode(infrastructureNode)
      
      return node
 
 
 
    
    def SetMineType(self,mineType):
      self.theMiningSystem.mineType = mineType
 
 
 
    def CaculateMineProductionAndValue(self,problemManager):
      """
      Determine after tax NPV for the mine
      """
      
      # Mining Model
      self.DetermineMiningSystem(problemManager)
      
      # Processing Model
      self.DetermineProcessingSystem(problemManager)
      
      # G&A Model
      self.CalculateGandAExpenses(problemManager)
      
      # Infrastructure Model
      self.CalculateInfrastructureCosts(problemManager)
      
      # Cash flow
      self.CalculateBeforeTaxCashFlow(problemManager)
      self.CalculateTaxes(problemManager)
      self.CalculateAfterTaxCashFlow(problemManager)
      
      # EconomicIndicators
      self.CalculateEconomicIndicators(problemManager)
      
      value = self.theEconomicDataManager.atNPV
    
      return value
      
           
      
    def SetMiningMethod(self,miningMethod):
      """
      Set the mining method (mine production and value need to be determined separately)
      """
      self.theMiningSystem.miningMethod = miningMethod
   
   
    def SetCoverDepth(self,cover):
      
      self.theOreBody.cover = cover
      
      return None
   
    def DetermineMiningSystem(self,problemManager):
      self.theMiningSystem.DetermineMiningSystem(problemManager,self)
      
      return self.theMiningSystem   
      
    def DetermineProcessingSystem(self,problemManager):
      self.theProcessingSystem.DetermineProcessingSystem(problemManager,self)
      
      return self.theProcessingSystem
      
    
      
    # G&A Model
    def CalculateGandAExpenses(self,problemManager):
      """
      General and administrative costs are estimated based on a fixed percentage of the overall mining and processing costs.
      """
      
      self.theEconomicDataManager.CalculateGandAExpenses(problemManager,self)
      
      return self.theEconomicDataManager.GandAOpex
    
      
      
    # Infrastructure Model
    def CalculateInfrastructureCosts(self,problemManager):
      self.theInfrastructureManager.DetermineDistanceToInfrastructure(problemManager,self)
      self.theInfrastructureManager.CalculateInfrastructureExpenses(problemManager,self)
      return None
      
    # Cash flow
    def CalculateBeforeTaxCashFlow(self,problemManager):
      self.theEconomicDataManager.CalculateBeforeTaxCashFlow(problemManager,self)
      return self.theEconomicDataManager.btNCF

    def CalculateTaxes(self,problemManager):
      self.theEconomicDataManager.CalculateTaxes(problemManager,self)
      return self.theEconomicDataManager.taxes

    def CalculateAfterTaxCashFlow(self,problemManager):
      self.theEconomicDataManager.CalculateAfterTaxCashFlow(problemManager,self)
      return self.theEconomicDataManager.atNCF
      
    # EconomicIndicators
    def CalculateEconomicIndicators(self,problemManager):
      
      self.theEconomicDataManager.CalculateBeforeTaxNPV(problemManager,self)
      self.theEconomicDataManager.CalculateAfterTaxNPV(problemManager,self)
      return self.theEconomicDataManager.atNPV
      
      
    def PlotResults(self):
      pl.plot(self.theEconomicDataManager.btNCF/1e6,"b-",label="btNCF")
      
      pl.plot(self.theEconomicDataManager.royalties/1e6,"r--",label="Royalties")
      pl.plot(self.theEconomicDataManager.taxes/1e6,"r-",label="Income tax")
      pl.plot(self.theEconomicDataManager.atNCF/1e6,'g-',label="atNCF")
      pl.legend()
      pl.figure()
      
      pl.plot(self.theMiningSystem.materialMined/1e6,"k-",label="total mined")
      pl.plot(self.theMiningSystem.oreMined/1e6,"b-",label="ore")
      pl.plot(self.theMiningSystem.wasteMined/1e6,"r-",label="waste")
      pl.plot(self.theProcessingSystem.oreProcessed/1e6,"g-",label="processed")
      
      
      pl.plot(self.theProcessingSystem.concentrateProduced/1e6,"g--",label="concentrate")
      pl.legend(loc="lower right")
      
      pl.figure()
      
      
      #print self.theMiningSystem.miningCapex
      
      inflation =  (1.0+self.theEconomicDataManager.inflation)**np.array( range( self.theMiningSystem.mineLife ) )
     
      
      pl.plot(self.theEconomicDataManager.btNCF/1e6,"b-",label="btNCF(MMAUD)")
      pl.plot(self.theEconomicDataManager.revenue/1e6,"g-",label="Revenue")
      pl.plot(inflation*self.theMiningSystem.miningCapex/1e6,"r--",label="Mining Startup")
      pl.plot(inflation*self.theMiningSystem.miningOpex/1e6,"r-",label="Mining Sustaining")
      pl.plot(inflation*self.theProcessingSystem.processingCapex/1e6,"m--",label="Processing Startup")
      pl.plot(inflation*self.theProcessingSystem.processingOpex/1e6,"m-",label="Processing Sustaining")
      pl.plot(inflation*self.theEconomicDataManager.GandAOpex/1e6,"r:",label="G&A")
      pl.plot(inflation*self.theInfrastructureManager.infrastructureCapex/1e6,"k--",label="Infrastructure Startup")
      pl.plot(inflation*self.theInfrastructureManager.infrastructureOpex/1e6,"k-",label="Infrastructure Sustaining")
      pl.legend(loc="lower right")
      
      
      pl.show()
      
    def RecordResults(self,problemManager):
      fid = open(problemManager.outputPrefix + ".txt","wt")
      fid.write("#BTNCF:\n")
      np.savetxt(fid,self.theEconomicDataManager.btNCF,newline=" ")
      fid.write("\n")
      fid.write("#Royalties:\n")
      np.savetxt(fid,self.theEconomicDataManager.royalties,newline=" ")
      fid.write("\n")
      fid.write("#taxes:\n")
      np.savetxt(fid,self.theEconomicDataManager.taxes,newline=" ")
      fid.write("\n")
      fid.write("#ATNCF:\n")
      np.savetxt(fid,self.theEconomicDataManager.atNCF,newline=" ")
      fid.write("\n")
      
      
      fid.write("#Total Mined:\n")
      np.savetxt(fid,self.theMiningSystem.materialMined,newline=" ")
      fid.write("\n")
      fid.write("#Ore Mined:\n")
      np.savetxt(fid,self.theMiningSystem.oreMined,newline=" ")
      fid.write("\n")
      fid.write("#Waste Mined:\n")
      np.savetxt(fid,self.theMiningSystem.wasteMined,newline=" ")
      fid.write("\n")
      fid.write("#Ore processed:\n")
      np.savetxt(fid,self.theProcessingSystem.oreProcessed,newline=" ")
      fid.write("\n")
      fid.write("#Concentrate:\n")
      np.savetxt(fid,self.theProcessingSystem.concentrateProduced,newline=" ")
      fid.write("\n")
      
      
      fid.write("#Revenue:\n")
      np.savetxt(fid,self.theEconomicDataManager.revenue,newline=" ")
      fid.write("\n")
      
      fid.write("#Mining Startup:\n")
      np.savetxt(fid,self.theMiningSystem.miningCapex,newline=" ")
      fid.write("\n")
      fid.write("#Mining Sustaining:\n")
      np.savetxt(fid,self.theMiningSystem.miningOpex,newline=" ")
      fid.write("\n")
      
      fid.write("#Processing Startup:\n")
      np.savetxt(fid,self.theProcessingSystem.processingCapex,newline=" ")
      fid.write("\n")
      fid.write("#Processing Sustaining:\n")
      np.savetxt(fid,self.theProcessingSystem.processingOpex,newline=" ")
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
      