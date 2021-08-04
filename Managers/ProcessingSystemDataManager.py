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



# IO
from IO.XML import HasChild,GetChild,AddChild
from IO.XML import HasAttribute
from IO.XML import GetAttributeString,GetAttributeStringOrDefault,GetAttributeValue,GetAttributeVectorString,SetAttributeString


#Units
from Units.UnitManager import UnitManager

#Functions
from Functions.FunctionManager import FunctionManager

# Managers

class ProcessingSystemDataManager():
    def __init__(self):
      """
      Create an empty processing system data manager and default variables. 
      """
      self.processingMethod = "Au"  
      self.processingLoss = 0.10
      self.refiningTake = 0.10
      self.userRefiningTake = np.copy(self.refiningTake)
      self.userProcessingLoss = np.copy(self.processingLoss)
      self.processingPower = []
      self.overwrite = False
      
      self.concentrateMetals = []   # produced in this processing system
      #self.externalProducts = []  # produced in other processing systems
      
      self.concentratePrimaryMetal = self.processingMethod  # 
      
      self.secondaryProcessingSystem = None # used to determine if REE extracted
      

    def ParseXMLNode(self, processingSystemNode):
      """
      Generate processing system data from xml tree node. 
      """
      
      # NB processing method, take etc for primary product is set in "Determine Processing system"
      if(HasAttribute(processingSystemNode,"refiningTake") or HasAttribute(processingSystemNode,"processingLoss")):
        self.overwrite = True
        self.userRefiningTake = GetAttributeValue(processingSystemNode,"refiningTake")
        self.userProcessingLoss = GetAttributeValue(processingSystemNode,"processingLoss")
      if(HasChild(processingSystemNode,"SecondarySystem")):
        # Secondary systems are a recursive instance of the ProcessingSystemDataManager
        self.secondaryProcessingSystem = ProcessingSystemDataManager()
        self.secondaryProcessingSystem.ParseSecondarySystemXMLNode(processingSystemNode,self) 
      
      return processingSystemNode
        

    def ParseSecondarySystemXMLNode(self, processingSystemNode,primaryProcessingSystem):
      """
      Generate processing system data from xml tree node. 
      """
      node = GetChild(processingSystemNode,"SecondarySystem")
      self.processingMethod = GetAttributeStringOrDefault(node,"method","")
      
      if(HasAttribute(node,"metals")):
        self.concentrateMetals = GetAttributeVectorString(node,"metals")
      
      if(self.processingMethod):
        self.concentrateMetals.append(self.processingMethod)
        
      self.concentrateMetals = set(self.concentrateMetals)
      
      self.concentratePrimaryMetal = self.processingMethod
      
      self.refiningTake = GetAttributeValue(node,"refiningTake")
      self.processingLoss = GetAttributeValue(node,"processingLoss")
      
      # For tertiary+ system
      if(HasChild(node,"SecondarySystem")):
        # Secondary systems are a recursive instance of the ProcessingSystemDataManager
        self.secondaryProcessingSystem = ProcessingSystemDataManager()
        self.secondaryProcessingSystem.ParseSecondarySystemXMLNode(node,self)
      
      return processingSystemNode

    def WriteXMLNode(self, node):
      """
      Write processing system data to xml node
      """
      # Processing data
      SetAttributeString(node,"method",self.processingMethod)
      SetAttributeString(node,"processingLoss",self.processingLoss)
      SetAttributeString(node,"refiningTake",self.refiningTake)
      
      if(self.secondaryProcessingSystem):
        secondaryNode = AddChild(node,"SecondarySystem")
        self.secondaryProcessingSystem.WriteXMLNode(secondaryNode) 
      
      return node
      
    def DetermineProcessingSystem(self, problemManager, mineDataManager):
      """
      Use available data to determine the most likely processing method/capacity/opex etc.
      """

      referenceMetalStr = mineDataManager.theOreBody.type[:2]  
      # first two letters of orebody type is assumed to be reference metal for determining processing grade
      # eg AuCu -> gold is reference metal,  
      self.processingMethod =  referenceMetalStr 
      
      #processing loss is fixed
      
      if(referenceMetalStr ==  "Au"):
        self.refiningTake = 0.01
      elif(referenceMetalStr ==  "Cu"):
        self.refiningTake = 0.10
      elif(referenceMetalStr ==  "Ni"):
        self.refiningTake = 0.30
      elif(referenceMetalStr ==  "Ag"):
        self.refiningTake = 0.05
      elif(referenceMetalStr ==  "Pb" or referenceMetalStr ==  "Zn" ):
        self.refiningTake = 0.17
      elif(mineDataManager.theOreBody.type ==  "P2O5"):   # phosphates - NB must be after Lead
        self.refiningTake = 0.0  # assumes sold at market price for concentrate 
        self.processingMethod =  "P2O5"
      elif(mineDataManager.theOreBody.type ==  "K2SO4"):   # Potash
        self.refiningTake = 0.0 # assumes sold at market price for concentrate 
        self.processingMethod =  "K2SO4"
      elif(referenceMetalStr ==  "RE"):
        self.refiningTake = 0.0 # assumes sold at market price for concentrate 
        self.processingMethod =  "REO"
      
      if self.overwrite:
        self.refiningTake = self.userRefiningTake
        self.processingLoss = self.userProcessingLoss
      
      allProcessedMetals = self.GetProcessedMetals()
      allContainedMetals =  list( mineDataManager.theOreBody.metalGrades) #python 2.0: list( mineDataManager.theOreBody.metalGrades.keys() )
      self.SetProcessedMetals(allProcessedMetals,allContainedMetals)
        
      self.CalculateProcessingCapacity(problemManager, mineDataManager)
      
      self.CalculateProcessingCapex(problemManager, mineDataManager)
      self.CalculateProcessingOpex(problemManager, mineDataManager)
      
      return self      
    
    def SetProcessedMetals(self,allProcessedMetals,allContainedMetals):
      
      #isFirst = True
      #if (isFirst):
      # everything that not processed in a secondary processing system is assumed processed in primary concentrate
      contained = set(allContainedMetals)
      processedOther = set(allProcessedMetals).difference(self.processingMethod)
        
      self.concentrateMetals = contained.difference(processedOther)
      
    
      """
      fixme - clean up: secondary circuits set metals and processing take explicitly now. 
      else:
        # secondary circuits only produce concentrate with the named commodity
        self.concentrateMetals = set([self.processingMethod])
    
      
      if(self.secondaryProcessingSystem):
        self.secondaryProcessingSystem.SetProcessedMetals(allProcessedMetals,allContainedMetals, False)
      
      """
      return self
    
    def GetProcessedMetals(self,isFirst = True):
      if(isFirst):
        alist = []
      else:  
        alist = list(self.concentrateMetals)
      
      if(self.secondaryProcessingSystem):
        alist.extend(self.secondaryProcessingSystem.GetProcessedMetals(False)) 
      return alist

    
    def CalculateProcessingCapacity(self,  problemManager, mineDataManager):
      """
      Calculate processing capacity, amount of ore processed each year and amount of concentrate produced
      """
      
      self.oreProcessed = np.zeros(len(mineDataManager.theMiningSystem.oreMined)) 
      self.processingPower = np.zeros(len(mineDataManager.theMiningSystem.oreMined))   
      self.processingCapacity = mineDataManager.theMiningSystem.mineOreProductionCapacity # ore is processed at a constant rate
      carryOver = 0.0
      for year in range( len(mineDataManager.theMiningSystem.oreMined )-1 ):
        processedOre = carryOver + mineDataManager.theMiningSystem.oreMined[year]
        
        if(processedOre > self.processingCapacity):
          carryOver = processedOre - self.processingCapacity
          processedOre = self.processingCapacity
        else:
          carryOver = 0.0
        self.oreProcessed[year] = processedOre
      
      self.oreProcessed[-1] = carryOver +  mineDataManager.theMiningSystem.oreMined[-1] # final year
      
      
      # fixme - need to make this an input function
      # convert tonnes processed each year to the number of Mwh based on powerlaw fit
      self.processingPower = 3.96*(self.oreProcessed  )**0.703  # in Mwh
      
      referenceMetalStr = self.processingMethod # mineDataManager.theOreBody.type[:2]  
      
      # first two letters of orebody type is assumed to be reference metal for determining processing grade
      # eg AuCu -> gold is reference metal - note that user must select correct method
      """
      if(referenceMetalStr[0] == "K"):
        referenceMetalStr = "K"  # fixme change to named type
      elif(referenceMetalStr[0] == "P" and referenceMetalStr[1] != "b"):
        referenceMetalStr = "P"  # fixme change to named type
      """
      
      
      referenceMetalOreConcentration = mineDataManager.theOreBody.metalGrades[self.processingMethod]

      self.concentrateMetalConcentration = 1.0
      
      # lookup concentrateMetalConcentrations based on reference metal type
      
      concentrateConcentrations = {"Au":0.75,"Ag":0.85,"Ni":0.1,"Cu":0.25,"Pb":0.5,"Zn":0.5,
                                   "P2O5":0.3,"K2SO4":1.0,"REO":0.2}  
                                   # P assumes concentrated to 30% (indicates P2O5 grade required for conversion to marketable product)
                                   # K assumes concentrated to 100% (indicates pure K2SO4) - alt could use 50% K grade.
                                   # alt could set to 100 and assume K2SO4 grade 
                                   # these will need to be updated with REOs etc
      
      # find the minimum amount of concentration needed to bring concentrate to market
      minConcentrationFactor = 1e64
      
      #for metal,metalOreGrade in mineDataManager.theOreBody.metalGrades.items():
      #  if metal in concentrateConcentrations:
      for metal in self.concentrateMetals:
        if metal in concentrateConcentrations:
          metalOreGrade = mineDataManager.theOreBody.metalGrades[metal]
          concentrateGrade = concentrateConcentrations[metal]
          concFactor = concentrateGrade/(metalOreGrade/(1.0+ mineDataManager.theMiningSystem.dilution) +1e-64)
          if concFactor < 1.0:
            concFactor = 1.0
            concentrateGrade = metalOreGrade
          #print "concFactor", metal, concFactor, metalOreGrade, concentrateGrade
          if(concFactor < minConcentrationFactor ):
            minConcentrationFactor = concFactor
            self.concentrateMetalConcentration = concentrateGrade
            self.concentratePrimaryMetal = metal
      
      # concentrate is calculated based on estimate of mineral content
      self.concentrateProduced = (1.0 - self.processingLoss) * np.array(self.oreProcessed)/minConcentrationFactor    
      
      if(self.secondaryProcessingSystem):  # secondary+ processing circuits
        secondaryOreProcessed = self.oreProcessed - self.concentrateProduced
        self.secondaryProcessingSystem.CalculateSecondaryProcessingCapacity(problemManager, mineDataManager,secondaryOreProcessed)
      
      return self.processingCapacity
    
    def CalculateSecondaryProcessingCapacity(self,  problemManager, mineDataManager,processedOre):
      """
      Calculate processing capacity, amount of ore processed each year and amount of concentrate produced
      """
      # TO DO: Include optional flag to optimise capacity to minimise costs
      self.oreProcessed = np.array(processedOre)
      self.processingCapacity = self.oreProcessed.max() # all ore is assumed processed in the same year for secondary circuits
      carryOver = 0.0  # should be no carry over. 
      
      # fixme - need to make this an input function
      # convert tonnes processed each year to the number of Mwh based on powerlaw fit
      self.processingPower = 3.96*(self.oreProcessed  )**0.703  # in Mwh
      
      referenceMetalStr = self.processingMethod 
      self.concentratePrimaryMetal = referenceMetalStr
      
      
      referenceMetalOreConcentration = mineDataManager.theOreBody.metalGrades[referenceMetalStr]  
      # NB. this assumes that secondary ore has the same concentration as the orebody itself

      self.concentrateMetalConcentration = 1.0
      
      # lookup concentrateMetalConcentrations based on reference metal type
      
      concentrateConcentrations = {"Au":0.75,"Ag":0.85,"Ni":0.1,"Cu":0.25,"Pb":0.5,"Zn":0.5,
                                   "P2O5":0.3,"K2SO4":0.5,"REO":0.2}  # these will need to be updated with REEs etc
      
      # find the minimum amount of concentration needed to bring secondary reference metal to market
      minConcentrationFactor = 1e64 #
      
      for metal in self.concentrateMetals:
        if referenceMetalStr in concentrateConcentrations:
            metalOreGrade = mineDataManager.theOreBody.metalGrades[metal]
            concentrateGrade = concentrateConcentrations[referenceMetalStr]
            concFactor = concentrateGrade/(metalOreGrade +1e-64)  # nb. no dilution assumed for secondary circuits
            if concFactor < 1.0:
              concFactor = 1.0
              concentrateGrade = metalOreGrade
            #print("concFactor", metal, concFactor, metalOreGrade, concentrateGrade)
            if(concFactor < minConcentrationFactor ):
              minConcentrationFactor = concFactor
              self.concentrateMetalConcentration = concentrateGrade
      
      # concentrate is calculated based on estimate of mineral content
      self.concentrateProduced = (1.0 - self.processingLoss) \
                                  * self.oreProcessed/minConcentrationFactor    
      
      if(self.secondaryProcessingSystem):  # tertiary+ processing circuits
        secondaryOreProcessed = self.oreProcessed - self.concentrateProduced
        self.secondaryProcessingSystem.CalculateSecondaryProcessingCapacity(problemManager, mineDataManager, secondaryOreProcessed)
      
      return self.processingCapacity
    
      
    def CalculateProcessingCapex(self, problemManager,  mineDataManager):
      
      self.processingCapex = np.zeros( mineDataManager.theMiningSystem.mineLife )
      
      theFunctionManager =  FunctionManager()
      
      processingCapexFunc = theFunctionManager.GetFunction("ProcessingCapex_" + self.processingMethod)
      
      theUnitManager = UnitManager()
      self.processingCapex[0] =  processingCapexFunc.f( [self.processingCapacity *  theUnitManager.ConvertTo("1e6 tonne")] )
      print("Processing Capex ", self.processingCapex[0])
      
      if (self.secondaryProcessingSystem):
        secondaryCapex = self.secondaryProcessingSystem.CalculateProcessingCapex(problemManager,  mineDataManager)
        self.processingCapex += secondaryCapex
      
      return self.processingCapex

    def ZeroProcessingCapex(self, problemManager,  mineDataManager):
      
      self.processingCapex = np.zeros( mineDataManager.theMiningSystem.mineLife )
      
      if (self.secondaryProcessingSystem):
        self.secondaryProcessingSystem.ZeroProcessingCapex(problemManager,  mineDataManager)
      
      return self.processingCapex

      
    def CalculateProcessingOpex(self,  problemManager, mineDataManager):
      
      self.processingOpex = np.zeros( mineDataManager.theMiningSystem.mineLife )
      
      theFunctionManager =  FunctionManager()
      theUnitManager = UnitManager()
      
      self.processingOpex = self.oreProcessed*theUnitManager.ConvertTo("tonne")   # tonnes processed per year
      
      processingOpexFunc = theFunctionManager.GetFunction("ProcessingOpex_" + self.processingMethod)
      opexPerTonne =  processingOpexFunc.f( [self.processingCapacity  *  theUnitManager.ConvertTo("1e6 tonne")] )
      
      print("Processing method: ", self.processingMethod)
      #print "Processing capacity (kg ore)", self.processingCapacity
      print("Processing capacity (Mt ore)", self.processingCapacity  *  theUnitManager.ConvertTo("1e6 tonne") )
      print("Processing Opex Per Tonne ", opexPerTonne )
      
      self.processingOpex *= opexPerTonne
      
      if (self.secondaryProcessingSystem):
        secondaryOpex = self.secondaryProcessingSystem.CalculateProcessingOpex(problemManager,  mineDataManager)
        self.processingOpex += secondaryOpex
      
      return self.processingOpex
            
    def GetTotalConcentrateProduced(self):
      rv = np.array(self.concentrateProduced)
      
      if(self.secondaryProcessingSystem):
        rv += self.secondaryProcessingSystem.GetTotalConcentrateProduced()
      return rv
      
      
    