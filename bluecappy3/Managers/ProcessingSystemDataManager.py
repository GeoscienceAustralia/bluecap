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
from IO.XML import GetAttributeString,GetAttributeValue, SetAttributeString


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
      self.processingPower = []
      

    def ParseXMLNode(self, processingSystemNode):
      """
      Generate processing system data from xml tree node. 
      """
      
      return processingSystemNode
        

    def WriteXMLNode(self, node):
      """
      Write processing system data to xml node
      """
      # Mine data
      SetAttributeString(node,"method",self.processingMethod)
      SetAttributeString(node,"processingLoss",self.processingLoss)
      SetAttributeString(node,"refiningTake",self.refiningTake)
      
      return node
      
    def DetermineProcessingSystem(self, problemManager, mineDataManager):
      """
      Use available data to determine the most likely processing method/capacity/opex etc.
      """

      self.CalculateProcessingCapacity(problemManager, mineDataManager)
      
      #Todo("determine processing method based on amount and type of ore mined")
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
      elif(referenceMetalStr ==  "Pb"):
        self.refiningTake = 0.17
        
      
      self.CalculateProcessingCapex(problemManager, mineDataManager)
      self.CalculateProcessingOpex(problemManager, mineDataManager)
      
      return self      
    
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
      
      
      # convert tonnes processed each year to the number of Mwh based on powerlaw fit
      self.processingPower = 3.96*(self.oreProcessed  )**0.703  # in Mwh
      
      referenceMetalStr = mineDataManager.theOreBody.type[:2]  
      # first two letters of orebody type is assumed to be reference metal for determining processing grade
      # eg AuCu -> gold is reference metal - note that user must select correct method
      
      
      referenceMetalOreConcentration = mineDataManager.theOreBody.metalGrades[referenceMetalStr]

      self.concentrateMetalConcentration = 1.0
      
      # lookup concentrateMetalConcentrations based on reference metal type
      
      concentrateConcentrations = {"Au":0.75,"Ag":0.85,"Ni":0.1,"Cu":0.25,"Pb":0.5}
      
      # find the minimum amount of concentration needed to bring concentrate to market
      minConcentrationFactor = 1e64
      
      for metal,metalOreGrade in mineDataManager.theOreBody.metalGrades.items():
        if metal in concentrateConcentrations:
          concentrateGrade = concentrateConcentrations[metal]
          concFactor = concentrateGrade/(metalOreGrade/(1.0+ mineDataManager.theMiningSystem.dilution) +1e-64)
          if concFactor < 1.0:
            concFactor = 1.0
          #print "concFactor", metal, concFactor, metalOreGrade, concentrateGrade
          if(concFactor < minConcentrationFactor ):
            minConcentrationFactor = concFactor
            self.concentrateMetalConcentration = concentrateGrade
      
      # concentrate is calculated based on estimate of mineral content
      self.concentrateProduced = (1.0 - self.processingLoss) \
                                  *np.array(mineDataManager.theMiningSystem.oreMined)/minConcentrationFactor    
      
      
      return self.processingCapacity
      
    def CalculateProcessingCapex(self, problemManager,  mineDataManager):
      
      self.processingCapex = np.zeros( mineDataManager.theMiningSystem.mineLife )
      
      theFunctionManager =  FunctionManager()
      
      processingCapexFunc = theFunctionManager.GetFunction("ProcessingCapex_" + self.processingMethod)
      
      theUnitManager = UnitManager()
      self.processingCapex[0] =  processingCapexFunc.f( [self.processingCapacity *  theUnitManager.ConvertTo("1e6 tonne")] )
      print("Processing Capex ", self.processingCapex[0])
      
      return self.processingCapex
      
    def CalculateProcessingOpex(self,  problemManager, mineDataManager):
      
      self.processingOpex = np.zeros( mineDataManager.theMiningSystem.mineLife )
      
      theFunctionManager =  FunctionManager()
      theUnitManager = UnitManager()
      
      self.processingOpex = self.oreProcessed*theUnitManager.ConvertTo("tonne")   # tonnes processed per year
      
      processingOpexFunc = theFunctionManager.GetFunction("ProcessingOpex_" + self.processingMethod)
      opexPerTonne =  processingOpexFunc.f( [self.processingCapacity  *  theUnitManager.ConvertTo("1e6 tonne")] )
      
      print("Processing method: ", self.processingMethod)
      print("Processing capacity (Mt ore)", self.processingCapacity  *  theUnitManager.ConvertTo("1e6 tonne"))
      print("Processing Opex Per Tonne ", opexPerTonne)
      
      self.processingOpex *= opexPerTonne
      
      return self.processingOpex
      