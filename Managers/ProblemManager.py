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

# IO
from IO.XML import NewXMLTree,GetXMLTreeRoot,SaveXMLFile,LoadXMLFile
from IO.XML import HasChild,GetChild,GetChildren,AddChild
from IO.XML import GetAttributeString,GetAttributeFileString,GetAttributeStringOrDefault,GetAttributeValueOrDefault,SetAttributeString
from IO.XML import PrettyXMLFormat

# Managers
from MineDataManager import MineDataManager
from ParameterManager import ParameterManager
from RegionalCalculationManager import RegionalCalculationManager
from Functions.FunctionManager import FunctionManager


# deep copy
from copy import deepcopy

class ProblemManager():
    """
      The Problem Manager  holds data and methods describing the problem.
    """
    
    def __init__(self):
      """
      Initialize an empty problem manager
      """
      self.theMineDataManager = MineDataManager()
      self.theRegionalCalculationManager = RegionalCalculationManager()
      
      
      self.outputPrefix = "output"
      self.outputType = ""
      self.recordRange = False
      
    @classmethod
    def FromXMLFile(cls,filename):
      """
      Initialize problem manager from xml file
      """
      rv = cls()
      
      xmlTreeRoot = LoadXMLFile(filename)
      rv.ParseXMLNode(xmlTreeRoot)
      
      return rv
    
    def ExportXMLFile(self,filename):
      """
      Write problem to xml file
      """
      problemTree = NewXMLTree()
      problemTreeRoot = GetXMLTreeRoot(problemTree) 
      self.WriteXMLNode(problemTreeRoot)
      PrettyXMLFormat(problemTreeRoot)
      SaveXMLFile(filename,problemTree)
      
    def ParseXMLNode(self, rootNode):
      """
      Generate Problem Manager data from xml tree node. 
      """
      
      # remove "Problem" branch if present
      if(HasChild(rootNode,"Problem")):
        root = GetChild(rootNode,"Problem")
        
      # load parameters not set from command line
      # done here to prevent recursion
      theParameterManager = ParameterManager()  
      if(HasChild(rootNode,"Parameters")):
        parameterNode = GetChild(rootNode,"Parameters")
        
        for child in GetChildren(parameterNode):
        
          name = GetAttributeString(child,"name")
          
          if(not theParameterManager.HasParameter(name) ):
            paramString = GetAttributeString(child,"value")
            theParameterManager.SetParameter(name,paramString)
            print "param: ", name, " value:",paramString
        
      # create functions before other classes (after UnitManager & Parameter Manager)
      theFunctionManager = FunctionManager()  
      if(HasChild(rootNode,"Functions")):
        functionNode = GetChild(rootNode,"Functions")
        theFunctionManager.ParseXMLNode(functionNode,self) 
      
      # Mine data
      if(HasChild(rootNode,"MineData")):
        mineDataNode = GetChild(rootNode,"MineData")
        self.theMineDataManager.ParseXMLNode(mineDataNode) 
        
       
      # Output - eventually will make into a manager
      if(HasChild(rootNode,"Output")):
        outputNode = GetChild(rootNode,"Output")
        self.outputPrefix = GetAttributeFileString(outputNode,"prefix")
        self.outputType = GetAttributeStringOrDefault(outputNode,"type","")
        
        self.recordRange = GetAttributeValueOrDefault(outputNode,"recordRange",False)
        
      # Regional calculation - optional
      if(HasChild(rootNode,"RegionalCalculation")):
        regionalCalculationNode = GetChild(rootNode,"RegionalCalculation")
        self.theRegionalCalculationManager.ParseXMLNode(regionalCalculationNode) 
    
    

    def WriteXMLNode(self, root):
      """
      Write problem to xml node
      """
      
      #functions
      theFunctionManager = FunctionManager()  
      funcMngrNode = AddChild(root,"Functions")
      theFunctionManager.WriteXMLNode(funcMngrNode) 
      
      # Mine data
      mineNode = AddChild(root,"MineData")
      self.theMineDataManager.WriteXMLNode(mineNode)
      
      # Output
      outputNode = AddChild(root,"Output")
      SetAttributeString(outputNode,"prefix",self.outputPrefix)
      
      return root
    
    def Initialize(self):
      """
      Set initial conditions/data for classes after initial data structure is loaded. 
      """
      
      return None  

    def EvaluateOCUGMines(self):
    
      # ugly - add dictionary once stable
      
      # Calculate mine value (OC vs UG)
      print "**** Underground calculation *****"
      self.theMineDataManager.SetMiningMethod("UG")
      UGValue = self.theMineDataManager.CaculateMineProductionAndValue(self)
      self.theUGMineDataManager = deepcopy(self.theMineDataManager)
      
      print "**** Open pit calculation *****"
      self.theMineDataManager.SetMiningMethod("OC")
      OCValue = self.theMineDataManager.CaculateMineProductionAndValue(self)
      self.theOCMineDataManager = deepcopy(self.theMineDataManager)
      
      
      print " OC vs UG"
      print "UGValue",UGValue
      print "OCValue",OCValue
      
      mineValue = UGValue
      if(UGValue < OCValue and OCValue > 0.0):
        self.theMineDataManager = self.theOCMineDataManager
        mineValue = OCValue
      else:
        self.theMineDataManager = self.theUGMineDataManager
      
      return mineValue
    
    
    def Run(self):
      """
      Run Problem
      """
      # Eventually will roll these into separate routines called from the input file
      
      if(self.theRegionalCalculationManager.type):
        self.theRegionalCalculationManager.Run(self)
      else:
        self.EvaluateOCUGMines()
        
      
      
      return None  
      
      

        