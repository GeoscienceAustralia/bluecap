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



# IO
from IO.XML import NewXMLTree,GetXMLTreeRoot,SaveXMLFile,LoadXMLFile
from IO.XML import HasChild,GetChild,GetChildren,AddChild
from IO.XML import GetAttributeString,GetAttributeFileString,GetAttributeStringOrDefault,GetAttributeValueOrDefault,SetAttributeString
from IO.XML import PrettyXMLFormat, WriteNotes

# Managers
from .MineDataManager import MineDataManager
from .ParameterManager import ParameterManager
from .RegionalCalculationManager import RegionalCalculationManager
from Functions.FunctionManager import FunctionManager
from .ActionManager import ActionManager
from .SolverManager import SolverManager


from .HydrogenDataManager import HydrogenDataManager

# actions
import Actions

# solvers
import Solvers

# Common
from Common.Version import GetGitSha

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
      self.theHydrogenDataManager = HydrogenDataManager()
      self.theRegionalCalculationManager = RegionalCalculationManager()
      self.theActionManager = None
      
      self.outputPrefix = "output"
      self.outputType = ""
      self.recordRange = False
      self.xmlTreeRoot = None
      self.inputFile = None
      
      
    @classmethod
    def FromXMLFile(cls,filename):
      """
      Initialize problem manager from xml file
      """
      rv = cls()
      
      rv.xmlTreeRoot = LoadXMLFile(filename)
      rv.inputFile = filename
      rv.ParseXMLNode(rv.xmlTreeRoot)
      
      return rv
    
    def ExportXMLFile(self,filename,recordSha):
      """
      Write problem to xml file
      """
      problemTree = NewXMLTree()
      problemTreeRoot = GetXMLTreeRoot(problemTree) 
      self.WriteXMLNode(problemTreeRoot,recordSha)
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
            print("param: ", name, " value:",paramString)
        
      # create functions before other classes (after UnitManager & Parameter Manager)
      theFunctionManager = FunctionManager()  
      if(HasChild(rootNode,"Functions")):
        functionNode = GetChild(rootNode,"Functions")
        theFunctionManager.ParseXMLNode(functionNode,self) 
      
      # Mine data
      if(HasChild(rootNode,"MineData")):
        mineDataNode = GetChild(rootNode,"MineData")
        self.theMineDataManager.ParseXMLNode(mineDataNode) 
      
      # Processing data
      if(HasChild(rootNode,"Processing")):
        processingNode = GetChild(rootNode,"Processing")
        self.theMineDataManager.theProcessingSystem.ParseXMLNode(processingNode)
       
      # Hydrogen data
      if(HasChild(rootNode,"HydrogenData")):
        hydrogenDataNode = GetChild(rootNode,"HydrogenData")
        self.theHydrogenDataManager.ParseXMLNode(hydrogenDataNode) 
        
      # Output - eventually will make into a manager
      if(HasChild(rootNode,"Output")):
        outputNode = GetChild(rootNode,"Output")
        self.outputPrefix = GetAttributeFileString(outputNode,"prefix")
        self.outputType = GetAttributeStringOrDefault(outputNode,"type","")
        
        self.recordRange = GetAttributeValueOrDefault(outputNode,"recordRange",False)
        
      # Regional calculation - optional - eventually will make into a solver
      if(HasChild(rootNode,"RegionalCalculation")):
        regionalCalculationNode = GetChild(rootNode,"RegionalCalculation")
        self.theRegionalCalculationManager.ParseXMLNode(regionalCalculationNode) 
    
      # Solver manager (after UnitManager & Parameter Manager & Function Manager)
      if(HasChild(rootNode,"Solvers")):
        solversNode = GetChild(rootNode,"Solvers")
        theSolverManager = SolverManager()
        theSolverManager.ParseXMLNode(self,solversNode) 
    
      # Create actions after other classes (after UnitManager & Parameter & Solver Manager)
      self.theActionManager = ActionManager()  
      if(HasChild(rootNode,"Actions")):
        actionNode = GetChild(rootNode,"Actions")
        self.theActionManager.ParseXMLNode(actionNode,self)   
      else:
      
        # Add actions for backwards compatibility
        # A) If regional manager present - run regional calculation
        # B) Otherwise add single site calculation
        
        # This will be depreciated (need to issue warning)
        if(self.theRegionalCalculationManager.type):
            self.theActionManager.actions.append( RegionalCalculationAction() )
        else:
            self.theActionManager.actions.append( SingleSiteEvaluationAction() )

        
        #print "Error: did not encounter Actions block in xml file."
        #exit(1)
        



    def WriteXMLNode(self, root,recordSha=True):
      """
      Write problem to xml node
      """
      if(recordSha):
        theGitSha = GetGitSha()
        SetAttributeString(root,"sha",theGitSha)
      
      #Functions
      theFunctionManager = FunctionManager()  
      funcMngrNode = AddChild(root,"Functions")
      theFunctionManager.WriteXMLNode(funcMngrNode) 
      
      # Mine data
      mineNode = AddChild(root,"MineData")
      self.theMineDataManager.WriteXMLNode(mineNode)
      
      # Output
      outputNode = AddChild(root,"Output")
      SetAttributeString(outputNode,"prefix",self.outputPrefix)
      
      # Actions
      #theActionManager = ActionManager()        
      actnMngrNode = AddChild(root,"Actions")
      self.theActionManager.WriteXMLNode(actnMngrNode) 
      
      return root
    
    def SaveNotes(self,filename):
      """Export notes included in input file (and added in the code) to a named text file."""
      
      fid = open(filename,"wt")
    
      WriteNotes(fid,self.xmlTreeRoot)
      
      fid.close()
    
    
    def Initialize(self):
      """
      Set initial conditions/data for classes after initial data structure is loaded. 
      """
      
      return None  

    def EvaluateOCUGMines(self):
      """Calculate value of open cut and underground mining projects."""
    
      if(self.theMineDataManager.theMiningSystem.miningMethod == "K2SO4"):
        mineValue = self.theMineDataManager.CaculateMineProductionAndValue(self)
      else:
          # ugly - add dictionary once stable
      
          # Calculate mine value (OC vs UG)
          print("**** Underground calculation *****")
          self.theMineDataManager.SetMiningMethod("UG")
          UGValue = self.theMineDataManager.CaculateMineProductionAndValue(self)
          self.theUGMineDataManager = deepcopy(self.theMineDataManager)
      
          if self.theMineDataManager.theOreBody.cover < 2000:
            print("**** Open pit calculation *****")
            self.theMineDataManager.SetMiningMethod("OC")
            OCValue = self.theMineDataManager.CaculateMineProductionAndValue(self)
            self.theOCMineDataManager = deepcopy(self.theMineDataManager)
          else:
            print("**** SKIPPING Open pit calculation *****")
            OCValue = UGValue - 1.
      
          print(" OC vs UG")
          print("UGValue",UGValue)
          print("OCValue",OCValue)
      
          mineValue = UGValue
          if(UGValue < OCValue and OCValue > 0.0):
            self.theMineDataManager = self.theOCMineDataManager
            mineValue = OCValue
          else:
            self.theMineDataManager = self.theUGMineDataManager
      
      return mineValue
    
    
    def EvaluateHydrogenProject(self):
      """Calculate value hydrogen projects."""
    
      projectValue = self.theHydrogenDataManager.CaculateHydrogenProductionAndValue(self)
      
      return projectValue
    
    def Run(self):
      """
      Run Problem
      """
      
      while True:
        self.theActionManager.RunNextAction(self)
        if self.theActionManager.HasCompleted():
          break
          
      return None  
      
      
    def RunUntil(self,count):
      """
      Run action manager until count is reached
      """
      
      self.theActionManager.RunUntil(count,self)
          
      return None  
        