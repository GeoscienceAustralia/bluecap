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
from IO.XML import GetAttributeString,HasAttribute
from IO.XML import SetAttributeString

# Managers
from Managers.ActionManager import ActionManager
from Managers.ActionManager import ActionFactory

from Managers.ParameterManager import ParameterManager

import pylab as pl
import numpy as np

####################

# Run regional calculation of a specific type.

class RegionalSensitivityAction():  
    typeStr = "RunRegionalSensitivity"

    def __init__(self):   
       self.type = "NPV"
       self.parameter = None
       self.minValue = None
       self.maxValue = None
       self.default = None
       
    def Run(self,problemManager):
        from Managers.ProblemManager import ProblemManager   # avoiding circular dependency
        
        # plotting
        doPlots = (problemManager.outputType == "plots")
        
        problemManager.theRegionalCalculationManager.type = self.type
        #problemManager.theRegionalCalculationManager.Run(problemManager)
        
        theParameterManager = ParameterManager()
        originalValue = theParameterManager.GetParameter(self.parameter)
        
        prevAction = problemManager.theActionManager.activeIndex - 1
        
        # reset params for max value
        theParameterManager.SetParameter(self.parameter,self.maxValue)
        
        # dummy manager based on run max value - nb may be able to use xml tree in prob manager rather than reloading...
        dummyProblemManager = ProblemManager.FromXMLFile(problemManager.inputFile)
        dummyProblemManager.Initialize()
        
        # rerun to this point
        dummyProblemManager.RunUntil(prevAction)
        
        # run regional calculation and record max value
        dummyProblemManager.theRegionalCalculationManager.type = self.type
        dummyProblemManager.theRegionalCalculationManager.saveResult = False
        maxValues = dummyProblemManager.theRegionalCalculationManager.Run(dummyProblemManager)
        
        # reset params for min value
        theParameterManager.SetParameter(self.parameter,self.minValue)
        dummyProblemManager = ProblemManager.FromXMLFile(problemManager.inputFile)
        
        # reset params
        dummyProblemManager.Initialize()
        
        # rerun to this point
        dummyProblemManager.RunUntil(prevAction)
        dummyProblemManager.theRegionalCalculationManager.type = self.type
        
        # run regional calculation and record min value
        minValues = dummyProblemManager.theRegionalCalculationManager.Run(dummyProblemManager)
        dummyProblemManager.theRegionalCalculationManager.saveResult = False
        
        diff = maxValues - minValues #np.abs(maxValues - minValues)
        if problemManager.theRegionalCalculationManager.type == 'NPV':
          diff *= 1e-6
        
        # fixme need to modify output control 
        filename = problemManager.outputPrefix+"_snstvty_"+self.parameter+"."+problemManager.outputType
        problemManager.theRegionalCalculationManager.type = self.parameter
        problemManager.theRegionalCalculationManager.SaveMap(filename, diff,problemManager.recordRange)
        problemManager.theRegionalCalculationManager.type = self.type
        
        theParameterManager.SetParameter(self.parameter,originalValue)
        
        if doPlots:
          pl.imshow(diff,origin="lower")
          pl.colorbar()
          pl.show()
           
    @classmethod
    def CreateFromXML(cls,actionNode,problemManager):
      # this creates the class from an xml node by initiating an empty class and parsing the xml node
      rv = cls()
      rv.ParseXMLNode(actionNode,problemManager)
      return rv
    
    def ParseXMLNode(self,actionNode,problemManager):
       """
       Generate named action from xml tree node. 
       """
       if(HasAttribute(actionNode,"type")):
         self.type = GetAttributeString(actionNode,"type")
         
       self.parameter = GetAttributeString(actionNode,"parameter")
       self.minValue = GetAttributeString(actionNode,"min")
       self.maxValue = GetAttributeString(actionNode,"max")
    
    def WriteXMLNode(self, node): 
       if(self.type):
         SetAttributeString(node,"type",self.type)
       return node
       
    # The factory that creates the RegionalSensitivityAction Object 
    class Factory:
      def CreateFromXML(self,xmlNode,problemManager): return RegionalSensitivityAction.CreateFromXML(xmlNode,problemManager)


# Factory Registrator
ActionFactory.AddFactory(RegionalSensitivityAction.typeStr,RegionalSensitivityAction.Factory())
