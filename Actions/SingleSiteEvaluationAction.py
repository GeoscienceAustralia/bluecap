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


####################

class SingleSiteEvaluationAction():  
    typeStr = "RunSiteEvaluation"

    def __init__(self):   
       self.dummy = ""
       
    def Run(self,problemManager):
        problemManager.EvaluateOCUGMines()
           
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
       # empty
       return
    
    def WriteXMLNode(self, node): 
       # empty
       return node
       
     # The factory that creates the ExitAction Object 
    class Factory:
      def CreateFromXML(self,xmlNode,problemManager): return SingleSiteEvaluationAction.CreateFromXML(xmlNode,problemManager)


# Factory Registrator
ActionFactory.AddFactory(SingleSiteEvaluationAction.typeStr,SingleSiteEvaluationAction.Factory())
