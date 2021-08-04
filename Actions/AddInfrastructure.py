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
from IO.XML import GetAttributeString,GetAttributeVector
from IO.XML import HasAttribute
from IO.XML import SetAttributeString

# Managers
from Managers.ActionManager import ActionManager
from Managers.ActionManager import ActionFactory


####################

# Action to add infrastructure to the regional calculation map

class AddInfrastructureAction():  
    typeStr = "AddInfrastructure"

    def __init__(self):   
       self.infrastructureType = ""
       self.coords = []
       
    def Run(self,problemManager):
        """Run the 'AddInfrastructure' action"""
        problemManager.theRegionalCalculationManager.AddInfrastructure(self.infrastructureType,self.coords)
        return
           
    @classmethod
    def CreateFromXML(cls,actionNode,problemManager):
      """Create the AddInfrastructureAction class from an xml node by initiating an empty class and parsing the xml node."""
      rv = cls()
      rv.ParseXMLNode(actionNode,problemManager)
      return rv
    
    def ParseXMLNode(self,actionNode,problemManager):
       """Generate AddInfrastructureAction action from xml tree node."""
       self.infrastructureType = GetAttributeString(actionNode,"type")
       self.coords = GetAttributeVector(actionNode,"coords")
    
    def WriteXMLNode(self, node): 
       """Export the AddInfrastructureAction action to an xml tree node."""
       SetAttributeString(node,"type",self.infrastructureType)
       SetAttributeVector(node,"coords",self.coords)
       return node
        
    class Factory:
      """The factory that creates the AddInfrastructureAction Object."""
      def CreateFromXML(self,xmlNode,problemManager): return AddInfrastructureAction.CreateFromXML(xmlNode,problemManager)


# Factory Registrator
ActionFactory.AddFactory(AddInfrastructureAction.typeStr,AddInfrastructureAction.Factory())