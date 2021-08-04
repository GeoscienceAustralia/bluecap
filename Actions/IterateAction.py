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
from IO.XML import GetChildren,AddChild
from IO.XML import GetAttributeValueOrDefault
from IO.XML import SetAttributeString
from IO.XML import GetXMLTag

# Managers
from Managers.ActionManager import ActionManager
from Managers.ActionManager import ActionFactory


####################

# Iterate over a loop.

class IterationAction():  
    typeStr = "Iterator"

    def __init__(self):   
       self.number = 1
       self.actions = []
       
    def Run(self,problemManager):
        """Run the iterator action."""
        for i in range(self.number):
          for action in self.actions:
            action.Run(problemManager)
           
    @classmethod
    def CreateFromXML(cls,actionNode,problemManager):
      """ This creates the class from an xml node by initiating an empty class and parsing the xml node."""
      rv = cls()
      rv.ParseXMLNode(actionNode,problemManager)
      return rv
    
    def ParseXMLNode(self,actionNode,problemManager):
        """
        Generate an Iterator action from the xml tree node. 
        """
      
        self.number = GetAttributeValueOrDefault(actionNode,"number",self.number)
        for child in GetChildren(actionNode):
          type = GetXMLTag(child)
          self.actions.append(  ActionFactory.CreateFromXML(type,child,problemManager) )
        return
    
    def WriteXMLNode(self, node): 
       """Write the iterator action to an xml tree node."""
       SetAttributeString(node,"number",self.number)
       for action in self.actions:
          type = action.typeStr
          actnNode = AddChild(node,type)
          action.WriteXMLNode(actnNode)
       return node
       
    class Factory:
      """The factory that creates the Iterator Action Object """
      def CreateFromXML(self,xmlNode,problemManager): return IterationAction.CreateFromXML(xmlNode,problemManager)


# Factory Registrator
ActionFactory.AddFactory(IterationAction.typeStr,IterationAction.Factory())
        