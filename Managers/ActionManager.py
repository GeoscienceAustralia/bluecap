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
from IO.XML import GetAttributeString,GetAttributeFileString
from IO.XML import HasAttribute,GetAttributeStringOrDefault,GetAttributeValueOrDefault
from IO.XML import SetAttributeString
from IO.XML import GetXMLTag
from IO.XML import PrettyXMLFormat,AddNote


from Common.Common import BluecapError

class ActionManager(object):
      """Manager controlling which actions are run in the problem manager."""
 
      def __init__(self):
            """ Create an empty ActionManager object and default variables. """
            self.actions = []
            self.activeIndex = 0
            
  
      def HasCompleted(self):
        """Returns true if all actions in this action manager have been completed."""
        rv = len(self.actions) <= self.activeIndex
        return rv
        
        
      def RunNextAction(self, problemManager):
            """Run the next actions in the list of actions."""
            if(self.activeIndex <  len(self.actions)):
              self.actions[self.activeIndex].Run(problemManager)
              self.activeIndex += 1

      def RunUntil(self, finalAction, problemManager):
            """Run actions until the final action is reached."""
            if(self.activeIndex <  len(self.actions) and self.activeIndex <= finalAction):
              self.actions[self.activeIndex].Run(problemManager)
              self.activeIndex += 1
            
      def ParseXMLNode(self, acManagerNode, problemManager):
        """Generate actions from the action manager xml tree node."""
      
        for child in GetChildren(acManagerNode):
          type = GetXMLTag(child)
          
          if(type != "note"):
            self.actions.append(  ActionFactory.CreateFromXML(type,child,problemManager) )
          

      def WriteXMLNode(self, node):
        """Write actions to xml node."""
      
        # functions
        for action in self.actions:
          type = action.typeStr
          funcNode = AddChild(node,type)
          action.WriteXMLNode(funcNode)
      
        return node

      
####################

## Action factory

class ActionFactory:
    """Factory class used to generate action objects for the action manager."""
    
    factories = {}
    
    @staticmethod
    def AddFactory(id, factory):
        """Add factory to the action factory manager."""
        ActionFactory.factories[id] = factory
    
    @staticmethod
    def CreateFromXML(id,xmlNode,problemManager):
        """Run factory to create Action 'id' based on input from xml node."""
        if not (id in ActionFactory.factories): 
            raise BluecapError("Error: Could not find " + id + " in ActionFactory")
            
        return ActionFactory.factories[id].CreateFromXML(xmlNode,problemManager)
