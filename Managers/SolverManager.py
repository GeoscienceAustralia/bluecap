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
from IO.XML import PrettyXMLFormat

# Managers
from Functions.FunctionManager import FunctionManager

# Common
from Common.Common import BluecapError

# deep copy
from copy import deepcopy



class SolverManager(object):
    """ 
    Object that controls system solvers.
    
    Usage:
      from Managers.SolverManager import SolverManager
      theSolverManager = SolverManager()
      theSolverManager.RunSolver(problemManager,solverName)
    """
 
    class __SolverManager:
      def __init__(self):
            self.solvers = {}
            
      def ParseXMLNode(self, problemManager,solvManagerNode):
        """
        Generate Solvers from xml tree node. 
        """
      
        for child in GetChildren(solvManagerNode):
          type = GetXMLTag(child)
          name = GetAttributeString(child,"name")
          self.solvers[name] = SolverFactory.CreateFromXML(type,child,problemManager) 
        
      def WriteXMLNode(self, node):
        """
        Write solvers to xml node.
        """
      
        # solvers
        for name,solver in self.solvers.items():
          type = solver.typeStr
          solvNode = AddChild(node,type)
          solvNode.SetAttributeString("name",name)
          solver.WriteXMLNode(solvNode)
      
        return node
          
          
      def RunSolver(self, solverName,problemManager,args):
        """Run the named solver."""
        self.solvers[solverName].Run(problemManager,args)
          
    instance = None
    
    def __new__(cls): # __new__ always a classmethod
        if not SolverManager.instance:
            SolverManager.instance = SolverManager.__SolverManager()
        return SolverManager.instance
        
    def __getattr__(self, name):
        return getattr(self.instance, name)
    def __setattr__(self, name):
        return setattr(self.instance, name)
      
      
####################

## Solver factory

class SolverFactory:
    """
    Factory class used to generate solver instances for the Solver Manager.
    
    Note: 
    To be available, Solvers must be registered with the SolverFactory using:
    SolverFactory.AddFactory(UpstreamImpactCalculator.typeStr,UpstreamImpactCalculator.Factory())
    
    As the "AddFactory" call is usually made within the script defining the solver, and 
    may not be referenced elsewhere, this implies that the script containing the solver 
    must be imported explicitly, for example in the "__init__.py" script. 
    eg:
    import Solvers.UpstreamImpactCalculator
    
    """
    factories = {}
    
    @staticmethod
    def AddFactory(id, factory):
        """Add factory to the solver factory manager"""
        SolverFactory.factories[id] = factory
    
    @staticmethod
    def CreateFromXML(id,xmlNode,problemManager):
        """Create new named solver from xml node data."""
        
        if not (id in SolverFactory.factories): 
            raise BluecapError("Error: Could not find " + id + " in SolverFactory")
            
        return SolverFactory.factories[id].CreateFromXML(xmlNode,problemManager)




