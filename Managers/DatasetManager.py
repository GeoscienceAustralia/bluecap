"""
Copyright (C) 2019-2024, Monash University, Geoscience Australia
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

import os 
import numpy as np

#Common
from Common.Common import BluecapError
from Common.BluecapBaseClass import BluecapClass


# Managers
from .FactoryManager import FactoryManager

# IO
from IO.XML import HasChild,GetChild,AddChild,GetChildren
from IO.XML import HasAttribute,GetAttributeString,GetAttributeStringOrDefault,GetAttributeFileString,GetAttributeFileStringOrDefault
from IO.XML import GetAttributeValue,GetAttributeVector, SetAttributeString,GetXMLTag



class DatasetManager(object):
    """Singleton object to store named datasets between calculations
    
    Usage:
      from Managers.DasetManager import DasetManager
      
    """
    
    class __DatasetManager:
      """ The unit manager is a singleton object."""
    
      def __init__(self):
            """Initialize the singleton object."""
            
            self.datasets = {}
  
      def SetDataset(self,name,dataset):
        """Set the named dataset to a given value."""
        self.datasets[name] =  dataset
        return self
  
      def GetDataset(self,name):
        """Get the named dataset."""
        return self.datasets[name]
    
      def GetDatasetData(self,name):
        """Convert from base units into the stated units."""
        return self.datasets[name].GetData()
  
      def ParseXMLNode(self, dsManagerNode,problemManager):
        """
        Generate Functions from xml tree node. 
        """
        
        theDatasetFactory = FactoryManager().GetFactoryGroup("DatasetFactory")
      
        for child in GetChildren(dsManagerNode):
          type = GetXMLTag(child)
          if(type != "note"):
            name = GetAttributeString(child,"name")
            self.datasets[name] = theDatasetFactory.CreateFromXML(type,child,problemManager)
        
  
      def WriteXMLNode(self, node):
        """
        Write problem to xml node
        """
      
        # functions
        for name,bluecapDataset in self.functions.items():
          type = dataset.typeStr
          dsNode = AddChild(node,type)
          bluecapDataset.WriteXMLNode(dsNode)
          
    instance = None    # this class variable points to the singleton object
    
    def __new__(cls): # __new__ always a classmethod
        if not DatasetManager.instance:
            DatasetManager.instance = DatasetManager.__DatasetManager()
        return DatasetManager.instance
    def __getattr__(self, name):
        return getattr(self.instance, name)
    def __setattr__(self, name):
        return setattr(self.instance, name)
        
        
####################

## Dataset factory

"""Factory class used to generate named datasets."""
DatasetFactory = FactoryManager().AddFactoryGroup("DatasetFactory")

# Dataset from file


class BluecapDataset(BluecapClass):

  def __init__(self):
     self.name = ""
     self.data = None
     self.filename = None
     
  @classmethod
  def GetXMLTypeStr(cls):
    return "Dataset"
  
  def ParseXMLNode(self, xmlnode,problemManager):
     self.filename = GetAttributeString(xmlnode,"file")
     self.data = np.load(self.filename)
     return self
 
  def WriteXMLNode(self, xmlnode):
     SetAttributeString(xmlnode,"file",self.filename)
     SetAttributeString(xmlnode,"name",self.name)
     return self
     
  def GetData(self):
     return self.data

# Factory Registrator
DatasetFactory.AddBluecapClassFactory(BluecapDataset)

      
        