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

class BluecapClass():
  # base class for bluecap classes - primarily used to define common factory methods
  
  @classmethod
  def GetXMLTypeStr(cls):
    # by default we set the xml string equal to the class name
    # however this can be redefined by the child class. 
    return cls.__name__
  
  @classmethod
  def CreateFromXML(cls,xmlNode,*args,**kwargs):
      # this creates the class from an xml node by initiating an empty class and parsing the xml node
      rv = cls()
      rv.ParseXMLNode(xmlNode,*args,**kwargs)
    
      return rv

  
  # a factory class used to generate instances of the class
  class InnerFactory(object):
      """The factory for creating Bluecap Instances from different sources
      
      Don't call this directly, instead use the Factory() function. 
      Eg. 
      aNewChildFactory = ChildClass.Factory()
      
      If necessary derived classes can override the Factory function with their own generators, however it is simpler to 
      define a class specific "ParseXMLNode" Function or "CreateFromXML" class method.
      """
      
      def __init__(self,cls):  
        self.FactoryOwnerClass = cls 
      
      
      def CreateFromXML(self,xmlNode,*args,**kwargs): 
         # doing this calls the parent function
         #return ProcessingNodeBase.CreateFromXML(xmlNode,problemManager)
         return self.FactoryOwnerClass.CreateFromXML(xmlNode,*args,**kwargs)

        
  @classmethod
  def Factory(cls):
    """Function used to call the inner factory class """
    rv = cls.InnerFactory(cls)
    return rv
    
