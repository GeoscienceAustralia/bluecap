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
from IO.XML import HasChild,GetChild,GetChildren,AddChild,GetXMLTag
from IO.XML import HasAttribute,GetAttributeString,SetAttributeString,GetAttributeValue
from IO.XML import GetAttributeFileString

import numpy as np
from scipy import interpolate


class FunctionManager(object):
    """
        Singleton object to contain user defined functions
        Usage:
    from Functions.FunctionManager import FunctionManager
    theFunctionManager = FunctionManager()
    taylorsRule = theFunctionManager.GetFunction("TaylorsRule_UG")
    value = taylorsRule.f(5000)
    """
 
    class __FunctionManager:
      def __init__(self):
            self.functions = {}
            # fixme load functions here
            
  
      def SetFunction(self,name,afunc):
            self.functions[name] = afunc
            
      def GetFunction(self,name):
            return self.functions[name]
            
      def HasFunction(self,name):
            return (name in self.functions)
            
      def ParseXMLNode(self, fnManagerNode,problemManager):
        """
        Generate Functions  from xml tree node. 
        """
      
        for child in GetChildren(fnManagerNode):
          type = GetXMLTag(child)
          name = GetAttributeString(child,"name")
          
          self.functions[name] = FunctionFactory.CreateFromXML(type,child,problemManager)
        
      def WriteXMLNode(self, node):
        """
        Write problem to xml node
        """
      
        # functions
        for name,usrFunc in self.functions.iteritems():
          type = usrFunc.typeStr
          funcNode = AddChild(node,type)
          usrFunc.WriteXMLNode(funcNode)
      
        return node
          
    instance = None
    
    def __new__(cls): # __new__ always a classmethod
        if not FunctionManager.instance:
            FunctionManager.instance = FunctionManager.__FunctionManager()
        return FunctionManager.instance
        
    def __getattr__(self, name):
        return getattr(self.instance, name)
    def __setattr__(self, name):
        return setattr(self.instance, name)

####################

## Function factory

class FunctionFactory:
    factories = {}
    
    @staticmethod
    def AddFactory(id, factory):
        """Add factory to the function factory manager"""
        FunctionFactory.factories[id] = factory
    #AddFactory = staticmethod(AddFactory)
    
    @staticmethod
    def CreateFromXML(id,xmlNode,problemManager):
        if not FunctionFactory.factories.has_key(id):
            print "Error!!! could not find " + id + " in FunctionFactory"
            
        return FunctionFactory.factories[id].CreateFromXML(xmlNode,problemManager)
    #Create = staticmethod(CreateFromXML)


####################




class UserDefinedFunction():  
    def __init__(self):   
       self.name = ""
       
    def f(self,args):
        return 0.0
    #   return  self.func(args) # maybe throw an error instead
       
    def GetName(self):
       return self.name
           
    @classmethod
    def CreateFromXML(cls,funcNode,problemManager):
      # this creates the class from an xml node by initiating an empty class and parsing the xml node
      # note that the class is the derived class - so functions inheriting from UserDefinedFunction only have to define their
      # own empty constructors and ParseXMLNode functions. 
      rv = cls()
      rv.ParseXMLNode(funcNode,problemManager)
      return rv
    
    def ParseXMLNode(self, functionObjectNode,problemManager):
       """
       Generate named function from xml tree node. 
       """
       self.name = GetAttributeString(functionObjectNode,"name")


 
class ConstantFunction(UserDefinedFunction):  

    typeStr = "Constant"

    def __init__(self):  
       UserDefinedFunction.__init__(self)
       self.fval = 0.0
    
    def ParseXMLNode(self, funcNode,problemManager):
       """
       Populate ConstantFunction object with data from xml tree node. 
       """
       UserDefinedFunction.ParseXMLNode(self,funcNode,problemManager)
       self.fval = GetAttributeValue(funcNode,"f")
     
     
    
    def WriteXMLNode(self, node): 
       SetAttributeString(node,"name",self.name)
       SetAttributeString(node,"f",self.fval)
       return node
          
    def f(self,args):
       """
       ConstantFunction function call
       """
       return self.fval

     # The factory that creates the ConstantFunction Object 
    class Factory:
      def CreateFromXML(self,xmlNode,problemManager): return ConstantFunction.CreateFromXML(xmlNode,problemManager)


# Factory Registrator
FunctionFactory.AddFactory(ConstantFunction.typeStr,ConstantFunction.Factory())

       

 
class PowerlawFunction(UserDefinedFunction):  


    typeStr = "PowerLaw"

    def __init__(self):  
       UserDefinedFunction.__init__(self)
       self.coeff = 0.0
       self.power = 0.0
    
    def ParseXMLNode(self, functionObjectNode,problemManager):
       """
       Generate named function from xml tree node. 
       """
       UserDefinedFunction.ParseXMLNode(self,functionObjectNode,problemManager)
       self.coeff = GetAttributeValue(functionObjectNode,"coeff")
       self.power = GetAttributeValue(functionObjectNode,"power")
       
    def f(self,args):
       rv = self.coeff*args[0]**self.power  
       return rv
    
    def WriteXMLNode(self, node): 
       SetAttributeString(node,"name",self.name)
       SetAttributeString(node,"coeff",self.coeff)
       SetAttributeString(node,"power",self.power)
       return node

     # The factory that creates the ConstantFunction Object 
    class Factory:
      def CreateFromXML(self,xmlNode,problemManager): return PowerlawFunction.CreateFromXML(xmlNode,problemManager)


# Factory Registrator
FunctionFactory.AddFactory(PowerlawFunction.typeStr,PowerlawFunction.Factory())




         
class Table2DFunction(UserDefinedFunction):  


    typeStr = "Table2D"

    def __init__(self):  
       UserDefinedFunction.__init__(self)
       self.xs = np.array([0])
       self.ys = np.array([0])
       self.values = np.array([0])
       self.xfile = ""
       self.yfile = ""
       self.valueFile = ""
       self.interpFunc = None
       self.transpose = False
    
    def ParseXMLNode(self, functionObjectNode,problemManager):
       """
       Generate named function from xml tree node. 
       """
       UserDefinedFunction.ParseXMLNode(self,functionObjectNode,problemManager)
       self.xfile = GetAttributeFileString(functionObjectNode,"xticks")
       self.yfile = GetAttributeFileString(functionObjectNode,"yticks")
       self.valueFile = GetAttributeFileString(functionObjectNode,"values")
       if(self.xfile[-4:] == ".npy"):
         self.xs = np.load(self.xfile)
       else:
         print "loading:", self.xfile
         self.xs = np.loadtxt(self.xfile)
         
       if(self.yfile[-4:] == ".npy"):
         self.ys = np.load(self.yfile)
       else:
         self.ys = np.loadtxt(self.yfile)
       
       if(self.valueFile[-4:] == ".npy"):
         self.values = np.load(self.valueFile)
       else:
         self.values = np.loadtxt(self.valueFile)
       
       if( HasAttribute(functionObjectNode,"transpose") ):
         self.transpose = GetAttributeValue(functionObjectNode,"transpose")
       
       if(self.transpose):
         self.values =  self.values.T
        
       self.values[np.isnan(self.values)]  = 2e3
       self.interpFunc = interpolate.interp2d(self.xs, self.ys, self.values, kind='linear')


    def f(self,args):
       rv = self.interpFunc(args[0],args[1])[0]
       return rv
    
    def WriteXMLNode(self, node): 
       SetAttributeString(node,"xticks",self.xfile)
       SetAttributeString(node,"yticks",self.yfile)
       SetAttributeString(node,"values",self.valueFile)
       return node

     # The factory that creates the ConstantFunction Object 
    class Factory:
      def CreateFromXML(self,xmlNode,problemManager): return Table2DFunction.CreateFromXML(xmlNode,problemManager)


# Factory Registrator
FunctionFactory.AddFactory(Table2DFunction.typeStr,Table2DFunction.Factory())

         