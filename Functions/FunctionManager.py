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
from IO.XML import HasChild,GetChild,GetChildren,AddChild,GetXMLTag
from IO.XML import HasAttribute,GetAttributeString,SetAttributeString
from IO.XML import GetAttributeValue, GetAttributeValueOrDefault
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
            
      
  
      def SetFunction(self,name,afunc):
            self.functions[name] = afunc
            
      def GetFunction(self,name):
            return self.functions[name]
            
      def HasFunction(self,name):
            return (name in self.functions)
            
      def PerturbRandomFunctions(self):
            for name,usrFunc in self.functions.items():
              usrFunc.UpdateRandomState()

      def ZeroPerturbations(self):
            for name,usrFunc in self.functions.items():
              usrFunc.ZeroPerturbation()
            
      def ParseXMLNode(self, fnManagerNode,problemManager):
        """
        Generate Functions from xml tree node. 
        """
      
        for child in GetChildren(fnManagerNode):
          type = GetXMLTag(child)
          if(type != "note"):
            name = GetAttributeString(child,"name")
          
            self.functions[name] = FunctionFactory.CreateFromXML(type,child,problemManager)
        
      def WriteXMLNode(self, node):
        """
        Write problem to xml node
        """
      
        # functions
        for name,usrFunc in self.functions.items():
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
        if not (id in FunctionFactory.factories): #.has_key(id):
            print("Error!!! could not find " + id + " in FunctionFactory")
            
        return FunctionFactory.factories[id].CreateFromXML(xmlNode,problemManager)
    #Create = staticmethod(CreateFromXML)


####################




class UserDefinedFunction():  
    def __init__(self):   
       self.name = ""
       self.doRandomSample = False
       
    def f(self,args):
        return 0.0
    #   return  self.func(args) # maybe throw an error instead
       
    def GetName(self):
       return self.name
       
    def UpdateRandomState(self):
       # empty
       return self

    def ZeroPerturbation(self):
       # empty
       return self
           
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
       
    # The factory that creates the UserDefinedFunction Object - attempt at making more general
    # @classmethod
    #def GetFactory(cls):
    #  
    #  class FactoryClass:
    #    def CreateFromXML(self,xmlNode,problemManager): return cls.CreateFromXML(xmlNode,problemManager)
    #    
    #  return FactoryClass()


 
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
    """
    Defines a powerlaw function of the form f(x) = A x^b
    """


    typeStr = "PowerLaw"

    def __init__(self):  
       UserDefinedFunction.__init__(self)
       self.coeff = 0.0
       self.power = 0.0
       self.unperturbed_coeff = 0.0
       self.stdErrorScale = None
    
    def ParseXMLNode(self, functionObjectNode,problemManager):
       """
       Generate a powerlaw function from xml tree node. 
       """
       UserDefinedFunction.ParseXMLNode(self,functionObjectNode,problemManager)
       self.coeff = GetAttributeValue(functionObjectNode,"coeff")
       self.unperturbed_coeff =  np.copy(self.coeff)
       self.power = GetAttributeValue(functionObjectNode,"power")
       self.stdErrorScale = GetAttributeValueOrDefault(functionObjectNode,"stdErrorScale",self.stdErrorScale)
       
    def f(self,args):
       """Return the value of the power-law function."""
       rv = self.coeff*args[0]**self.power 
       return rv

    def UpdateRandomState(self):
      self.coeff = np.copy(self.unperturbed_coeff)
      if(self.stdErrorScale):
        self.coeff *= self.stdErrorScale** np.random.normal()
      return self

    def ZeroPerturbedState(self):
      self.coeff = np.copy(self.unperturbed_coeff)
      return self
    
    def GetStdErrorScale(self):
       return self.stdErrorScale
    
    def WriteXMLNode(self, node): 
       """Write the power law function to the xml node."""
       SetAttributeString(node,"name",self.name)
       SetAttributeString(node,"coeff",self.unperturbed_coeff)
       SetAttributeString(node,"power",self.power)
       if(self.stdErrorScale is not None):
          SetAttributeString(node,"stdErrorScale",self.stdErrorScale)
       return node

     # The factory that creates the PowerlawFunction Object 
    class Factory:
      def CreateFromXML(self,xmlNode,problemManager): return PowerlawFunction.CreateFromXML(xmlNode,problemManager)


# Factory Registrator
FunctionFactory.AddFactory(PowerlawFunction.typeStr,PowerlawFunction.Factory())


class WrapperFunction(UserDefinedFunction):
    """
    Function used to change the name of another function. It can be used with parameters to change active function call. 
    """
    
    typeStr = "WrapperFunction"

    def __init__(self):  
       UserDefinedFunction.__init__(self)
       self.wrappedFunctionName = ""
       self.wrappedFunction = None
    
    def ParseXMLNode(self, functionObjectNode,problemManager):
       """
       Generate a wrapper function from an xml tree node. 
       """
       UserDefinedFunction.ParseXMLNode(self,functionObjectNode,problemManager)
       self.wrappedFunctionName = GetAttributeString(functionObjectNode,"function")
       theFunctionManager = FunctionManager()
       if(theFunctionManager.HasFunction(self.wrappedFunctionName)):
         self.wrappedFunction = theFunctionManager.GetFunction(self.wrappedFunctionName)
       else:
         print("Error - the function '"+ self.wrappedFunctionName + "' must be defined before the wrapper function is called.")
         exit(1)
       
    def f(self,args):
       """Return the value of the wrapped function."""
       rv = self.wrappedFunction.f(args)
       return rv

    def UpdateRandomState(self):
      """Updates the random state of the wrapped function."""
      self.wrappedFunction.UpdateRandomState()  # need to be a little careful here
      return self

    def ZeroPerturbedState(self):
      """Zeros the perturbed state of the wrapped function."""
      self.wrappedFunction.ZeroPerturbedState() 
      return self
    
    def GetStdErrorScale(self):
       """Returns the standard error scale of the wrapped function."""
       self.wrappedFunction.GetStdErrorScale() 
       return self.wrappedFunction.stdErrorScale
    
    def WriteXMLNode(self, node): 
       """Writes the wrapper function (but not the underlying function) to the xml node."""
       SetAttributeString(node,"name",self.name)
       SetAttributeString(node,"function",self.wrappedFunctionName)
       return node

     # The factory that creates the WrapperFunction Object 
    class Factory:
      def CreateFromXML(self,xmlNode,problemManager): return WrapperFunction.CreateFromXML(xmlNode,problemManager)


# Factory Registrator
FunctionFactory.AddFactory(WrapperFunction.typeStr,WrapperFunction.Factory())
    
    

         
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
       Generate the two dimensional table function from an xml tree node. 
       """
       UserDefinedFunction.ParseXMLNode(self,functionObjectNode,problemManager)
       self.xfile = GetAttributeFileString(functionObjectNode,"xticks")
       self.yfile = GetAttributeFileString(functionObjectNode,"yticks")
       self.valueFile = GetAttributeFileString(functionObjectNode,"values")
       if(self.xfile[-4:] == ".npy"):
         self.xs = np.load(self.xfile)
       else:
         #print("loading:", self.xfile)
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
       """Lookup the value of the 2D table - NB requires two arguments"""
       rv = self.interpFunc(args[0],args[1])[0]
       return rv
    
    def WriteXMLNode(self, node): 
       """Write the 2D table function to the xml node"""
       SetAttributeString(node,"xticks",self.xfile)
       SetAttributeString(node,"yticks",self.yfile)
       SetAttributeString(node,"values",self.valueFile)
       return node

     # The factory that creates the Table2DFunction Object 
    class Factory:
      def CreateFromXML(self,xmlNode,problemManager): return Table2DFunction.CreateFromXML(xmlNode,problemManager)


# Factory Registrator
FunctionFactory.AddFactory(Table2DFunction.typeStr,Table2DFunction.Factory())

         