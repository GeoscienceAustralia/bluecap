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



from string import Template


# makes # the delimiter
class HashTemplate(Template):
    delimiter = '#'
    


class ParameterManager(object):
    """
        Singleton object to contain user defined parameters
        
    """
 
    class __ParameterManager:
      def __init__(self):
            self.params = {}
            
  
      def SetParameter(self,name,paramVal):
            self.params[name] = paramVal
            
      def GetParameter(self,name):
            return self.params[name]
            
      def HasParameter(self,name):
            return (name in self.params)
            
      def ReplaceParams(self,aString):
        s = HashTemplate(aString)
        rv = s.substitute(self.params)
        
        return rv
          
            
      """  
      def WriteXMLNode(self, node): - must be done outside to avoid recursion
      
        # parameters 
        for name,paramVal in self.params.iteritems():
          paramNode = AddChild(node,"parameter")
          SetAttributeString(paramNode,"name",name)
          SetAttributeString(paramNode,"value",paramVal)
      
        return node
      """
    
    instance = None
    
    def __new__(cls): # __new__ always a classmethod
        if not ParameterManager.instance:
            ParameterManager.instance = ParameterManager.__ParameterManager()
        return ParameterManager.instance
        
    def __getattr__(self, name):
        return getattr(self.instance, name)
    def __setattr__(self, name):
        return setattr(self.instance, name)




####################


         