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

import xml.etree.ElementTree as ElementTree

from Units.UnitManager import UnitManager
from Managers.ParameterManager import ParameterManager

import os.path


## Load and parse an xml file and return the root node
def LoadXMLFile(filename):
  tree = ElementTree.parse(filename)
  ## fixme - we may want to do additional parsing here to add included files etc. 
  root = tree.getroot()
  
  activePath = os.path.abspath(os.path.dirname(filename))
  
  ResolveIncludes(root,activePath)
  
  return root
  

def NewXMLTree(rootType = "Problem"):
  tree = ElementTree.ElementTree(ElementTree.XML("<"+rootType+"/>"))
  return tree

def GetXMLTreeRoot(tree):
  root = tree.getroot()
  return root

def SaveXMLFile(filename,tree):  
  tree.write(filename)
 

def GetXMLTag(node):
  return node.tag

def HasAttribute(node,attribute):
  return attribute in node.attrib
  
def GetAttributeValue(node,attribute):
  theString = GetAttributeString(node,attribute)
  theUnitManager = UnitManager()
  return theUnitManager.ConvertToBaseUnits(theString)
  

def GetAttributeFileString(node,attribute):
  theFileString = GetAttributeString(node,attribute)
  if( HasAttribute(node,"ActiveFilePath") ):
    theFileString = os.path.join(GetAttributeString(node,"ActiveFilePath"), theFileString)
  return theFileString

def GetAttributeString(node,attribute):
  rv = node.attrib[attribute]
  theParameterManager = ParameterManager()
  rv = theParameterManager.ReplaceParams(rv)
  return rv

def GetAttributeStringOrDefault(node,attribute,default):
  if( HasAttribute(node,attribute) ):
    return GetAttributeString(node,attribute)
  else:
    return default

def SetAttributeString(node,attribute,value):
  node.attrib[attribute] = str(value)


def GetAttributeVectorString(node,attribute):
  #print node.tag
  theString = GetAttributeString(node,attribute)
  theUnitManager = UnitManager()
  rv = [x.strip() for x in theString.split(",")]
  return rv


def GetAttributeVector(node,attribute):
  #print node.tag
  theString = GetAttributeString(node,attribute)
  theUnitManager = UnitManager()
  rv = [theUnitManager.ConvertToBaseUnits(x) for x in theString.split(",")]
  return rv

def SetAttributeVector(node,attribute,vectorValue):
  #print node.tag
  node.attrib[attribute] = ",".join( [str(i) for i in vectorValue] )

def GetAttributeValueOrDefault(node,attribute,default):
  rv = default # fixme convert string to default value?
  if(HasAttribute(node,attribute)):
    rv = GetAttributeValue(node,attribute)
  return rv

def GetXMLTag(node):
  return node.tag
  
def HasChild(node,type):
  rv = False
  if( not (node.find(type) is None) ):
    rv = True
  return rv
  
def GetChildren(node,type=""):
  if(type):
    rv = node.findall(type)
  else:
    rv= node.getchildren() 
  return rv
  
def GetChild(node,type):
  rv = node.find(type)  # first child only
  return rv
  
def AddChild(node,type):
  rv = ElementTree.SubElement(node,type)  
  return rv
  
# in-place prettyprint formatter

def PrettyXMLFormat(node, level=0):
    theStr = "\n" + level*"    "
    if len(node):
        if not node.text or not node.text.strip():
            node.text = theStr + "    "
        if not node.tail or not node.tail.strip():
            node.tail = theStr
        for i,child in enumerate(node):
            PrettyXMLFormat(child, level+1)
        if not child.tail or not child.tail.strip():
            child.tail = theStr
    else:
        if level and (not node.tail or not node.tail.strip()):
            node.tail = theStr
            


           
def ResolveIncludes(node,activePath):

  SetAttributeString(node,"ActiveFilePath",activePath)

  for child in GetChildren(node):
    if GetXMLTag(child) == "Include":
      filename = GetAttributeString(child,"file")
      
      
      filename = os.path.join(activePath, filename)
      newActivePath = os.path.abspath(os.path.dirname(filename))
      includeRoot =  LoadXMLFile(filename)
      
      for includedChild in GetChildren(includeRoot):
        ResolveIncludes(includedChild,newActivePath)
        node.append(includedChild)
      #node.remove(child)
    else:
      ResolveIncludes(child,activePath)
      
  # clean up - bit ugly as removing a child node seems to break for loop
  hasIncludes = True
  while(hasIncludes):  
    hasIncludes = False   
    for child in GetChildren(node):
      if GetXMLTag(child) == "Include":
        hasIncludes = True
        node.remove(child)
        break





