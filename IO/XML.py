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

import xml.etree.ElementTree as ElementTree

from Units.UnitManager import UnitManager
from Managers.ParameterManager import ParameterManager

import os.path



def LoadXMLFile(filename):
  """Load and parse an xml file and return the root node."""
  tree = ElementTree.parse(filename)
   
  root = tree.getroot()
  
  activePath = os.path.abspath(os.path.dirname(filename))
  
  ResolveIncludes(root,activePath)
  
  return root

def NewXMLTree(rootType = "Problem"):
  """Create a new empty xml tree."""
  tree = ElementTree.ElementTree(ElementTree.XML("<"+rootType+"/>"))
  return tree

def GetXMLTreeRoot(tree):
  """Get the root node of an xml tree."""
  root = tree.getroot()
  return root

def SaveXMLFile(filename,tree):  
  """Save the xml tree structure to an xml file."""
  tree.write(filename)
 

def GetXMLTag(node):
  """Return the tag associated with an xml node."""
  return node.tag

def HasAttribute(node,attribute):
  """Return true if the node has a named attribute."""
  return attribute in node.attrib
  
def GetAttributes(node):
  """Return the list of attributes associated with the node. """
  attributeList = list(node.attrib.keys() )
  return attributeList
  
def GetAttributeValue(node,attribute):
  """Get the attribute value."""
  theString = GetAttributeString(node,attribute)
  theUnitManager = UnitManager()
  return theUnitManager.ConvertToBaseUnits(theString)

def GetAttributeFileString(node,attribute):
  """ Gets the attribute file string."""
  theFileString = GetAttributeString(node,attribute)
  if( HasAttribute(node,"ActiveFilePath") ):
    theFileString = os.path.join(GetAttributeString(node,"ActiveFilePath"), theFileString)
  return theFileString

def GetAttributeFileStringOrDefault(node,attribute,default):
  """ Gets the attribute file string or a default value if not defined."""
  if( HasAttribute(node,attribute) ):
    return GetAttributeFileString(node,attribute)
  else:
    return default

def GetAttributeString(node,attribute):
  """ Gets the attribute string."""
  rv = node.attrib[attribute]
  theParameterManager = ParameterManager()
  rv = theParameterManager.ReplaceParams(rv)
  return rv

def GetAttributeStringOrDefault(node,attribute,default):
  """ Gets the attribute string or a default value if not defined."""
  if( HasAttribute(node,attribute) ):
    return GetAttributeString(node,attribute)
  else:
    return default

def SetAttributeString(node,attribute,value):
  """ Sets the attribute string to a particular value."""
  node.attrib[attribute] = str(value)

def GetAttributeVectorString(node,attribute):
  """ Returns a vector of strings associated with the attribute (containing comma separated values)."""
  #print node.tag
  theString = GetAttributeString(node,attribute)
  theUnitManager = UnitManager()
  rv = [x.strip() for x in theString.split(",")]
  return rv

def GetAttributeVector(node,attribute):
  """ Returns vector valued attribute (containing comma separated values)."""
  
  theString = GetAttributeString(node,attribute)
  theUnitManager = UnitManager()
  rv = []
  for x in theString.split(","): # Allow paths to be passed to sensitivity calculations
    if ('.npy' in x) or ('.txt' in x):
      rv.append(x)
    else:
      rv.append(theUnitManager.ConvertToBaseUnits(x))
  # rv = [theUnitManager.ConvertToBaseUnits(x) for x in theString.split(",")]
  return rv

def SetAttributeVector(node,attribute,vectorValue):
  """ Sets node attribute to comma separated vector values."""
  node.attrib[attribute] = ",".join( [str(i) for i in vectorValue] )

def GetAttributeValueOrDefault(node,attribute,default):
  """ Returns the attribute value or a default value if not provided."""
  rv = default
  if(HasAttribute(node,attribute)):
    rv = GetAttributeValue(node,attribute)
  return rv

def GetXMLTag(node):
  """ Get XML tag of current node."""
  return node.tag

def HasChild(node,type):
  """ Return true if node has a child of type 'type'."""
  rv = False
  if( not (node.find(type) is None) ):
    rv = True
  return rv

def GetChildren(node,type=""):
  """ Return all children of the current node."""
  if(type):
    rv = node.findall(type)
  else:
    try: # TODO: getchildren is now depreciated
      rv = node.getchildren()
    except AttributeError:
      rv = node
  return rv

def GetChild(node,type):
  """ Return the first child of 'type' type."""
  rv = node.find(type)  # first child only
  return rv

def AddChild(node,type):
  """ Add a child of 'type' to present node."""
  rv = ElementTree.SubElement(node,type)  
  return rv

def PrettyXMLFormat(node, level=0):
    """ Add whitespace to nodes to make text output more visually appealing."""
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
  """ Traverse the xml tags - load and resolve 'Include' tags. """

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

# Notes
#######

def WriteNotes(fid,node,lvl=0):
    """ Traverse the xml tags - write any notes to file id (fid). """
    for child in GetChildren(node):
      if GetXMLTag(child) == "Note":
        fid.write(child.text)
      else:
        WriteNotes(fid,child,lvl+1)

def AddNote(node,noteString):
  """ Add a note to an xml node. """
  noteNode = AddChild(node,"Note")
  noteNode.text =  noteString
