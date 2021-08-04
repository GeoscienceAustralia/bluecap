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

# Test definitions

import os, sys

import pytest

import numpy as np

# Paths 
###############

bluecapPath =os.environ.get("BLUECAP_PATH",None)

if(bluecapPath is None):
   raise Exception('The environment variable "BLUECAP_PATH" has not been set. ')

testsPath = os.environ.get("BLUECAP_TESTS_PATH", bluecapPath + "/Tests")

baselinePath = os.environ.get("BLUECAP_BASELINE_PATH", bluecapPath + "/../bluecap_test_baseline" )  # baseline tests

# normalize paths
bluecapPath = os.path.abspath(bluecapPath)
testsPath = os.path.abspath(testsPath)
baselinePath = os.path.abspath(baselinePath)

# Files
###############

xmlSchema = os.environ.get("BLUECAP_XML_SCHEMA",bluecapPath + "/IO/bluecap_schema.xsd")

# Executables and args
########################

bluecapExecutable = os.environ.get("BLUECAP_EXECUTABLE",bluecapPath +"/main.py")

xmlValidatorExecutable = os.environ.get("XML_VALIDATOR","/usr/bin/xmllint")
xmlCommandLineArgs = os.environ.get("XML_VALIDATOR_ARGS","--noout --schema " + xmlSchema)

diffExecutable = os.environ.get("DIFF_EXECUTABLE","/usr/bin/diff")
diffCommandLineArgs = os.environ.get("DIFF_EXECUTABLE_ARGS","-u") # output diffs in unified mode

makeExecutable = os.environ.get("MAKE_EXECUTABLE","/usr/bin/make")

# Test functions
#################

def DoFileDiff(fileA,fileB):
  """Compare two files using the diff executable.
  
  Returns:
  	0 if successfull, diff command line exit code otherwise. 
  """
  
  diffString = diffExecutable + " " + diffCommandLineArgs + " " + fileA + " " +  fileB 
  rv = os.system(diffString)
  return rv
  
def DoNpyFileDiff(fileA,fileB):
  """Compare two numpy files using the numpy assert equal executable.
  
  Returns:
  	0 if successfull, 1 otherwise (keeping with command line error code convention). 
  """
  dataA = np.load(fileA)
  dataB = np.load(fileB)
  
  #areEqual = np.array_equal(dataA,dataB)  # doesn't handle nans
  rv = 0
  try:
    np.testing.assert_equal(dataA,dataB)
  except AssertionError:
    np.save(fileA+".diff.npy",dataA-dataB)
    rv = 1
    
  return rv


def DoXmlSchemaCheck(xmlFile):
  """Checks that the XML file complies with the bluecap XML schema.
  
  Returns:
  	0 if successfull, xml validator exit code otherwise. 
  """
  diffString = xmlValidatorExecutable + " " + xmlCommandLineArgs + " " + xmlFile 
  rv = os.system(diffString)
  return rv

def RunBluecap(filename="",commandLineArgs=""):
  """Runs bluecap using the input file provided.
  
  Returns:
  	0 if successfull, xml validator exit code otherwise. 
  """
  if(filename):
    runString = bluecapExecutable + " " + commandLineArgs + " -i " + filename 
  else:
    runString = bluecapExecutable + " " + commandLineArgs
  
  rv = os.system(runString)
  return rv


def GetCurrentTestPath():
  """Returns path of folder containing the current test."""

  currentTestPath = os.path.dirname(os.getenv('PYTEST_CURRENT_TEST').split(":")[0])
  
  #currentTestPath = testsPath + "/" + 
  currentTestPath =  currentTestPath + "/"
  return currentTestPath

def GetCurrentBaselinePath():
  """Returns path of folder containing baseline file corresponding to the current test."""

  currentTestPath = os.path.dirname(os.getenv('PYTEST_CURRENT_TEST').split(":")[0])
  
  currentBaselinePath = baselinePath + "/" +  currentTestPath + "/"
  return currentBaselinePath


# Test Clases
#############


  

class XmlSchema_Base():
  """Base class used to create XML schema tests."""
  
  # We can't use init to store these variables and still use pytest
  # which leads to some ugly code below. 
  _xmlFile = ""
  
  def test_XmlSchema_Base(self):
    fullXMLFile = GetCurrentTestPath() + self.__class__._xmlFile
    assert DoXmlSchemaCheck(fullXMLFile) == 0


class XmlRunDiff_Base(XmlSchema_Base):
  """Base class used to create tests that check the file against the XML schema, 
  run bluecap on the file and then compare the result to the output (text) file 
  using diff.
  
  Usage:

  class Test_ExampleXmlRunDiff(XmlRunDiff_Base):
  
    _xmlFile = "one_input.xml"
    _outputFile = "one.xml"
  """
  
  # We can't use init to store these variables and still use pytest
  # which leads to some ugly code below. 
  _xmlFile = ""
  _outputFile = ""
  _commandLineArgs=""
  
  
  def test_RunAndDiff(self):
    fullXMLFile = GetCurrentTestPath() + self.__class__._xmlFile
    rv = RunBluecap(fullXMLFile,self.__class__._commandLineArgs)
    if(rv is not 0):
      print("Error: " + self.__class__._xmlFile + " did not run.")
      assert 0
    
    fullOutputFile = GetCurrentTestPath() + self.__class__._outputFile
    baselineOutputFile = GetCurrentBaselinePath()  + self.__class__._outputFile
    
    assert DoFileDiff(fullOutputFile,baselineOutputFile) == 0



class XmlRunNpyDiff_Base(XmlSchema_Base):
  """Base class used to create tests that check the file against the XML schema, 
  run bluecap on the file and then compare the result to an output numpy (npy) 
  file.
  
  Usage:

  class Test_ExampleXmlRunNpyDiff(XmlRunNpyDiff_Base):
  
    _xmlFile = "one_input.xml"
    _outputFile = "one.npy"
  """
  
  # We can't use init to store these variables and still use pytest
  # which leads to some ugly code below. 
  _xmlFile = ""
  _outputFile = ""
  _commandLineArgs=""
  
  
  def test_RunAndNpyDiff(self):
    fullXMLFile = GetCurrentTestPath() + self.__class__._xmlFile
    rv = RunBluecap(fullXMLFile,self.__class__._commandLineArgs)
    if(rv is not 0):
      print("Error: " + self.__class__._xmlFile + " did not run.")
      assert 0
    
    
    fullOutputFile = GetCurrentTestPath() + self.__class__._outputFile
    baselineOutputFile = GetCurrentBaselinePath() + self.__class__._outputFile
    
    assert DoNpyFileDiff(fullOutputFile,baselineOutputFile) == 0




