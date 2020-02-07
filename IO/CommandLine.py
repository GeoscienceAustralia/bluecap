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

import argparse
import string

from Managers.ParameterManager import ParameterManager

def ParseCommandLineArgs():


  theParameterManager = ParameterManager()  
      
  # set input options
  parser = argparse.ArgumentParser(description='Command line input.')
  parser.add_argument("-i","--input",action="store", help = "Load XML input file")
  parser.add_argument("-p","--param",action="append", help = "Add user-defined parameter")
  
  rv = parser.parse_args()
  
  # strip out params
  if (rv.param):
	  for param in rv.param:
		splt = string.split(param,"=",1)
		if (len(splt) < 2 ): 
		  print "error - input params should be in the form -p key=value "
		  exit()
		name = splt[0].strip()
		paramString = splt[1]
		theParameterManager.SetParameter(name,paramString)
  
  
  # sanity check
  if (not rv.input):
    print ("Warning - failed to provide input file")
  
  
  return rv