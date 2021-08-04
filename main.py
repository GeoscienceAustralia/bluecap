#!/usr/bin/env python
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

import numpy as np
import matplotlib.pyplot as pl


from Managers.ProblemManager import ProblemManager
from IO.CommandLine import ParseCommandLineArgs

# Initialisation
rv = ParseCommandLineArgs()

# from xml
theProblemManager = ProblemManager.FromXMLFile(rv.input)

theProblemManager.Initialize()

theProblemManager.Run()

# Record final state in xml - backwards compatibility
if(not theProblemManager.theRegionalCalculationManager.type):

  if(not theProblemManager.outputType):
    finalOutputFile = theProblemManager.outputPrefix + "_final.xml"
    theProblemManager.ExportXMLFile(finalOutputFile)
    theProblemManager.theMineDataManager.PlotResults()
  elif (theProblemManager.outputType == "txt"):
    theProblemManager.theMineDataManager.RecordResults(theProblemManager)
    

exit()