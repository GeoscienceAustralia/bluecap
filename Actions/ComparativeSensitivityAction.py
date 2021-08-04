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
from IO.XML import GetAttributeString,HasAttribute
from IO.XML import SetAttributeString
from IO.XML import GetAttributeVectorString, GetAttributeVector

# Managers
from Managers.ActionManager import ActionManager
from Managers.ActionManager import ActionFactory

from Managers.ParameterManager import ParameterManager


import pylab as pl
import numpy as np

####################

# Compare sensitivity of regional calculation to variations in input parameters.

class ComparativeSensitivityAction():  
    typeStr = "CompareSensitivities"

    def __init__(self):   
       self.type = "NPV"
       self.parameters = ["","",""]
       self.minValues = [0,0,0]
       self.maxValues = [0,0,0]
       self.defaults = [0,0,0]
       self.abs = True
       # param list helps with plotting
       self.paramlist   = ['param1','param2','param3']
       self.sens_keys   = []
       
    def Run(self,problemManager):
        """Run the ComparativeSensitivityAction action."""
        
        from Managers.ProblemManager import ProblemManager   # avoiding circular dependency
        
        problemManager.theRegionalCalculationManager.type = self.type
        #problemManager.theRegionalCalculationManager.Run(problemManager)
        
        diffs = []
        for i in range(3):
            theParameterManager = ParameterManager()
            originalValue = theParameterManager.GetParameter(self.parameters[i])
        
            prevAction = problemManager.theActionManager.activeIndex - 1
        
            # reset params for max value
            theParameterManager.SetParameter(self.parameters[i],self.maxValues[i])
        
            # dummy manager based on run max value - nb may be able to use xml tree in prob manager rather than reloading...
            dummyProblemManager = ProblemManager.FromXMLFile(problemManager.inputFile)
            dummyProblemManager.Initialize()

            # rerun to this point
            dummyProblemManager.RunUntil(prevAction)
        
            # run regional calculation and record max value
            dummyProblemManager.theRegionalCalculationManager.type = self.type
            dummyProblemManager.theRegionalCalculationManager.saveResult = False
            maxValues = dummyProblemManager.theRegionalCalculationManager.Run(dummyProblemManager)
        
            # reset params for min value
            theParameterManager.SetParameter(self.parameters[i],self.minValues[i])
            dummyProblemManager = ProblemManager.FromXMLFile(problemManager.inputFile)

            # reset params
            dummyProblemManager.Initialize()

            # rerun to this point
            dummyProblemManager.RunUntil(prevAction)
            dummyProblemManager.theRegionalCalculationManager.type = self.type
        
            # run regional calculation and record min value
            minValues = dummyProblemManager.theRegionalCalculationManager.Run(dummyProblemManager)
            dummyProblemManager.theRegionalCalculationManager.saveResult = False
        
            diffs.append( maxValues - minValues )
        
            ## fixme need to modify output control
            self.sens_keys.append(self.parameters[i])
            filename = problemManager.outputPrefix+"_snstvty_"+self.parameters[i]+"."+problemManager.outputType
            if problemManager.theRegionalCalculationManager.type == 'NPV':
              data = (maxValues - minValues)*1e-6
            else:
              data = (maxValues - minValues)
            problemManager.theRegionalCalculationManager.type = self.paramlist[self.sens_keys.index(self.parameters[i])]
            problemManager.theRegionalCalculationManager.SaveMap(filename, data, problemManager.recordRange)
            problemManager.theRegionalCalculationManager.type = self.type
        
            theParameterManager.SetParameter(self.parameters[i],originalValue)
        
        diffs = np.array(diffs)
        
        # reshape to image configuration [3,M,N] -> [M,N,3]
        diffs = np.moveaxis(diffs,0,-1)
        
        # pl.imshow(diffs,origin="lower")
        
        # pl.figure()
        # pl.imshow(diffs[:,:,0],origin="lower")
        # pl.colorbar()
        
        # pl.figure()
        # pl.imshow(diffs[:,:,1],origin="lower")
        # pl.colorbar()
        
        # pl.figure()
        # pl.imshow(diffs[:,:,2],origin="lower")
        # pl.colorbar()
        
        # pl.figure()
        # # trying hsv
        
        # np.save("diffs.npy",diffs)
        
        # from matplotlib.colors import hsv_to_rgb
        # x = 0.5*(diffs[:,:,1] - diffs[:,:,2])
        # y =  diffs[:,:,0] - 0.8660254037844386*(diffs[:,:,1] + diffs[:,:,2])
        # h = np.arctan2(y,x)/(2*np.pi)+0.5
        # # print "max h", np.nanmax(h)
        # s = np.sqrt(x**2 +  y**2)
        # v = np.sum(diffs,2)/3.
        # hsv = np.array([h,s,v])
        # hsv = np.moveaxis(hsv,0,-1)
        # hh = hsv_to_rgb(hsv)
        # pl.imshow(hh,origin="lower")
        
        # pl.show()
        
        if self.abs:
          diffs = np.abs(diffs)
        diffScale = np.nanmax(diffs) +1e-64 # np.nanmean(diffs)+2*np.nanstd(diffs)+1e-64  # np.nanmax(diffs)
        #print "diffScale", diffScale
        diffs = diffs/diffScale
        
        filename = problemManager.outputPrefix+"_snstvty_RGB."+problemManager.outputType
        # Process per-pixel diff values for better hue contrast
        diffs = diffs/diffs.sum(2).reshape((diffs.shape[0],diffs.shape[1],1))
        diffs = diffs/diffs.mean(2).reshape((diffs.shape[0],diffs.shape[1],1))
        diffs[diffs>1.] = 1.
        # Add transparency layer
        trans = (np.isnan(diffs[:,:,0]) & np.isnan(diffs[:,:,1]) & np.isnan(diffs[:,:,2])) == False
        diffs = np.concatenate((diffs, trans.reshape((trans.shape[0],trans.shape[1],1))), 2)
        # Adjust calculation manager type for appropriate plotting
        problemManager.theRegionalCalculationManager.type = "sensitivity"
        problemManager.theRegionalCalculationManager.SaveMap(filename, diffs, problemManager.recordRange)
        problemManager.theRegionalCalculationManager.type = self.type
        # Save parameter order
        filename = 'parameters.txt'
        with open(filename, 'w') as F:
            F.write(' '.join(self.parameters))
           
    @classmethod
    def CreateFromXML(cls,actionNode,problemManager):
      """This creates a ComparativeSensitivityAction object from an xml node by initiating an empty class and parsing the xml node."""
      rv = cls()
      rv.ParseXMLNode(actionNode,problemManager)
      return rv
    
    def ParseXMLNode(self,actionNode,problemManager):
       """Generate comparitive sensitivity action from xml tree node."""
       if(HasAttribute(actionNode,"type")):
         self.type = GetAttributeString(actionNode,"type")
         
       self.parameters = GetAttributeVectorString(actionNode,"parameters")
       self.minValues = GetAttributeVector(actionNode,"minima")
       self.maxValues = GetAttributeVector(actionNode,"maxima")
    
    def WriteXMLNode(self, node): 
       """Write the ComparativeSensitivityAction action to an XML tree node."""
       if(self.type):
         SetAttributeString(node,"type",self.type)
       return node
       
    class Factory:
      """The factory that creates the ComparativeSensitivityAction Object."""
      def CreateFromXML(self,xmlNode,problemManager): return ComparativeSensitivityAction.CreateFromXML(xmlNode,problemManager)


# Factory Registrator
ActionFactory.AddFactory(ComparativeSensitivityAction.typeStr,ComparativeSensitivityAction.Factory())
