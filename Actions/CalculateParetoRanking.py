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

# IO
from IO.XML import GetAttributes, GetAttributeValueOrDefault
from IO.XML import GetAttributeString,GetAttributeStringOrDefault,HasAttribute,SetAttributeString

# Managers
from Managers.ActionManager import ActionManager
from Managers.ActionManager import ActionFactory



from Managers.DatasetManager import DatasetManager

from Functions.ParetoRank import FindParetoRanks

import numpy as np
import pylab as pl

####################

class CalculateParetoRanking():  
    typeStr = "CalculateParetoRanking"

    def __init__(self):   
       self.datasetNames = []
       self.types = []
       self.outputFile = "pareto.npy"
       self.stride = 1
       
    def Run(self,problemManager):
       """
       Calculate the pareto ranking of the named datasets
       """
       theDatasetManager = DatasetManager()
       
       pts = []
       cc = 0
       for dsString in self.datasetNames:
          print("loading", dsString)
          if(dsString[-4:] ==".npy"):
            ds =  np.load(dsString)
          elif(dsString[-4:] ==".tif"):
            ds = np.array(pl.imread(dsString))
            if(ds.ndim == 3):
              ds = ds[:,:,0]
            #dataA = np.array(Image.open(dsString) )
          else:
            ds = theDatasetManager.GetDatasetData(dsString)
          
          if(self.stride > 1):
            ds = ds[::self.stride,::self.stride]
          
          if(self.types[cc] == "negative"):
             ds = -1*ds
          
          
          originalShape = ds.shape
          pts.append(ds.ravel())
          
          cc += 1
          
       # remove nans 
       isNan = False   
       for ds in pts:
         isNan = np.logical_or( np.isnan( ds ),isNan )
         
         
       isNotNan = np.logical_not(isNan)
       
       pts = [ds[isNotNan] for ds in pts]
       pts = [ (ds - np.min(ds))/(np.max(ds)-np.min(ds)) for ds in pts]  # normalize points
       pts = np.column_stack(pts)
       
       ranksMap = np.zeros(originalShape)
       
       ptsUnique, indices = np.unique(pts, axis=0, return_inverse=True)
       
       print("running pareto ranks")
       #ranks,parents = FindParetoRanks(pts)
       ranksUnique,parents = FindParetoRanks(ptsUnique)
       ranks = np.zeros(pts.shape[0])
       ranks = ranksUnique[indices]
       
       ranksMap[:] =  np.max(ranks)+1
       ranksMap.ravel()[isNotNan] = ranks
       
       print("finished pareto ranks")
       
       if(self.outputFile[-4:] ==".npy"):
         np.save(self.outputFile,ranksMap)
       elif(self.outputFile[-4:] ==".txt"):
         np.savetxt(self.outputFile,ranksMap)
       elif(self.outputFile[-4:] ==".png"):
         pl.imsave(self.outputFile,ranksMap)
         
       
       
       pl.imshow(ranksMap,cmap="inferno_r")
       pl.colorbar()
       pl.show()
       
           
    @classmethod
    def CreateFromXML(cls,actionNode,problemManager):
      """ This creates the class from an xml node by initiating an empty class and parsing the xml node."""
      rv = cls()
      rv.ParseXMLNode(actionNode,problemManager)
      return rv
    
    def ParseXMLNode(self,actionNode,problemManager):
       """
       Generate the solver call from the xml input
       """
       datasetNames = GetAttributeString(actionNode,"datasets")
       
       self.datasetNames =  [ name.strip() for name in datasetNames.split(",") ]
       typeNames = GetAttributeString(actionNode,"types")
       self.types =  [ name.strip() for name in typeNames.split(",") ]
       self.outputFile = GetAttributeStringOrDefault(actionNode,"outputFile",self.outputFile)
       self.stride = int( GetAttributeValueOrDefault(actionNode,"stride",self.stride) )
       
       return
    
    def WriteXMLNode(self, node): 
       """
       Write the call to the xml tree
       """
       datasetsString = ",".join(self.datasetNames)
       SetAttributeString(node,"datasets",datasetsString)
       return node
       
    
    class Factory:
      """The factory that creates the Calculate Pareto Ranking Object"""
      def CreateFromXML(self,xmlNode,problemManager): 
         return CalculateParetoRanking.CreateFromXML(xmlNode,problemManager)


# Factory Registrator
ActionFactory.AddFactory(CalculateParetoRanking.typeStr,CalculateParetoRanking.Factory())
