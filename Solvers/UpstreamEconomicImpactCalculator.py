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


# Managers
from Managers.SolverManager import SolverManager
from Managers.SolverManager import SolverFactory


# IO
from IO.XML import HasChild,GetChild,AddChild
from IO.XML import GetAttributeString,GetAttributeValue,GetAttributeValueOrDefault, SetAttributeString, GetAttributeStringOrDefault
from IO.XML import GetAttributeFileString,GetAttributeFileStringOrDefault

# Units
from Units.UnitManager import UnitManager



class UpstreamImpactCalculator():
    typeStr = "UpstreamImpactCalculator"

    def __init__(self):
      """
      Create Empty Upstream Impact Calculator. 
      """
      
      self.ioMatrixPrefix = []
      self.ioMatrixRows = [] 
      #self.ioMatrixColumns = []  # rows and columns are assumed to correspond
      self.ioMatrix = None
      
      
      self.ioMatrixExcludedRows = []   # rows to exclude when calculating input contribution
      self._ioMatrixReducedRows = [] # rows included when calculating input contribution
      
      #self.ioMatrixRegions = []  # indicate regions with "_N" where N is the region number
      
      self.industryColStringPrefix = "Non-ferrous metal ore mining"  # name of column/prefix of column (row) used to indicate the industry in question
      
      self.totalRequirementsMatrix = None
      
      self.secondaryMatrixPrefix = []
      self.secondaryMatrixRows = [] 
      self.secondaryMatrix = None
      
      self.stateBasedIOMatrix = False # Use state based regional io matrix
      self.varyOutputByState = True 
        # Change output industry depending on state (state based calculations only)
        # i.e. if true then output reflects contribution to industry in each state eg contribution to NT farming in NT, QLD farming in QLD
        #      if false then output is the same i.e. reflects contribution to NT farming from all states. 
      

      
    

      
    def Run(self,theProblemManager,args):
    
      # calculate NPC for the project 
      # We will determine the project benefits from the costs rather than from the profits
      # Costs will assume to be distributed over the following quantities for the iomatrix
      #
      
      netCost = -theProblemManager.theRegionalCalculationManager.CalculateNetRealCosts(theProblemManager)
      
      output = args["output"]
      
      ### Io matrix calculation
      
      # determine excluded indexes
      excludedIndxs = self.GetExcludedIndexes()
      
      
      # calculate total requirements matrix
      totalRequirementsMatrix = self.CalculateTotalRequirementsMatrix(excludedIndxs)
      
      if(self.stateBasedIOMatrix):
      
        outputMap = np.zeros_like(netCost)
        
        stateStrs = theProblemManager.theRegionalCalculationManager.stateStrs
        stateIds = theProblemManager.theRegionalCalculationManager.stateIds
        stateIdsMap =  theProblemManager.theRegionalCalculationManager.LoadMap(theProblemManager.theRegionalCalculationManager.stateIdsMapFile)
      
        for i in range(len(stateStrs)):
          
          # calculate demand
          demandVector = self.CalculateDemandVector(self.ioMatrix,excludedIndxs,"_"+stateStrs[i])
          
          # normalize demand vector (to reflect that we are considering costs not outputs)
          demandVector /= np.sum(demandVector)
          
          # calculate totalInputsVector =  totalRequirements * demand
          totalInputsVector = np.dot(totalRequirementsMatrix,demandVector)
          
          # determine scaling factors for output
          if(self.varyOutputByState):
            outputScale = self.CalculateOutputScale(output+"_"+stateStrs[i],totalInputsVector,excludedIndxs)
          else:
            outputScale = self.CalculateOutputScale(output,totalInputsVector,excludedIndxs)
          
          # scale map
          indx = (stateIdsMap == stateIds[i])
          outputMap[ indx ] = outputScale*netCost[indx]
          
      else:
      
          # calculate demand
          demandVector = self.CalculateDemandVector(self.ioMatrix,excludedIndxs)
      
          # normalize demand vector (to reflect that we are considering costs not outputs)
          demandVector /= np.sum(demandVector)
      
          # calculate totalInputsVector =  totalRequirements * demand
          totalInputsVector = np.dot(totalRequirementsMatrix,demandVector)
      
          # determine scaling factors for output
          outputScale = self.CalculateOutputScale(output,totalInputsVector,excludedIndxs)
      
      
          # scale map
          outputMap = outputScale*netCost
      
      
      outputFilename = "Indirect_" + output + "." + theProblemManager.outputType
      
      theProblemManager.theRegionalCalculationManager.SaveMap(outputFilename, outputMap,theProblemManager.recordRange)
      
      # return result
      return outputMap
      
    
    def GetExcludedIndexes(self):
      """ Identify indexes of excluded entries from io matrix """
      indxs = np.where( [x in self.ioMatrixExcludedRows for x in self.ioMatrixRows ] )[0]
      self._ioMatrixReducedRows = list(self.ioMatrixRows)
      for i in indxs:
        del self._ioMatrixReducedRows[i]

      return  indxs
      
    def CalculateDemandVector(self,ioMatrix,excludedIndxs,regionSuffix = ""):
      """ Determing upstream sectors supporting project output """
      
      rv = None
      # find industry associated with output 
      rowIndxStr = self.industryColStringPrefix+regionSuffix
      
      
      # return column of inputs associated with demand
      if(rowIndxStr in self.ioMatrixRows):
        industryIndx = np.where( [x == rowIndxStr for x in self.ioMatrixRows ] )[0][0]
        rv = np.array(ioMatrix[:,industryIndx])
        rv = np.delete(rv, excludedIndxs)
        
      else:
        print("Error: Output column "+ rowIndxStr + " was not found." )
        exit(0)
        
      return rv
      
      
      
    def CalculateTotalRequirementsMatrix(self,excludedIndxs):
      """ Leontief inverse of IO matrix = (I-A)^{-1} ~= A + A^2 + A^3 + ... """
      
      niters = 1000 # fixme - this is overkill
      
      AA = np.array(self.ioMatrix)
      
      # remove excluded rows/columns
      AA =np.delete(AA,excludedIndxs,0)
      AA =np.delete(AA,excludedIndxs,1)
      
      reducedIOMatrix = np.array(AA)   # = reduced IO matrix "A"
      self.totalRequirementsMatrix = np.array(AA)  # = A + A^2 + A^3 etc

      for i in range(niters):
        AA = np.dot(AA,reducedIOMatrix)  # AA = A, A^2 etc. 
        self.totalRequirementsMatrix += AA
      
      return self.totalRequirementsMatrix
    
    def CalculateOutputScale(self,field,totalInputsVector,excludedIndxs):
      
      rv = 0.0
      
      ioIndx = np.where( [x == field for x in self._ioMatrixReducedRows ] )[0]
      if(ioIndx):
        ioIndx = ioIndx[0]
        rv = totalInputsVector[ioIndx]
      elif(self.secondaryMatrix):  # check if member of secondary matrix
        secondaryIndx = np.where( [x == field for x in self.secondaryMatrixRows ] )[0]
        if(secondaryIndx):
          secondaryIndx=secondaryIndx[0]
          SS = np.array(self.secondaryMatrix)
          SS = np.delete(SS,excludedIndxs,1)
          rv = np.dot(SS[secondaryIndx,:],totalInputsVector)
        else:
          print("Error: Output field "+ field + " was not found.")
          exit()
      else:
          print("Error: Output field "+ field + " was not found.")
          exit()
          
      return rv
      

    @classmethod
    def CreateFromXML(cls,solverNode,problemManager):
      """ Creates an UpstreamImpactCalculator instance from an xml node by initiating an empty instance and parsing the xml node """
      rv = cls()
      rv.ParseXMLNode(solverNode,problemManager)
      return rv
        
    def ParseXMLNode(self, node,problemManager):
      """
      Generate UpstreamImpactCalculator from xml tree node. 
      """
      self.ioMatrixPrefix = GetAttributeFileString(node,"ioMatrixPrefix")
      self.secondaryMatrixPrefix = GetAttributeFileStringOrDefault(node,"secondaryMatrixPrefix",self.secondaryMatrixPrefix)
      
      
      self.industryColStringPrefix = GetAttributeStringOrDefault(node,"industryColumnString",self.industryColStringPrefix)
      
      ioMatrixFilename = self.ioMatrixPrefix + ".npy"
      self.ioMatrix = np.load(ioMatrixFilename)
      
      ioMatrixRowsFilename = self.ioMatrixPrefix + "_rows.txt"
      fid = open(ioMatrixRowsFilename)
      self.ioMatrixRows = [x.strip() for x in fid.readlines()]
      
      if(self.secondaryMatrixPrefix):
        secondaryMatrixFilename = self.econdaryMatrixPrefix + ".npy"
        self.secondaryMatrix = np.load(secondaryMatrixFilename)
      
        secondaryMatrixRowsFilename = self.secondaryMatrixPrefix + "_rows.txt"
        self.secondaryMatrixRows = np.loadtxt(secondaryMatrixRowsFilename)        
      
      self.stateBasedIOMatrix =  GetAttributeFileStringOrDefault(node,"stateBasedIOMatrix",self.stateBasedIOMatrix) 
      self.varyOutputByState =  GetAttributeValueOrDefault(node,"varyOutputByState",self.varyOutputByState) 



    def WriteXMLNode(self, node):
      """
      Write data to xml node
      """
      
      SetAttributeString(node,"ioMatrixPrefix",self.ioMatrixPrefix)
      if(self.secondaryMatrixPrefix):
        SetAttributeString(node,"secondaryMatrixPrefix",self.secondaryMatrixPrefix)
      
      if(self.stateBasedIOMatrix):
        SetAttributeString(node,"stateBasedIOMatrix",self.stateBasedIOMatrix)
      
      return node


    class Factory:
      """ The factory that creates the UpstreamImpactCalculator Solver. """
      def CreateFromXML(self,xmlNode,problemManager): return UpstreamImpactCalculator.CreateFromXML(xmlNode,problemManager)



# Factory Registrator
SolverFactory.AddFactory(UpstreamImpactCalculator.typeStr,UpstreamImpactCalculator.Factory())



