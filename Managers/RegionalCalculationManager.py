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

from scipy import interpolate  # for interp1d, interpn
from scipy import ndimage

from scipy.sparse import load_npz  # may not be required (could load coords instead)

from skimage.draw import line_aa  # for mapping line onto array # todo should move to another file (to reduce dependencies)

import time  

# common
from Common.GeoUtils import LatLonDegreeToKMCoords, UsingGDAL
from Common.Common import BluecapError

# IO
from IO.XML import HasChild,GetChild,AddChild,GetChildren
from IO.XML import HasAttribute,GetAttributeString,GetAttributeStringOrDefault,GetAttributeFileString,GetAttributeFileStringOrDefault
from IO.XML import GetAttributeValue,GetAttributeVector, SetAttributeString,GetXMLTag

from IO.outputs import array_to_raster  # nb requires osgeo
from IO.outputs import GetMapOutputStyle

# Managers
from Units.UnitManager import UnitManager

from Functions.FunctionManager import FunctionManager

# Beta Distribution (for uncertainty calculations)
from Functions.BetaDistribution import ScaledBetaDistribution, UncertainInterp1d, GetExpectations, GetVariances





class RegionalCalculationManager():
    def __init__(self):
      """
      Create an empty regional calculation data manager and default variables. 
      """
      self.type = "NPV"   # NPV, benefit_cost_ratio, breakeven_grade
      
      
      # Mining
      ########
      
      self.coverDepthMapFile = ""
      self.coverDepthStdMapFile = "" # std deviation of cover depth
      self.waterDistanceMapFile = ""
      self.paleovalleyDistanceMapFile = ""
      self.powerDistanceMapFile = ""
      self.gasDistanceMapFile = ""
      self.co2DistanceMapFile = ""
      self.railDistanceMapFile = ""
      self.roadDistanceMapFile = ""
      self.railDistanceToPortMapFile = ""
      self.roadDistanceToPortMapFile = ""
      
      self.roadNetworkFile = ""
      self.roadNetworkOrigin = []
      self.railNetworkFile = ""
      self.railNetworkOrigin = []
      
      self.roadClosureMapFile = None
      self.railClosureMapFile = None
      
      # Hydrogen
      ##########
      
      self.pvCapacityFactorFile = ""
      self.cspCapacityFactorFile = ""
      self.windCapacityFactorFile = ""
      self.steamMethaneCapacityFactorFile = ""
      self.brownCoalCapacityFactorFile = ""
      self.blackCoalCapacityFactorFile = ""
      
      
      self.hybridCapacityFactorPrefix = None
      
      self.h2pipelineDistanceMapFile  = ""
      
      self.brownCoalTransportationDistanceMapFile = ""
      self.blackCoalTransportationDistanceMapFile = ""
      
      
      self.distanceToHydrogenStorageMapFile = ""
      self.railDistanceToHydrogenStorageMapFile = ""
      self.roadDistanceToHydrogenStorageMapFile = ""
      
      self.connectToPowerGrid = True # False for hydrogen as storage or local trasport
      self.doHydrogenStorageCalculation = False
      self.hydrogenStorageOnly = False
      
      # States
      ########
      self.stateStrs = ["NSW","NT","QLD","SA","WA","VIC","ACT","TAS"]
      self.stateIds = [0,1,2,3,4,5,6,7]

      
      # Added infrastructure
      ######################
      
      self.addedInfrastructure = {}
      self.addedInfrastructure["powerlines"] = []
      self.addedInfrastructure["waterlines"] = []
      self.addedInfrastructure["gaslines"] = []
      self.addedInfrastructure["roads"] = []
      self.addedInfrastructure["raillines"] = []
      self.addedInfrastructure["ports"] = []
      self.addedInfrastructure["refineries"] = []
      
      self.portAndRefiningCenterCoordinates = []
      
      # Uncertainty
      #############
      self.numRealizations = 1000
      
      self.coverOffset = 0.0
      self.coverMin    = 1e-06 #0.0
      self.coverMax    = 30000.0
      
      self.dataBoundingBox = []
      self.outputBoundingBox = []
      
      self.stride = 1

      self.stateIdsMapFile = ""
      
      self.saveResult = True  # only false for sensitivity calculations
      
      # Output flags
      self.doEmploymentCalculation = False
      self.recordTaxationCalculation = False
      self.recordWaterUse = False
      
      # Default data extents
      self.DATA_EXTENT = [112.0, -44.0, 154.0, -10.0]
      self.RESOLUTION  = [0.01, 0.01]
      self.LONGITUDES  = np.arange(self.DATA_EXTENT[0], self.DATA_EXTENT[2]+0.5*self.RESOLUTION[0], self.RESOLUTION[0])
      self.LATITUDES   = np.arange(self.DATA_EXTENT[1], self.DATA_EXTENT[3]+0.5*self.RESOLUTION[1], self.RESOLUTION[1])
      

    def ParseXMLNode(self, regionalDataNode):
      """
      Generate regional calculation data from xml tree node. 
      """
      
      self.type = GetAttributeStringOrDefault(regionalDataNode,"type",self.type)
          
      self.coverDepthMapFile = GetAttributeFileString(regionalDataNode,"coverDepth")
      if(HasAttribute(regionalDataNode,"coverDepthSTD")):
        self.coverDepthStdMapFile = GetAttributeFileString(regionalDataNode,"coverDepthSTD")
      self.waterDistanceMapFile = GetAttributeFileString(regionalDataNode,"waterDistance") 
      self.powerDistanceMapFile = GetAttributeFileString(regionalDataNode,"powerDistance")
      self.railDistanceMapFile = GetAttributeFileString(regionalDataNode,"railDistance")
      self.roadDistanceMapFile = GetAttributeFileString(regionalDataNode,"roadDistance")
      self.railDistanceToPortMapFile = GetAttributeFileString(regionalDataNode,"railTransportationDistance")
      self.roadDistanceToPortMapFile = GetAttributeFileString(regionalDataNode,"roadTransportationDistance")
      
      self.stateIdsMapFile = GetAttributeFileString(regionalDataNode,"states")
            
      self.gasDistanceMapFile = GetAttributeFileStringOrDefault(regionalDataNode,"gasDistance",self.gasDistanceMapFile)  # optional
      self.paleovalleyDistanceMapFile = GetAttributeFileStringOrDefault(regionalDataNode,"paleovalleyDistance",self.paleovalleyDistanceMapFile)  # optional
      
      self.roadClosureMapFile = GetAttributeFileStringOrDefault(regionalDataNode,"roadClosureDays",self.roadClosureMapFile)  # optional
      self.railClosureMapFile = GetAttributeFileStringOrDefault(regionalDataNode,"railClosureDays",self.railClosureMapFile)  # optional
      
      # uncertainty
      self.numRealizations = GetAttributeStringOrDefault(regionalDataNode,"samples",self.numRealizations)
      
      # hydrogen
      self.pvCapacityFactorFile = GetAttributeFileStringOrDefault(regionalDataNode,"pvCapacityFactor",self.pvCapacityFactorFile)  # optional
      self.cspCapacityFactorFile = GetAttributeFileStringOrDefault(regionalDataNode,"cspCapacityFactor",self.cspCapacityFactorFile)  # optional
      self.windCapacityFactorFile = GetAttributeFileStringOrDefault(regionalDataNode,"windCapacityFactor",self.windCapacityFactorFile)  # optional
      self.hybridCapacityFactorPrefix = GetAttributeFileStringOrDefault(regionalDataNode,"hybridCapacityPrefix",self.pvCapacityFactorFile)  # optional
      
      self.h2pipelineDistanceMapFile = GetAttributeFileStringOrDefault(regionalDataNode,"hydrogenPipelineDistance",self.h2pipelineDistanceMapFile)  # optional

      
      self.steamMethaneCapacityFactorFile = GetAttributeFileStringOrDefault(regionalDataNode,"steamMethaneCapacityFactor",self.steamMethaneCapacityFactorFile)  # optional
      self.brownCoalCapacityFactorFile = GetAttributeFileStringOrDefault(regionalDataNode,"blackCoalCapacityFactor",self.brownCoalCapacityFactorFile)  # optional
      self.blackCoalCapacityFactorFile = GetAttributeFileStringOrDefault(regionalDataNode,"brownCoalCapacityFactor",self.blackCoalCapacityFactorFile)  # optional
      
      self.co2DistanceMapFile = GetAttributeFileStringOrDefault(regionalDataNode,"co2Distance",self.co2DistanceMapFile)  # optional
      
      self.brownCoalTransportationDistanceMapFile = GetAttributeFileStringOrDefault(regionalDataNode,"brownCoalTransportationDistance",self.brownCoalTransportationDistanceMapFile) # optional
      self.blackCoalTransportationDistanceMapFile = GetAttributeFileStringOrDefault(regionalDataNode,"blackCoalTransportationDistance",self.blackCoalTransportationDistanceMapFile) # optional
      
      
      self.distanceToHydrogenStorageMapFile = GetAttributeFileStringOrDefault(regionalDataNode,"pipelineDistanceToHydrogenStorage",self.distanceToHydrogenStorageMapFile) # optional
      self.railDistanceToHydrogenStorageMapFile = GetAttributeFileStringOrDefault(regionalDataNode,"railDistanceToHydrogenStorage",self.railDistanceToHydrogenStorageMapFile) # optional
      self.roadDistanceToHydrogenStorageMapFile = GetAttributeFileStringOrDefault(regionalDataNode,"roadDistanceToHydrogenStorage",self.roadDistanceToHydrogenStorageMapFile) # optional
      

      
      if(HasAttribute(regionalDataNode,"dataBoundingBox")):
        self.dataBoundingBox = GetAttributeVector(regionalDataNode,"dataBoundingBox")
      if(HasAttribute(regionalDataNode,"outputBoundingBox")):
        self.outputBoundingBox = GetAttributeVector(regionalDataNode,"outputBoundingBox")
      if(HasAttribute(regionalDataNode,"coverOffset")):
        self.coverOffset = GetAttributeValue(regionalDataNode,"coverOffset")
      
      # Re-calculate input data resolution
      if(len(self.dataBoundingBox) == 4):
        dummy = np.load(self.powerDistanceMapFile)
        #print(self.powerDistanceMapFile)
        self.DATA_EXTENT = self.dataBoundingBox
        self.RESOLUTION  = [(self.DATA_EXTENT[2]-self.DATA_EXTENT[0])/dummy.shape[1], (self.DATA_EXTENT[3]-self.DATA_EXTENT[1])/dummy.shape[0]]
        self.LONGITUDES  = np.arange(self.DATA_EXTENT[0], self.DATA_EXTENT[2]+0.5*self.RESOLUTION[0], self.RESOLUTION[0])
        self.LATITUDES   = np.arange(self.DATA_EXTENT[1], self.DATA_EXTENT[3]+0.5*self.RESOLUTION[1], self.RESOLUTION[1])

      if(HasAttribute(regionalDataNode,"stride")):
        self.stride = GetAttributeValue(regionalDataNode,"stride")
        self.LONGITUDES  = self.LONGITUDES[self.stride//2::self.stride]
        self.LATITUDES   = self.LATITUDES[self.stride//2::self.stride]
        self.RESOLUTION  = [self.LONGITUDES[1]-self.LONGITUDES[0], self.LATITUDES[1]-self.LATITUDES[0]]
        
      if(HasAttribute(regionalDataNode,"estimateEmployment")):
        self.doEmploymentCalculation = GetAttributeValue(regionalDataNode,"estimateEmployment")
        
      if(HasAttribute(regionalDataNode,"estimateTax")):
        self.recordTaxationCalculation = GetAttributeValue(regionalDataNode,"estimateTax")
      
      if(HasAttribute(regionalDataNode,"estimateWaterUse")):
        self.recordWaterUse = GetAttributeValue(regionalDataNode,"estimateWaterUse")
        
      if(HasAttribute(regionalDataNode,"sellToHydrogenStorage")):
        self.doHydrogenStorageCalculation = GetAttributeValue(regionalDataNode,"sellToHydrogenStorage")
        if(HasAttribute(regionalDataNode,"hydrogenStorageOnly")):
          self.hydrogenStorageOnly = GetAttributeValue(regionalDataNode,"hydrogenStorageOnly")
        
      if(HasAttribute(regionalDataNode,"portCoordinates")):
        portFile = GetAttributeFileString(regionalDataNode,"portCoordinates")
        self.portAndRefiningCenterCoordinates = list(np.loadtxt(portFile))
        
      if(HasAttribute(regionalDataNode,"roadNetwork")):
        self.roadNetworkFile = GetAttributeFileString(regionalDataNode,"roadNetwork")
        networkOriginFile = GetAttributeFileString(regionalDataNode,"roadNetworkOrigin")
        self.roadNetworkOrigin = np.loadtxt(networkOriginFile)

      if(HasAttribute(regionalDataNode,"railNetwork")):
        self.railNetworkFile = GetAttributeFileString(regionalDataNode,"railNetwork")
        networkOriginFile = GetAttributeFileString(regionalDataNode,"railNetworkOrigin")
        self.railNetworkOrigin = np.loadtxt(networkOriginFile)
        
    def WriteXMLNode(self, node):
      """
      Write regional calculation to xml node
      """
      SetAttributeString(node,"type",self.type)
      
      return node
    
    
    def LoadMap(self,filename,offset=0.,Min=None,Max=None):
      if filename[-3:] == "npy":
        rv = np.load(filename)
      else:
        rv = np.loadtxt(filename)

      rv += offset
      # bound map
      if Min:
        rv[rv < Min] = Min
      if Max:
        rv[rv > Max] = Max
      
      if(self.stride > 1):
        rv = rv[self.stride//2::self.stride,self.stride//2::self.stride]
      return rv


    def LoadMapWUncertainty(self,meanfilename,stdDevFilename,offset=0.,Min=1e-06,Max=None):
      if meanfilename[-3:] == "npy":
        mu = np.load(meanfilename)
      else:
        mu = np.loadtxt(meanfilename)
      if stdDevFilename[-3:] == "npy":
        stddev = np.load(stdDevFilename)
      else:
        stddev = np.loadtxt(stdDevFilename)
      
      mu += offset
      # bound map
      if Min:
        mu[mu < Min] = Min
        # stddev[stddev < Min] = 0.
      if Max:
        mu[mu > Max] = Max
      
      stddev[stddev <= 1e-06] = 1e-06

      if(self.stride > 1):
        mu = mu[self.stride//2::self.stride,self.stride//2::self.stride]
        stddev = stddev[self.stride//2::self.stride,self.stride//2::self.stride]

      rv = ScaledBetaDistribution.FromMeanAndStd(mu,stddev,min=(mu-3*stddev), max=(mu+3*stddev))
      return rv    
    
    
    def LoadInfrastructureMap(self,theProblemManager,type,filename):
      rv = self.LoadMap(filename)
      
      if( (type in self.addedInfrastructure) and self.addedInfrastructure[type]):
      
        if(type in ["powerlines","waterlines","gaslines"]):
          rv = self.UpdateDistanceMap(self,self.addedInfrastructure[type],rv)
        
        else:
          raise BluecapError("Error unsupported infrastructure type"+ type)
        
      return rv
    
    
    def LoadWaterDistance(self,theProblemManager):
    
      rv = self.LoadInfrastructureMap(theProblemManager,"waterlines",self.waterDistanceMapFile)
      
      if(self.paleovalleyDistanceMapFile):
        rvb = self.LoadMap(self.paleovalleyDistanceMapFile)
        rv = np.minimum(rv,rvb)
        
      return rv
    
    
    def LoadInfrastructureDistanceAndTravelDistanceMaps(self,theProblemManager,type,distanceFilename, travelDistanceFilename,nLat,nLon):
      
      doGenerateNewMaps = False
      if ( (type in self.addedInfrastructure) and self.addedInfrastructure[type]):
        #i.e. new roads/raillines added
        doGenerateNewMaps = True
      elif( ("ports"  in self.addedInfrastructure) and self.addedInfrastructure["ports"]): 
        #i.e. new port(s) added
        doGenerateNewMaps = True
      elif( ("refineries" in self.addedInfrastructure) and self.addedInfrastructure["refineries"]): 
        #i.e. new refineries added
        doGenerateNewMaps = True
      
      if(doGenerateNewMaps):
        # adding additional routes to transport network or checking the existing routes
        
        networkData = []
        networkOrigin =[]
        
        if(type is "roads"):

          if(self.roadNetworkFile[-3:] == "npz"):
            sparseData = load_npz(self.roadNetworkFile)
            networkData = sparseData.toarray()
          else:
            networkData = np.load(self.roadNetworkFile)
          networkOrigin = self.roadNetworkOrigin
          
        elif(type is "raillines"):
        
          if(self.railNetworkFile[-3:] == "npz"):
            sparseData = load_npz(self.railNetworkFile)
            networkData = sparseData.toarray()
          else:
            networkData = np.load(self.railNetworkFile)
          networkOrigin = self.railNetworkOrigin
        else:
          raise BluecapError("Error unsupported/unrecognized infrastructure type" + type)
          
        distanceMap,travelMap = self.UpdateDistanceAndTravelMap(theProblemManager,self.addedInfrastructure[type],networkData,networkOrigin,nLat,nLon)
      else:
        distanceMap = self.LoadMap(distanceFilename)
        travelMap = self.LoadMap(travelDistanceFilename)
        
      return distanceMap,travelMap
    
    
    def SaveMap(self,filename, data, reportRange):
          
      if (not self.saveResult):   # set to false for sensitivity calculations
          return 0
      else: 

          if(filename[-3:] in ["jpg","tif","png"]  or  filename[-4:] == "tiff"):
            # geotiff
            if ((filename[-6:] == "geotif") or (filename[-7:] == "geotiff")) and UsingGDAL():
          
              stateIdsMap =  self.LoadMap(self.stateIdsMapFile)
              nLat,nLon = stateIdsMap.shape
          
              dataMinLon,dataMinLat,dataMaxLon,dataMaxLat = self.dataBoundingBox[:]
              outputMinLon,outputMinLat,outputMaxLon,outputMaxLat = self.outputBoundingBox[:]
              maxLatIndx = self.LATITUDES > outputMaxLat
              maxLatIndx = int( np.ceil( nLat* (outputMaxLat-dataMinLat)/(dataMaxLat-dataMinLat) ) ) - 1
              minLonIndx = int( np.floor( nLon* (outputMinLon-dataMinLon)/(dataMaxLon-dataMinLon) ) )
              minLon = self.LONGITUDES[minLonIndx] - 0.5*self.RESOLUTION[0]
              maxLat = self.LATITUDES[maxLatIndx] + 0.5*self.RESOLUTION[1]
      
              if (self.type == "NPV") or (self.type == "uncertainNPV") or (self.type == "tax") or (self.type == "HYDROGEN") or (self.type == "HYDROGEN CSP"):
                array_to_raster(filename.replace('geotif','tif'), data[::-1]*1e-06, xpix=self.RESOLUTION[0], ypix=self.RESOLUTION[1], xmin=minLon, ymax=maxLat)
              elif self.type == "water":
                array_to_raster(filename.replace('geotif','tif'), data[::-1]*1e-06, xpix=self.RESOLUTION[0], ypix=self.RESOLUTION[1], xmin=minLon, ymax=maxLat)
              elif self.type == "breakeven_grade":
                array_to_raster(filename.replace('geotif','tif'), data[::-1], xpix=self.RESOLUTION[0], ypix=self.RESOLUTION[1], xmin=minLon, ymax=maxLat)
              else:
                array_to_raster(filename.replace('geotif','tif'), data[::-1], xpix=self.RESOLUTION[0], ypix=self.RESOLUTION[1], xmin=minLon, ymax=maxLat)
              filename = filename.replace("tiff", "tif")
              filename = filename.replace("geotif",  "png")
            
            # regular output (and geotiff also)
            
            data = data.copy(order='C') # fixme - SW: not clear/documented why this is needed, little concerned that we are changing data when writing output. 
            
            if(self.type == "benefit_cost_ratio"):
              logData = np.log(data)
              cmap,data_max, data_min, label_min, label_max = GetMapOutputStyle(logData,self.type,style="GA")
              pl.imsave(filename,logData,origin="lower",cmap=pl.get_cmap(cmap),vmin=data_min,vmax=data_max)
            else:
              cmap,data_max, data_min, label_min, label_max = GetMapOutputStyle(data,self.type,style="GA")
              pl.imsave(filename,data,origin="lower",cmap=pl.get_cmap(cmap),vmin=data_min,vmax=data_max)
              
            
          elif( filename[-3:] == "npy"):
            np.save(filename,data)
          elif( filename[-3:] == "txt"):
            np.savetxt(filename,data)
          else:
            raise BluecapError("Error unrecognized output type: "+ filename)
            #print("Error unrecognized output type: ", filename)
      
          if(reportRange):
                np.savetxt(filename+"_range.txt",[np.round(label_min,2), np.round(label_max,2)])
      
      return 0
 
 
    def Run(self,problemManager):
      """
      Run regional calculation
      """
      
      if(self.type == "NPV" ):
        rv = self.RunNPVcalculation(problemManager)
      elif(self.type == "uncertainNPV"):
        rv = self.RunUncertainNPVcalculation(problemManager)
      elif(self.type == "benefit_cost_ratio"):
        rv = self.RunBenefitCostRatioCalculation(problemManager)
      elif(self.type == "breakeven_grade"):
        rv = self.RunBreakevenGradeCalculation(problemManager)
      elif(self.type == "employment"):
        rv = self.RunEmploymentCalculation(problemManager)
      elif(self.type == "HydrogenNPV" or self.type == "HydrogenPositiveNPV"):
        rv = self.RunHydrogenNPVcalculation(problemManager)
        if(self.type == "HydrogenPositiveNPV" ):
          rv[rv < 0] = 0.0
      else:
      
        raise BluecapError("Error: Unrecognized regional calculation: " +  self.type)
    
      return rv
 
    def RunHydrogenCalculation(self,problemManager):
      """
      Run regional calculations for hydrogen projects
      """
      
      if(self.type == "NPV" or self.type == "HydrogenNPV" ):
        rv = self.RunHydrogenNPVcalculation(problemManager)
        #elif(self.type == "employment"):
        #  rv = self.RunEmploymentCalculation(problemManager)
      else:
        raise BluecapError("Error: Unrecognized regional hydrogen calculation: " +  self.type)
    
      return rv
 
    
    def ApplyOutputBoundingBoxMask(self,stateIdsMap,replaceValue = -1):   
      """Masks off regions outside the outputBoundingBox."""
      if(len(self.dataBoundingBox) == 4 and len(self.outputBoundingBox) == 4):
          dataMinLon,dataMinLat,dataMaxLon,dataMaxLat = self.dataBoundingBox[:]
          outputMinLon,outputMinLat,outputMaxLon,outputMaxLat = self.outputBoundingBox[:]
          
          nLat,nLon = stateIdsMap.shape
          minLatIndx = int( np.floor( nLat* (outputMinLat-dataMinLat)/(dataMaxLat-dataMinLat) ) )
          maxLatIndx = int( np.ceil( nLat* (outputMaxLat-dataMinLat)/(dataMaxLat-dataMinLat) ) )
          
          minLonIndx = int( np.floor( nLon* (outputMinLon-dataMinLon)/(dataMaxLon-dataMinLon) ) )
          maxLonIndx = int( np.ceil( nLon* (outputMaxLon-dataMinLon)/(dataMaxLon-dataMinLon)  ) )
          
          
          if(minLatIndx > 0 ):
            stateIdsMap[:minLatIndx,:] = replaceValue
          
          if(outputMaxLat < dataMaxLat):
            stateIdsMap[maxLatIndx:,:] = replaceValue
        
          
          if(minLonIndx > 0 ):
            stateIdsMap[:,:minLonIndx] = replaceValue
          
          if(outputMaxLon < dataMaxLon):
            stateIdsMap[:,maxLonIndx:] = replaceValue
    
    def GetCoverDepths(self):
        """Returns a list of trial cover depths on which to base interpolation functions."""
        if self.coverMax > 2000:
          coverSteps = np.linspace(0.,self.coverMax-200.,20)
          coverSteps += coverSteps[1]
          coverSteps += 200
          coverDepths = np.append ( np.linspace(0.,200,21), coverSteps )
        else:
          coverDepths = np.linspace(0.,self.coverMax,21)
        return coverDepths
      
    
    def ClipToOutputBounds(self,data,outputBoundingBox = None):  
      """Clip data to rectangular region specified as fractions of the original box dimension."""
    
      if( outputBoundingBox is None):
         outputBoundingBox = self.outputBoundingBox  # can't set self.outputBoundingBox as default
         
      # clips the output to within the bounding box
      if(len(self.dataBoundingBox) == 4 and len(outputBoundingBox) == 4):
          dataMinLon,dataMinLat,dataMaxLon,dataMaxLat = self.dataBoundingBox[:]
          outputMinLon,outputMinLat,outputMaxLon,outputMaxLat = self.outputBoundingBox[:]
          
          nLat,nLon = data.shape
          minLatIndx = int( np.floor( nLat* (outputMinLat-dataMinLat)/(dataMaxLat-dataMinLat) ) )
          maxLatIndx = int( np.ceil( nLat* (outputMaxLat-dataMinLat)/(dataMaxLat-dataMinLat) ) )
          
          minLonIndx = int( np.floor( nLon* (outputMinLon-dataMinLon)/(dataMaxLon-dataMinLon) ) )
          maxLonIndx = int( np.ceil( nLon* (outputMaxLon-dataMinLon)/(dataMaxLon-dataMinLon)  ) )
          
          data = data[minLatIndx:maxLatIndx,minLonIndx:maxLonIndx]
      return data
    
    def ClipToOutputLatLon(self,data,outputBoundingBox = None):  
      """Clip data to rectangular region specified in lat,lon coordinates."""
      outputBoundingBox = np.array(outputBoundingBox)
         
      # clips the output to within the bounding box
      if(len(self.dataBoundingBox) == 4 and len(outputBoundingBox) == 4):
          # dataMinLat,dataMinLon,dataMaxLat,dataMaxLon = self.dataBoundingBox[:]
          dataMinLon,dataMinLat,dataMaxLon,dataMaxLat = self.dataBoundingBox[:]
          outputMinLat,outputMinLon,outputMaxLat,outputMaxLon = outputBoundingBox[:]
          
          nLat,nLon = data.shape
          minLatIndx = int( np.floor( nLat* (outputMinLat-dataMinLat)/(dataMaxLat-dataMinLat) ) )
          maxLatIndx = int( np.ceil( nLat* (outputMaxLat-dataMinLat)/(dataMaxLat-dataMinLat) ) )
          
          minLonIndx = int( np.floor( nLon* (outputMinLon-dataMinLon)/(dataMaxLon-dataMinLon) ) )
          maxLonIndx = int( np.ceil( nLon* (outputMaxLon-dataMinLon)/(dataMaxLon-dataMinLon)  ) )
          
          data = data[minLatIndx:maxLatIndx,minLonIndx:maxLonIndx]
      return data
      
    def GetEffectiveDiscountRate(self,theProjectManager,calculationType="atNPV"):
       """Return the discount rate employed in the calcuation.
       
       Note: 
           This function wraps the discount rate used in the regional calculations to simplify swapping between company and government outputs.
           The discount rate does not affect the discount rate recorded internally by the economic data manager - it only affects the discount rate applied in the regional calculation manager.
       """
       
       rv = theProjectManager.theEconomicDataManager.discountRate
       if (calculationType == "NetRealCost"):
         rv = 0.0
       
       return rv
    
    def GetEffectiveTaxRate(self,theProjectManager,calculationType="atNPV"):
       """Return the effective tax rate employed in the calcuation.
       
       Note: 
           This function wraps the discount rate used in the regional calculations to simplify swapping between company and government outputs.
           The effective tax rate does not affect the tax rate recorded internally by the economic data manager - it only affects the tax rate applied in the regional calculation manager.
       """
       
       rv = theProjectManager.theEconomicDataManager.incomeTaxRate 
       if (calculationType == "NetRealCost"):
         rv = 0.0
       
       return rv 
      
            
    #########################
    ## NPV calculation
    
    def RunNPVcalculation(self,theProblemManager,valueType="atNPV"):
        """Calculate the regional Net Present Value for mining projects."""
    
        doPlots = (theProblemManager.outputType == "plots")

        doAll = (theProblemManager.outputType == "all")
        
        # profiling
        doTimer = False
        timesList = []

        if doTimer:
          timesList.append(["start", time.time()])
        
          
        theProblemManager.theMineDataManager.theInfrastructureManager.ZeroDistanceToInfrastructure()

        # calculate mine cost, life etc as function of depth of cover (in meters)
        
        coverDepths = self.GetCoverDepths()
  
        discountRate = self.GetEffectiveDiscountRate(theProblemManager.theMineDataManager,valueType)

        mineValues = []
        mineTypes =[]
        mineYears = []
        concentrateDiscountedTotalMasses = []
        concentrateCapacities = []
        
        mineCapExs = []  # for employment calculation
        stateAndFederalTaxes = {}  # for tax calculation
        waterUse = []  # for water use calculation

        stateMineValues = {}
        
        for state in self.stateStrs:
          stateMineValues[state] = []
          if(self.recordTaxationCalculation):
            stateAndFederalTaxes[state] = []

        for depthOfCover in coverDepths:
        
          theProblemManager.theMineDataManager.SetCoverDepth(depthOfCover)
  
          value = theProblemManager.EvaluateOCUGMines()
          mineValues.append(value)
  
          mineLife = theProblemManager.theMineDataManager.theMiningSystem.mineLife
          mineYears.append(mineLife)
          
          if(self.doEmploymentCalculation):
            mineAndProcCapex = theProblemManager.theMineDataManager.theEconomicDataManager.netCapex[0]
            mineCapExs.append(mineAndProcCapex)
          
          if(self.recordWaterUse):
            waterUsePerYear = theProblemManager.theMineDataManager.theInfrastructureManager.CalculateWaterUsePerYear(theProblemManager)
            waterUse.append(waterUsePerYear)
          
            
  
          mineTypeStr = theProblemManager.theMineDataManager.theMiningSystem.miningMethod  # 1 if UG, 0.5 mixed, 0 OC 
  
          mineType = 0.0   # open cut
          if(mineTypeStr[:2] == "UG"):
            mineType = 1.0
          elif(mineTypeStr == "OCUG"):
            mineType = 0.5
    
          for state in self.stateStrs:
            theProblemManager.theMineDataManager.theEconomicDataManager.SetState(state)
            theProblemManager.theMineDataManager.CalculateTaxes(theProblemManager)
            theProblemManager.theMineDataManager.CalculateAfterTaxCashFlow(theProblemManager)
            
            value = theProblemManager.theMineDataManager.CalculateEconomicIndicator(theProblemManager,valueType)
            
            stateMineValues[state].append(value)
            
            if(self.recordTaxationCalculation):
              taxNPV = np.sum(theProblemManager.theMineDataManager.theEconomicDataManager.btNCF - theProblemManager.theMineDataManager.theEconomicDataManager.atNCF)
              #taxNPV = theProblemManager.theMineDataManager.theEconomicDataManager.btNPV - theProblemManager.theMineDataManager.theEconomicDataManager.atNPV 
              stateAndFederalTaxes[state].append(taxNPV)
    
  
          mineTypes.append(mineType)
  
          # discounted total mass
          discountFactors = (1+theProblemManager.theMineDataManager.theEconomicDataManager.discountRate)**np.array(range( int( np.ceil(mineLife) ) ) )
          concentrateMass = np.sum(theProblemManager.theMineDataManager.theProcessingSystem.GetTotalConcentrateProduced()/discountFactors)
          concentrateDiscountedTotalMasses.append(concentrateMass)
  
          concentrateCapacity = np.max(theProblemManager.theMineDataManager.theProcessingSystem.GetTotalConcentrateProduced())
          concentrateCapacities.append(concentrateCapacity)
  
        if doTimer:
          timesList.append(["coverFunc", time.time()])
          
        #mineValues = np.array(mineValues,dtype=np.float)
        #mineTypes = np.array(mineTypes,dtype=np.float)
        #mineYears = np.array(mineYears,dtype=np.float)
        #concentrateDiscountedTotalMasses = np.array(concentrateDiscountedTotalMasses,dtype=np.float)
        #concentrateCapacities = np.array(concentrateCapacities,dtype=np.float)

        mineValueFunc = interpolate.interp1d(coverDepths,mineValues)
        mineTypeFunc = interpolate.interp1d(coverDepths,mineTypes)
        mineLifeFunc = interpolate.interp1d(coverDepths,mineYears)
        concentrateDiscountedTotalMassFunc = interpolate.interp1d(coverDepths,concentrateDiscountedTotalMasses)
        concentrateCapacityFunc = interpolate.interp1d(coverDepths,concentrateCapacities)
        
        if(self.doEmploymentCalculation):
          minecapexFunc = interpolate.interp1d(coverDepths,mineCapExs)
        
            
        if(doPlots):

            for state in self.stateStrs:
              pl.plot(coverDepths,stateMineValues[state],label=state)
          
            pl.legend()
            pl.plot(coverDepths,mineValues)
            pl.figure()
            pl.plot(coverDepths,mineTypes)

        ## regional mine cost calculation
        
        if doTimer:
          #startMapsTime = time.time()
          timesList.append(["startMaps", time.time()])

        stateIdsMap =  self.LoadMap(self.stateIdsMapFile)
        self.ApplyOutputBoundingBoxMask(stateIdsMap) 
        


        ###################################
        
        coverMap =  self.LoadMap(self.coverDepthMapFile,offset=self.coverOffset,Min=self.coverMin,Max=self.coverMax)  
          
        countryIndxs = stateIdsMap >= 0.0

        mineValueMap = np.zeros(coverMap.shape)
        
        taxRate = self.GetEffectiveTaxRate(theProblemManager.theMineDataManager,valueType)
        taxRelief = 1.0 - taxRate  # = tax relief from federal income tax (for costs)
        
        
        if(self.recordTaxationCalculation):
            taxEstimate = np.zeros(coverMap.shape)
          

        ####
        for i in range(len(self.stateStrs)):
          stateIndxs = stateIdsMap == self.stateIds[i]
          mineValueFunc = interpolate.interp1d(coverDepths,stateMineValues[self.stateStrs[i]])
          mineValueMap[stateIndxs] = mineValueFunc(coverMap[stateIndxs])

          if(self.recordTaxationCalculation):
            stateAndFederalTaxFunc = interpolate.interp1d(coverDepths,stateAndFederalTaxes[self.stateStrs[i]])
            taxEstimate[stateIndxs] = stateAndFederalTaxFunc(coverMap[stateIndxs])
        
        
        if(doAll):
          filename = theProblemManager.outputPrefix+"_mineValue.npy"  
          self.SaveMap(filename, mineValueMap,False)

        if doTimer:
          timesList.append(["mineValue", time.time()])

        # Water Expenses
        distanceToWater =  self.LoadWaterDistance(theProblemManager)
        mineValueMap[countryIndxs] -=  taxRelief* theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalWaterExpenses(theProblemManager,distanceToWater[countryIndxs])
        
        if(self.recordTaxationCalculation):
          taxEstimate[countryIndxs] -= taxRate*theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalWaterExpenses(theProblemManager,distanceToWater[countryIndxs])


        if(doAll):
          filename = theProblemManager.outputPrefix+"_waterCost.npy"  
          dd = np.zeros(coverMap.shape)
          dd[countryIndxs] = taxRelief* theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalWaterExpenses(theProblemManager,distanceToWater[countryIndxs])
          #print "Saving: ", filename
          self.SaveMap(filename, dd,False)

        if doTimer:
          timesList.append(["water", time.time()])

        # Power Expenses
        distanceToPower =  self.LoadInfrastructureMap(theProblemManager,"powerlines",self.powerDistanceMapFile)
        
        if(theProblemManager.theMineDataManager.theInfrastructureManager.calculateGas):
          powerCost = theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalPowerExpenses(theProblemManager,distanceToPower[countryIndxs])
          
          distanceToGas =  self.LoadInfrastructureMap(theProblemManager,"gaslines",self.gasDistanceMapFile)
          gasCost =  theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalGasExpenses(theProblemManager,distanceToGas[countryIndxs])

          mineValueMap[countryIndxs] -=  taxRelief*  np.minimum(powerCost,gasCost)
          
          if(self.recordTaxationCalculation):
            taxEstimate[countryIndxs]  -= taxRate*np.minimum(powerCost,gasCost)

        else:
          mineValueMap[countryIndxs] -=  taxRelief* theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalPowerExpenses(theProblemManager,distanceToPower[countryIndxs])

          if(self.recordTaxationCalculation):
            taxEstimate[countryIndxs]  -= taxRate*theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalPowerExpenses(theProblemManager,distanceToPower[countryIndxs])




        if(doAll):
          filename = theProblemManager.outputPrefix+"_powerCost.npy"  
          dd = np.zeros(coverMap.shape)
          dd[countryIndxs] = taxRelief* theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalPowerExpenses(theProblemManager,distanceToPower[countryIndxs])
          #print "Saving: ", filename
          self.SaveMap(filename, dd,False)
          
          if(theProblemManager.theMineDataManager.theInfrastructureManager.calculateGas):
              filename = theProblemManager.outputPrefix+"_gasCost.npy"  
              dd = np.zeros(coverMap.shape)
              dd[countryIndxs] = taxRelief* theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalGasExpenses(theProblemManager,distanceToGas[countryIndxs])
              #print "Saving: ", filename
              self.SaveMap(filename, dd,False)

        if doTimer:
          timesList.append(["power", time.time()])

        # Transportation Expenses
        nLat, nLon = coverMap.shape 
        distanceToRoad,roadTransportationDistance = self.LoadInfrastructureDistanceAndTravelDistanceMaps(theProblemManager,"roads",self.roadDistanceMapFile,self.roadDistanceToPortMapFile,nLat, nLon)
        distanceToRail,railTransporationDistance = self.LoadInfrastructureDistanceAndTravelDistanceMaps(theProblemManager,"raillines",self.railDistanceMapFile,self.railDistanceToPortMapFile,nLat, nLon)
        
        
        if(theProblemManager.theMineDataManager.theInfrastructureManager.doTransportClosureCostEstimate):
        
            roadClosureDays = self.LoadMap(self.roadClosureMapFile)
            if(self.railClosureMapFile): 
              railClosureDays = self.LoadMap(self.railClosureMapFile)
            else:
              railClosureDays = np.zeros(coverMap.shape )  # just include effect of road closure
              
            mineValueMap[countryIndxs] -= taxRelief* theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalTransportationExpenses( \
                                                              distanceToRoad[countryIndxs], roadTransportationDistance[countryIndxs], \
                                                              distanceToRail[countryIndxs], railTransporationDistance[countryIndxs], \
                                                              coverMap[countryIndxs],concentrateCapacityFunc,concentrateDiscountedTotalMassFunc,
                                                              roadClosureDays[countryIndxs],railClosureDays[countryIndxs])
        else:
            mineValueMap[countryIndxs] -= taxRelief* theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalTransportationExpenses( \
                                                              distanceToRoad[countryIndxs], roadTransportationDistance[countryIndxs], \
                                                              distanceToRail[countryIndxs], railTransporationDistance[countryIndxs], \
                                                              coverMap[countryIndxs],concentrateCapacityFunc,concentrateDiscountedTotalMassFunc)


        if(self.recordTaxationCalculation):
            taxEstimate[countryIndxs]  -= taxRate*theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalTransportationExpenses( \
                                                              distanceToRoad[countryIndxs], roadTransportationDistance[countryIndxs], \
                                                              distanceToRail[countryIndxs], railTransporationDistance[countryIndxs], \
                                                              coverMap[countryIndxs],concentrateCapacityFunc,concentrateDiscountedTotalMassFunc)


        if(doAll):
          filename = theProblemManager.outputPrefix+"_transportCost.npy"  
          dd = np.zeros(coverMap.shape)
          dd[countryIndxs] = taxRelief* theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalTransportationExpenses( \
                                                              distanceToRoad[countryIndxs], roadTransportationDistance[countryIndxs], \
                                                              distanceToRail[countryIndxs], railTransporationDistance[countryIndxs], \
                                                              coverMap[countryIndxs],concentrateCapacityFunc,concentrateDiscountedTotalMassFunc)
          #print "Saving: ", filename
          self.SaveMap(filename, dd,False)


        if doTimer:
          timesList.append(["transport", time.time()])


        # no tax if no profit
        if(self.recordTaxationCalculation):
          taxEstimate[mineValueMap <= 0.0] = 0 
          
        
        # clean up map
        mineValueMap[stateIdsMap< 0.0] = np.nan  
        
        if(self.recordTaxationCalculation):
          
          filename = theProblemManager.outputPrefix+"_tax.npy"  
          
          mineLife = mineLifeFunc(coverMap)
          taxEstimate[countryIndxs==False] = np.nan

        if(self.doEmploymentCalculation):
          minecapexFunc = interpolate.interp1d(coverDepths,mineCapExs)
          
          initialCapex = np.zeros(coverMap.shape)
          positiveNPVindxs = np.nan_to_num(mineValueMap) >0
          initialCapex[positiveNPVindxs] = minecapexFunc(coverMap[positiveNPVindxs])
          
          employmentEstimate = np.zeros(coverMap.shape)
          employmentEstimate = theProblemManager.theMineDataManager.theEconomicDataManager.EstimateDirectEmployment(initialCapex)
          employmentEstimate[countryIndxs==False] = np.nan
          
          # clip to bounding box
          if(len(self.dataBoundingBox) == 4 and len(self.outputBoundingBox) == 4):
            employmentEstimate = self.ClipToOutputBounds(employmentEstimate)
            


        if(self.recordWaterUse):
          waterUseFunc = interpolate.interp1d(coverDepths,waterUse)
          
          waterUseEstimate = np.zeros(coverMap.shape)
          positiveNPVindxs = np.nan_to_num(mineValueMap) >0
          waterUseEstimate = waterUseFunc(coverMap)
          # waterUseEstimate[positiveNPVindxs] = waterUseFunc(coverMap[positiveNPVindxs])
          waterUseEstimate[countryIndxs==False] = np.nan
          
          # clip to bounding box
          if(len(self.dataBoundingBox) == 4 and len(self.outputBoundingBox) == 4):
            waterUseEstimate = self.ClipToOutputBounds(waterUseEstimate)
          
        
        # clip to bounding box
        if(len(self.dataBoundingBox) == 4 and len(self.outputBoundingBox) == 4):
          mineValueMap = self.ClipToOutputBounds(mineValueMap)



        if doTimer:
          timesList.append(["end", time.time()])

        ######################################################################


        if(doPlots):
          pl.figure()
          
          pl.imshow(mineValueMap,origin="lower")

          pl.colorbar()

        ######################################################################

        if doTimer:
            pl.figure()
            dts = np.array( [tt[1] for tt in timesList] )
            dts -=dts[0]
            labels =  [tt[0] for tt in timesList] 

            pl.bar(range(len(dts)),dts,align="center")
            pl.xticks(range(len(dts)), labels)


            pl.bar(range(1,len(dts)),dts[1:]-dts[:-1],color="r",align="center")

            for i in range(1,len(dts)):
              print(labels[i], dts[i], dts[i]-dts[i-1])


        
        if(doPlots):
          pl.show()    
        else:
          filename = theProblemManager.outputPrefix+"."+theProblemManager.outputType
          
          if(doAll):
            filename = theProblemManager.outputPrefix+".npy"
            
          # Save the maps
          self.SaveMap(filename, mineValueMap,theProblemManager.recordRange)
          if(self.doEmploymentCalculation):
            filename = theProblemManager.outputPrefix+"_employment."+theProblemManager.outputType
            self.type = "employment"
            self.SaveMap(filename, employmentEstimate, theProblemManager.recordRange)
          if(self.recordTaxationCalculation):
            filename = theProblemManager.outputPrefix+"_taxes."+theProblemManager.outputType
            self.type = "tax"
            self.SaveMap(filename, taxEstimate, theProblemManager.recordRange)
          if(self.recordWaterUse):
            filename = theProblemManager.outputPrefix+"_waterUse."+theProblemManager.outputType
            self.type = "water"
            self.SaveMap(filename, waterUseEstimate, theProblemManager.recordRange)
        
        return mineValueMap
    
    
    #########################
    ## Benefit to cost ratio
    
    def RunBenefitCostRatioCalculation(self,theProblemManager):
        """Calculate ratio of discounted returns to discounted costs for the project over a region."""
        # profiling
        doTimer = False
        timesList = []
        
        doPlots = (theProblemManager.outputType == "plots")
        
        theProblemManager.theMineDataManager.theInfrastructureManager.ZeroDistanceToInfrastructure()

        if doTimer:
          timesList.append(["start", time.time()])


        # calculate mine cost, life etc as function of depth of cover (in meters)
        
        coverDepths = self.GetCoverDepths()

        discountRate = self.GetEffectiveDiscountRate(theProblemManager.theMineDataManager)

        mineValues = []
        mineYears = []
        concentrateDiscountedTotalMasses = []
        concentrateCapacities = []
        discountedRevenues = []


        stateDiscountedCosts = {}
        for state in self.stateStrs:
          stateDiscountedCosts[state] = []

        for depthOfCover in coverDepths:

          theProblemManager.theMineDataManager.SetCoverDepth(depthOfCover)
  
          value = theProblemManager.EvaluateOCUGMines()
          mineValues.append(value)
  
          mineLife = theProblemManager.theMineDataManager.theMiningSystem.mineLife
          mineYears.append(mineLife)
  
          discountFactors = (1+discountRate)**np.array(range( int( np.ceil(mineLife) ) ) )

  
          discountedRevenue = np.sum( theProblemManager.theMineDataManager.theEconomicDataManager.revenue*discountFactors)
  
          discountedRevenues.append(discountedRevenue)
    
          for state in self.stateStrs:
            theProblemManager.theMineDataManager.theEconomicDataManager.SetState(state)
            theProblemManager.theMineDataManager.CalculateTaxes(theProblemManager)
            theProblemManager.theMineDataManager.CalculateAfterTaxCashFlow(theProblemManager)
            theProblemManager.theMineDataManager.CalculateEconomicIndicators(theProblemManager)
    
            discountedCost = np.sum( ( theProblemManager.theMineDataManager.theEconomicDataManager.revenue \
                                       - theProblemManager.theMineDataManager.theEconomicDataManager.btNCF \
                                       + theProblemManager.theMineDataManager.theEconomicDataManager.royalties )*discountFactors)
            stateDiscountedCosts[state].append(discountedCost)
  
          # discounted total mass
          concentrateMass = np.sum(theProblemManager.theMineDataManager.theProcessingSystem.GetTotalConcentrateProduced()*discountFactors)
          concentrateDiscountedTotalMasses.append(concentrateMass)
  
          concentrateCapacity = np.max(theProblemManager.theMineDataManager.theProcessingSystem.GetTotalConcentrateProduced())
          concentrateCapacities.append(concentrateCapacity)
  

  
        if doTimer:
          timesList.append(["coverFunc", time.time()])
        
        
        #mineYears = np.array(mineYears,dtype=np.float)
        #concentrateDiscountedTotalMasses = np.array(concentrateDiscountedTotalMasses,dtype=np.float)
        #concentrateCapacities = np.array(concentrateCapacities,dtype=np.float)

        mineLifeFunc = interpolate.interp1d(coverDepths,mineYears)
        concentrateDiscountedTotalMassFunc = interpolate.interp1d(coverDepths,concentrateDiscountedTotalMasses)
        concentrateCapacityFunc = interpolate.interp1d(coverDepths,concentrateCapacities)

        mineRevenueFunc = interpolate.interp1d(coverDepths,discountedRevenues)


        
        ## regional mine cost calculation
        
        if doTimer:
          timesList.append(["startMaps", time.time()])

        stateIdsMap =  self.LoadMap(self.stateIdsMapFile)
        self.ApplyOutputBoundingBoxMask(stateIdsMap) 

        countryIndxs = stateIdsMap >= 0.0

        ###################################
        
        coverMap =  self.LoadMap(self.coverDepthMapFile,offset=self.coverOffset,Min=self.coverMin,Max=self.coverMax)

        mineCostMap = np.zeros(coverMap.shape)


        ####
        for i in range(len(self.stateStrs)):
          stateIndxs = stateIdsMap == self.stateIds[i]
          mineCostFunc = interpolate.interp1d(coverDepths,stateDiscountedCosts[self.stateStrs[i]])
          mineCostMap[stateIndxs] = mineCostFunc(coverMap[stateIndxs])

        if doTimer:
          timesList.append(["mineValue", time.time()])

        # Water Expenses
        distanceToWater =  self.LoadMap(self.waterDistanceMapFile)
        mineCostMap[countryIndxs] +=  theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalWaterExpenses(theProblemManager,distanceToWater[countryIndxs])

        if doTimer:
          timesList.append(["water", time.time()])

        # Power Expenses
        distanceToPower =  self.LoadMap(self.powerDistanceMapFile)
        mineCostMap[countryIndxs] +=  theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalPowerExpenses(theProblemManager,distanceToPower[countryIndxs])


        if doTimer:
          timesList.append(["power", time.time()])

        # Transportation Expenses
        distanceToRoad =  self.LoadMap(self.roadDistanceMapFile)
        distanceToRail =  self.LoadMap(self.railDistanceMapFile)
        roadTransportationDistance =  self.LoadMap(self.roadDistanceToPortMapFile)
        railTransporationDistance =  self.LoadMap(self.railDistanceToPortMapFile)


        mineCostMap[countryIndxs] += theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalTransportationExpenses( \
                                                              distanceToRoad[countryIndxs], roadTransportationDistance[countryIndxs], \
                                                              distanceToRail[countryIndxs], railTransporationDistance[countryIndxs], \
                                                              coverMap[countryIndxs],concentrateCapacityFunc,concentrateDiscountedTotalMassFunc)


        if doTimer:
          timesList.append(["transport", time.time()])


        revenueCostRatioMap = np.zeros(mineCostMap.shape)
        revenueCostRatioMap[countryIndxs] =  mineRevenueFunc( coverMap[countryIndxs] )/ mineCostMap[countryIndxs]
        
        revenueCostRatioMap[stateIdsMap< 0.0] = np.nan  

        # clip to bounding box
        if(len(self.dataBoundingBox) == 4 and len(self.outputBoundingBox) == 4):
          revenueCostRatioMap = self.ClipToOutputBounds(revenueCostRatioMap)


        if doTimer:
          timesList.append(["end", time.time()])

        ######################################################################

        
        if(doPlots):
        
          pl.figure()
          pl.imshow(revenueCostRatioMap,origin="lower")
          pl.colorbar()
        
          pl.show()    
        else:
          filename = theProblemManager.outputPrefix+"."+theProblemManager.outputType
          self.SaveMap(filename, revenueCostRatioMap,theProblemManager.recordRange)

        ######################################################################

        if doTimer:
            pl.figure()
            dts = np.array( [tt[1] for tt in timesList] )
            dts -=dts[0]
            labels =  [tt[0] for tt in timesList] 

            pl.bar(range(len(dts)),dts,align="center")
            pl.xticks(range(len(dts)), labels)


            pl.bar(range(1,len(dts)),dts[1:]-dts[:-1],color="r",align="center")

            for i in range(1,len(dts)):
              print(labels[i], dts[i], dts[i]-dts[i-1])

    
        return revenueCostRatioMap
    
    
    #########################
    ## Breakeven calculation
    
    def RunBreakevenGradeCalculation(self,theProblemManager,valueType="atNPV"):
        """Estimate the Regional Breakeven grade"""
    
        assert len(theProblemManager.theMineDataManager.theOreBody.metalGrades.keys()) == 1
        
        theProblemManager.theMineDataManager.theInfrastructureManager.ZeroDistanceToInfrastructure()
        
        # plotting
        doPlots = (theProblemManager.outputType == "plots")

        # profiling
        doTimer = False
        timesList = []

        if doTimer:
          timesList.append(["start", time.time()])

        # calculate mine cost, life etc as function of depth of cover (in meters)
          
        coverDepths = self.GetCoverDepths()
  
        discountRate = self.GetEffectiveDiscountRate(theProblemManager.theMineDataManager,valueType)

        mineValues = []
        mineTypes =[]
        mineYears = []
        concentrateDiscountedTotalMasses = []
        concentrateCapacities = []


        stateMineValues = {}
        stateMineValuesDoubleGrade = {}
        stateMineValuesHalfGrade = {}


        ### check that current grade gives economic results
        ### if not increase grade by a factor (gradeAdjustment) until an economic grade is encountered at depth=0
        gradeAdjustment = 1.0
        theProblemManager.theMineDataManager.SetCoverDepth(coverDepths[10])
        
        
        theProblemManager.EvaluateOCUGMines()
        value = theProblemManager.theMineDataManager.theEconomicDataManager.atNPV
          
        theProblemManager.theMineDataManager.theOreBody.ScaleCommodityGrades(1./gradeAdjustment)
        
        adjustUp = True
        if(value >= 0.0):
            gradeAdjustment *= 0.75
            adjustUp =False
        else:
            gradeAdjustment *= 1.5
        
        
        for i in range(20):
          
          theProblemManager.theMineDataManager.theOreBody.ScaleCommodityGrades(gradeAdjustment)
        
          theProblemManager.EvaluateOCUGMines()
          value = theProblemManager.theMineDataManager.theEconomicDataManager.atNPV
          
          theProblemManager.theMineDataManager.theOreBody.ScaleCommodityGrades(1./gradeAdjustment)
          
          if(value >= 0.0 and not adjustUp):
            gradeAdjustment *= 0.75
          elif(value <= 0.0 and adjustUp):
            gradeAdjustment *= 1.5
          else:
            #print "gradeAdjustment", gradeAdjustment, value
            break
        

        ### original mine values

        for state in self.stateStrs:
          stateMineValues[state] = []
          stateMineValuesHalfGrade[state] = []
          stateMineValuesDoubleGrade[state] = []

        for depthOfCover in coverDepths:

          theProblemManager.theMineDataManager.SetCoverDepth(depthOfCover)
  
          # original grades (multiplied by grade adjustment)
          theProblemManager.theMineDataManager.theOreBody.ScaleCommodityGrades(gradeAdjustment)
  
          value = theProblemManager.EvaluateOCUGMines()
          mineValues.append(value)
  
          mineLife = theProblemManager.theMineDataManager.theMiningSystem.mineLife
          mineYears.append(mineLife)
  
          mineTypeStr = theProblemManager.theMineDataManager.theMiningSystem.miningMethod  # 1 if UG, 0.5 mixed, 0 OC 
  
          mineType = 0.0
          if(mineTypeStr[:2] == "UG"):
            mineType = 1.0
          elif(mineTypeStr == "OCUG"):
            mineType = 0.5
  
          mineTypes.append(mineType)
    
          for state in self.stateStrs:
            theProblemManager.theMineDataManager.theEconomicDataManager.SetState(state)
            theProblemManager.theMineDataManager.CalculateTaxes(theProblemManager)
            theProblemManager.theMineDataManager.CalculateAfterTaxCashFlow(theProblemManager)
            value = theProblemManager.theMineDataManager.CalculateEconomicIndicator(theProblemManager,valueType)
            
            stateMineValues[state].append(value)
    
          # discounted total mass - these will need to be scaled according to changes in grade in transportation cost calculation
          discountFactors = (1+theProblemManager.theMineDataManager.theEconomicDataManager.discountRate)**np.array(range( int( np.ceil(mineLife) ) ) )
          concentrateMass = np.sum(theProblemManager.theMineDataManager.theProcessingSystem.GetTotalConcentrateProduced()/discountFactors)
          concentrateDiscountedTotalMasses.append(concentrateMass)
  
          concentrateCapacity = np.max(theProblemManager.theMineDataManager.theProcessingSystem.GetTotalConcentrateProduced())
          concentrateCapacities.append(concentrateCapacity)
  
          ### half grade ###
          theProblemManager.theMineDataManager.theOreBody.ScaleCommodityGrades(0.5)
  
          for state in self.stateStrs:
            theProblemManager.theMineDataManager.theEconomicDataManager.CalculateBeforeTaxCashFlow(theProblemManager, theProblemManager.theMineDataManager)
            theProblemManager.theMineDataManager.theEconomicDataManager.SetState(state)
            theProblemManager.theMineDataManager.CalculateTaxes(theProblemManager)
            theProblemManager.theMineDataManager.CalculateAfterTaxCashFlow(theProblemManager)
            halfGradeValue = theProblemManager.theMineDataManager.CalculateEconomicIndicator(theProblemManager,valueType)
            
            stateMineValuesHalfGrade[state].append(halfGradeValue)
  
          ### double grade ###
  
          theProblemManager.theMineDataManager.theOreBody.ScaleCommodityGrades(4)
  
          for state in self.stateStrs:
            theProblemManager.theMineDataManager.theEconomicDataManager.CalculateBeforeTaxCashFlow(theProblemManager, theProblemManager.theMineDataManager)
            theProblemManager.theMineDataManager.theEconomicDataManager.SetState(state)
            theProblemManager.theMineDataManager.CalculateTaxes(theProblemManager)
            theProblemManager.theMineDataManager.CalculateAfterTaxCashFlow(theProblemManager)
            doubleGradeValue = theProblemManager.theMineDataManager.CalculateEconomicIndicator(theProblemManager,valueType)
            
            stateMineValuesDoubleGrade[state].append(doubleGradeValue)
  
          # revert to original grade
  
          theProblemManager.theMineDataManager.theOreBody.ScaleCommodityGrades(0.5/gradeAdjustment)
  
        if doTimer:
          timesList.append(["coverFunc", time.time()])


        #mineValues = np.array(mineValues,dtype=np.float)
        #mineTypes = np.array(mineTypes,dtype=np.float)
        #mineYears = np.array(mineYears,dtype=np.float)
        concentrateDiscountedTotalMasses = np.array(concentrateDiscountedTotalMasses,dtype=np.float) # needs to be array - as we double and halve concentration
        concentrateCapacities = np.array(concentrateCapacities,dtype=np.float)  # needs to be array


        mineValueFunc = interpolate.interp1d(coverDepths,mineValues)
        mineTypeFunc =interpolate.interp1d(coverDepths,mineTypes)
        mineLifeFunc = interpolate.interp1d(coverDepths,mineYears)
        concentrateDiscountedTotalMassFunc = interpolate.interp1d(coverDepths,concentrateDiscountedTotalMasses)
        concentrateCapacityFunc = interpolate.interp1d(coverDepths,concentrateCapacities)
        halfConcentrateDiscountedTotalMassFunc = interpolate.interp1d(coverDepths,0.5*concentrateDiscountedTotalMasses)
        halfConcentrateCapacityFunc = interpolate.interp1d(coverDepths,0.5*concentrateCapacities)
        doubleConcentrateDiscountedTotalMassFunc = interpolate.interp1d(coverDepths,2*concentrateDiscountedTotalMasses)
        doubleConcentrateCapacityFunc = interpolate.interp1d(coverDepths,2*concentrateCapacities)

        if(doPlots):
          for state in self.stateStrs:
            pl.plot(coverDepths,stateMineValues[state],label=state)
  
          pl.legend()
          pl.plot(coverDepths,mineValues)


        ## regional mine cost calculation

        if doTimer:
          timesList.append(["startMaps", time.time()])

        stateIdsMap =  self.LoadMap(self.stateIdsMapFile)
        self.ApplyOutputBoundingBoxMask(stateIdsMap) 

        countryIndxs = stateIdsMap >= 0.0

        ###################################
        
        coverMap =  self.LoadMap(self.coverDepthMapFile,self.coverOffset,Min=self.coverMin,Max=self.coverMax)

        mineValueMap = np.zeros(coverMap.shape)
        halfGradeMineValueMap = np.zeros(coverMap.shape)
        doubleGradeMineValueMap = np.zeros(coverMap.shape)

        taxRate = self.GetEffectiveTaxRate(theProblemManager.theMineDataManager,valueType)
        taxRelief = 1.0 - taxRate  #  tax relief from incometax

        ####
        
        for i in range(len(self.stateStrs)):
          stateIndxs = stateIdsMap == self.stateIds[i]
          mineValueFunc = interpolate.interp1d(coverDepths,stateMineValues[self.stateStrs[i]])
          mineValueMap[stateIndxs] = mineValueFunc(coverMap[stateIndxs])
          halfGradeMineValueFunc = interpolate.interp1d(coverDepths,stateMineValuesHalfGrade[self.stateStrs[i]])
          halfGradeMineValueMap[stateIndxs] = halfGradeMineValueFunc(coverMap[stateIndxs])
          doubleGradeMineValueFunc = interpolate.interp1d(coverDepths,stateMineValuesDoubleGrade[self.stateStrs[i]])
          doubleGradeMineValueMap[stateIndxs] = doubleGradeMineValueFunc(coverMap[stateIndxs])

        if doTimer:
          timesList.append(["mineValue", time.time()])

        # Water Expenses
        distanceToWater =  self.LoadMap(self.waterDistanceMapFile)
        distanceToWaterCost =  taxRelief* theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalWaterExpenses(theProblemManager,distanceToWater[countryIndxs])
        mineValueMap[countryIndxs] -=  distanceToWaterCost
        halfGradeMineValueMap[countryIndxs] -= distanceToWaterCost
        doubleGradeMineValueMap[countryIndxs] -=  distanceToWaterCost
        del distanceToWater
        del distanceToWaterCost

        if doTimer:
          timesList.append(["water", time.time()])

        # Power Expenses
        distanceToPower =  self.LoadMap(self.powerDistanceMapFile)
        distanceToPowerCost = taxRelief* theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalPowerExpenses(theProblemManager,distanceToPower[countryIndxs])
        mineValueMap[countryIndxs] -=  distanceToPowerCost
        halfGradeMineValueMap[countryIndxs] -= distanceToPowerCost
        doubleGradeMineValueMap[countryIndxs] -=  distanceToPowerCost
        del distanceToPower
        del distanceToPowerCost


        if doTimer:
          timesList.append(["power", time.time()])

        # Transportation Expenses
        distanceToRoad =  self.LoadMap(self.roadDistanceMapFile)
        distanceToRail =  self.LoadMap(self.railDistanceMapFile)
        roadTransportationDistance =  self.LoadMap(self.roadDistanceToPortMapFile)
        railTransporationDistance =  self.LoadMap(self.railDistanceToPortMapFile)


        mineValueMap[countryIndxs] -= taxRelief* theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalTransportationExpenses( \
                                                              distanceToRoad[countryIndxs], roadTransportationDistance[countryIndxs], \
                                                              distanceToRail[countryIndxs], railTransporationDistance[countryIndxs], \
                                                              coverMap[countryIndxs],concentrateCapacityFunc,concentrateDiscountedTotalMassFunc)

        halfGradeMineValueMap[countryIndxs] -= taxRelief* theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalTransportationExpenses( \
                                                              distanceToRoad[countryIndxs], roadTransportationDistance[countryIndxs], \
                                                              distanceToRail[countryIndxs], railTransporationDistance[countryIndxs], \
                                                              coverMap[countryIndxs],halfConcentrateCapacityFunc,halfConcentrateDiscountedTotalMassFunc)

        doubleGradeMineValueMap[countryIndxs] -= taxRelief* theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalTransportationExpenses( \
                                                              distanceToRoad[countryIndxs], roadTransportationDistance[countryIndxs], \
                                                              distanceToRail[countryIndxs], railTransporationDistance[countryIndxs], \
                                                              coverMap[countryIndxs],doubleConcentrateCapacityFunc,doubleConcentrateDiscountedTotalMassFunc)

        if doTimer:
          timesList.append(["transport", time.time()])



        mineValueMap[stateIdsMap< 0.0] = np.nan  

        # interpolate to find the grade scale that gives value = 0

        c = mineValueMap[countryIndxs]
        a = ( doubleGradeMineValueMap[countryIndxs]+ 2*halfGradeMineValueMap[countryIndxs]-3*c )/1.5
        b = doubleGradeMineValueMap[countryIndxs] - a - c
        # from matplotlib import pyplot as plt
        # aaa = ( doubleGradeMineValueMap+ 2*halfGradeMineValueMap-3*mineValueMap )/1.5
        # bbb = (doubleGradeMineValueMap - aaa - mineValueMap)
        # plt.subplot(131)
        # plt.imshow(np.sqrt( bbb*bbb-4*aaa*mineValueMap )[::-1])
        # plt.colorbar()
        # plt.subplot(132)
        # ttt = np.sqrt( bbb*bbb-4*aaa*mineValueMap )
        # ggg = -np.sqrt( 4*aaa*mineValueMap-bbb*bbb )
        # eee = np.zeros(coverMap.shape)
        # eee[np.isnan(ttt)==False] = ttt[np.isnan(ttt)==False]
        # eee[np.isnan(ggg)==False] = ggg[np.isnan(ggg)==False]
        # plt.imshow( eee[::-1] )
        # plt.colorbar()
        # plt.subplot(133)
        # plt.imshow((4*aaa*mineValueMap)[::-1])
        # plt.colorbar()
        # plt.show()

        #
        breakEvenFactor = np.zeros(coverMap.shape)

        # r1 = np.sqrt( b*b-4*a*c )
        # r2 = 1./(-np.sqrt( 4*a*c-b*b ))
        # rrr = np.zeros(coverMap.shape)[countryIndxs]
        # rrr[np.isnan(r1)==False] = r1[np.isnan(r1)==False]
        # rrr[np.isnan(r2)==False] = r2[np.isnan(r2)==False]
        breakEvenFactor[countryIndxs] = 1+(-b + np.sqrt( b*b-4*a*c ) )/(2*a+1e-64)#1+(-b + rrr )/(2*a+1e-64)#
        # breakEvenFactor[countryIndxs][np.isnan(r2)==False] = r2[np.isnan(r2)==False]

        breakEvenFactor[stateIdsMap< 0.0] = np.nan  
        
        # correct for break even factor
        if(gradeAdjustment != 1.0):
          #print "gradeAdjustment",gradeAdjustment
          breakEvenFactor *= gradeAdjustment
        
        
        # clip to bounding box
        if(len(self.dataBoundingBox) == 4 and len(self.outputBoundingBox) == 4):
          breakEvenFactor = self.ClipToOutputBounds(breakEvenFactor)

        if doTimer:
          timesList.append(["end", time.time()])

        ######################################################################



        if(doPlots):
        
          pl.figure()
          pl.plot(breakEvenFactor[500,:])
          
          pl.figure()
          pl.imshow(breakEvenFactor,origin="lower")
          pl.colorbar()
          pl.show()    
        else:
          filename = theProblemManager.outputPrefix+"."+theProblemManager.outputType
          Commodity = list(theProblemManager.theMineDataManager.theOreBody.metalGrades.keys())[0]
          Grade = theProblemManager.theMineDataManager.theOreBody.metalGrades[Commodity]
          if Commodity in ['Cu','Zn','Pb','Ni','P2O5','K2SO4']:
            Grade *= 100.
          elif Commodity in ['Co','PGE','REO','Au','Ag']:
            theUnitManager =  UnitManager()
            gramsPerTonne  = theUnitManager.ConvertToBaseUnits("g/t")
            Grade /= gramsPerTonne
          breakEvenGrade = Grade * breakEvenFactor
          self.SaveMap(filename, breakEvenGrade,theProblemManager.recordRange)

        ######################################################################

        if doTimer:
            pl.figure()
            dts = np.array( [tt[1] for tt in timesList] )
            dts -=dts[0]
            labels =  [tt[0] for tt in timesList] 

            pl.bar(range(len(dts)),dts,align="center")
            pl.xticks(range(len(dts)), labels)


            pl.bar(range(1,len(dts)),dts[1:]-dts[:-1],color="r",align="center")

            for i in range(1,len(dts)):
              print(labels[i], dts[i], dts[i]-dts[i-1])

    
        return breakEvenFactor

        ######################################################################

    def CalculateNetRealCosts(self,theProblemManager):
      """Calculates the real (undiscounted) amount spent on the mining project over its lifetime. Used to determine indirect from IO Matrix. 
      
       When Net Real Cost is passed to the NPV calculation, the regional calculation manager sets the discount rate, federal and state taxes and revenues to 0. 
       
       """
      
      # store current values 
      original_saveResult = self.saveResult
      original_outputType = theProblemManager.outputType
      
      # Replace values 
      self.saveResult = False
      theProblemManager.outputType = ""
      
      # run NPV calculation
      realNetExpediture = self.RunNPVcalculation(theProblemManager, "NetRealCost")
      
      self.saveResult  = original_saveResult
      theProblemManager.outputType  = original_outputType
      
      return realNetExpediture


    #########################
    ## Uncertain NPV calculation
    
    def RunUncertainNPVcalculation(self,theProblemManager,valueType="atNPV"):
        """Calculate the regional Net Present Value and associated uncertainty in mining projects"""
    
        doPlots = (theProblemManager.outputType == "plots")

        doAll = (theProblemManager.outputType == "all")
        
        # profiling
        doTimer = True
        timesList = []

        if doTimer:
          timesList.append(["start", time.time()])
        
          
        theProblemManager.theMineDataManager.theInfrastructureManager.ZeroDistanceToInfrastructure()

        # calculate mine cost, life etc as function of depth of cover (in meters)
        
        coverDepths = self.GetCoverDepths()
  
        discountRate = self.GetEffectiveDiscountRate(theProblemManager.theMineDataManager,valueType)

        mineValues = []
        mineTypes =[]
        mineYears = []
        concentrateDiscountedTotalMasses = []
        concentrateCapacities = []
        
        mineCapExs = []  # for employment calculation
        stateAndFederalTaxes = {}  # for tax calculation

        stateMineValues = {}
        
        
        
        for state in self.stateStrs:
          stateMineValues[state] = []
          if(self.recordTaxationCalculation):
            stateAndFederalTaxes[state] = []
        
        theFunctionManager = FunctionManager()
        
        for cc in range(self.numRealizations):
            print("realization: ", cc)
            
            # switch on uncertainty in function evaluations
            theFunctionManager.PerturbRandomFunctions()
            theProblemManager.theMineDataManager.theEconomicDataManager.GenerateRandomWalkCommodityPrices(200)
            for depthOfCover in coverDepths:
        
              theProblemManager.theMineDataManager.SetCoverDepth(depthOfCover)
  
              value = theProblemManager.EvaluateOCUGMines()
              mineValues.append(value)
  
              mineLife = theProblemManager.theMineDataManager.theMiningSystem.mineLife
              mineYears.append(mineLife)
          
              if(self.doEmploymentCalculation):
                mineAndProcCapex = theProblemManager.theMineDataManager.theEconomicDataManager.netCapex[0]
                mineCapExs.append(mineAndProcCapex)
            
  
              mineTypeStr = theProblemManager.theMineDataManager.theMiningSystem.miningMethod  # 1 if UG, 0.5 mixed, 0 OC 
  
              mineType = 0.0   # open cut or K
              if(mineTypeStr[:2] == "UG"):
                mineType = 1.0
              elif(mineTypeStr == "OCUG"):
                mineType = 0.5
    
              for state in self.stateStrs:
                theProblemManager.theMineDataManager.theEconomicDataManager.SetState(state)
                theProblemManager.theMineDataManager.CalculateTaxes(theProblemManager)
                theProblemManager.theMineDataManager.CalculateAfterTaxCashFlow(theProblemManager)
                
                value = theProblemManager.theMineDataManager.CalculateEconomicIndicator(theProblemManager,valueType)
            
                stateMineValues[state].append(value)
            
    
  
              mineTypes.append(mineType)
  
              # discounted total mass
              discountFactors = (1+theProblemManager.theMineDataManager.theEconomicDataManager.discountRate)**np.array(range( int( np.ceil(mineLife) ) ) )
              concentrateMass = np.sum(theProblemManager.theMineDataManager.theProcessingSystem.GetTotalConcentrateProduced()/discountFactors)
              concentrateDiscountedTotalMasses.append(concentrateMass)
  
              concentrateCapacity = np.max(theProblemManager.theMineDataManager.theProcessingSystem.GetTotalConcentrateProduced())
              concentrateCapacities.append(concentrateCapacity)
        
        # deactivate random sampler   
        theFunctionManager.ZeroPerturbations()
        theProblemManager.theMineDataManager.theEconomicDataManager.ResetCommodityPrices()
        
        
        # record realization stats
        #np.savetxt("mineValues.txt",mineValues)
        #np.savetxt("concentrateDiscountedTotalMasses.txt",concentrateDiscountedTotalMasses)
        #np.savetxt("concentrateCapacities.txt",concentrateCapacities)
        
        def convertToBetaDistributions(values):
          values = np.array(values,dtype=np.float).reshape(self.numRealizations,len(coverDepths))
          meanValues = np.mean(values,0)
          stdValues = np.std(values,0)
          minValues = np.percentile(values,1,0)
          maxValues = np.percentile(values,99,0)
          
          rv = ScaledBetaDistribution.FromMeanAndStd(meanValues,stdValues,minValues,maxValues)
          
          return rv
        
        
        concentrateDiscountedTotalMasses = convertToBetaDistributions(concentrateDiscountedTotalMasses)
        
        concentrateCapacities = convertToBetaDistributions(concentrateCapacities)
        
        
        for state in self.stateStrs:
            stateMineValues[state] = convertToBetaDistributions(stateMineValues[state])
        

        if doTimer:
          timesList.append(["coverFunc", time.time()])
          

        concentrateDiscountedTotalMassFunc = UncertainInterp1d(coverDepths,concentrateDiscountedTotalMasses)
        concentrateCapacityFunc = UncertainInterp1d(coverDepths,concentrateCapacities)
        
        
        if doTimer:
          timesList.append(["startMaps", time.time()])
        
        stateIdsMap =  self.LoadMap(self.stateIdsMapFile)
        self.ApplyOutputBoundingBoxMask(stateIdsMap) 
        


        ###################################
        
        if(self.coverDepthStdMapFile and self.coverDepthMapFile):
          coverMap = self.LoadMapWUncertainty(self.coverDepthMapFile,self.coverDepthStdMapFile,offset=self.coverOffset,Min=5e-6,Max=self.coverMax)
        else:
          coverMap =  self.LoadMap(self.coverDepthMapFile,offset=self.coverOffset,Min=self.coverMin,Max=self.coverMax)
        
        countryIndxs = stateIdsMap >= 0.0
        zeroDistribution = ScaledBetaDistribution(1.0,1.0,0.0,1e-64)
        mineValueMap = np.full(coverMap.shape,zeroDistribution)
        
        taxRate = self.GetEffectiveTaxRate(theProblemManager.theMineDataManager,valueType)
        taxRelief = 1.0 - taxRate  # = 30% tax relief from federal income tax (for costs)
        
        ####
        for i in np.unique(stateIdsMap.astype(int)):
          if i < 0: pass
          else:
            stateIndxs = stateIdsMap == self.stateIds[i]
            mineValueFunc = UncertainInterp1d(coverDepths,stateMineValues[self.stateStrs[i]])
            mineValueMap[stateIndxs] = mineValueFunc(coverMap[stateIndxs])
        

        if doTimer:
          timesList.append(["mineValue", time.time()])
          

        # Water Expenses
        # BACCARINI, D. & LOVE, P. 2014. Statistical Characteristics of Cost Contingency in Water Infrastructure Projects. Journal of Construction Engineering and Management, 140.
        waterCostOverruns = ScaledBetaDistribution.FromMeanAndStd(1.0512 ,0.2595 ,min=0.5608,max= 2.654) #,min=-0.4392,max= 1.654)
        
        distanceToWater =  self.LoadWaterDistance(theProblemManager)
        mineValueMap[countryIndxs] -=  (waterCostOverruns*taxRelief) * theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalWaterExpenses(theProblemManager,distanceToWater[countryIndxs])
        
        # attempt to trigger garbage cleanup
        distanceToWater = None
        del distanceToWater
        

        if doTimer:
          timesList.append(["water", time.time()])
        
        
        # transmission cost overrun - based on data from Sovacool et al "Construction Cost Overruns and Electricity Infrastructure: An Unavoidable Risk?", The Electricity Journal 2014
        # ignoring outlier
        transmean,transstd,transmin,transmax = [1.0539, 0.1971, 0.6666, 1.7143]
        transmissionCostOverruns = ScaledBetaDistribution.FromMeanAndStd(transmean,transstd,transmin,transmax)
        
        # Power Expenses
        distanceToPower =  self.LoadInfrastructureMap(theProblemManager,"powerlines",self.powerDistanceMapFile)
        
        mineValueMap[countryIndxs] -=  (transmissionCostOverruns*taxRelief) * theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalPowerExpenses(theProblemManager,distanceToPower[countryIndxs])
        distanceToPower =  None
        del distanceToPower


        if doTimer:
          timesList.append(["power", time.time()])
        

        # Transportation Expenses
       
        nLat, nLon = coverMap.shape 
        distanceToRoad,roadTransportationDistance =  self.LoadInfrastructureDistanceAndTravelDistanceMaps(theProblemManager,"roads",self.roadDistanceMapFile,self.roadDistanceToPortMapFile,nLat, nLon)
        distanceToRail,railTransporationDistance = self.LoadInfrastructureDistanceAndTravelDistanceMaps(theProblemManager,"raillines",self.railDistanceMapFile,self.railDistanceToPortMapFile,nLat, nLon)
        


        mineValueMap[countryIndxs] -= taxRelief* theProblemManager.theMineDataManager.theInfrastructureManager.CalculateUncertainRegionalTransportationExpenses( \
                                                              distanceToRoad[countryIndxs], roadTransportationDistance[countryIndxs], \
                                                              distanceToRail[countryIndxs], railTransporationDistance[countryIndxs], \
                                                              coverMap[countryIndxs],concentrateCapacityFunc,concentrateDiscountedTotalMassFunc)
        
        # free memory
        distanceToRoad = []
        roadTransportationDistance = []
        distanceToRail = []
        railTransporationDistance = []
        del distanceToRoad
        del roadTransportationDistance
        del distanceToRail
        del railTransporationDistance
        
        
        if doTimer:
          timesList.append(["end", time.time()])

        ######################################################################
        
        if(doPlots):
            print("Done - plotting now")
            
            expectations = GetExpectations(mineValueMap)
            vars = GetVariances(mineValueMap)
            stds = vars**0.5
            
            pl.figure()
            pl.imshow(expectations,origin="lower")
            
            pl.colorbar()
        
            pl.figure()
        
            pl.imshow(stds,origin="lower",clim=[0,1e9])
            print("max std:", np.max(stds))
            pl.colorbar()
        
        ######################################################################

        if doTimer:
            pl.figure()
            dts = np.array( [tt[1] for tt in timesList] )
            dts -=dts[0]
            labels =  [tt[0] for tt in timesList] 

            pl.bar(range(len(dts)),dts,align="center")
            pl.xticks(range(len(dts)), labels)


            pl.bar(range(1,len(dts)),dts[1:]-dts[:-1],color="r",align="center")

            for i in range(1,len(dts)):
              print(labels[i], dts[i], dts[i]-dts[i-1])
        
        if(doPlots):
          pl.show()    
        else:
          E = GetExpectations(mineValueMap)
          E[countryIndxs==False] = np.nan
          V = np.sqrt(GetVariances(mineValueMap))  # V = std dev (not variance)
          V[countryIndxs==False] = np.nan
          # clip to bounding box
          if(len(self.dataBoundingBox) == 4 and len(self.outputBoundingBox) == 4):
            E = self.ClipToOutputBounds(E)
            V = self.ClipToOutputBounds(V)
          filename = theProblemManager.outputPrefix+"_mean."+theProblemManager.outputType
          if(doAll):
            filename = theProblemManager.outputPrefix+"_mean.npy"
          self.SaveMap(filename, E, theProblemManager.recordRange)
          self.type = "tax"
          filename = theProblemManager.outputPrefix+"_std."+theProblemManager.outputType
          self.SaveMap(filename, V, theProblemManager.recordRange)
          self.type = "uncertainNPV"
        
        return mineValueMap


    ######################################################################
    # Hydrogen 
    ######################################################################

    #########################
    ## GenerateHybridCapacityFactorMap

    def GenerateHybridCapacityFactorMap(self,theProblemManager,thePowerPlant):
      """Combine two or more hydrogen capacity factor maps to produce a capacity factor map for the whole plant."""
      
      if(thePowerPlant.useHybridPlantWithCurtailment):  # ugly
        lowerBound = 0.1* np.floor(thePowerPlant.powerFraction*10)  # to lower fraction
        upperBound = lowerBound +0.1
        
        lowerCapacityFactorMapFile = self.hybridCapacityFactorPrefix+"_"+str(lowerBound)[:3] + ".npy"
        upperCapacityFactorMapFile = self.hybridCapacityFactorPrefix+"_"+str(upperBound)[:3] + ".npy"
        
        lw = (upperBound - thePowerPlant.powerFraction)/(upperBound-lowerBound)   # lower weight
        
        rv = lw *np.load(lowerCapacityFactorMapFile) + (1-lw)*np.load(upperCapacityFactorMapFile)
        
      else:
          # original method - a straight weighted average of the two sources for hybrid plants
          if(thePowerPlant.type == "Photovoltaic"):
              capacityFactorMapFile = self.pvCapacityFactorFile
          elif(thePowerPlant.type == "CSP"):
              capacityFactorMapFile = self.cspCapacityFactorFile
          elif(thePowerPlant.type == "Wind"):
              capacityFactorMapFile = self.windCapacityFactorFile
      
          rv = thePowerPlant.powerFraction * np.load(capacityFactorMapFile)
          
          #print("thePowerPlant.powerFraction",thePowerPlant.powerFraction)
      
          if thePowerPlant.secondaryPowerSource:
            rv += self.GenerateHybridCapacityFactorMap(theProblemManager,thePowerPlant.secondaryPowerSource)
      
      return rv


    #########################
    ## Hydrogen Regional NPV calculation
    
    def RunHydrogenNPVcalculation(self,theProblemManager,valueType="atNPV"):
        """Calculate the regional Net Present Value (NPV) for hydrogen projects"""
        
        # Hydrogen storage - alternative point of sale
        if(self.doHydrogenStorageCalculation):
          NPVwStorage = self.RunHydrogenStorageNPVcalculation(theProblemManager)
          
    
        doPlots = (theProblemManager.outputType == "plots")

        doAll = (theProblemManager.outputType == "all")
        
        # profiling
        doTimer = False
        timesList = []
        
        if doTimer:
          timesList.append(["start", time.time()])
        
        theProblemManager.theHydrogenDataManager.theInfrastructureManager.ZeroDistanceToInfrastructure()

        # calculate project cost as a function of the capacity factor
        # nb: will get better results if we can bound this based on map data. 
        capacityFactorMin = 0.01
        capacityFactorMax = 0.99
        originalCapacityFactor = theProblemManager.theHydrogenDataManager.GetCapacityFactor() # to restore

        capacityFactors = np.linspace(capacityFactorMin,capacityFactorMax,50)
  
        discountRate = self.GetEffectiveDiscountRate(theProblemManager.theHydrogenDataManager,valueType)

        projectValues = []
        projectYears = []
        hydrogenDiscountedTotalMasses = []
        powerPlantCapacities = []
        hydrogenPlantCapacity = theProblemManager.theHydrogenDataManager.theHydrogenPlant.hydrogenProductionCapacity
        
        projectCapExs = []  # for employment calculation
        stateAndFederalTaxes = {}  # for tax calculation
        waterUse = []  # for water use calculation

        stateProjectValues = {}
        
        for state in self.stateStrs:
          stateProjectValues[state] = []
          if(self.recordTaxationCalculation):
            stateAndFederalTaxes[state] = []

        for capacityFactor in capacityFactors:
        
          theProblemManager.theHydrogenDataManager.SetCapacityFactor(capacityFactor)
  
          value = theProblemManager.EvaluateHydrogenProject()
          projectValues.append(value)
  
          projectLife = theProblemManager.theHydrogenDataManager.theHydrogenPlant.projectLife
          projectYears.append(projectLife)
          
          if(self.doEmploymentCalculation):
            capex = theProblemManager.theHydrogenDataManager.theEconomicDataManager.netCapex[0]
            projectCapExs.append(capex)
          
          if(self.recordWaterUse):
            waterUsePerYear = theProblemManager.theHydrogenDataManager.theInfrastructureManager.CalculateWaterUsePerYear(theProblemManager)
            waterUse.append(waterUsePerYear)
          
            
          for state in self.stateStrs:
            theProblemManager.theHydrogenDataManager.theEconomicDataManager.SetState(state)
            theProblemManager.theHydrogenDataManager.CalculateTaxes(theProblemManager)
            theProblemManager.theHydrogenDataManager.CalculateAfterTaxCashFlow(theProblemManager)
            
            value = theProblemManager.theHydrogenDataManager.CalculateEconomicIndicator(theProblemManager,valueType)
            stateProjectValues[state].append(value)
            
            if(self.recordTaxationCalculation):
              taxNPV = np.sum(theProblemManager.theHydrogenDataManager.theEconomicDataManager.btNCF - theProblemManager.theHydrogenDataManager.theEconomicDataManager.atNCF)
              #taxNPV = theProblemManager.theHydrogenDataManager.theEconomicDataManager.btNPV - theProblemManager.theHydrogenDataManager.theEconomicDataManager.atNPV 
              stateAndFederalTaxes[state].append(taxNPV)
    
  
  
          # discounted total mass
          discountFactors = (1+theProblemManager.theHydrogenDataManager.theEconomicDataManager.discountRate)**np.array(range( int( np.ceil(projectLife) ) ) )
          hydrogenDiscountedMass = np.sum(theProblemManager.theHydrogenDataManager.theHydrogenPlant.hydrogenProduced/discountFactors)
          hydrogenDiscountedTotalMasses.append(hydrogenDiscountedMass)
  
          powerPlantCapacity = theProblemManager.theHydrogenDataManager.thePowerPlant.plantCapacity
          powerPlantCapacities.append(powerPlantCapacity)
          

        # restore capacity factor
        theProblemManager.theHydrogenDataManager.SetCapacityFactor(originalCapacityFactor)
  
        if doTimer:
          timesList.append(["coverFunc", time.time()])
          
        projectValues = np.array(projectValues,dtype=np.float)
        projectYears = np.array(projectYears,dtype=np.float)
        
        hydrogenDiscountedTotalMasses = np.array(hydrogenDiscountedTotalMasses,dtype=np.float)
        powerPlantCapacities = np.array(powerPlantCapacities,dtype=np.float)


        projectValueFunc = interpolate.interp1d(capacityFactors,projectValues)
        projectLifeFunc = interpolate.interp1d(capacityFactors,projectYears)
        hydrogenDiscountedTotalMassFunc = interpolate.interp1d(capacityFactors,hydrogenDiscountedTotalMasses)
        powerPlantCapacityFunc = interpolate.interp1d(capacityFactors,powerPlantCapacities)
        
        if(self.doEmploymentCalculation):
          projectcapexFunc = interpolate.interp1d(capacityFactors,projectCapExs)
        
            
        if(doPlots):

            for state in self.stateStrs:
              pl.plot(capacityFactors,stateProjectValues[state],label=state)
          
            pl.legend()
            pl.plot(capacityFactors,projectValues)

        
        ## regional project cost calculation
        connectToGas = False
        connectToCO2 = False
        calculateCoalTransport = False
        
        if(theProblemManager.theHydrogenDataManager.thePowerPlant.type == "Photovoltaic"):
          capacityFactorMapFile = self.pvCapacityFactorFile
        elif(theProblemManager.theHydrogenDataManager.thePowerPlant.type == "CSP"):
          capacityFactorMapFile = self.cspCapacityFactorFile
        elif(theProblemManager.theHydrogenDataManager.thePowerPlant.type == "Wind"):
          capacityFactorMapFile = self.windCapacityFactorFile
        elif(theProblemManager.theHydrogenDataManager.thePowerPlant.type == "SteamMethane"):
          capacityFactorMapFile = self.steamMethaneCapacityFactorFile
          connectToGas = True
          connectToCO2 = True
        elif(theProblemManager.theHydrogenDataManager.thePowerPlant.type == "BrownCoalGasification"):
          capacityFactorMapFile = self.brownCoalCapacityFactorFile
          roadDistanceToCoalMapFile = self.brownCoalTransportationDistanceMapFile
          calculateCoalTransport = True
          connectToCO2 = True
        elif(theProblemManager.theHydrogenDataManager.thePowerPlant.type == "BlackCoalGasification"):
          capacityFactorMapFile = self.blackCoalCapacityFactorFile
          roadDistanceToCoalMapFile = self.blackCoalTransportationDistanceMapFile
          calculateCoalTransport = True
          connectToCO2 = True
        else:
        
          raise BluecapError("Error Regional Calculation: Unrecognized power plant type: " + theProblemManager.theHydrogenDataManager.thePowerPlant.type)
        
        
        if doTimer:
          timesList.append(["startMaps", time.time()])

        stateIdsMap = self.LoadMap(self.stateIdsMapFile)
        self.ApplyOutputBoundingBoxMask(stateIdsMap) 
        
        
        ###################################
        
        if(theProblemManager.theHydrogenDataManager.thePowerPlant.secondaryPowerSource):
          capacityFactorMap = self.GenerateHybridCapacityFactorMap(theProblemManager,theProblemManager.theHydrogenDataManager.thePowerPlant)
        else:
          capacityFactorMap =  self.LoadMap(capacityFactorMapFile)
          
        # bound capacityFactor map
        capacityFactorMap[capacityFactorMap > capacityFactorMax] = capacityFactorMax
        capacityFactorMap[capacityFactorMap < capacityFactorMin] = capacityFactorMin
          
          
          
        countryIndxs = stateIdsMap >= 0.0

        hydrogenValueMap = np.zeros(capacityFactorMap.shape)
        
        taxRate = self.GetEffectiveTaxRate(theProblemManager.theHydrogenDataManager,valueType)
        taxRelief = 1.0 - taxRate  # = 30% tax relief from federal income tax (for costs)
        
        
        if(self.recordTaxationCalculation):
            taxEstimate = np.zeros(capacityFactorMap.shape)
          

        #### Run project value calculations for each state
        for i in range(len(self.stateStrs)):
          stateIndxs = stateIdsMap == self.stateIds[i]
          projectValueFunc = interpolate.interp1d(capacityFactors,stateProjectValues[self.stateStrs[i]])
          hydrogenValueMap[stateIndxs] = projectValueFunc(capacityFactorMap[stateIndxs])

          if(self.recordTaxationCalculation):
            stateAndFederalTaxFunc = interpolate.interp1d(capacityFactors,stateAndFederalTaxes[self.stateStrs[i]])
            taxEstimate[stateIndxs] = stateAndFederalTaxFunc(capacityFactorMap[stateIndxs])
        
        
        if(doAll):
          filename = theProblemManager.outputPrefix+"_projectValue.npy"  
          self.SaveMap(filename, hydrogenValueMap,False)

        if doTimer:
          timesList.append(["projectValue", time.time()])

        # Water Expenses
        distanceToWater =  self.LoadWaterDistance(theProblemManager)
        hydrogenValueMap[countryIndxs] -=  taxRelief* theProblemManager.theHydrogenDataManager.theInfrastructureManager.CalculateRegionalWaterExpenses(theProblemManager,distanceToWater[countryIndxs])
        
        if(self.recordTaxationCalculation):
          taxEstimate[countryIndxs] -= taxRate*theProblemManager.theHydrogenDataManager.theInfrastructureManager.CalculateRegionalWaterExpenses(theProblemManager,distanceToWater[countryIndxs])


        if(doAll):
          filename = theProblemManager.outputPrefix+"_waterCost.npy"  
          dd = np.zeros(capacityFactorMap.shape)
          dd[countryIndxs] = taxRelief* theProblemManager.theHydrogenDataManager.theInfrastructureManager.CalculateRegionalWaterExpenses(theProblemManager,distanceToWater[countryIndxs])
          #print "Saving: ", filename
          self.SaveMap(filename, dd,False)

        if doTimer:
          timesList.append(["water", time.time()])

        if(self.connectToPowerGrid):
            # Power Expenses
            distanceToPower =  self.LoadInfrastructureMap(theProblemManager,"powerlines",self.powerDistanceMapFile)

            hydrogenValueMap[countryIndxs] -=  taxRelief* theProblemManager.theHydrogenDataManager.theInfrastructureManager.CalculateRegionalPowerExpenses(theProblemManager,distanceToPower[countryIndxs])

            if(self.recordTaxationCalculation):
              taxEstimate[countryIndxs]  -= taxRate*theProblemManager.theHydrogenDataManager.theInfrastructureManager.CalculateRegionalPowerExpenses(theProblemManager,distanceToPower[countryIndxs])
        
        if(connectToGas):
          distanceToGas =  self.LoadInfrastructureMap(theProblemManager,"gaslines",self.gasDistanceMapFile)
          gasCost =  theProblemManager.theHydrogenDataManager.theInfrastructureManager.CalculateRegionalGasExpenses(theProblemManager,distanceToGas[countryIndxs])
          hydrogenValueMap[countryIndxs] -=  taxRelief* gasCost
          
          if(self.recordTaxationCalculation):
              taxEstimate[countryIndxs]  -= taxRate*gasCost

        if(connectToCO2):
          distanceToCO2 =  self.LoadMap(self.co2DistanceMapFile)
          flowRate = theProblemManager.theHydrogenDataManager.theHydrogenPlant.hydrogenProductionCapacity # annual production in kg. 
          flowRate *= theProblemManager.theHydrogenDataManager.theHydrogenPlant.co2PerKgH2
            
          CO2Cost =  theProblemManager.theHydrogenDataManager.theInfrastructureManager.CalculateRegionalCO2Expenses(theProblemManager,distanceToCO2[countryIndxs], flowRate)
          hydrogenValueMap[countryIndxs] -=  taxRelief* CO2Cost
          
          if(self.recordTaxationCalculation):
              taxEstimate[countryIndxs]  -= taxRate*CO2Cost


        if(calculateCoalTransport):
          coalTransportationDistance =  self.LoadMap(roadDistanceToCoalMapFile)
            
          coalTransportCost =  theProblemManager.theHydrogenDataManager.theInfrastructureManager.CalculateRegionalCoalTransportationExpenses(theProblemManager, coalTransportationDistance[countryIndxs],\
                                                                                                                                            theProblemManager.theHydrogenDataManager.theHydrogenPlant.coalPerKgH2,\
                                                                                                                                            capacityFactorMap[countryIndxs],hydrogenPlantCapacity,hydrogenDiscountedTotalMassFunc)
          hydrogenValueMap[countryIndxs] -=  taxRelief* coalTransportCost
          
          if(self.recordTaxationCalculation):
              taxEstimate[countryIndxs]  -= taxRate*coalTransportCost


        if(doAll):
          filename = theProblemManager.outputPrefix+"_powerCost.npy"  
          dd = np.zeros(capacityFactorMap.shape)
          dd[countryIndxs] = taxRelief* theProblemManager.theHydrogenDataManager.theInfrastructureManager.CalculateRegionalPowerExpenses(theProblemManager,distanceToPower[countryIndxs])
          
          self.SaveMap(filename, dd,False)
          
          if(theProblemManager.theHydrogenDataManager.theInfrastructureManager.calculateGas):
              filename = theProblemManager.outputPrefix+"_gasCost.npy"  
              dd = np.zeros(capacityFactorMap.shape)
              dd[countryIndxs] = taxRelief* theProblemManager.theHydrogenDataManager.theInfrastructureManager.CalculateRegionalGasExpenses(theProblemManager,distanceToGas[countryIndxs])
              
              self.SaveMap(filename, dd,False)

        if doTimer:
          timesList.append(["power", time.time()])

        # Transportation Expenses
        nLat, nLon = capacityFactorMap.shape 
        distanceToRoad,roadTransportationDistance =  self.LoadInfrastructureDistanceAndTravelDistanceMaps(theProblemManager,"roads",self.roadDistanceMapFile,self.roadDistanceToPortMapFile,nLat, nLon)
        distanceToRail,railTransporationDistance = self.LoadInfrastructureDistanceAndTravelDistanceMaps(theProblemManager,"raillines",self.railDistanceMapFile,self.railDistanceToPortMapFile,nLat, nLon)
                
        if(self.h2pipelineDistanceMapFile):
          h2pipelineDistanceToPort = self.LoadInfrastructureMap(theProblemManager,"h2pipelines",self.h2pipelineDistanceMapFile)
          hydrogenValueMap[countryIndxs] -= taxRelief* theProblemManager.theHydrogenDataManager.theInfrastructureManager.CalculateRegionalTransportationExpenses( \
                                                              distanceToRoad[countryIndxs], roadTransportationDistance[countryIndxs], \
                                                              distanceToRail[countryIndxs], railTransporationDistance[countryIndxs], \
                                                              capacityFactorMap[countryIndxs],hydrogenPlantCapacity,hydrogenDiscountedTotalMassFunc,
                                                              h2pipelineDistanceToPort[countryIndxs])
        else:
          hydrogenValueMap[countryIndxs] -= taxRelief* theProblemManager.theHydrogenDataManager.theInfrastructureManager.CalculateRegionalTransportationExpenses( \
                                                              distanceToRoad[countryIndxs], roadTransportationDistance[countryIndxs], \
                                                              distanceToRail[countryIndxs], railTransporationDistance[countryIndxs], \
                                                              capacityFactorMap[countryIndxs],hydrogenPlantCapacity,hydrogenDiscountedTotalMassFunc)
          


        if(self.recordTaxationCalculation):
            taxEstimate[countryIndxs]  -= taxRate*theProblemManager.theHydrogenDataManager.theInfrastructureManager.CalculateRegionalTransportationExpenses( \
                                                              distanceToRoad[countryIndxs], roadTransportationDistance[countryIndxs], \
                                                              distanceToRail[countryIndxs], railTransporationDistance[countryIndxs], \
                                                              capacityFactorMap[countryIndxs],hydrogenPlantCapacity,hydrogenDiscountedTotalMassFunc)


        if(doAll):
          filename = theProblemManager.outputPrefix+"_transportCost.npy"  
          dd = np.zeros(capacityFactorMap.shape)
          dd[countryIndxs] = taxRelief* theProblemManager.theHydrogenDataManager.theInfrastructureManager.CalculateRegionalTransportationExpenses( \
                                                              distanceToRoad[countryIndxs], roadTransportationDistance[countryIndxs], \
                                                              distanceToRail[countryIndxs], railTransporationDistance[countryIndxs], \
                                                              capacityFactorMap[countryIndxs],hydrogenPlantCapacity,hydrogenDiscountedTotalMassFunc)
          #print "Saving: ", filename
          self.SaveMap(filename, dd,False)


        if doTimer:
          timesList.append(["transport", time.time()])

        # no tax if no profit
        if(self.recordTaxationCalculation):
          taxEstimate[hydrogenValueMap <= 0.0] = 0 
          
        
        # clean up map
        hydrogenValueMap[stateIdsMap< 0.0] = np.nan  
        
        if(self.recordTaxationCalculation):
          
          filename = theProblemManager.outputPrefix+"_tax.npy"  
          
          #print "Saving: ", filename
          self.SaveMap(filename, taxEstimate,False)
          
          filename = theProblemManager.outputPrefix+"_projectLife.npy"  
          projectLife = projectLifeFunc(capacityFactorMap)
          #print "Saving: ", filename
          self.SaveMap(filename, projectLife,False)


        if(self.doEmploymentCalculation):
          projectcapexFunc = interpolate.interp1d(capacityFactors,projectCapExs)
          
          initialCapex = np.zeros(capacityFactorMap.shape)
          positiveNPVindxs = np.nan_to_num(hydrogenValueMap) >0
          initialCapex[positiveNPVindxs] = projectcapexFunc(capacityFactorMap[positiveNPVindxs])
          
          employmentEstimate = np.zeros(capacityFactorMap.shape)
          employmentEstimate = theProblemManager.theHydrogenDataManager.theEconomicDataManager.EstimateDirectEmployment(initialCapex)
          
          # clip to bounding box
          if(len(self.dataBoundingBox) == 4 and len(self.outputBoundingBox) == 4):
            employmentEstimate = self.ClipToOutputBounds(employmentEstimate)
            
          filename = theProblemManager.outputPrefix+"_employment.npy"  
          
          #print "Saving: ", filename
          self.SaveMap(filename, employmentEstimate,False)


        if(self.recordWaterUse):
          waterUseFunc = interpolate.interp1d(capacityFactors,waterUse)
          
          waterUseEstimate = np.zeros(capacityFactorMap.shape)
          positiveNPVindxs = np.nan_to_num(hydrogenValueMap) >0
          waterUseEstimate[positiveNPVindxs] = waterUseFunc(capacityFactorMap[positiveNPVindxs])
          
          # clip to bounding box
          if(len(self.dataBoundingBox) == 4 and len(self.outputBoundingBox) == 4):
            waterUseEstimate = self.ClipToOutputBounds(waterUseEstimate)
            
          filename = theProblemManager.outputPrefix+"_waterUse.npy"  
          
          #print "Saving: ", filename
          self.SaveMap(filename, waterUseEstimate,False)
          
          
        
        # clip to bounding box
        if(len(self.dataBoundingBox) == 4 and len(self.outputBoundingBox) == 4):
          hydrogenValueMap = self.ClipToOutputBounds(hydrogenValueMap)
          
        
        # choose better of hydrogen value and hydrogen value with storage
        if(self.doHydrogenStorageCalculation):
          if(self.hydrogenStorageOnly):
            hydrogenValueMap = NPVwStorage
          else:
            hydrogenValueMap = np.maximum(hydrogenValueMap,NPVwStorage)



        if doTimer:
          timesList.append(["end", time.time()])

        ######################################################################



        if(doPlots):
          pl.figure()

          hydrogenValueMap[hydrogenValueMap< 0.0] = 0.0 
          pl.imshow(hydrogenValueMap,origin="lower")

          pl.colorbar()

        ######################################################################

        if doTimer:
            pl.figure()
            dts = np.array( [tt[1] for tt in timesList] )
            dts -= dts[0]
            labels =  [tt[0] for tt in timesList] 

            pl.bar(range(len(dts)),dts,align="center")
            pl.xticks(range(len(dts)), labels)


            pl.bar(range(1,len(dts)),dts[1:]-dts[:-1],color="r",align="center")

            for i in range(1,len(dts)):
              print(labels[i], dts[i], dts[i]-dts[i-1])




        if(doPlots):
          pl.show()    
        else:
          filename = theProblemManager.outputPrefix+"."+theProblemManager.outputType
          
          if(doAll):
            filename = theProblemManager.outputPrefix+".npy"
            
          if (theProblemManager.theHydrogenDataManager.thePowerPlant.type == 'CSP'):
            self.type = 'HYDROGEN CSP'
          else:
            self.type = 'HYDROGEN'
          self.SaveMap(filename, hydrogenValueMap,theProblemManager.recordRange)
          self.type = "NPV"
        
        return hydrogenValueMap



    def RunHydrogenStorageNPVcalculation(self,theProblemManager):
      """Calculates the NPV for with hydrogen sent to a storage facility rather than to port.
      
       Notes:
       - Price paid for hydrogen sold to storage facility may differ from hydrogen sold at port (price is set as "StoredH2" rather than "H2" in commodity prices).
       - Calculation is performed by replacing transportation distance maps to port with those for storage
       - Original map filenames are restored after the calculation is completed. 
       """
      
      # store current values 
      
      original_doHydrogenStorageCalculation = self.doHydrogenStorageCalculation
      
      original_h2pipelineDistanceMapFile = self.h2pipelineDistanceMapFile
      original_roadDistanceToPortMapFile = self.roadDistanceToPortMapFile
      original_railDistanceToPortMapFile = self.railDistanceToPortMapFile
      
      original_hydrogenPrice = theProblemManager.theHydrogenDataManager.theEconomicDataManager.commodityPrices["H2"]
      
      if( "StoredH2" in theProblemManager.theHydrogenDataManager.theEconomicDataManager.commodityPrices ):
        stored_hydrogenPrice = theProblemManager.theHydrogenDataManager.theEconomicDataManager.commodityPrices["StoredH2"]
      else:
        stored_hydrogenPrice = original_hydrogenPrice
      
      # output options
      original_outputType = theProblemManager.outputType
      original_saveResult = self.saveResult
      
      # Replace values with storage maps and prices before NPV calculation
      
      self.doHydrogenStorageCalculation = False  # little ugly - have to set to false so we don't end up back here. 
      
      self.distanceToHydrogenStorageMapFile = self.distanceToHydrogenStorageMapFile
      self.roadDistanceToPortMapFile = self.roadDistanceToHydrogenStorageMapFile
      self.railDistanceToPortMapFile = self.railDistanceToHydrogenStorageMapFile
      
      theProblemManager.theHydrogenDataManager.theEconomicDataManager.commodityPrices["H2"] = stored_hydrogenPrice
      
      # suppress output
      theProblemManager.outputType = ""
      self.saveResult = False
      
      # run NPV calculation
      hydrogenValueMap = self.RunHydrogenNPVcalculation(theProblemManager)
      
      # reset previous values before returning
      
      self.doHydrogenStorageCalculation = original_doHydrogenStorageCalculation
      
      self.h2pipelineDistanceMapFile = original_h2pipelineDistanceMapFile
      self.roadDistanceToPortMapFile = original_roadDistanceToPortMapFile
      self.railDistanceToPortMapFile = original_railDistanceToPortMapFile
      
      theProblemManager.theHydrogenDataManager.theEconomicDataManager.commodityPrices["H2"] = original_hydrogenPrice
      theProblemManager.outputType = original_outputType
      self.saveResult  = original_saveResult
      
      return hydrogenValueMap


    ######################################################################

    def AddInfrastructure(self,type,coords):
      """Records added roads, rail and power water and gas transmission lines.
      
      Args:
          type: The type of infrastructure to be added. 
          coords: List of lat,lon coordinates denoting the locations or path of the added infrastructure. 
          
      Note:
          Adding infrastructure does not automatically update the associated maps. 
          UpdateDistancMap and UpdateDistanceAndTravelMap must be called explicitly. 
      """
    
      numPoints = int(len(coords)/2)
      
      isLineType = False
      if(type in ["roads","raillines","powerlines","gaslines","waterlines"]):
        isLineType =  True
      
      if(isLineType):  
        for i in range(numPoints-1):  #check - this may not be needed (i.e. just store lines as lists of coords)
          self.addedInfrastructure[type].append( [ [coords[2*i],coords[2*i+1]], [coords[2*i+2],coords[2*i+3]] ]  )
      else:
        for i in range(numPoints):
          self.addedInfrastructure[type].append( [coords[2*i],coords[2*i+1]]  )
          
      return


    def ResetInfrastructure(self,type):
      """Removes added roads, rail and power water and gas transmission lines.
      
      Args: 
        type: the type of infrastructure to be reset.
      """
      self.addedInfrastructure[type] = []
      return 


    def UpdateDistanceMap(self,theProblemManager,lines,distanceData):
        """ Update distance to infrastructure/locations.
        
        Args:
          lines: New road and rail segments to be added.
          distanceData: Distance map updated by the function (data overwritten).
        
        Note:
          Radial distance function is used to update distances to new features.  
        """
            
        dataMinLon,dataMinLat,dataMaxLon,dataMaxLat = self.dataBoundingBox[:]

        # find km extent of lines
        maxY,maxX = LatLonDegreeToKMCoords(dataMaxLat,dataMaxLon)
        minY,minX = LatLonDegreeToKMCoords(dataMinLat,dataMinLon)

        for line in lines:
          xy = np.array(line)

          xy[:,1], xy[:,0] = LatLonDegreeToKMCoords(xy[:,1], xy[:,0])

          mX = np.min(xy[:,0])
          mY = np.min(xy[:,1])
          mxX = np.max(xy[:,0])
          mxY = np.max(xy[:,1])

          if(minX > mX):
            minX = mX

          if(minY > mY):
            minY = mY

          if(maxX < mxX):
            maxX = mxX

          if(maxY < mxY):
            maxY = mxY


        minX = np.floor(minX)
        minY = np.floor(minY)

        maxX = np.ceil(maxX)
        maxY = np.ceil(maxY)

        dX = (maxX-minX)    # ~ 1 pix per km
        dY = (maxY-minY)     #  ~ 1 pix per km
        
        # image of newly added transmission lines
        img = np.zeros([int(dY+1),int(dX+1)])

        for line in lines:
          xy = np.array(line)
          xy[:,1], xy[:,0] = LatLonDegreeToKMCoords(xy[:,1], xy[:,0])
          xy[:,0] -= minX
          xy[:,1] -= minY
          xy = np.round( xy )
          for i in range(len(xy)-1 ):
            x0,y0 = xy[i,:]
            x1,y1 = xy[i+1,:]
            rr, cc, val = line_aa(int(y0),int(x0), int(y1),int(x1) )
            img[rr, cc] = val 
        
        # get distances to lines
        distimg, distIndx = ndimage.distance_transform_edt(img == 0, return_indices=True)

        nLat,nLon = distanceData.shape
        
        lats = np.linspace(dataMinLat, dataMaxLat,nLat)
        lons = np.linspace(dataMinLon, dataMaxLon,nLon)

        latGrid,lonGrid = np.meshgrid( lats,lons, indexing="ij")
        yCoords,xCoords = LatLonDegreeToKMCoords(latGrid,lonGrid)
        
        minX = np.floor(minX)
        minY = np.floor(minY)

        xs = np.array(range(int(dX+1))) + minX
        ys = np.array(range(int(dY+1))) + minY

        finalDistances = np.zeros(latGrid.shape )
        finalDistances.ravel()[:] = interpolate.interpn([ys,xs],distimg,np.array([yCoords.ravel(),xCoords.ravel()]).T, bounds_error=False ,  fill_value=None ) # non interpolates outside grid
        
        distanceData = np.minimum(finalDistances,distanceData)  #nb distance data is overwritten
        
        return distanceData


    def UpdateDistanceAndTravelMap(self,theProblemManager,lines,networkData,networkOrigin,nLat,nLon):
        """ Update road and rail network distance and transmission distances to port.
        
        Args:
          lines: New road and rail segments to be added.
          networkData: contains the map of the entire road/rail network (not just in the region of interest) at a 1km scale. 
          networkOrigin: contains is the origin of the map at the 1km scale resolution 
          nLat,nLon: is the output resolution for the lat/lon grid.
        
        Notes:
          Port coordinates should be connected to/consistent with the road/rail network 
        """

        # Add new ports/refineries to networkData
        #########################################
                
        portLatLongData =  self.portAndRefiningCenterCoordinates
        if(self.addedInfrastructure["ports"]):
          for newCoord in self.addedInfrastructure["ports"]:
            if not (np.array(portLatLongData)==np.array(newCoord[::-1])).any(): # do not duplicate points
              portLatLongData.append(newCoord[::-1])
        
        if(self.addedInfrastructure["refineries"]):
          for newCoord in self.addedInfrastructure["refineries"]:
            if not (np.array(portLatLongData)==np.array(newCoord[::-1])).any(): # do not duplicate points
              portLatLongData.append(newCoord[::-1])
    
        dataMinLon,dataMinLat,dataMaxLon,dataMaxLat = self.dataBoundingBox[:] # i.e. output data dimesions

        # add new transport lines to networkData
        ########################################
        minX = networkOrigin[1]
        minY = networkOrigin[0]

        # add labels to each road node
        isNetworkNode = (networkData > 0).astype(int)
        orgDistimg, orgDistIndx = ndimage.distance_transform_edt(networkData == 0, return_indices=True)
        
        for line in lines:
          xy = np.array(line)
          xy[:,1], xy[:,0] = LatLonDegreeToKMCoords(xy[:,1], xy[:,0])
          xy[:,0] -= minX
          xy[:,1] -= minY
          xy = np.round( xy )
          for i in range(len(xy)-1 ):
            x0,y0 = xy[i,:]
            x1,y1 = xy[i+1,:]
            rr, cc, val = line_aa(int(y0),int(x0), int(y1),int(x1) )
            networkData[rr, cc] = val 

        # Check for connectivity
        isNewNetworkNode = (networkData > 0).astype(int)
        distimg, distIndx = ndimage.distance_transform_edt(networkData == 0, return_indices=True)
        
        minDist = 0.
        isNewNetworkBool = (isNewNetworkNode-isNetworkNode).astype(bool)
        if isNewNetworkBool.any():
          minDist = orgDistimg[isNewNetworkBool].min()
        
        # Add line to achieve connectivity
        if minDist > 1:
          # Find nearest node in existing network
          newDistimg, newDistIndx = ndimage.distance_transform_edt(isNewNetworkNode-isNetworkNode == 0, return_indices=True)
          newDistimg[isNetworkNode==False] = np.inf
          y1,x1 = newDistIndx[:,newDistimg==newDistimg.min()].reshape((2,))
          # Find nearest node in new feature
          orgDistimg[(isNewNetworkNode-isNetworkNode).astype(bool)==False] = np.inf
          y0,x0 = orgDistIndx[:,orgDistimg==orgDistimg.min()].reshape((2,))
          # Add connection between features
          rr, cc, val = line_aa(int(y0),int(x0), int(y1),int(x1) )
          networkData[rr, cc] = val 

        # get distances to lines
        distimg, distIndx = ndimage.distance_transform_edt(networkData == 0, return_indices=True)
        
        # calculate travel distances on transportation network
        #######################################################
        
        isNetworkNode = networkData > 0 

        # get the network node coords
        coords  = np.where(isNetworkNode)
        nodeCoords = np.column_stack([coords[0],coords[1]])

        numNodes = len(nodeCoords)


        # add labels to each road node
        nodeIds = -np.ones(networkData.shape,dtype=np.int)

        for i in range(numNodes):
          nodeIds[nodeCoords[i,0],nodeCoords[i,1]] = i

        # get the node neighbours
        dirs = np.array( [[0,1],[0,-1],[1,0],[-1,0],[-1,-1],[-1,1],[1,-1],[-1,-1] ])

        nodeNbrs = -np.ones([numNodes,9],dtype=np.int)
        for i in range(numNodes):
          nn = nodeIds[nodeCoords[i,0]-1:nodeCoords[i,0]+2 ,nodeCoords[i,1]-1:nodeCoords[i,1]+2].ravel()
          if(len(nn)):
            nodeNbrs[i,:len(nn)] = nn
  
        #
        parentNodes = -np.ones(numNodes,dtype=int)  # needed?
        distances = np.zeros(numNodes,dtype=float)  # assumes no disconnected components


        nodesToCheck = []

        # loop over ports 
        #print("Seeding ports")
        for portLatLong in portLatLongData:
          pY,pX = LatLonDegreeToKMCoords(portLatLong[0],portLatLong[1])
          pY -= minY
          pX -= minX
  
          pY = np.int(pY)
          pX = np.int(pX)
  
          root = nodeIds[pY,pX] 

          parentNodes[root] = root
          nodesToCheck.append(root)

        while (len(nodesToCheck) > 0):
            nextRank = []
            for node in nodesToCheck:
              newNbrs = nodeNbrs[node,:]
              nodeKM_y,nodeKM_x = nodeCoords[node,0],nodeCoords[node,1]
      
              for id in newNbrs:
                if ( id > 0 and id != node):
        
                  nbrKM_y,nbrKM_x = nodeCoords[id,0],nodeCoords[id,1]
          
                  if(parentNodes[id] < 0):  # unassigned
                    parentNodes[id] = node 
                    nextRank.append(id)
            
            
                    distances[id] = distances[node] + ( ( nodeKM_x - nbrKM_x )**2  + (  nodeKM_y - nbrKM_y   )**2 )**0.5
                  else:
          
                    dd = distances[node] + ( ( nodeKM_x - nbrKM_x )**2  + (  nodeKM_y - nbrKM_y   )**2 )**0.5
                    if (dd < distances[id] ):
                      parentNodes[id] = node 
                      distances[id] =  dd
                      nextRank.append(id)

            nodesToCheck = list(nextRank)         



        ###  put distance to port on road map

        distanceToPortMap = np.zeros(networkData.shape,dtype=np.float)

        for i in range(numNodes):
          distanceToPortMap[nodeCoords[i,0],nodeCoords[i,1]] = distances[i]
          #print distances[i]

        
        # map grid from lat lon back to distimage to make lookup table
        ###############################################################
        
        lats = np.linspace(dataMinLat, dataMaxLat,nLat)
        lons = np.linspace(dataMinLon, dataMaxLon,nLon)

        latGrid,lonGrid = np.meshgrid( lats,lons, indexing="ij")
        yCoords,xCoords = LatLonDegreeToKMCoords(latGrid,lonGrid)
        
        dY,dX = distimg.shape

        xs = np.array(range(int(dX))) + minX
        ys = np.array(range(int(dY))) + minY
        
        
        # map travel distance to output (lat-lon grid) shape
        ############################################
        
        distanceData = np.zeros(latGrid.shape )
        distanceData.ravel()[:] = interpolate.interpn([ys,xs],distimg,np.array([yCoords.ravel(),xCoords.ravel()]).T, bounds_error=False ,  fill_value=None ) # none interpolates outside grid
        
        
        # generate total travel distance km grid (bit inefficient) and then map result to output (lat-lon grid) shape
        #####################################################################################
        
        totalTravelDistanceOnKmGrid = np.array(distimg)
        totalTravelDistanceOnKmGrid.ravel()[:] += distanceToPortMap[(distIndx[0].ravel(),distIndx[1].ravel() )]

        travelDistanceData = np.zeros(latGrid.shape )
        travelDistanceData.ravel()[:] = interpolate.interpn([ys,xs],totalTravelDistanceOnKmGrid,np.array([yCoords.ravel(),xCoords.ravel()]).T, bounds_error=False ,  fill_value=None ) # none interpolates outside grid
        
        
        
        return distanceData, travelDistanceData
  
  