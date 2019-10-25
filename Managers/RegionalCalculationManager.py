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

import numpy as np
import pylab as pl
from scipy import interpolate  # for interp1d

# GeoTIFF functionality
gdal_flag = False
try:
    from osgeo import gdal,osr
    gdal_flag = True
except ImportError:
    pass

import time  

# Common
from Common.Common import Todo,Fixme

from IO.XML import HasChild,GetChild,AddChild,GetChildren
from IO.XML import HasAttribute,GetAttributeString,GetAttributeFileString,GetAttributeFileStringOrDefault
from IO.XML import GetAttributeValue,GetAttributeVector, SetAttributeString,GetXMLTag

# Managers

from Units.UnitManager import UnitManager
#from ProblemManager import ProblemManager
from IO.CommandLine import ParseCommandLineArgs


def array_to_raster(name, array, xpix=0.010002603, ypix=0.010002603, xmin=111.994998809, ymax=-9.99499852898, epsg=4326):
    """Array > Raster
    Save a raster from an array.

    :name array:  str
    :param array: ndarray
    :xpix array:  float
    :ypix array:  float
    :xmin array:  float
    :ymax array:  float
    :epsg array:  int
    """
    
    ny,nx = array.shape
    
    proj = osr.SpatialReference()
    proj.ImportFromEPSG(epsg)
    driver = gdal.GetDriverByName('GTiff')
    
    dataset = driver.Create(
        name,
        nx,
        ny,
        1,
        gdal.GDT_Float32 )

    dataset.SetGeoTransform((
        xmin,
        xpix,
        0,
        ymax,
        0,
        -ypix))

    dataset.SetProjection(proj.ExportToWkt())
    dataset.GetRasterBand(1).WriteArray(array)
    dataset.GetRasterBand(1).SetNoDataValue(-99999)
    dataset.FlushCache()  # Write to disk.
    dataset = None

class RegionalCalculationManager():
    def __init__(self):
      """
      Create an empty regional calculation data manager and default variables. 
      """
      self.type = ""   # NPV, benefit_cost_ratio, breakeven_grade
      
      self.coverDepthMapFile = ""
      self.waterDistanceMapFile = ""
      self.powerDistanceMapFile = ""
      self.gasDistanceMapFile = ""
      self.railDistanceMapFile = ""
      self.roadDistanceMapFile = ""
      self.railDistanceToPortMapFile = ""
      self.roadDistanceToPortMapFile = ""
      
      self.coverOffset = 0.0
      
      self.dataBoundingBox = []
      self.outputBoundingBox = []
      
      self.stride = 1

      self.stateIdsMapFile = ""
      
      self.doEmploymentCalculation = False

    def ParseXMLNode(self, regionalDataNode):
      """
      Generate regional calculation data from xml tree node. 
      """
      
      self.type = GetAttributeString(regionalDataNode,"type")
          
      self.coverDepthMapFile = GetAttributeFileString(regionalDataNode,"coverDepth")
      self.waterDistanceMapFile = GetAttributeFileString(regionalDataNode,"waterDistance") 
      self.powerDistanceMapFile = GetAttributeFileString(regionalDataNode,"powerDistance")
      self.railDistanceMapFile = GetAttributeFileString(regionalDataNode,"railDistance")
      self.roadDistanceMapFile = GetAttributeFileString(regionalDataNode,"roadDistance")
      self.railDistanceToPortMapFile = GetAttributeFileString(regionalDataNode,"railTransportationDistance")
      self.roadDistanceToPortMapFile = GetAttributeFileString(regionalDataNode,"roadTransportationDistance")
      
      self.gasDistanceMapFile = GetAttributeFileStringOrDefault(regionalDataNode,"gasDistance",self.gasDistanceMapFile)  # optional

      self.stateIdsMapFile = GetAttributeFileString(regionalDataNode,"states")
      
      if(HasAttribute(regionalDataNode,"dataBoundingBox")):
        self.dataBoundingBox = GetAttributeVector(regionalDataNode,"dataBoundingBox")
      if(HasAttribute(regionalDataNode,"outputBoundingBox")):
        self.outputBoundingBox = GetAttributeVector(regionalDataNode,"outputBoundingBox")
      if(HasAttribute(regionalDataNode,"coverOffset")):
        self.coverOffset = GetAttributeValue(regionalDataNode,"coverOffset")
  
      if(HasAttribute(regionalDataNode,"stride")):
        self.stride = GetAttributeValue(regionalDataNode,"stride")
        
      if(HasAttribute(regionalDataNode,"estimateEmployment")):
        self.doEmploymentCalculation = GetAttributeValue(regionalDataNode,"estimateEmployment")
        
    def WriteXMLNode(self, node):
      """
      Write regional calculation to xml node
      """
      SetAttributeString(node,"type",self.type)
      
      return node
    
    
    def LoadMap(self,filename):
      if filename[-3:] == "npy":
        rv = np.load(filename)
      else:
        rv = np.loadtxt(filename)

      if(self.stride > 1):
        rv = rv[::self.stride,::self.stride]
      return rv
    
    
    def SaveMap(self,filename, data, reportRange):
        
      data_max = np.nanmax(data)
      data_min = np.nanmin(data)
      label_min = data_min
      label_max = data_max
      cmap = 'viridis'

      if(filename[-3:] in ["jpg","tif","png"]  or  filename[-4:] == "tiff"):
        if ((filename[-6:] == "geotif") or (filename[-7:] == "geotiff")) and gdal_flag:
          if self.type == "NPV":
            array_to_raster(filename.replace('geotif','tif'), data[::-1]*1e-06)
          else:
            array_to_raster(filename.replace('geotif','tif'), data[::-1])
          filename = filename.replace("tiff", "tif")
          filename = filename.replace("geotif",  "png")
        if self.type == "NPV":
          if abs(data_min) > data_max:
            data_max = -1.*data_min
          if data_max > abs(data_min):
            data_min = -1.*data_max
          cmap = 'seismic'
          label_min = data_min * 1e-06
          label_max = data_max * 1e-06
        elif self.type == "benefit_cost_ratio":
          cmap = 'PiYG'
          data_min = np.log(data_min)
          data_max = np.log(data_max)
          if abs(data_min) > data_max:
            data_max = -1.*data_min
          if data_max > abs(data_min):
            data_min = -1.*data_max
          label_min = np.exp(data_min)
          label_max = np.exp(data_max)
          data = np.log(data)
        elif self.type == "breakeven_grade":
          cmap = 'viridis'
        elif self.type == "employment":
          cmap = 'viridis'
        else:
          cmap = 'viridis'
        #if not gdal_flag:
        pl.imsave(filename,data,origin="lower",cmap=pl.get_cmap(cmap),vmin=data_min,vmax=data_max)
      elif( filename[-3:] == "npy"):
        np.save(filename,data)
      elif( filename[-3:] == "txt"):
        np.savetxt(filename,data)
      else:
        print "Error unrecognized output type: ", filename
      
      if(reportRange):
        np.savetxt(filename+"_range.txt",[np.round(label_min,2), np.round(label_max,2)])
      
      return 0
 
    def Run(self,problemManager):
      """
      Run regional calculation
      """
      
      if(self.type == "NPV" ):
        self.RunNPVcalculation(problemManager)
      elif(self.type == "benefit_cost_ratio"):
        self.RunBenefitCostRatioCalculation(problemManager)
      elif(self.type == "breakeven_grade"):
        self.RunBreakevenGradeCalculation(problemManager)
      elif(self.type == "employment"):
        self.RunEmploymentCalculation(problemManager)
      else:
        print "Error: Unrecognized regional calculation: " +  self.type
        exit()
    
      return 
    
    def ApplyOutputBoundingBoxMask(self,stateIdsMap,replaceValue = -1):   
      # masks off regions outside the outputBoundingBox
      if(len(self.dataBoundingBox) == 4 and len(self.outputBoundingBox) == 4):
          dataMinLat,dataMinLon,dataMaxLat,dataMaxLon = self.dataBoundingBox[:]
          outputMinLat,outputMinLon,outputMaxLat,outputMaxLon = self.outputBoundingBox[:]
          
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
    
    
    def ClipToOutputBounds(self,data):   
      # clips the output to within the bounding box
      if(len(self.dataBoundingBox) == 4 and len(self.outputBoundingBox) == 4):
          dataMinLat,dataMinLon,dataMaxLat,dataMaxLon = self.dataBoundingBox[:]
          outputMinLat,outputMinLon,outputMaxLat,outputMaxLon = self.outputBoundingBox[:]
          
          nLat,nLon = data.shape
          minLatIndx = int( np.floor( nLat* (outputMinLat-dataMinLat)/(dataMaxLat-dataMinLat) ) )
          maxLatIndx = int( np.ceil( nLat* (outputMaxLat-dataMinLat)/(dataMaxLat-dataMinLat) ) )
          
          minLonIndx = int( np.floor( nLon* (outputMinLon-dataMinLon)/(dataMaxLon-dataMinLon) ) )
          maxLonIndx = int( np.ceil( nLon* (outputMaxLon-dataMinLon)/(dataMaxLon-dataMinLon)  ) )
          
          data = data[minLatIndx:maxLatIndx,minLonIndx:maxLonIndx]
      return data
            
    #########################
    ## NPV calculation
    
    def RunNPVcalculation(self,theProblemManager):
    
        doPlots = (theProblemManager.outputType == "")

        doAll = (theProblemManager.outputType == "all")
        # profiling
        doTimer = False
        timesList = []
        
        # self.doEmploymentCalculation = True

        if doTimer:
          #startTime = time.time()
          timesList.append(["start", time.time()])
        
          
        theProblemManager.theMineDataManager.theInfrastructureManager.ZeroDistanceToInfrastructure()

        # calculate mine cost, life etc as function of depth of cover (in meters)
        coverMin = 0
        coverMax = 2000.

        if coverMax > 200:
          coverDepths = np.append ( np.linspace(coverMin,200,20), np.linspace(200,coverMax,20) )
        else:
          coverDepths = np.linspace(coverMin,coverMax,20)
  
        discountRate = theProblemManager.theMineDataManager.theEconomicDataManager.discountRate

        mineValues = []
        mineTypes =[]
        mineYears = []
        concentrateDiscountedTotalMasses = []
        concentrateCapacities = []
        
        mineCapExs = []  # for employment calculation


        stateStrs = ["NSW","NT","QLD","SA","WA","VIC","ACT","TAS"]
        stateIds = [0,1,2,3,4,5,6,7]
        stateMineValues = {}
        for state in stateStrs:
          stateMineValues[state] = []

        for depthOfCover in coverDepths:

          # fixme need to set distance to infrastructure to 0
          theProblemManager.theMineDataManager.SetCoverDepth(depthOfCover)
  
          value = theProblemManager.EvaluateOCUGMines()
          mineValues.append(value)
  
          mineLife = theProblemManager.theMineDataManager.theMiningSystem.mineLife
          mineYears.append(mineLife)
          
          if(self.doEmploymentCalculation):
            mineAndProcCapex = theProblemManager.theMineDataManager.theEconomicDataManager.netCapex[0]
            mineCapExs.append(mineAndProcCapex)
  
          mineTypeStr = theProblemManager.theMineDataManager.theMiningSystem.miningMethod  # 1 if UG, 0.5 mixed, 0 OC 
  
          mineType = 0.0   # open cut
          if(mineTypeStr[:2] == "UG"):
            mineType = 1.0
          elif(mineTypeStr == "OCUG"):
            mineType = 0.5
    
          for state in stateStrs:
            theProblemManager.theMineDataManager.theEconomicDataManager.SetState(state)
            theProblemManager.theMineDataManager.CalculateTaxes(theProblemManager)
            theProblemManager.theMineDataManager.CalculateAfterTaxCashFlow(theProblemManager)
            theProblemManager.theMineDataManager.CalculateEconomicIndicators(theProblemManager)
            value = theProblemManager.theMineDataManager.theEconomicDataManager.atNPV
            stateMineValues[state].append(value)
    
  
          mineTypes.append(mineType)
  
          # discounted total mass
          discountFactors = (1+theProblemManager.theMineDataManager.theEconomicDataManager.discountRate)**np.array(range( int( np.ceil(mineLife) ) ) )
          concentrateMass = np.sum(theProblemManager.theMineDataManager.theProcessingSystem.concentrateProduced/discountFactors)
          concentrateDiscountedTotalMasses.append(concentrateMass)
  
          concentrateCapacity = np.max(theProblemManager.theMineDataManager.theProcessingSystem.concentrateProduced)
          concentrateCapacities.append(concentrateCapacity)

        #print coverDepths
        #print mineValues
  
        if doTimer:
          #coverFuncTime = time.time()
          timesList.append(["coverFunc", time.time()])
          


        mineValues = np.array(mineValues,dtype=np.float)
        mineTypes = np.array(mineTypes,dtype=np.float)
        mineYears = np.array(mineYears,dtype=np.float)
        concentrateDiscountedTotalMasses = np.array(concentrateDiscountedTotalMasses,dtype=np.float)
        concentrateCapacities = np.array(concentrateCapacities,dtype=np.float)


        mineValueFunc = interpolate.interp1d(coverDepths,mineValues)
        mineTypeFunc =interpolate.interp1d(coverDepths,mineTypes)
        mineLifeFunc = interpolate.interp1d(coverDepths,mineYears)
        concentrateDiscountedTotalMassFunc = interpolate.interp1d(coverDepths,concentrateDiscountedTotalMasses)
        concentrateCapacityFunc = interpolate.interp1d(coverDepths,concentrateCapacities)
        
        if(self.doEmploymentCalculation):
          minecapexFunc = interpolate.interp1d(coverDepths,mineCapExs)



        if(doPlots):

			for state in stateStrs:
			  pl.plot(coverDepths,stateMineValues[state],label=state)
		  
			pl.legend()
			pl.plot(coverDepths,mineValues)
			pl.figure()
			pl.plot(coverDepths,mineTypes)


        #pl.show()


        ## regional mine cost calculation

        coverDepthMapFile = self.coverDepthMapFile
        waterDistanceMapFile = self.waterDistanceMapFile
        powerDistanceMapFile = self.powerDistanceMapFile
        gasDistanceMapFile = self.gasDistanceMapFile
        railDistanceMapFile = self.railDistanceMapFile
        roadDistanceMapFile = self.roadDistanceMapFile
        railDistanceToPortMapFile = self.railDistanceToPortMapFile
        roadDistanceToPortMapFile = self.roadDistanceToPortMapFile
        
        stateIdsMapFile = self.stateIdsMapFile

        if doTimer:
          #startMapsTime = time.time()
          timesList.append(["startMaps", time.time()])

        stateIdsMap =  self.LoadMap(stateIdsMapFile)
        self.ApplyOutputBoundingBoxMask(stateIdsMap) 
        


        ###################################
        
        if(coverDepthMapFile):
          coverMap =  self.LoadMap(coverDepthMapFile)
          coverMap += self.coverOffset
          
          # for MD basin only - remove
          #cc = coverMap==0
          #cc[5000::,:] = 0 
          #cc[:,5000::] = 0 
          #stateIdsMap[cc] = -1
          
          # bound cover map
          coverMap[coverMap > coverMax] = coverMax
          coverMap[coverMap < coverMin] = coverMin
  
        else:
          # nonsense values for testing
          coverMap =  self.LoadMap(waterDistanceMapFile) # fixme
          XX,YY = np.mgrid[0:1:coverMap.shape[0]*1j,0:1:coverMap.shape[1]*1j]   
          coverMap[:] = 25 # + 1.*XX - 2*YY

          coverMap = 25  + 1.*XX - 2*YY
          
          
        countryIndxs = stateIdsMap >= 0.0

        mineValueMap = np.zeros(coverMap.shape)

        taxRelief = 1.0 - 0.3  # 30% tax relief from incometax

        ####
        for i in range(len(stateStrs)):
          stateIndxs = stateIdsMap == stateIds[i]
          mineValueFunc = interpolate.interp1d(coverDepths,stateMineValues[stateStrs[i]])
          mineValueMap[stateIndxs] = mineValueFunc(coverMap[stateIndxs])

        
        
        
        if(doAll):
          filename = theProblemManager.outputPrefix+"_mineValue.npy"  
          #print "Saving: ", filename
          self.SaveMap(filename, mineValueMap,False)
          
          filename = theProblemManager.outputPrefix+"_mineValueNT.npy"  
          NTValueFunc = interpolate.interp1d(coverDepths,stateMineValues["NT"])
          
          dd = np.zeros(coverMap.shape)
          dd[countryIndxs] = NTValueFunc(coverMap[countryIndxs])
          #print "Saving: ", filename
          self.SaveMap(filename, dd,False)

        if doTimer:
          #mineValueFuncTime = time.time()
          timesList.append(["mineValue", time.time()])

        # Water Expenses
        distanceToWater =  self.LoadMap(waterDistanceMapFile)
        mineValueMap[countryIndxs] -=  taxRelief* theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalWaterExpenses(theProblemManager,distanceToWater[countryIndxs])
        


        if(doAll):
          filename = theProblemManager.outputPrefix+"_waterCost.npy"  
          dd = np.zeros(coverMap.shape)
          dd[countryIndxs] = taxRelief* theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalWaterExpenses(theProblemManager,distanceToWater[countryIndxs])
          #print "Saving: ", filename
          self.SaveMap(filename, dd,False)

        if doTimer:
          #waterTime = time.time()
          timesList.append(["water", time.time()])

        # Power Expenses
        #Todo("Pre-calculate power supply costs")
        distanceToPower =  self.LoadMap(powerDistanceMapFile)
        
        if(theProblemManager.theMineDataManager.theInfrastructureManager.calculateGas):
          powerCost = theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalPowerExpenses(theProblemManager,distanceToPower[countryIndxs])
          
          distanceToGas =  self.LoadMap(gasDistanceMapFile)
          gasCost =  theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalGasExpenses(theProblemManager,distanceToGas[countryIndxs])

          mineValueMap[countryIndxs] -=  taxRelief*  np.minimum(powerCost,gasCost)
        else:
          mineValueMap[countryIndxs] -=  taxRelief* theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalPowerExpenses(theProblemManager,distanceToPower[countryIndxs])





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
          #powerTime = time.time()
          timesList.append(["power", time.time()])

        # Transportation Expenses
        distanceToRoad =  self.LoadMap(roadDistanceMapFile)
        distanceToRail =  self.LoadMap(railDistanceMapFile)
        roadTransportationDistance =  self.LoadMap(roadDistanceToPortMapFile)
        railTransporationDistance =  self.LoadMap(railDistanceToPortMapFile)


        mineValueMap[countryIndxs] -= taxRelief* theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalTransportationExpenses( \
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
          #transportTime = time.time()
          timesList.append(["transport", time.time()])



        mineValueMap[stateIdsMap< 0.0] = np.nan  


        if(self.doEmploymentCalculation):
          minecapexFunc = interpolate.interp1d(coverDepths,mineCapExs)
          
          initialCapex = np.zeros(coverMap.shape)
          positiveNPVindxs = np.nan_to_num(mineValueMap) >0
          initialCapex[positiveNPVindxs] = minecapexFunc(coverMap[positiveNPVindxs])
          
          employmentEstimate = np.zeros(coverMap.shape)
          employmentEstimate = theProblemManager.theMineDataManager.theEconomicDataManager.EstimateDirectEmployment(initialCapex)
          
          # clip to bounding box
          if(len(self.dataBoundingBox) == 4 and len(self.outputBoundingBox) == 4):
            employmentEstimate = self.ClipToOutputBounds(employmentEstimate)
            
          filename = theProblemManager.outputPrefix+"_employment.npy"  
          
          #print "Saving: ", filename
          self.SaveMap(filename, employmentEstimate,False)

        
        # clip to bounding box
        if(len(self.dataBoundingBox) == 4 and len(self.outputBoundingBox) == 4):
          mineValueMap = self.ClipToOutputBounds(mineValueMap)



        if doTimer:
          #transportTime = time.time()
          timesList.append(["end", time.time()])

        ######################################################################



        if(doPlots):
          pl.figure()

          #mineValueMap[mineValueMap< 0.0] = 0.0 
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
              print labels[i], dts[i], dts[i]-dts[i-1]




        #pl.figure()
        #pl.imshow(coverMap,origin="lower",cmap="pink")
        
        if(doPlots):
          pl.show()    
        else:
          filename = theProblemManager.outputPrefix+"."+theProblemManager.outputType
          
          if(doAll):
            filename = theProblemManager.outputPrefix+".npy"
            
          #print "Saving: ", filename
          self.SaveMap(filename, mineValueMap,theProblemManager.recordRange)
        
        return 
    
    
    #########################
    ## Benefit to cost ratio
    
    def RunBenefitCostRatioCalculation(self,theProblemManager):
        # profiling
        doTimer = False
        timesList = []
        
        doPlots = (theProblemManager.outputType == "")
        
          
        theProblemManager.theMineDataManager.theInfrastructureManager.ZeroDistanceToInfrastructure()

        if doTimer:
          #startTime = time.time()
          timesList.append(["start", time.time()])


        # calculate mine cost, life etc as function of depth of cover (in meters)
        coverMin = 0
        coverMax = 2000.

        if coverMax > 200:
          coverDepths = np.append ( np.linspace(coverMin,200,20), np.linspace(200,coverMax,20) )
        else:
          coverDepths = np.linspace(coverMin,coverMax,20)

        discountRate = theProblemManager.theMineDataManager.theEconomicDataManager.discountRate

        mineValues = []
        mineTypes =[]
        mineYears = []
        concentrateDiscountedTotalMasses = []
        concentrateCapacities = []
        discountedRevenues = []


        stateStrs = ["NSW","NT","QLD","SA","WA","VIC","ACT","TAS"]
        stateIds = [0,1,2,3,4,5,6,7]
        stateDiscountedCosts = {}
        for state in stateStrs:
          stateDiscountedCosts[state] = []

        for depthOfCover in coverDepths:

          # fixme need to set distance to infrastructure to 0
          theProblemManager.theMineDataManager.SetCoverDepth(depthOfCover)
  
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
  
          discountFactors = (1+theProblemManager.theMineDataManager.theEconomicDataManager.discountRate)**np.array(range( int( np.ceil(mineLife) ) ) )
  
          discountedRevenue = np.sum( theProblemManager.theMineDataManager.theEconomicDataManager.revenue*discountFactors)
  
          discountedRevenues.append(discountedRevenue)
    
          for state in stateStrs:
            theProblemManager.theMineDataManager.theEconomicDataManager.SetState(state)
            theProblemManager.theMineDataManager.CalculateTaxes(theProblemManager)
            theProblemManager.theMineDataManager.CalculateAfterTaxCashFlow(theProblemManager)
            theProblemManager.theMineDataManager.CalculateEconomicIndicators(theProblemManager)
    
            discountedCost = np.sum( ( theProblemManager.theMineDataManager.theEconomicDataManager.revenue \
                                       - theProblemManager.theMineDataManager.theEconomicDataManager.btNCF \
                                       + theProblemManager.theMineDataManager.theEconomicDataManager.royalties )*discountFactors)
            stateDiscountedCosts[state].append(discountedCost)
    
  
          mineTypes.append(mineType)
  
          # discounted total mass
          concentrateMass = np.sum(theProblemManager.theMineDataManager.theProcessingSystem.concentrateProduced*discountFactors)
          concentrateDiscountedTotalMasses.append(concentrateMass)
  
          concentrateCapacity = np.max(theProblemManager.theMineDataManager.theProcessingSystem.concentrateProduced)
          concentrateCapacities.append(concentrateCapacity)
  

        #print coverDepths
        #print mineValues
  
        if doTimer:
          #coverFuncTime = time.time()
          timesList.append(["coverFunc", time.time()])


        mineTypes = np.array(mineTypes,dtype=np.float)
        mineYears = np.array(mineYears,dtype=np.float)
        concentrateDiscountedTotalMasses = np.array(concentrateDiscountedTotalMasses,dtype=np.float)
        concentrateCapacities = np.array(concentrateCapacities,dtype=np.float)


        mineTypeFunc =interpolate.interp1d(coverDepths,mineTypes)
        mineLifeFunc = interpolate.interp1d(coverDepths,mineYears)
        concentrateDiscountedTotalMassFunc = interpolate.interp1d(coverDepths,concentrateDiscountedTotalMasses)
        concentrateCapacityFunc = interpolate.interp1d(coverDepths,concentrateCapacities)

        mineRevenueFunc = interpolate.interp1d(coverDepths,discountedRevenues)


        if(doPlots):
          pl.plot(coverDepths,mineValues)
          pl.figure()
          pl.plot(coverDepths,mineTypes)

          #pl.show()


        ## regional mine cost calculation

        coverDepthMapFile = self.coverDepthMapFile
        waterDistanceMapFile = self.waterDistanceMapFile
        powerDistanceMapFile = self.powerDistanceMapFile
        railDistanceMapFile = self.railDistanceMapFile
        roadDistanceMapFile = self.roadDistanceMapFile
        railDistanceToPortMapFile = self.railDistanceToPortMapFile
        roadDistanceToPortMapFile = self.roadDistanceToPortMapFile

        stateIdsMapFile  = self.stateIdsMapFile

        if doTimer:
          #startMapsTime = time.time()
          timesList.append(["startMaps", time.time()])

        stateIdsMap =  self.LoadMap(stateIdsMapFile)
        self.ApplyOutputBoundingBoxMask(stateIdsMap) 

        countryIndxs = stateIdsMap >= 0.0

        ###################################
        if(coverDepthMapFile):
          coverMap =  self.LoadMap(coverDepthMapFile)
          coverMap += self.coverOffset
          # bound cover map
          coverMap[coverMap > coverMax] = coverMax
          coverMap[coverMap < coverMin] = coverMin
  
        else:
          # nonsense values
          coverMap =  self.LoadMap(waterDistanceMapFile) # fixme
          XX,YY = np.mgrid[0:1:coverMap.shape[0]*1j,0:1:coverMap.shape[1]*1j]    #random.rand(coverMap.shape[0],coverMap.shape[1]) # fixme
          coverMap[:] = 25 # + 1.*XX - 2*YY

          coverMap = 25  + 1.*XX - 2*YY

        mineCostMap = np.zeros(coverMap.shape)


        ####
        for i in range(len(stateStrs)):
          stateIndxs = stateIdsMap == stateIds[i]
          mineCostFunc = interpolate.interp1d(coverDepths,stateDiscountedCosts[stateStrs[i]])
          mineCostMap[stateIndxs] = mineCostFunc(coverMap[stateIndxs])

        if doTimer:
          #mineValueFuncTime = time.time()
          timesList.append(["mineValue", time.time()])

        # Water Expenses
        distanceToWater =  self.LoadMap(waterDistanceMapFile)
        mineCostMap[countryIndxs] +=  theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalWaterExpenses(theProblemManager,distanceToWater[countryIndxs])

        if doTimer:
          #waterTime = time.time()
          timesList.append(["water", time.time()])

        # Power Expenses
        #Todo("Pre-calculate power supply costs")
        distanceToPower =  self.LoadMap(powerDistanceMapFile)
        mineCostMap[countryIndxs] +=  theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalPowerExpenses(theProblemManager,distanceToPower[countryIndxs])


        if doTimer:
          #powerTime = time.time()
          timesList.append(["power", time.time()])

        # Transportation Expenses
        distanceToRoad =  self.LoadMap(roadDistanceMapFile)
        distanceToRail =  self.LoadMap(railDistanceMapFile)
        roadTransportationDistance =  self.LoadMap(roadDistanceToPortMapFile)
        railTransporationDistance =  self.LoadMap(railDistanceToPortMapFile)


        mineCostMap[countryIndxs] += theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalTransportationExpenses( \
                                                              distanceToRoad[countryIndxs], roadTransportationDistance[countryIndxs], \
                                                              distanceToRail[countryIndxs], railTransporationDistance[countryIndxs], \
                                                              coverMap[countryIndxs],concentrateCapacityFunc,concentrateDiscountedTotalMassFunc)


        if doTimer:
          #transportTime = time.time()
          timesList.append(["transport", time.time()])




        revenueCostRatioMap = np.zeros(mineCostMap.shape)
        revenueCostRatioMap[countryIndxs] =  mineRevenueFunc( coverMap[countryIndxs] )/ mineCostMap[countryIndxs]
        
        revenueCostRatioMap[stateIdsMap< 0.0] = np.nan  

        # clip to bounding box
        if(len(self.dataBoundingBox) == 4 and len(self.outputBoundingBox) == 4):
          revenueCostRatioMap = self.ClipToOutputBounds(revenueCostRatioMap)


        if doTimer:
          #transportTime = time.time()
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
              print labels[i], dts[i], dts[i]-dts[i-1]

        """
        """

    
        return 
    
    
    #########################
    ## Breakeven calculation
    
    def RunBreakevenGradeCalculation(self,theProblemManager):
    
        
        theProblemManager.theMineDataManager.theInfrastructureManager.ZeroDistanceToInfrastructure()
        
        # plotting
        doPlots = (theProblemManager.outputType == "")

        # profiling
        doTimer = False
        timesList = []

        if doTimer:
          #startTime = time.time()
          timesList.append(["start", time.time()])

        # calculate mine cost, life etc as function of depth of cover (in meters)
        coverMin = 0
        coverMax = 2000.

        if coverMax > 200:
          coverDepths = np.append ( np.linspace(coverMin,200,20), np.linspace(200,coverMax,20) )
        else:
          coverDepths = np.linspace(coverMin,coverMax,20)
  
        discountRate = theProblemManager.theMineDataManager.theEconomicDataManager.discountRate

        mineValues = []
        mineTypes =[]
        mineYears = []
        concentrateDiscountedTotalMasses = []
        concentrateCapacities = []


        stateStrs = ["NSW","NT","QLD","SA","WA","VIC","ACT","TAS"]
        stateIds = [0,1,2,3,4,5,6,7]
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

        for state in stateStrs:
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
    
          for state in stateStrs:
            theProblemManager.theMineDataManager.theEconomicDataManager.SetState(state)
            theProblemManager.theMineDataManager.CalculateTaxes(theProblemManager)
            theProblemManager.theMineDataManager.CalculateAfterTaxCashFlow(theProblemManager)
            theProblemManager.theMineDataManager.CalculateEconomicIndicators(theProblemManager)
            value = theProblemManager.theMineDataManager.theEconomicDataManager.atNPV
            stateMineValues[state].append(value)
    
          # discounted total mass - these will need to be scaled according to changes in grade in transportation cost calculation
          discountFactors = (1+theProblemManager.theMineDataManager.theEconomicDataManager.discountRate)**np.array(range( int( np.ceil(mineLife) ) ) )
          concentrateMass = np.sum(theProblemManager.theMineDataManager.theProcessingSystem.concentrateProduced/discountFactors)
          concentrateDiscountedTotalMasses.append(concentrateMass)
  
          concentrateCapacity = np.max(theProblemManager.theMineDataManager.theProcessingSystem.concentrateProduced)
          concentrateCapacities.append(concentrateCapacity)
  
          ### half grade ###
          theProblemManager.theMineDataManager.theOreBody.ScaleCommodityGrades(0.5)
  
          for state in stateStrs:
            theProblemManager.theMineDataManager.theEconomicDataManager.CalculateBeforeTaxCashFlow(theProblemManager, theProblemManager.theMineDataManager)
            theProblemManager.theMineDataManager.theEconomicDataManager.SetState(state)
            theProblemManager.theMineDataManager.CalculateTaxes(theProblemManager)
            theProblemManager.theMineDataManager.CalculateAfterTaxCashFlow(theProblemManager)
            theProblemManager.theMineDataManager.CalculateEconomicIndicators(theProblemManager)
            halfGradeValue = theProblemManager.theMineDataManager.theEconomicDataManager.atNPV
            stateMineValuesHalfGrade[state].append(halfGradeValue)
  
          ### double grade ###
  
          theProblemManager.theMineDataManager.theOreBody.ScaleCommodityGrades(4)
  
          for state in stateStrs:
            theProblemManager.theMineDataManager.theEconomicDataManager.CalculateBeforeTaxCashFlow(theProblemManager, theProblemManager.theMineDataManager)
            theProblemManager.theMineDataManager.theEconomicDataManager.SetState(state)
            theProblemManager.theMineDataManager.CalculateTaxes(theProblemManager)
            theProblemManager.theMineDataManager.CalculateAfterTaxCashFlow(theProblemManager)
            theProblemManager.theMineDataManager.CalculateEconomicIndicators(theProblemManager)
            doubleGradeValue = theProblemManager.theMineDataManager.theEconomicDataManager.atNPV
            stateMineValuesDoubleGrade[state].append(doubleGradeValue)
  
          # revert to original grade
  
          theProblemManager.theMineDataManager.theOreBody.ScaleCommodityGrades(0.5/gradeAdjustment)
  
  
  
  

        #print coverDepths
        #print mineValues
  
        if doTimer:
          #coverFuncTime = time.time()
          timesList.append(["coverFunc", time.time()])




        mineValues = np.array(mineValues,dtype=np.float)
        mineTypes = np.array(mineTypes,dtype=np.float)
        mineYears = np.array(mineYears,dtype=np.float)
        concentrateDiscountedTotalMasses = np.array(concentrateDiscountedTotalMasses,dtype=np.float)
        concentrateCapacities = np.array(concentrateCapacities,dtype=np.float)


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
          for state in stateStrs:
            pl.plot(coverDepths,stateMineValues[state],label=state)
  
          pl.legend()
          pl.plot(coverDepths,mineValues)
          pl.figure()
          pl.plot(coverDepths,mineTypes)


        #pl.show()


        ## regional mine cost calculation

        coverDepthMapFile = self.coverDepthMapFile
        waterDistanceMapFile = self.waterDistanceMapFile
        powerDistanceMapFile = self.powerDistanceMapFile
        railDistanceMapFile = self.railDistanceMapFile
        roadDistanceMapFile = self.roadDistanceMapFile
        railDistanceToPortMapFile = self.railDistanceToPortMapFile
        roadDistanceToPortMapFile = self.roadDistanceToPortMapFile
  
        stateIdsMapFile  = self.stateIdsMapFile

        if doTimer:
          #startMapsTime = time.time()
          timesList.append(["startMaps", time.time()])

        stateIdsMap =  self.LoadMap(stateIdsMapFile)
        self.ApplyOutputBoundingBoxMask(stateIdsMap) 

        countryIndxs = stateIdsMap >= 0.0

        ###################################
        if(coverDepthMapFile):
          coverMap =  self.LoadMap(coverDepthMapFile)
          coverMap += self.coverOffset
          # bound cover map
          coverMap[coverMap > coverMax] = coverMax
          coverMap[coverMap < coverMin] = coverMin
  
        else:
          # nonsense values
          coverMap =  self.LoadMap(waterDistanceMapFile) 
          XX,YY = np.mgrid[0:1:coverMap.shape[0]*1j,0:1:coverMap.shape[1]*1j]    #random.rand(coverMap.shape[0],coverMap.shape[1]) # fixme
          coverMap[:] = 25 # + 1.*XX - 2*YY

          coverMap = 25  + 1.*XX - 2*YY

        mineValueMap = np.zeros(coverMap.shape)
        halfGradeMineValueMap = np.zeros(coverMap.shape)
        doubleGradeMineValueMap = np.zeros(coverMap.shape)

        taxRelief = 1.0 - 0.3  # 30% tax relief from incometax

        ####
        for i in range(len(stateStrs)):
          stateIndxs = stateIdsMap == stateIds[i]
          mineValueFunc = interpolate.interp1d(coverDepths,stateMineValues[stateStrs[i]])
          mineValueMap[stateIndxs] = mineValueFunc(coverMap[stateIndxs])
          halfGradeMineValueFunc = interpolate.interp1d(coverDepths,stateMineValuesHalfGrade[stateStrs[i]])
          halfGradeMineValueMap[stateIndxs] = halfGradeMineValueFunc(coverMap[stateIndxs])
          doubleGradeMineValueFunc = interpolate.interp1d(coverDepths,stateMineValuesDoubleGrade[stateStrs[i]])
          doubleGradeMineValueMap[stateIndxs] = doubleGradeMineValueFunc(coverMap[stateIndxs])

        if doTimer:
          #mineValueFuncTime = time.time()
          timesList.append(["mineValue", time.time()])

        # Water Expenses
        distanceToWater =  self.LoadMap(waterDistanceMapFile)
        distanceToWaterCost =  taxRelief* theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalWaterExpenses(theProblemManager,distanceToWater[countryIndxs])
        mineValueMap[countryIndxs] -=  distanceToWaterCost
        halfGradeMineValueMap[countryIndxs] -= distanceToWaterCost
        doubleGradeMineValueMap[countryIndxs] -=  distanceToWaterCost
        del distanceToWater
        del distanceToWaterCost

        if doTimer:
          #waterTime = time.time()
          timesList.append(["water", time.time()])

        # Power Expenses
        # ** here we could Pre-calculate power supply costs to increase speed slightly**
        distanceToPower =  self.LoadMap(powerDistanceMapFile)
        distanceToPowerCost = taxRelief* theProblemManager.theMineDataManager.theInfrastructureManager.CalculateRegionalPowerExpenses(theProblemManager,distanceToPower[countryIndxs])
        mineValueMap[countryIndxs] -=  distanceToPowerCost
        halfGradeMineValueMap[countryIndxs] -= distanceToPowerCost
        doubleGradeMineValueMap[countryIndxs] -=  distanceToPowerCost
        del distanceToPower
        del distanceToPowerCost


        if doTimer:
          #powerTime = time.time()
          timesList.append(["power", time.time()])

        # Transportation Expenses
        distanceToRoad =  self.LoadMap(roadDistanceMapFile)
        distanceToRail =  self.LoadMap(railDistanceMapFile)
        roadTransportationDistance =  self.LoadMap(roadDistanceToPortMapFile)
        railTransporationDistance =  self.LoadMap(railDistanceToPortMapFile)


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
          #transportTime = time.time()
          timesList.append(["transport", time.time()])



        mineValueMap[stateIdsMap< 0.0] = np.nan  

        # interpolate to find the grade scale that gives value = 0

        c = mineValueMap[countryIndxs]
        a = ( doubleGradeMineValueMap[countryIndxs]+ 2*halfGradeMineValueMap[countryIndxs]-3*c )/1.5
        b = doubleGradeMineValueMap[countryIndxs] - a - c

        #
        breakEvenFactor = np.zeros(coverMap.shape)

        breakEvenFactor[countryIndxs] = 1+(-b + np.sqrt( b*b-4*a*c ) )/(2*a+1e-64)

        breakEvenFactor[stateIdsMap< 0.0] = np.nan  
        
        # correct for break even factor
        if(gradeAdjustment != 1.0):
          #print "gradeAdjustment",gradeAdjustment
          breakEvenFactor *= gradeAdjustment
        
        
        # clip to bounding box
        if(len(self.dataBoundingBox) == 4 and len(self.outputBoundingBox) == 4):
          breakEvenFactor = self.ClipToOutputBounds(breakEvenFactor)

        if doTimer:
          #transportTime = time.time()
          timesList.append(["end", time.time()])

        ######################################################################

        #pl.figure()
        #pl.imshow(mineValueMap,origin="lower")
        #pl.colorbar()



        if(doPlots):
        
          pl.figure()
          pl.plot(breakEvenFactor[500,:])
          
          pl.figure()
          pl.imshow(breakEvenFactor,origin="lower")
          pl.colorbar()
          pl.show()    
        else:
          filename = theProblemManager.outputPrefix+"."+theProblemManager.outputType
          self.SaveMap(filename, breakEvenFactor,theProblemManager.recordRange)

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
              print labels[i], dts[i], dts[i]-dts[i-1]

    
        return 
