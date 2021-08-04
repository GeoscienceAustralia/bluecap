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

import pylab as pl  # used to show output
import numpy as np

from Functions.Astar import AStar

from scipy.optimize import minimize, differential_evolution

# Managers
from Managers.SolverManager import SolverManager
from Managers.SolverManager import SolverFactory


# IO
from IO.XML import HasChild,GetChild,AddChild
from IO.XML import GetAttributeString,GetAttributeValue,GetAttributeValueOrDefault, SetAttributeString, GetAttributeStringOrDefault
from IO.outputs import write_point_geojson

# Units
from Units.UnitManager import UnitManager


# common
from Common.GeoUtils import LatLonDegreeToKMCoords


"""
def latLonDegreeToKMCoords(lat,lon):
  avlat = lat*np.pi/180.
  yCoord = 1e-3*( 111132.954*(lat+28) - 559.822 * np.sin( 2 * avlat )*180/(np.pi*2.0) + 1.175 * np.sin( 4 * avlat)*180/(np.pi*4.0) - 0.0023*np.sin( 6 * avlat)*180/(np.pi*6.0) );
  xCoord = 1e-3*( (111412.84 * np.cos( avlat ) - 93.5 * np.cos(3*avlat) + 0.118 * np.cos(5*avlat) ) * (lon - 134 ) )
  return yCoord,xCoord

from osgeo import ogr, osr

def latLonDegreeToKMCoords(lat,lon):
  source = osr.SpatialReference()
  source.ImportFromEPSG(4326)
  
  target = osr.SpatialReference()
  target.ImportFromEPSG(3112)
  
  transform = osr.CoordinateTransformation(source, target)
  
  point = ogr.CreateGeometryFromWkt("POINT (%f %f)" % (lon, lat))
  point.Transform(transform)
  
  return 1e-3*point.GetY(), 1e-3*point.GetX()
"""

# helper class for Processing Facility Locator
class DepositDataStruct():

    def __init__(self,name,equivalentAnnuity):
      self.name = name
      self.equivalentAnnuity = equivalentAnnuity
      self.roadTransportCostPerKm = None
      self.roadCapitalCostPerKm = None
      self.startTime  = None
      self.mineLife = None
      
    def CalculateTransportCosts(self,problemManager,mineDataManager):
    
      discountRate = mineDataManager.theEconomicDataManager.discountRate
      self.mineLife =  mineDataManager.theMiningSystem.mineLife
      
      # discountedOreMass must be in base units
      discountFactors = (1.0+discountRate)**np.array(range( int( np.ceil(self.mineLife) ) ) )
      discountedOreMass = np.sum(mineDataManager.theMiningSystem.oreMined/discountFactors)
      
      theUnitManager = UnitManager()
      oneKm = theUnitManager.ConvertToBaseUnits("1 km")
      
      self.roadTransportCostPerKm = mineDataManager.theInfrastructureManager.RoadTransportOpex(oneKm,discountedOreMass)
    
      self.roadCapitalCostPerKm = mineDataManager.theInfrastructureManager.RoadConstructionCapex(oneKm)
      
      # discount for late start date
      
      startDiscountFactor = 1.0/((1.0+ discountRate)**self.startTime)
      self.roadTransportCostPerKm *= startDiscountFactor
      self.roadCapitalCostPerKm *= startDiscountFactor
      
      return self
      
      

class ProcessingFacilityLocator():
    typeStr = "ProcessingFacilityLocator"

    def __init__(self):
      """
      Create Empty Processing Facility Locator. 
      """
      
      ### Dummy resolution - nX, nY are updated base on output region  ###
      self.nX = 100
      self.nY = 100
      self.gridXScale = 1.0
      self.gridYScale = 1.0
      self.initialCostData = np.zeros([100,100])  #2.00*np.random.rand(100,100)+10.   #road data
      self.costData = None   #road data
      self.roadTransportCostPerKm = 0.5 # dummy value
      self.roadCapitalCostPerKm = 4.5  # dummy value
      self.heuristicCostPerDistance = self.roadTransportCostPerKm  
      
      self.targetX = 0.0
      self.targetY = 0.0
      self.activeOrebody = 0
      
      self.depositList = []
      self.plotting = False

    def CostFunction(self,x,y,nbrx,nbry):
      """Function used to determine the cost of a segment of road from x,y to nbrx,nbry."""
      cc = 0.0
      if( self.costData[int(x),int(y)] <= 0 ):  # road construction per km
        cc += self.roadCapitalCostPerKm
      
      if( self.costData[int(nbrx),int(nbry)] <= 0): # road construction per km
        cc += self.roadCapitalCostPerKm
      
      cc = cc * 0.5 # average
      
      cc += self.roadTransportCostPerKm ## contribution from transport per distance on road cells
      
      
      #cc = 0.5*(self.costData[int(x),int(y)] + self.costData[int(nbrx),int(nbry)] )
      #print cc
      
      dx = (x - nbrx)*self.gridXScale
      dy = (y - nbry)*self.gridYScale
      return cc*(dx**2 + dy**2)**0.5      
    
    

      
      
    def Run(self,problemManager,args):
        """Run the solver."""
    
        mineDataManager = problemManager.theMineDataManager
      
        # loop over orebodies and calculated NPV of each (less transport costs)
        # rank according to equivalent annuity  
      
        print("Calculating production capacity")
        print("#####")
        print("")
        totalMass = 0.0
        for name,orebody in mineDataManager.theOrebodies.items():
          
          if name is not "Active":
             print(name, orebody.latLong)
             totalMass += orebody.CalculateDepositMass()
        # convert to tonne 
        theUnitManager =  UnitManager()
        totalMass = totalMass*theUnitManager.ConvertTo("tonne")
        
        print("Total mass: ",totalMass)
        # setting mine type to stoping for the purpose of calculating the orebody life
        mineDataManager.theMiningSystem.miningMethod = "UG_ST"
        productionCapacity = mineDataManager.theMiningSystem.CalculateMineCapacityFromMass(totalMass)
        
        for name,orebody in mineDataManager.theOrebodies.items(): 
          if name is not "Active":
             mineDataManager.SetActiveOrebody(name)  # also sets active mine
             mine = mineDataManager.theMiningSystem
             mine.mineOreProductionCapacity =  productionCapacity
             mine.mineCapacity = productionCapacity    # NB this will be adjusted later
             
             # set mining method for individual deposits     
             mine.DetermineDepositMineType(problemManager,mineDataManager)
             
             # set mine life for individual deposits - NB production rate based on total from all deposits 
             mine.DetermineMineLifeProductionAndCosts(problemManager,mineDataManager)
             
        
        # problemManager.theMineDataManager.SetActiveOrebody(self.orebodyName)
        print("Calculating equivalent annuity")
        print("#####")
        print("")
        
        # Get npv and minelife and calculate equivalent annuity for each deposit from BTNCF
        
        for name,orebody in mineDataManager.theOrebodies.items():
          if name is not "Active":
          
             # Calculate npv (nb includes processing system capex which must be zeroed out)
             ######
             
             mineDataManager.SetActiveOrebody(name)
             
             # Processing Model - this is repeditive, but need to capture changing mine life
             mineDataManager.DetermineProcessingSystem(problemManager)
             # need to remove processing facility capex
             mineDataManager.ZeroProcessingCapex(problemManager)
      
             # G&A Model
             mineDataManager.CalculateGandAExpenses(problemManager)
      
             # Infrastructure Model - not called but must be zeroed 
             mineDataManager.ZeroInfrastructureCosts(problemManager) 
      
             # Calculate BTNCF and annuity
             ######
             mineDataManager.theEconomicDataManager.CalculateBeforeTaxCashFlow(problemManager,mineDataManager)
             npv = mineDataManager.theEconomicDataManager.CalculateBeforeTaxNPV(problemManager,mineDataManager)
             numYears = mineDataManager.theMiningSystem.mineLife
             
             equivalentAnnuity = mineDataManager.theEconomicDataManager.CalculateEquivalentAnnuity(npv,numYears)
             
             dd = DepositDataStruct(name, equivalentAnnuity)
             self.depositList.append(dd)
        
        # order orebodies according to equivalent annuities from largest to smallest
        self.depositList.sort(key=lambda xx: -xx.equivalentAnnuity) 
        #print orebodyEquivalentAnnuities
        
        # get start times and costs
        startTime = 0.0
        for deposit in self.depositList:
          deposit.startTime = startTime
          
          deposit.CalculateTransportCosts(problemManager,mineDataManager)
          
          print("deposit.roadTransportCostPerKm",deposit.roadTransportCostPerKm)
          print("deposit.roadCapitalCostPerKm",deposit.roadCapitalCostPerKm)
          startTime += deposit.mineLife
        
      
        ######
        print("Initializing orebody map")
        print("#####")
        print("")
        globalCoords = []

        minLatLong = np.array( [1000.,1000.])  # unlikely values
        maxLatLong = np.array( [-1000.,-1000.]) # unlikely values
        
        for deposit in self.depositList:
          name = deposit.name
          orebody = mineDataManager.theOrebodies[name]
          
          if name is not "Active":
            print(name, orebody.latLong)
            loc = orebody.latLong
            [yy,xx] = LatLonDegreeToKMCoords(loc[0],loc[1])
            minLatLong = np.minimum(minLatLong,loc)
            maxLatLong = np.maximum(maxLatLong,loc)
            globalCoords.append([xx,yy])
          
        
        globalCoords = np.array(globalCoords)
        
        # ~0.05 deg buffer around outer deposits
        minLatLong -= 0.5 #0.05
        maxLatLong += 0.5 #0.05
        
        ### set road locations 
        ## nb. currently done from distance to road (coud be cleaned up)
        regionalCalculationManager = problemManager.theRegionalCalculationManager
        distanceToRoad =  regionalCalculationManager.LoadMap(regionalCalculationManager.roadDistanceMapFile)
        distanceToRoad = regionalCalculationManager.ClipToOutputLatLon(distanceToRoad,[minLatLong[0],minLatLong[1],maxLatLong[0],maxLatLong[1]])
        [minY,minX] = LatLonDegreeToKMCoords(minLatLong[0],minLatLong[1])
        [maxY,maxX] = LatLonDegreeToKMCoords(maxLatLong[0],maxLatLong[1])
        
        rangeXY = np.array([maxX - minX,maxY - minY])
        
        """
        minX,minY = np.min(globalCoords,0)
        
        if (rangeXY[0] <self.nX):
          minX -= (self.nX-rangeXY[0])/2
        else:
          self.nX = int( np.ceil( rangeXY[0] ) )+10
          minX -= 5
        
        if (rangeXY[1] <self.nY):
          minY -= (self.nY-rangeXY[1])/2
        else:
          self.nY = int( np.ceil( rangeXY[1] ))+10
          minY -= 5
        """
        #self.nX = int( np.ceil( rangeXY[0] ))
        #self.nY = int( np.ceil( rangeXY[1] ))
        
        self.nY,self.nX = distanceToRoad.shape
        # fixme initialCostData (and cost data) needs to be renamed - indicates the presence (1 or more) or absence (0) of roads
        self.initialCostData = np.zeros([self.nX,self.nY],dtype=np.float)  # ie no roads
        
        
        minDistToRoad = 0.5*np.minimum(self.gridXScale,self.gridYScale)   #(self.gridXScale**2 + self.gridYScale**2)**0.5
        #print "minDistToRoad",minDistToRoad
        
        indxA = np.zeros(distanceToRoad.T.shape,dtype=np.bool)
        indxA[:,1:-1] = np.logical_and(distanceToRoad.T[:,1:-1] <= distanceToRoad.T[:,2:], distanceToRoad.T[:,1:-1] <= distanceToRoad.T[:,:-2])
        
        indxB = np.zeros(distanceToRoad.T.shape,dtype=np.bool)
        indxB[1:-1,:] = np.logical_and(distanceToRoad.T[1:-1,:] <= distanceToRoad.T[2:,:], distanceToRoad.T[1:-1,:] <= distanceToRoad.T[:-2,:])
        
        indx = np.logical_or(indxA,indxB)
        indx = np.logical_and(indx, distanceToRoad.T < minDistToRoad)
        
        self.initialCostData[indx] = 1.0  # indicates the presence of road
        #### fixme need to rename
        
        """
        pl.imshow(distanceToRoad)
        pl.figure()
        
        pl.imshow(self.initialCostData.T) 
        pl.show()
        pl.exit()
        """

        self.costData = np.array(self.initialCostData)
        #### fixme need to rename
        
        
        # convert to km scale
        depositCoords = np.array(globalCoords) #= np.array([[35,5],[20,10],[16,70],[19,35]])
        depositCoords[:,0] -= minX
        depositCoords[:,1] -= minY
        
        # convert to grid scale
        self.gridXScale = (maxX-minX)/(self.nX-1.0)
        self.gridYScale = (maxY-minY)/(self.nY-1.0)
        depositCoords[:,0] /= self.gridXScale
        depositCoords[:,1] /= self.gridYScale

        #processorCoords = np.array([56,60.])

        processorCoords = np.mean(depositCoords,0)  # initial guess

        #joinDepositsAndProcessingCenter(costData, depositCoords,processorCoords)

        ###
        
        print("Calculating processor location")
        print("#####")
        print("")
        
        def dummyFunc(pC,dC,dummy):
          return self.ProjectCostfunc(pC,dC)

        # res = minimize(dummyFunc,processorCoords,args=(depositCoords,None), method = "COBYLA",options={'rhobeg': 1.0})  
 
        bounds = [(1,self.nX-1), (1,self.nY-1)]

        print("Running differential_evolution")
        # fixme should be faster if initial guess is closer to first deposit and if the cutoff is closer to 1 cell. 
        res = differential_evolution(dummyFunc, bounds, args=(depositCoords,None))#, updating='immediate', polish=False) 
  
        #print "done"
        #print "self.gridXScale, self.gridYScale", self.gridXScale, self.gridYScale 
        ####
        self.costData = np.array(self.initialCostData)
        processorCoords = res.x
        print(processorCoords, res.message, res.success)
        self.JoinDepositsAndProcessingCenter(depositCoords,processorCoords)
        
        filename = problemManager.outputPrefix
        #write_point_geojson(filename+'_deposits', np.array(deposits), False)
        resultLonLat = np.array([[res.x[0],res.x[1]]])
        resultLonLat *= np.array([self.gridXScale,self.gridYScale])
        resultLonLat += np.array([minX,minY])
        write_point_geojson(filename, 'Processing Optimise', [1e3*globalCoords,1e3*resultLonLat], [{'id':'deposits','name':'Deposit(s)'},{'id':'processing_plant','name':'Processing plant'}], XY=True, in_epsg=3112)
        
        ## plotting
        if self.plotting:
          X,Y = np.meshgrid(np.linspace(0.,self.nX-1,self.nX),np.linspace(0.,self.nY-1,self.nY))
          X = X.reshape(np.prod(X.shape))
          Y = Y.reshape(np.prod(Y.shape))
          # print(self.nX,X)
          Z = np.array([dummyFunc([X[i],Y[i]],depositCoords,None) for i in range(len(X))]).reshape((self.nY,self.nX))

          pl.imshow(Z,origin="lower")
          pl.colorbar()
		  
          for depositX,depositY in depositCoords:
            pl.plot(depositX,depositY,"go", markersize=5)
		  
          pl.plot(res.x[0],res.x[1],"rd", markersize=5)
		  
          pl.axis("equal")
          
          pl.show()
      
   
    def HeuristicCostFunction(self,x,y):
      """Heuristic function used to estimate the cost of connecting x,y, to the target. It should this should always be less than or equal to the actual cost for best performance. """
      
      return self.heuristicCostPerDistance*( ((self.targetX-x)*self.gridXScale)**2. + ((self.targetY-y)*self.gridYScale)**2. )**0.5
   
      
    def FindPathAndUpdateCostData(self,initialX,initialY,targetX,targetY):
      """Finds the least cost path between an initial location and a target location. Returns the cost of the path."""
      
      # update self target
      self.targetX = targetX
      self.targetY = targetY
      
      # field size and rounded targets
      nX,nY = self.nX, self.nY
      tx,ty = int(np.round(self.targetX)),int(np.round(self.targetY))
        
      if(tx<0): tx = 0
      if(ty<0): ty = 0
      if(tx>nX-1): tx = nX-1
      if(ty>nY-1): ty = nY-1
          
      
      def heuristicFunc(x,y):
        return self.HeuristicCostFunction(x,y)

      def costFunc(x,y,nbrX,nbrY):
        return self.CostFunction(x,y,nbrX,nbrY)

      #path,pathCost = AStar(nX,nY,initialX,initialY,targetX,targetY, self.HeuristicCostFunction, self.CostFunction)
      path,pathCost = AStar(nX,nY,initialX,initialY,targetX,targetY, heuristicFunc, costFunc)  ## fixme -visited for debug
      


      for pp in path:
        if(self.costData[pp[0],pp[1]] == 0):
          self.costData[pp[0],pp[1]] =  2 #self.activeOrebody+1  ## <----<<<< updated road (labeled by when added)
      
      
      return pathCost


    def ResetCostData(self):
      """Revert the cost data back to its initial value"""
      self.costData =  np.array(self.initialCostData)
    
    def UpdateRoadCosts(self,roadTransportCostPerKm,roadCapitalCostPerKm ): 
      """Set the road transportation, capital costs and heuristic cost per distance."""
      self.roadTransportCostPerKm = roadTransportCostPerKm
      self.roadCapitalCostPerKm= roadCapitalCostPerKm
      self.heuristicCostPerDistance = roadTransportCostPerKm
      return self 

    def JoinDepositsAndProcessingCenter(self,depositCoords,processorCoords):
      """Connect deposits to the processing center and return the total cost."""
    
      # fixme - need to clean this up - merge deposit coords with deposit data struct
      self.ResetCostData()
      processorX,processorY = processorCoords
      totalCosts = 0.0
      depositId = 0 
      
      for depositX,depositY in depositCoords:
      
          self.UpdateRoadCosts(self.depositList[depositId].roadTransportCostPerKm,\
                               self.depositList[depositId].roadCapitalCostPerKm )
          
          pathCost = self.FindPathAndUpdateCostData(depositX,depositY,processorX,processorY)
          totalCosts += pathCost
          depositId += 1
          
      return totalCosts

    
    def ProjectCostfunc(self,processorCoords, depositCoords):
      """Determine the cost of connecting the processing center to each of the deposits. 
         Currently the same as JoinDepositsAndProcessingCenter, but with additional scope for new costs later."""
      totalCosts = self.JoinDepositsAndProcessingCenter(depositCoords,processorCoords)
      # add additional costs here
      print(processorCoords,totalCosts)
      return totalCosts

    @classmethod
    def CreateFromXML(cls,solverNode,problemManager):
      # this creates the class from an xml node by initiating an empty class and parsing the xml node
      rv = cls()
      rv.ParseXMLNode(solverNode,problemManager)
      return rv
        
    def ParseXMLNode(self, node,problemManager):
      """
      Generate Processing Facility Locator from xml tree node. 
      """
      self.nX = GetAttributeValueOrDefault(node,"nX",self.nX)
      self.nY = GetAttributeValueOrDefault(node,"nY",self.nY)
      


    def WriteXMLNode(self, node):
      """
      Write data to xml node
      """
      
      SetAttributeString(node,"nX",self.nX)
      SetAttributeString(node,"nY",self.nY)
      
      return node


    class Factory:
      """ The factory that creates the ProcessingFacilityLocator Solver. """
      def CreateFromXML(self,xmlNode,problemManager): return ProcessingFacilityLocator.CreateFromXML(xmlNode,problemManager)



# Factory Registrator
SolverFactory.AddFactory(ProcessingFacilityLocator.typeStr,ProcessingFacilityLocator.Factory())



