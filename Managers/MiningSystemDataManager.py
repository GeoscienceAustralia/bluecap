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

from scipy.optimize import brentq

# Common
from Common.Common import Todo,Fixme

# IO
from IO.XML import HasChild,GetChild,AddChild
from IO.XML import GetAttributeString,GetAttributeValue, SetAttributeString


#Units
from Units.UnitManager import UnitManager

#Functions
from Functions.FunctionManager import FunctionManager


class MiningSystemDataManager():
    """
        This class holds information representing the mining system. 
        Note it is a subset of the Mine data manager, which holds information about the mining operation, processing/milling and infrastructure.
    """
    def __init__(self):
      """
      Create an empty mining system data manager and default variables. 
      """
      self.miningMethod = ""
      self.mineLife = 0
      self.mineCapacity = 0.0
      self.mineOreProductionCapacity = 0.0
      self.orebodyMass = 0.0
      self.oreMined = np.array([0.0])
      self.dilution = 0.05
      self.wasteMined = np.array([0.0])
      self.miningCapex = np.array([0.0])
      self.miningOpex = np.array([0.0])
      self.depths = np.array([0.0])
      self.rampGrade = 0.1
      self.alpha = 40*np.pi/180   # pit slope angle
      self.mineStartupTime = 0  # in years
      self.actualMineStartupTime = 0.0
      self.miningMethod = "OC"
      
      theUnitManager =  UnitManager()
      self.ugHaulCostPerDepth = theUnitManager.ConvertToBaseUnits("0.0073 AUD/tonne/m")

    def ParseXMLNode(self, miningSystemNode):
      """
      Generate mining system data from xml tree node.
       - Note this may not be needed in many cases as the mining system is derived from the orebody shape  
      """
    
      if(HasAttribute(miningSystemNode,"alpha")):  
        self.alpha = GetAttributeValue(miningSystemNode,"alpha")
      
      return miningSystemNode
    
    
      
    def UpdateOpenPitDepth(self,w,l,d,excavatedVol,alpha):
      cotAlpha = 1.0/np.tan(alpha)
      func = lambda newD : w*l*(newD-d) + (w+l)*cotAlpha * (newD **2 - d**2) + np.pi/3 * cotAlpha**2 *(newD**3 - d**3) - excavatedVol
  
      dmax = d+excavatedVol/float(w*l)   # might be able to improve bound by using 2nd order terms
  
      #print func(d), func(dmax), cotAlpha
  
      theNewD = brentq(func,d,dmax)  # using current d as lower bound
  
      return theNewD


    def UpdateSlopedOpenPitDepth(self,w,l,d,excavatedVol,alpha,beta):
      cotAlpha = 1.0/np.tan(alpha)  # just store grade instead?
      cotBeta = 1.0/np.tan(beta)
      func = lambda newD : w*l*newD - (w-(cotBeta-cotAlpha)*(newD-d))*l*d\
                         + (w+0.5*l)*cotAlpha * newD **2 - (0.5*l+w-(cotBeta-cotAlpha)*(newD-d))*cotAlpha*d**2 \
                         + np.pi/6 * cotAlpha**2 *(newD**3 - d**3) - excavatedVol
  
      dmax = d+excavatedVol/float(w*l)   # might be able to improve by using 2nd order
  
      #print func(d), func(dmax), cotAlpha
  
      theNewD = brentq(func,d,dmax)
  
      return theNewD
  

    def UpdateOpenPitDepthAndOreFraction(self,w,l,d,excavatedVol,alpha):
      dd = self.UpdateOpenPitDepth(w,l,d,excavatedVol,alpha)
      oreFraction =  (dd-d)*w*l/excavatedVol
      return dd,oreFraction
  
    def UpdateSlopedOpenPitDepthAndOreFraction(self,w,l,d,excavatedVol,alpha,beta):
      dd = self.UpdateSlopedOpenPitDepth(w,l,d,excavatedVol,alpha,beta)
      oreFraction =  (dd-d)*w*l/excavatedVol
      return dd,oreFraction   
    
    ###
    
    def OpenPitExcavatedVolume(self,w,l,d,newD,alpha):
      """
      Excavated volume between depths d and newD for an open pit with an orebody of width w, length l and pit slope alpha 
      """
      cotAlpha = 1.0/np.tan(alpha)
      excavatedVolume = w*l*(newD-d) + (w+l)*cotAlpha * (newD **2 - d**2) + np.pi/3 * cotAlpha**2 *(newD**3 - d**3)
      return excavatedVolume

    def OpenPitOreFraction(self,w,l,d,newD,alpha):
      """
      Fraction of ore excavated between depths d and newD for an open pit with an orebody of width w, length l and pit slope alpha 
      """
      orefrac = w*l*(newD-d)/(self.OpenPitExcavatedVolume(w,l,d,newD,alpha)+1e-64)
      return orefrac


    def OpenPitMarginalOreFraction(self,w,l,newD,alpha):
      """
      Fraction of ore excavated at depth newD for an open pit with an orebody of width w, length l and pit slope alpha 
      """
      cotAlpha = 1.0/np.tan(alpha)
      borderWidth = newD*cotAlpha
      excavatedArea = w*l + (w+l)*borderWidth  + np.pi * borderWidth**2
      orefrac = w*l/excavatedArea
      return orefrac



    def RequiredTotalCapacity(self,oreProductionRate, w,l,d,newD,alpha):
      """
      Mining capacity (total material moved) to maintain a given ore production rate between depths d and newd for an open pit with an orebody of width w, length l and pit slope alpha 
      """
      oreFrac = self.OpenPitOreFraction(w,l,d,newD,alpha)
      totalCapacity = oreProductionRate/(oreFrac + 1e-64)
      return totalCapacity
  
    def MarginalCostAtDepthPerUnitOre(self,costAtDepthPerUnitMaterial, w,l,newD,alpha):
      """
      Marginal cost per unit ore at depth newd for an open pit with an orebody of width w, length l and pit slope alpha 
      """
      marginalOreFrac = self.OpenPitMarginalOreFraction(w,l,newD,alpha)
      costAtDepthPerUnitOre = costAtDepthPerUnitMaterial/(marginalOreFrac + 1e-64)
      return costAtDepthPerUnitOre


    def FindOPMarginalOreCostAtDepth(self,oreProductionRate, w,l,d,newD,alpha,materialCostFunc ):
      """
      Find marginal cost at depth per unit ore from a material cost function based on the mine capacity in terms of total material mined
      """
      totalCapacity = self.RequiredTotalCapacity(oreProductionRate, w,l,d,newD,alpha)
      theUnitManager =  UnitManager()  
      costAtDepthPerUnitMaterial = materialCostFunc.f([totalCapacity *theUnitManager.ConvertTo("1e6 tonne") ])
      marginalcostAtDepth = self.MarginalCostAtDepthPerUnitOre(costAtDepthPerUnitMaterial, w,l,newD,alpha)
      return marginalcostAtDepth

    def FindOpUgSwitchingDepth(self,oreProductionRate, w,l,d,maxD,alpha,materialCostFunc,ugMaterialCost, ugMaterialCostPerTonMDepth ):
      """
      Find switching depth between an open pit and underground mine
      """
      func = lambda newD : ugMaterialCost + ugMaterialCostPerTonMDepth*(newD-60)  - self.FindOPMarginalOreCostAtDepth(oreProductionRate, w,l,d,newD,alpha,materialCostFunc ) 
  
      theSwitchingDepth = maxD
      if func(theSwitchingDepth) < 0:
        theSwitchingDepth = brentq(func,d,maxD)
  
      return theSwitchingDepth
  
    
    
    
    ###
    
    ###
    
    
    
     
    def DetermineMiningSystem(self, problemManager, mineManager):
      """
      Use available data to determine the most likely mining method/capacity/opex etc.
      """
    
      
      flatPlunge = 20 * np.pi/180.
      steepPlunge = 55 * np.pi/180.
      
      theUnitManager =  UnitManager()
      narrowWidth = theUnitManager.ConvertToBaseUnits("10m")
      thickWidth = theUnitManager.ConvertToBaseUnits("30m")
      
      # Thick > 30 m
      # Intermediate 10-30m
      # Narrow < 10 m
      
      # Plunge
      # Flat < 20 degrees
      # Intermediate  20-55 degrees
      # Steep > 55 degrees  
      
      doUnderground = self.miningMethod == "UG"
      
      if (doUnderground):
          if (mineManager.theOreBody.dip < flatPlunge ):  # see selection process for hard rock mining by Carter
            self.miningMethod = "UG_RP"
            # Room and pilar
            # "Typically flat and tabular"
          elif (mineManager.theOreBody.width < narrowWidth ):
            # Cut and fill
            self.miningMethod = "UG_CF"
      
          elif (mineManager.theOreBody.dip > steepPlunge and mineManager.theOreBody.width > thickWidth and mineManager.theOreBody.height > mineManager.theOreBody.width ):
            # Block cave
            self.miningMethod = "UG_BC"
          else:
            # Stoping
            self.miningMethod = "UG_ST"
      
      
      self.CalculateMineCapacity(problemManager, mineManager)
      self.CalculateMineLife(mineManager)
      self.CalculateOreMined(mineManager)
      
      
      self.CalculateMiningCapex(mineManager)
      self.CalculateMiningOpex(mineManager)
      
      return self

      
    def CalculateMineCapacity(self,problemManager, mineManager):
      """
      Determine the maximum annual extraction for the mine using Taylor's rule
      """
      theUnitManager =  UnitManager()
      
      orebodyMass =  mineManager.theOreBody.CalculateDepositMass()*theUnitManager.ConvertTo("tonne")
      
      
      #self.mineCapacity =  problemManager.theFunctionManager.Evaluate("TaylorsRule_" + self.miningMethod, [orebodyMass])
      
      theFunctionManager =  FunctionManager()
      
      if(self.miningMethod == "OCUG"):
        taylorsRuleFunc = theFunctionManager.GetFunction("TaylorsRule_" + self.miningMethod[:2])
      else:
        taylorsRuleFunc = theFunctionManager.GetFunction("TaylorsRule_" + self.miningMethod)
      daysPerYear = 350 
      #Todo("Check operating days per year - assuming 350 days")
      self.mineCapacity = daysPerYear * taylorsRuleFunc.f( [orebodyMass] )* theUnitManager.ConvertToBaseUnits("tonne")
      self.mineOreProductionCapacity = self.mineCapacity
      
      print "Mine ore capacity from Taylor's rule in Mt/year", self.mineCapacity* theUnitManager.ConvertTo("1e6 tonne")
      #print "Mine ore capacity from Taylor's rule in Mt/day", self.mineCapacity* theUnitManager.ConvertTo("1e6 tonne")/350.
      
      
      return self.mineCapacity
      
    
    def CalculateMineLife(self, mineManager):
      """
      Determine the life of the mine
      """
      # fixme - taylor's rule may need to be tweaked for open cut 
      #       - original tracks amount of ore produced (based on processing cost)
      #       - therefore assume that mine material moved is sufficient to meet this on average
      
      # current plan - calculate break even SR based on Open cut capacity
      #              - determine depth where break even SR is reached
      #              - mine to depth at open cut capacity based on average SR to reach that depth
      #              - continue mining at processing capacity
      
      self.orebodyMass =  mineManager.theOreBody.CalculateDepositMass()
      
      self.mineLife = self.orebodyMass/(self.mineOreProductionCapacity+1e-64) # rampup/rampdown are not accounted for
      self.mineLife = int( np.ceil(self.mineLife) )  # round up to years
      
    
      theUnitManager =  UnitManager()
      
      orebodyMass =  mineManager.theOreBody.CalculateDepositMass()*theUnitManager.ConvertTo("tonne")
      print "orebodyMass in 1e6 tonne", orebodyMass/1e6
      
      
      # undergound
      rampLength = mineManager.theOreBody.cover *( 1. + 1./self.rampGrade**2)**0.5
      rampVolume = 25*theUnitManager.ConvertToBaseUnits("m^2")*rampLength
      
      
      # opencut
      if(self.miningMethod[:2] == "OC"):
        w = mineManager.theOreBody.width 
        l = mineManager.theOreBody.length
        d = mineManager.theOreBody.cover
        rampVolume = self.OpenPitExcavatedVolume(w,l,0,d,self.alpha)
        
      
      overburdenDensity = mineManager.theOreBody.specificDensity*1000*theUnitManager.ConvertToBaseUnits("kg/m^3")
      self.actualMineStartupTime = (overburdenDensity*rampVolume)/(self.mineCapacity+1e-64)
      
      self.mineStartupTime = np.max([int( np.ceil(self.actualMineStartupTime) ),1])  # round up to years
    
      

      self.mineLife += self.mineStartupTime
      
      self.miningCapex = np.zeros(self.mineLife)
      self.miningOpex = np.zeros(self.mineLife)
      self.depths = np.zeros(self.mineLife)
      
      
      print "mineLife", self.mineLife
      print "mineStartupTime", self.mineStartupTime
      
      return self.mineLife
      
    def CalculateOreMined(self, mineManager):
      
      self.materialMined = np.zeros(self.mineLife)
      self.oreMined = np.zeros(self.mineLife)
      self.wasteMined = np.zeros(self.mineLife)
      self.depths = np.zeros(self.mineLife)
      
      self.materialMined[0] =self.mineCapacity*self.actualMineStartupTime
      
      self.materialMined[self.mineStartupTime:-1] = self.mineCapacity   # constant in all years but last
      self.materialMined[-1] =  self.orebodyMass - self.mineCapacity*(self.mineLife-self.mineStartupTime-1)  # remainder in last year
  
      # assuming all material mined is ore in underground mines
      
      if(self.miningMethod[:2] == "OC" ):  # open cut
    
        self.oreMined = np.array(self.materialMined)
        self.oreMined[:self.mineStartupTime] = 0.0
        #self.wasteMined =  (1 - oreFraction)*self.oreMined
        #self.oreMined *=  oreFraction
        
        theUnitManager =  UnitManager()
        theFunctionManager =  FunctionManager()
        
        ugMiningOpexFunc = theFunctionManager.GetFunction("MiningOpex_UG_ST")  # underground stoping assumed as alternative mining method
        ugOpexPerTonne =  ugMiningOpexFunc.f( [self.mineCapacity *theUnitManager.ConvertTo("1e6 tonne")] ) 
        ugOpexPerTonneMDepth = self.ugHaulCostPerDepth*theUnitManager.ConvertTo("AUD/tonne/m")
        
        
        ocMiningOpexFunc = theFunctionManager.GetFunction("MiningOpex_OC")  
        
        oreProductionRate = self.mineOreProductionCapacity # taylor's rule gives average ore produced not total material
        
        minD = mineManager.theOreBody.cover
        maxD = mineManager.theOreBody.cover + mineManager.theOreBody.height
        w = mineManager.theOreBody.length
        l = mineManager.theOreBody.width
        alpha = self.alpha  
        
        switchingDepth = self.FindOpUgSwitchingDepth(oreProductionRate, w,l,minD,maxD,alpha,ocMiningOpexFunc,ugOpexPerTonne,ugOpexPerTonneMDepth )+1e-3 # add a mm to prevent roundoff error
        
        
        print "OC/UG Switching Depth",switchingDepth
        print "minD",minD
        print "maxD",maxD
        
        if (switchingDepth > minD):
          
          maxOCdepth = np.minimum(switchingDepth,maxD)
          
          ocTotalMaterialProductionRate = self.RequiredTotalCapacity(oreProductionRate, w,l,minD,maxOCdepth,alpha)  # total material produced (mass Rate)
          self.mineCapacity = ocTotalMaterialProductionRate
          rho = mineManager.theOreBody.orebodyMass/mineManager.theOreBody.orebodyVolume
          ocTotalMaterialProductionVolume = ocTotalMaterialProductionRate/ rho # assumes ore body density = material density
          
          
          if (switchingDepth >= maxD):
          
            print "Open cut only"
            # all open cut
            self.materialMined[self.mineStartupTime:] = ocTotalMaterialProductionRate
            
            # ore mined varies over the years (but on average is equal to the mine ore production capacity) - this may not be needed. 
            d = minD
            
            self.depths[self.mineStartupTime-1] = minD
            for year in range(self.mineStartupTime,self.mineLife-1):
              dd = self.UpdateOpenPitDepth(w,l,d,ocTotalMaterialProductionVolume,alpha)
              self.oreMined[year] = ( (dd-d)/mineManager.theOreBody.height )*mineManager.theOreBody.orebodyMass
              d = dd
              self.depths[year] = dd
            
            self.oreMined[-1] =  mineManager.theOreBody.orebodyMass - np.sum( self.oreMined[:-1] )  # add any remaining ore to the last year ore
            self.depths[-1] = maxD
            
            
            fractionOre = self.OpenPitOreFraction(w,l,dd,maxD,self.alpha)
            
            
            self.materialMined[-1] =  self.oreMined[-1]/fractionOre
            print "self.materialMined[-1]",self.materialMined[-1],dd,maxD

            
          else:
          
            print "Mixed opencut and underground"
            
            self.miningMethod = "OCUG"
            self.ugMaterialMined = np.zeros(self.mineLife)
            self.ugOpexPerTonne = ugOpexPerTonne   # ugly - should store both above and below for all mining methods
            
            switchingTime = ( (switchingDepth-minD)/mineManager.theOreBody.height)*mineManager.theOreBody.orebodyMass/oreProductionRate  
            switchingTime += self.mineStartupTime
            theSwitchingYear = int( np.floor(switchingTime) )
          
            print "Op/UG Switching Time,year",switchingTime,theSwitchingYear
          
            print "OC opex",ocMiningOpexFunc.f([self.mineCapacity*theUnitManager.ConvertTo("1e6 tonne")])
            print "UG opex",ugOpexPerTonne
            print "OC capacity (Mt)",self.mineCapacity*theUnitManager.ConvertTo("1e6 tonne")
            print "UG capacity (Mt)",self.mineOreProductionCapacity*theUnitManager.ConvertTo("1e6 tonne")
            
            firstUGYear = theSwitchingYear+1
            if(firstUGYear < self.mineLife):
              self.materialMined[firstUGYear:] = self.oreMined[firstUGYear:] # already true
              self.ugMaterialMined[firstUGYear:] = self.oreMined[firstUGYear:] 
              
              self.depths[theSwitchingYear:] = np.linspace(switchingDepth,maxD,self.mineLife-theSwitchingYear)
            
            if(theSwitchingYear < self.mineLife):
              self.materialMined[theSwitchingYear] = ocTotalMaterialProductionRate*(switchingTime-theSwitchingYear) +  self.mineOreProductionCapacity*(1 -  switchingTime+theSwitchingYear)
              #self.oreMined[theSwitchingYear] = self.mineOreProductionCapacity*(switchingTime-theSwitchingYear)   # we will add OC contribution later
            self.materialMined[:theSwitchingYear] = ocTotalMaterialProductionRate  # oc material production up to time of switch
          
            # ore mined varies over the years (but on average is equal to the mine ore production capacity) - this may not be needed. 
            d = minD
            self.depths[self.mineStartupTime-1] = minD
            for year in range(self.mineStartupTime,theSwitchingYear):
              dd = self.UpdateOpenPitDepth(w,l,d,ocTotalMaterialProductionVolume,alpha)
              self.oreMined[year] = ( (dd-d)/mineManager.theOreBody.height )*mineManager.theOreBody.orebodyMass
              #print "dd",dd,ocTotalMaterialProductionVolume, (dd-d)*w*l, alpha
              d = dd
              self.depths[year] = dd
            
            self.oreMined[theSwitchingYear] =  self.mineOreProductionCapacity*(switchingTime-self.mineStartupTime) - np.sum( self.oreMined[:theSwitchingYear] )  # add any remaining oc ore to the switching year ore
            
            self.ugMaterialMined[theSwitchingYear] = self.mineOreProductionCapacity*(1-switchingTime+theSwitchingYear)   # ug component
            self.oreMined[theSwitchingYear] +=  self.ugMaterialMined[theSwitchingYear]  # ug component
            self.depths[theSwitchingYear] = switchingDepth  # nqr - should be depth at end of year.
            self.depths[-1] = maxD
          # account for 
          
        else:
          print "Open cut mining is less economic than underground mining"
          self.materialMined[:] = 0.0   # ug is better at all depths
          self.oreMined[:] = 0.0
        
      else:
        # assuming all material mined is ore in underground mines
        
        minD = mineManager.theOreBody.cover
        maxD = mineManager.theOreBody.cover + mineManager.theOreBody.height
        
        self.oreMined[self.mineStartupTime:] = self.materialMined[self.mineStartupTime:]
        print self.mineStartupTime-1, self.mineLife+1-self.mineStartupTime,self.mineLife
        print len(self.depths[self.mineStartupTime-1:]),len(self.depths)
        self.depths[self.mineStartupTime-1:] = np.linspace(minD,maxD,self.mineLife+1-self.mineStartupTime )
    
      self.wasteMined =  self.materialMined - self.oreMined

      return self.oreMined
      
    def CalculateMiningCapex(self, mineManager):

      self.miningCapex = np.zeros(self.mineLife)
      
      theFunctionManager =  FunctionManager()
      theUnitManager =  UnitManager()
      
      if(self.miningMethod == "OCUG"):
        miningCapexFunc = theFunctionManager.GetFunction("MiningCapex_" + self.miningMethod[:2])
      else:
        miningCapexFunc = theFunctionManager.GetFunction("MiningCapex_" + self.miningMethod)
        
      
      self.miningCapex[0] =  miningCapexFunc.f( [self.mineCapacity*theUnitManager.ConvertTo("1e6 tonne")] )
      
      #print "self.miningCapex", self.miningCapex
      
      return self.miningCapex
      
    def CalculateMiningOpex(self, mineManager):
      #  no distinction between opex and sustaining capital (adjust tax calculation)
      
      theFunctionManager =  FunctionManager()
      theUnitManager = UnitManager()
      
      
      miningOpexFunc = theFunctionManager.GetFunction("MiningOpex_" + self.miningMethod[:2])  #  NB. just opencut or underground
      opexPerTonne =  miningOpexFunc.f( [self.mineCapacity *theUnitManager.ConvertTo("1e6 tonne")] )
      ugOpexPerTonnePerMDepth = self.ugHaulCostPerDepth*theUnitManager.ConvertTo("AUD/tonne/m")
      
      
      print "Mining method: ", self.miningMethod
      print "Mining capacity (Mt material)", self.mineCapacity  *  theUnitManager.ConvertTo("1e6 tonne")
      print "Mining Opex Per Tonne ", opexPerTonne
      

      if(self.miningMethod == "OCUG"):
      
        #print "self.ugOpexPerTonne",self.ugOpexPerTonne
        self.miningOpex = (self.materialMined - self.ugMaterialMined)*theUnitManager.ConvertTo("tonne")   # tonnes mined per year
        self.miningOpex *= opexPerTonne
        
        self.miningOpex += self.ugMaterialMined*theUnitManager.ConvertTo("tonne")*( self.ugOpexPerTonne+  ugOpexPerTonnePerMDepth*(self.depths-60) )
     
      else:
        self.miningOpex = self.materialMined*theUnitManager.ConvertTo("tonne")   # tonnes mined per year
        if(self.miningMethod[:2] == "UG"):
           self.miningOpex *= opexPerTonne +  ugOpexPerTonnePerMDepth*(self.depths-60)
        else:
           self.miningOpex *= opexPerTonne
           
      # pre-production expenses are capitalized
      self.miningCapex[:self.mineStartupTime] = self.miningOpex[:self.mineStartupTime] 
      self.miningOpex[:self.mineStartupTime] = 0.0

      return self.miningOpex

    def WriteXMLNode(self, node):
      """
      Write mining system data to xml node
      """
      # Mine data
      SetAttributeString(node,"method",self.miningMethod)
      SetAttributeString(node,"capacity",self.mineCapacity)
      SetAttributeString(node,"oreProductionCapacity",self.mineOreProductionCapacity)
      SetAttributeString(node,"life",self.mineLife)
      SetAttributeString(node,"rampGrade",self.rampGrade)
      
      return node
      
      
      