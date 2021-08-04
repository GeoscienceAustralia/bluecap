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



# IO
from IO.XML import HasChild,GetChild,GetChildren,AddChild
from IO.XML import GetAttributeString,GetAttributeValue, SetAttributeString, HasAttribute

# Functions
from Functions.Royalties import AustralianRoyalties

# Managers

class EconomicDataManager():

    def __init__(self):
      """
      Create an empty economic data manager and default variables. 
      """
      
      self.GandAOpex = np.array([0.0])
      self.GandAFraction = 0.14   # measures proportion of expenses dedicated to G&A
      self.commodityPrices = {}      # commodity prices used in calculation
      self.referenceCommodityPrices = {} # reference prices 
      self.commodityPriceSigmas = {}  # std deviation in the log ratios of metal prices
      
      # NCF
      self.revenue = np.array([0.0])
      self.netOpex = np.array([0.0])
      self.netCapex = np.array([0.0])
      
      self.btNCF = np.array([0.0])
      
      self.state = "WA"
      self.royalties =  np.array([0.0])
      
      self.taxes =  np.array([0.0])
      self.incomeTaxRate = 0.30
      
      self.depreciation =  np.array([0.0])
      self.atNCF = np.array([0.0])
      
      # inflation
      self.inflation = 0.02
      
      # sustaining capex fraction of total sustaining costs
      self.sustainingCapexFraction = 0.2
      
      # NPV
      self.discountRate = 0.0
      self.btNPV = 0.0
      self.atNPV = 0.0

    def ParseXMLNode(self, economicDataNode):
      """
      Generate the Economic Data Manager data from xml tree node. 
      """
      
      self.discountRate = GetAttributeValue(economicDataNode,"discountRate")
      
      for child in GetChildren(economicDataNode,"Commodity"):
        name = GetAttributeString(child,"name")
        price = GetAttributeValue(child,"price")
        self.commodityPrices[name] = price
        self.referenceCommodityPrices[name] = price
        self.commodityPriceSigmas[name] = 0.0
        if(HasAttribute(child,"sigma")):
          self.commodityPriceSigmas[name] = GetAttributeValue(child,"sigma")
        
      if(HasAttribute(economicDataNode,"GandAFraction")):
        self.GandAFraction = GetAttributeValue(economicDataNode,"GandAFraction")
        
      if(HasAttribute(economicDataNode,"state")):
        self.state = GetAttributeString(economicDataNode,"state")
        
      if(HasAttribute(economicDataNode,"inflation")):
        self.inflation = GetAttributeValue(economicDataNode,"inflation")



    def WriteXMLNode(self, node):
      """
      Write the Economic Data Manager data to an xml node
      """
      
      # discount rate
      SetAttributeString(node,"discountRate",self.discountRate)
      SetAttributeString(node,"incomeTaxRate",self.incomeTaxRate)
      SetAttributeString(node,"GandAFraction",self.GandAFraction)
      SetAttributeString(node,"state",self.state)
      
      # price data
      
      for name,price in self.referenceCommodityPrices.items():
        commNode = AddChild(node,"Commodity")
        SetAttributeString(commNode,"name",name)
        SetAttributeString(commNode,"price",price)
        SetAttributeString(commNode,"sigma",self.commodityPriceSigmas[name])
      
      
      # inflation
      SetAttributeString(node,"inflation",self.inflation)
      
      
      SetAttributeString(node,"btNPV",self.btNPV)
      SetAttributeString(node,"atNPV",self.atNPV)
      
      return node
    
    def SetState(self, stateStr):
      """Set the state of the economic data manager."""
      self.state = stateStr
    
    def ResetCommodityPrices(self):
      """Reset commodity prices to their reference values (used in random walk simulations)."""
      for name,price in self.referenceCommodityPrices.items():
        self.commodityPrices[name] = price
    
    def GenerateRandomWalkCommodityPrices(self,years):
      """Create a commodity price series for each commodity using a random walk simulator."""
      for name,price in self.referenceCommodityPrices.items():
        randomPrices = np.zeros(years)
        sigma = self.commodityPriceSigmas[name]
        randomPrices[0] = price
        for i in range(1,years):
           ratio = np.exp( sigma * np.random.normal() )
           randomPrices[i] = randomPrices[i-1]*ratio
        self.commodityPrices[name] = randomPrices
    
    def CalculateGandAExpenses(self, problemManager, mineDataManager):
      """
      Calculate the general and administrative expense component of the mine costs 
      Assumed equal to a fixed percentage of the overall mining, G&A, and processing costs.
      """
      self.GandAOpex = mineDataManager.theMiningSystem.miningOpex + mineDataManager.theMiningSystem.miningCapex \
                      + mineDataManager.theProcessingSystem.processingOpex + mineDataManager.theProcessingSystem.processingCapex
      self.GandAOpex *= self.GandAFraction/(1.0 - self.GandAFraction)

    def CalculateGandAExpensesFromManagerList(self, problemManager, managerList):
      """
      Calculate the general and administrative expense component of the mine costs 
      Assumed equal to a fixed percentage of the startup and sustaining costs of each manager in the list.
      """
      self.GandAOpex = 0
      for manager in CalculateGandAExpensesFromManagerList:
        self.GandAOpex += manager.opex + manager.capex 
      self.GandAOpex *= self.GandAFraction/(1.0 - self.GandAFraction)
      
    def CalculateRevenue(self,problemManager, mineDataManager):
      """
      Calculate revenue from sale of concentrate 
      """
      pricePerUnitMassOre = 0.0
      mineYears = len(mineDataManager.theProcessingSystem.oreProcessed)
      
      processingSystem = mineDataManager.theProcessingSystem
      
      self.revenue = np.zeros(mineYears)
      
      # loop through each processing system
      while(processingSystem):
          #for metal,grade in mineDataManager.theOreBody.metalGrades.items():
          for metal in processingSystem.concentrateMetals:
              grade = mineDataManager.theOreBody.metalGrades[metal]
              prefactor = 1.0
              if(processingSystem.processingMethod == "Ni" and metal == "Co"):  # cobalt price is halved in ni concentrates
                prefactor = 0.5 
              if(processingSystem.processingMethod == "P2O5"):
                prefactor = 1.0/0.3  # Phosphate rock price is based on 30% P2O5
              #print metal, grade,  self.commodityPrices[metal],grade*self.commodityPrices[metal]
              if( np.isscalar( self.commodityPrices[metal] ) ):
                pricePerUnitMassOre += prefactor*( grade/(1.0+ mineDataManager.theMiningSystem.dilution) )*self.commodityPrices[metal]
              elif(len( self.commodityPrices[metal] ) >=  mineYears  ):
                pricePerUnitMassOre += prefactor*( grade/(1.0+ mineDataManager.theMiningSystem.dilution) )*self.commodityPrices[metal][:mineYears]
              else:
                pricePerUnitMassOre *= np.ones(mineYears)
                priceYears = len(self.commodityPrices[metal])
                pricePerUnitMassOre[:priceYears] += prefactor*( grade/(1.0+ mineDataManager.theMiningSystem.dilution) )*self.commodityPrices[metal]
                pricePerUnitMassOre[priceYears:] += prefactor*( grade/(1.0+ mineDataManager.theMiningSystem.dilution) )*self.commodityPrices[metal][-1]
        
          pricePerUnitMassOre *= (1.0 - processingSystem.processingLoss)*(1.0 - processingSystem.refiningTake)
          self.revenue += processingSystem.oreProcessed  * pricePerUnitMassOre
          
          processingSystem = processingSystem.secondaryProcessingSystem # secondary products
          
      #print "revenue", self.revenue
      #print "pricePerUnitMassOre",pricePerUnitMassOre
      #print "processingSystem.oreProcessed", processingSystem.oreProcessed
      return self.revenue
     
    def CalculateBreakevenFactor(self,problemManager, mineDataManager,cost): 
      """Estimate the breakeven factor for all commodities(i.e. the multiplier for all commodity prices to return NPV=0)."""
      pricePerUnitMassOre = 0.0
      
      for metal,grade in mineDataManager.theOreBody.metalGrades.items():
         #print metal, grade,  self.commodityPrices[metal],grade*self.commodityPrices[metal]
         pricePerUnitMassOre += grade*self.commodityPrices[metal]
      
      pricePerUnitMassOre *= (1.0 - mineDataManager.theProcessingSystem.processingLoss)*(1.0 - mineDataManager.theProcessingSystem.refiningTake)
      #print "pricePerUnitMassOre ", pricePerUnitMassOre
            
      self.revenue = mineDataManager.theProcessingSystem.oreProcessed  * pricePerUnitMassOre

      breakevenFactor = cost/self.revenue
    
      return breakevenFactor    
          
    def CalculateBeforeTaxCashFlow(self, problemManager, mineDataManager):
      """
      Calculate the before-tax net cash flow for the project - in real terms
      NB - depreciated - will replace with the "CalculateBeforeTaxCashFlowFromManager"
      below but need to implement start up and sustaining cost calculations in mine manager
      and test - in the meantime non-mine projects (eg. hydrogen) should  use version below
      in their own economic data manager
      """
      
      self.CalculateRevenue(problemManager, mineDataManager)
      
      
      # calculate net capex/startup costs
      self.netCapex = mineDataManager.theMiningSystem.miningCapex \
                      + mineDataManager.theProcessingSystem.processingCapex \
                      + mineDataManager.theInfrastructureManager.infrastructureCapex
      
      # calculate net opex/sustaining costs
      self.netOpex =  mineDataManager.theMiningSystem.miningOpex \
                      + mineDataManager.theProcessingSystem.processingOpex  \
                      + mineDataManager.theInfrastructureManager.infrastructureOpex \
                      + self.GandAOpex
      
      if(mineDataManager.theRehabilitationManager):
        self.netOpex +=  mineDataManager.theRehabilitationManager.rehabilitationCosts
                      
      # add sustaining capital cost contribution to capex
      self.netCapex += self.sustainingCapexFraction*self.netOpex
      
      # remove sustaining capital cost contribution from opex
      self.netOpex *= 1.0-self.sustainingCapexFraction
      
      # inflate revenue opex and (sustaining) capex
      inflation =  (1.0+self.inflation)**np.array( range( mineDataManager.theMiningSystem.mineLife ) )
      self.revenue *= inflation
      self.netOpex *= inflation
      self.netCapex *= inflation
      
      # before tax net cash flow
      self.btNCF =  self.revenue - self.netOpex - self.netCapex 
      
      
      return self.btNCF    
    
                
    def CalculateTaxes(self, problemManager, mineDataManager):
      """
      Calculate the royalties and income tax for the project
      Note that the state must be set prior to this calculation
      """
      
      # linear depreciation over remainder of mine life
      yearsRemaining = np.array( range( mineDataManager.theMiningSystem.mineLife ,0,-1 ) ) 
      self.depreciatedCapex = self.netCapex/yearsRemaining  # annual depreciation rate from each year's capex
      
      for i in range(1, mineDataManager.theMiningSystem.mineLife ): # spread annual depreciation over remaining years
        self.depreciatedCapex[i] += self.depreciatedCapex[i-1]
    
    
      # depreciation deductions for mine mouth value
      mineMouthCapex = mineDataManager.theProcessingSystem.processingCapex \
                      + mineDataManager.theInfrastructureManager.infrastructureCapex 
      mineMouthDepreciatedCapex = mineMouthCapex/yearsRemaining
      for i in range(1, mineDataManager.theMiningSystem.mineLife ): # spread annual depreciation over remaining years
        mineMouthDepreciatedCapex[i] += mineMouthDepreciatedCapex[i-1]
      
      # Note that the state needs to be set prior to the calculation
      state = self.state
      processedState = self.state
      commodity = mineDataManager.theOreBody.type[:2]
      if(commodity == "P2" or commodity =="K2" or commodity == "RE"):
        commodity = mineDataManager.theOreBody.type
      
      commodityPrice = self.commodityPrices[commodity]
      profit = self.revenue - self.netOpex - self.depreciatedCapex
      value = self.revenue
      # https://resourcesandgeoscience.nsw.gov.au/__data/assets/pdf_file/0009/691713/Deductions_NonCoal.pdf
      # mine mouth value - calculated by subtracting total allowable deductions from the value of the minearl recovered
      mineMouthValue = self.revenue - mineDataManager.theProcessingSystem.processingOpex \
                        - mineDataManager.theInfrastructureManager.infrastructureOpex \
                        - (1.0/3.0)*self.GandAOpex \
                        - mineMouthDepreciatedCapex
      if type(commodityPrice) == float:
        self.royalties = AustralianRoyalties(state,commodity,profit,value,mineMouthValue,commodityPrice,processedState)
      else:
        self.royalties = AustralianRoyalties(state,commodity,profit,value,mineMouthValue,commodityPrice[:len(value)],processedState)
      
      ##
    
      #print self.netCapex
      #print self.depreciatedCapex
      #print sum(self.netCapex)
      #print sum(self.depreciatedCapex)
        
      # Taxes
      self.taxes = self.incomeTaxRate * ( self.revenue - self.netOpex - self.depreciatedCapex  - self.royalties) 
      
      # account for loss carry forward
      for i in range(mineDataManager.theMiningSystem.mineLife-1):
        if (self.taxes[i] < 0.0):  # i.e. loss in year 0.0
          self.taxes[i+1] += self.taxes[i] # carry loss forward
          self.taxes[i] = 0.0  
      
      if (self.taxes[-1] < 0.0):  # tax loss in final year 
         self.taxes[-1] = 0.0
         
      return self.taxes  
      
      
    def CalculateAfterTaxCashFlow(self, problemManager, mineDataManager):
      """
      Calculate the after-tax net cash flow for the project.
      """
      self.atNCF =  self.btNCF - self.taxes - self.royalties
      return self.atNCF  
      
      
      
    # Economic Indicators  
      
    def CalculateNPV(self,ncf):
      """
      Calculate the discounted net present value of a cash flow. 
      """
      numYears = len(ncf)
      years = np.array(range(numYears)) + 1.0 # NB end of year discounting (and no year 0)
      discountFactors = (1.0+self.discountRate)**years
      return np.sum(ncf/discountFactors)      
      
      
    def CalculateBeforeTaxNPV(self, problemManager, mineDataManager):
      """
      Calculate the before tax net present value for the project.
      """
      self.btNPV  = self.CalculateNPV(self.btNCF)
      
      # additional NPC due to rehabilitation bond
      if(mineDataManager.theRehabilitationManager):
        self.btNPV -= mineDataManager.theRehabilitationManager.rehabilitationNPC
      
      return self.btNPV 
      
    
      
    def CalculateAfterTaxNPV(self, problemManager, projectManager):
      """
      Calculate the after tax net present value for the project.
      """
      self.atNPV  = self.CalculateNPV(self.atNCF)
      
      # additional NPC due to rehabilitation bond
      if(projectManager.theRehabilitationManager):
        self.atNPV -= projectManager.theRehabilitationManager.rehabilitationNPC
      
      return self.atNPV 
     
 
    def CalculateEquivalentAnnuity(self,npv,numYears):
      """
      Calculate the annuity that would deliver the same NPV .
      """
      #years = np.array(range(numYears)) + 1.0 # NB end of year discounting (and no year 0)
      denom = (1-1./(1+self.discountRate)**numYears)*(1+ self.discountRate)/self.discountRate
      ea = npv/(denom+1e-64)
      return ea   
 
     
    def EstimateDirectEmployment(self, totalStartupCosts):
      """
      Estimate employment for the project.
      """
      self.employment  = 10.24 * (totalStartupCosts/1e6)**0.53  
      # jobs = 10.24 x (M$AUD 2018)^0.5314 
      #(Based on operating employment data from Resources and Energy Major Projects List: 2015, 2018) 
      
      return self.employment  