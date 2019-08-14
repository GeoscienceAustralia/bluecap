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


# Common
from Common.Common import Todo,Fixme

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
      self.metalPrices = {}
      
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
      Generate Economic Data Manager data from xml tree node. 
      """
      
      self.discountRate = GetAttributeValue(economicDataNode,"discountRate")
      
      for child in GetChildren(economicDataNode,"Commodity"):
      	name = GetAttributeString(child,"name")
      	price = GetAttributeValue(child,"price")
      	self.metalPrices[name] = price
      	
      if(HasAttribute(economicDataNode,"GandAFraction")):
        self.GandAFraction = GetAttributeValue(economicDataNode,"GandAFraction")
        
      if(HasAttribute(economicDataNode,"state")):
        self.state = GetAttributeString(economicDataNode,"state")
        
      if(HasAttribute(economicDataNode,"inflation")):
        self.inflation = GetAttributeValue(economicDataNode,"inflation")



    def WriteXMLNode(self, node):
      """
      Write data to xml node
      """
      
      # discount rate
      SetAttributeString(node,"discountRate",self.discountRate)
      SetAttributeString(node,"incomeTaxRate",self.incomeTaxRate)
      SetAttributeString(node,"GandAFraction",self.GandAFraction)
      SetAttributeString(node,"state",self.state)
      
      # price data
      for name,price in self.metalPrices.iteritems():
        commNode = AddChild(node,"Commodity")
        SetAttributeString(commNode,"name",name)
        SetAttributeString(commNode,"price",price)
      
      
      # inflation
      SetAttributeString(node,"inflation",self.inflation)
      
      
      SetAttributeString(node,"btNPV",self.btNPV)
      SetAttributeString(node,"atNPV",self.atNPV)
      
      return node
    
    def SetState(self, stateStr):
      self.state = stateStr
      
    def CalculateGandAExpenses(self, problemManager, mineDataManager):
      """
      Calculate the general and administrative expense component of the mine costs 
      Assumed equal to a fixed percentage of the overall mining, G&A, and processing costs.
      """
      self.GandAOpex = mineDataManager.theMiningSystem.miningOpex + mineDataManager.theMiningSystem.miningCapex \
                      + mineDataManager.theProcessingSystem.processingOpex + mineDataManager.theProcessingSystem.processingCapex
      self.GandAOpex *= self.GandAFraction/(1.0 - self.GandAFraction)
      
    def CalculateRevenue(self,problemManager, mineDataManager):
      """
      Calculate revenue from sale of concentrate 
      """
      pricePerUnitMassOre = 0.0
      
      for metal,grade in mineDataManager.theOreBody.metalGrades.iteritems():
         #print metal, grade,  self.metalPrices[metal],grade*self.metalPrices[metal]
         pricePerUnitMassOre += ( grade/(1.0+ mineDataManager.theMiningSystem.dilution) )*self.metalPrices[metal]
      
      pricePerUnitMassOre *= (1.0 - mineDataManager.theProcessingSystem.processingLoss)*(1.0 - mineDataManager.theProcessingSystem.refiningTake)
      #print "pricePerUnitMassOre ", pricePerUnitMassOre
            
      self.revenue = mineDataManager.theProcessingSystem.oreProcessed  * pricePerUnitMassOre

      
      return self.revenue
     
    def CalculateBreakevenFactor(self,problemManager, mineDataManager,cost): 
    
      pricePerUnitMassOre = 0.0
      
      for metal,grade in mineDataManager.theOreBody.metalGrades.iteritems():
         #print metal, grade,  self.metalPrices[metal],grade*self.metalPrices[metal]
         pricePerUnitMassOre += grade*self.metalPrices[metal]
      
      pricePerUnitMassOre *= (1.0 - mineDataManager.theProcessingSystem.processingLoss)*(1.0 - mineDataManager.theProcessingSystem.refiningTake)
      #print "pricePerUnitMassOre ", pricePerUnitMassOre
            
      self.revenue = mineDataManager.theProcessingSystem.oreProcessed  * pricePerUnitMassOre

      breakevenFactor = cost/self.revenue
    
      return breakevenFactor    
          
    def CalculateBeforeTaxCashFlow(self, problemManager, mineDataManager):
      """
      Calculate the before-tax net cash flow for the project - in real terms
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
                      
      # add sustaining capital cost contribution to capex
      self.netCapex += self.sustainingCapexFraction*self.netOpex
      
      # remove from opex
      self.netOpex *= 1.0-self.sustainingCapexFraction
      
      # inflate revenue opex and capex
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
      
      # Note that the state needs to be set prior to the
      state = self.state
      processedState = self.state
      commodity = mineDataManager.theOreBody.type[:2]  
      commodityPrice = self.metalPrices[commodity]
      profit = self.revenue - self.netOpex - self.depreciatedCapex
      value = self.revenue
      # https://resourcesandgeoscience.nsw.gov.au/__data/assets/pdf_file/0009/691713/Deductions_NonCoal.pdf
      # mine mouth value - calculated by subtracting total allowable deductions from the value of the minearl recovered
      mineMouthValue = self.revenue - mineDataManager.theProcessingSystem.processingOpex \
                        - mineDataManager.theInfrastructureManager.infrastructureOpex \
                        - (1.0/3.0)*self.GandAOpex \
                        - mineMouthDepreciatedCapex
      self.royalties = AustralianRoyalties(state,commodity,profit,value,mineMouthValue,commodityPrice,processedState)
      
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
      Calculate the after-tax net cash flow for the project
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
      Calculate the before tax net present value for the project
      """
      self.btNPV  = self.CalculateNPV(self.btNCF)
      
      return self.btNPV 
      
    
      
    def CalculateAfterTaxNPV(self, problemManager, mineDataManager):
      """
      Calculate the after tax net present value for the project
      """
      self.atNPV  = self.CalculateNPV(self.atNCF)
      
      return self.atNPV 
     
     
    def EstimateDirectEmployment(self, totalStartupCosts):
      """
      Estimate employment for the project
      """
      self.employment  = 10.24 * (totalStartupCosts/1e6)**0.53  
      # jobs = 10.24 x (M$AUD 2018)^0.5314 
      #(Based on operating employment data from Resources and Energy Major Projects List: 2015, 2018) 
      
      return self.employment  