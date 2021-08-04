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

from .EconomicDataManager import EconomicDataManager

# Managers

class HydrogenEconomicDataManager(EconomicDataManager):

    def __init__(self):
      EconomicDataManager.__init__(self)
      """
      Create an empty economic data manager and default variables. 
      """
      
      self.GandAOpex = np.array([0.0])
      self.GandAFraction = 0.14   # measures proportion of expenses dedicated to G&A
      
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
      
      # calculate contribution of commodity costs (i.e. not included in cost model)
      self.calculateGasCosts = False
      self.calculateBlackCoalCosts = False
      self.calculateBrownCoalCosts = False
      
      # NPV
      self.discountRate = 0.0
      self.btNPV = 0.0
      self.atNPV = 0.0

    def ParseXMLNode(self, economicDataNode):
      """
      Generate Economic Data Manager data from xml tree node. 
      """
      
      EconomicDataManager.ParseXMLNode(self, economicDataNode)
      
      self.discountRate = GetAttributeValue(economicDataNode,"discountRate")
      
      for child in GetChildren(economicDataNode,"Commodity"):
      	name = GetAttributeString(child,"name")
      	price = GetAttributeValue(child,"price")
      	self.commodityPrices[name] = price
      	self.referenceCommodityPrices[name] = price
      	self.commodityPriceSigmas[name] = 0.0
      	if(HasAttribute(economicDataNode,"sigma")):
      	  self.commodityPriceSigmas[name] = GetAttributeValue(child,"sigma")
      	
      if(HasAttribute(economicDataNode,"GandAFraction")):
        self.GandAFraction = GetAttributeValue(economicDataNode,"GandAFraction")
        
      if(HasAttribute(economicDataNode,"state")):
        self.state = GetAttributeString(economicDataNode,"state")
        
      if(HasAttribute(economicDataNode,"inflation")):
        self.inflation = GetAttributeValue(economicDataNode,"inflation")


      if(HasAttribute(economicDataNode,"calculateGasCosts")):
        self.calculateGasCosts = GetAttributeValue(economicDataNode,"calculateGasCosts")
      
      if(HasAttribute(economicDataNode,"calculateBlackCoalCosts")):
        self.calculateBlackCoalCosts = GetAttributeValue(economicDataNode,"calculateBlackCoalCosts")

      if(HasAttribute(economicDataNode,"calculateBrownCoalCosts")):
        ## NB: This model uses CSIRO brown coal cost estimate (good) but employs black coal sensitivity to coal price per GJ due to lack of data.
        self.calculateBrownCoalCosts = GetAttributeValue(economicDataNode,"calculateBrownCoalCosts")

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
    
    
          
    
    def CalculateGandAExpenses(self, problemManager, hydrogenDataManager):
      """
      Calculate the general and administrative expense component of the project costs 
      Assumed equal to a fixed percentage of the overall plant, G&A, and energy costs.
      """
      self.GandAOpex = hydrogenDataManager.theHydrogenPlant.opex + hydrogenDataManager.theHydrogenPlant.capex \
                      + hydrogenDataManager.thePowerPlant.opex + hydrogenDataManager.thePowerPlant.capex
      self.GandAOpex *= self.GandAFraction/(1.0 - self.GandAFraction)


      
    def CalculateRevenue(self,problemManager, hydrogenManager):
      """
      Calculate revenue from sale of hydrogen 
      """
      hydrogenPrice = self.commodityPrices["H2"]
      projectYears = hydrogenManager.GetProjectDuration()
      
      self.revenue = np.zeros(projectYears)
      
      pricePerUnitMassHydrogen = 0.0
      
      if( np.isscalar( self.commodityPrices["H2"] ) ):
        pricePerUnitMassHydrogen += self.commodityPrices["H2"]
      elif(len( self.commodityPrices["H2"] ) >=  projectYears  ):
        pricePerUnitMassHydrogen += self.commodityPrices["H2"][:projectYears]
      else:
        pricePerUnitMassHydrogen *= np.ones(projectYears)  # converts scalar to array (only required if multiple commodities)
        priceYears = len(self.commodityPrices["H2"])
        pricePerUnitMassHydrogen[:priceYears] += self.commodityPrices["H2"]
        pricePerUnitMassHydrogen[priceYears:] += self.commodityPrices["H2"][-1]


      #pricePerUnitMassOre *= (1.0 - processingSystem.processingLoss)*(1.0 - processingSystem.refiningTake)
      
      self.revenue += hydrogenManager.theHydrogenPlant.hydrogenProduced  * pricePerUnitMassHydrogen
      
      #processingSystem = processingSystem.secondaryProcessingSystem # secondary products? eg energy?
          
      #print "revenue", self.revenue
      #print "pricePerUnitMassHydrogen",pricePerUnitMassHydrogen
      return self.revenue
     
    def CalculateBreakevenFactor(self,problemManager, hydrogenDataManager,cost): 
      """Not yet implemented."""
    
      return None    
          
    def CalculateBeforeTaxCashFlow(self, problemManager, projectManager):
      """
      Calculate the before-tax net cash flow for the project - in nominal terms.
      """
      
      self.CalculateRevenue(problemManager, projectManager)
      
      # calculate net capex/startup costs
      self.netCapex = projectManager.CalculateStartupCosts()
      
      # calculate net opex/sustaining costs
      self.netOpex =  projectManager.CalculateSustainingCosts()
                      
      # add sustaining capital cost contribution to capex
      self.netCapex += self.sustainingCapexFraction*self.netOpex
      
      # remove sustaining capital cost contribution from opex
      self.netOpex *= 1.0-self.sustainingCapexFraction
      
      # inflate revenue opex and (sustaining) capex
      inflation =  (1.0+self.inflation)**np.array( range( len(self.netOpex) ) )
      self.revenue *= inflation
      self.netOpex *= inflation
      self.netCapex *= inflation
      
      # before tax net cash flow
      self.btNCF =  self.revenue - self.netOpex - self.netCapex 
      
      return self.btNCF  
      
    
                
    def CalculateTaxes(self, problemManager, hydrogenDataManager):
      """
      Calculate the royalties and income tax for the project.
      Note that the state must be set prior to this calculation.
      """
      
      # linear depreciation over remainder of project life
      yearsRemaining = np.array( range( hydrogenDataManager.theHydrogenPlant.projectLife ,0,-1 ) ) 
      self.depreciatedCapex = self.netCapex/yearsRemaining  # annual depreciation rate from each year's capex
      
      for i in range(1, hydrogenDataManager.theHydrogenPlant.projectLife ): # spread annual depreciation over remaining years
        self.depreciatedCapex[i] += self.depreciatedCapex[i-1]
    
      
      # Note that the state needs to be set prior to the calculation
      state = self.state
      processedState = self.state
      commodity = "H2"
      
      commodityPrice = self.commodityPrices[commodity]
      profit = self.revenue - self.netOpex - self.depreciatedCapex
      value = self.revenue
      
      # self.royalties = AustralianRoyalties(state,commodity,profit,value,mineMouthValue,commodityPrice,processedState)
      self.royalties = np.zeros(profit.shape)  # no royalties for time being
      self.rebates = np.zeros(profit.shape) # no rebates for time being nb - should be adjusted for each year
      
        
      # Taxes
      # income taxes assumed applied after profits are adjusted for royalties and rebates 
      self.taxes = self.incomeTaxRate * ( self.revenue - self.netOpex - self.depreciatedCapex  - self.royalties + self.rebates) 
      
      # account for loss carry forward
      for i in range(hydrogenDataManager.theHydrogenPlant.projectLife-1):
        if (self.taxes[i] < 0.0):  # i.e. loss in year 0.0
          self.taxes[i+1] += self.taxes[i] # carry loss forward
          self.taxes[i] = 0.0  
      
      if (self.taxes[-1] < 0.0):  # tax loss in final year 
         self.taxes[-1] = 0.0
         
      return self.taxes  
      
      
    def CalculateAfterTaxCashFlow(self, problemManager, hydrogenDataManager):
      """
      Calculate the after-tax net cash flow for the project.
      """
      self.atNCF =  self.btNCF - self.taxes - self.royalties + self.rebates
      return self.atNCF  
      
      
    # Economic Indicators      
      
     
    def EstimateDirectEmployment(self, totalStartupCosts):
      """
      Estimate employment for the project.
      """
      #self.employment  = 10.24 * (totalStartupCosts/1e6)**0.53    - fixme based on mining projects
      # jobs = 10.24 x (M$AUD 2018)^0.5314 
      #(Based on operating employment data from Resources and Energy Major Projects List: 2015, 2018) 
      
      self.employment = None
      
      return self.employment  