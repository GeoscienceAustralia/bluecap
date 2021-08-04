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
from IO.XML import HasChild,GetChild,AddChild
from IO.XML import GetAttributeString,GetAttributeValue, SetAttributeString


#Units
from Units.UnitManager import UnitManager

#Functions
from Functions.FunctionManager import FunctionManager

# Managers

class RehabilitationDataManager():
    def __init__(self):
      """
      Create an empty rehabilitation data manager and default variables. 
      """
      self.rehabilitationCosts = []  # treated as opex costs
      self.rehabSecurityCost = 0.0
      self.rehabilitationNPC = 0.0
      

    def ParseXMLNode(self, rehabilitationSystemNode):
      """
      Generate rehabilitation system data from xml tree node. 
      """
      
      return rehabilitationSystemNode
        

    def WriteXMLNode(self, node):
      """
      Write rehabilitation system data to xml node.
      """
      
      return node  
    
    
    def ZeroRehabilitationExpenses(self, problemManager,  mineDataManager):
      """
      Set rehabilitation costs to zero (used for regional calculations).
      """
      
      self.rehabilitationCosts = np.zeros( mineDataManager.theMiningSystem.mineLife )
      self.rehabSecurityCost = 0.0
      
      return self.rehabilitationCosts

      
    def CalculateRehabilitationExpenses(self,  problemManager, mineDataManager):
      """
      Calculate rehabilitation costs.
      """
      
      self.rehabilitationCosts = np.zeros( mineDataManager.theMiningSystem.mineLife )
      
      theFunctionManager =  FunctionManager()
      theUnitManager = UnitManager()
      
      processingCapacity =  mineDataManager.theProcessingSystem.processingCapacity  *  theUnitManager.ConvertTo("1e6 tonne")
      rehabSecurityFunc = theFunctionManager.GetFunction("RehabilitationSecurity_" + mineDataManager.theProcessingSystem.processingMethod)
      
      
      self.rehabSecurityCost = rehabSecurityFunc.f([processingCapacity])
      
      # NT model: 1% of rehabSecurity paid each year of operation
      self.rehabilitationCosts = np.ones( mineDataManager.theMiningSystem.mineLife )*self.rehabSecurityCost*0.01
      
      self.rehabilitationCosts[0] = self.rehabSecurityCost
      
      self.rehabilitationNPC = self.CalculateRehabilitationBondNPC(problemManager, mineDataManager)
      
      return self.rehabilitationCosts
      
    def CalculateRehabilitationBondNPC(self,  problemManager, mineDataManager):   
      """ Calculate net present cost of rehabilitation bond. NB this is a crude calculation as it assumes mine pays all security at start and is reimbursed all at end, with the bond value discounted over that period.""" 
      
      life = mineDataManager.theMiningSystem.mineLife
      discountRate = mineDataManager.theEconomicDataManager.discountRate
      npc = self.rehabSecurityCost*(1-1/(1.0+discountRate)**life)
    
      return npc
    