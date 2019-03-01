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

#Units
from Units.UnitManager import UnitManager



def AustralianRoyalties(state, commodity,profit,value,minemouth,commodityPrice,processedState):
  """ Australian royalty calculation"""
  rv =0.0
  
  if (state == "NSW" or state == "ACT"):   # ACT royalties assumed equal to NSW
    if (commodity == "Au" or commodity == "Cu" or commodity == "Ni" or commodity =="Zn" or commodity =="Pb"):
      rv = minemouth*0.04   # 4% of mine mouth
  
  elif(state == "VIC"):
  
    if (commodity == "Au"):  # Vic has no royalties on gold
      rv = 0.0
    elif (commodity == "Cu" or commodity == "Ni" or commodity =="Zn" or commodity =="Pb"):
      rv = 0.0275*value    # 2.75% of market value
    
  elif(state == "NT"):
    rv = 0.20*profit   # 20% of profit
    
  elif(state == "WA"):
    if (commodity == "Au"):
      theUnitManager =  UnitManager()
      troyOz = theUnitManager.ConvertToBaseUnits("toz")
      rv = ( value  - commodityPrice * 2500*troyOz ) * 0.025    # first 2500 troy oz are exempt
      rv = np.maximum( 0.0,rv)
      
    elif (commodity == "Ni" ):
      rv = value * 0.025  # actually 2.5% of market value of the contained nickel (not value of the concentrate)

    elif (commodity =="Cu" or commodity =="Zn" or commodity =="Pb" ):
      rv = value * 0.025  # if sold as pure metal, otherwise 5% if sold as concentrate

  
  elif(state == "QLD"):
    inStateAdjustment = 1.0
    theUnitManager =  UnitManager()
    if (commodity =="Au"):
      AUDperTroyOz = theUnitManager.ConvertToBaseUnits("AUD/toz")
      lowerThresholdPrice = 600*AUDperTroyOz
      upperThresholdPrice = 890*AUDperTroyOz
      
    elif (commodity =="Cu" ):
      AUDperTonne = theUnitManager.ConvertToBaseUnits("AUD/tonne")
      lowerThresholdPrice = 3600*AUDperTonne
      upperThresholdPrice = 9200*AUDperTonne
      # Where copper is processed in Queensland to a metal content of at least 95%, the royalty payable is discounted by 20%.
      inStateAdjustment=0.8
    elif (commodity =="Ni" ):
      AUDperTonne = theUnitManager.ConvertToBaseUnits("AUD/tonne")
    
      lowerThresholdPrice = 12500*AUDperTonne
      upperThresholdPrice = 38100*AUDperTonne
      # Where copper is processed in Queensland to a metal content of at least 95%, the royalty payable is discounted by 20%.
      inStateAdjustment=0.8
    elif (commodity =="Pb"):
      AUDperTonne = theUnitManager.ConvertToBaseUnits("AUD/tonne")
      lowerThresholdPrice = 1100*AUDperTonne
      upperThresholdPrice = 2500*AUDperTonne
      # Where zinc is processed in Queensland to a metal content of at least 95%, the royalty payable is discounted by 35%.
      inStateAdjustment=0.65
    elif (commodity =="Zn" ):
      AUDperTonne = theUnitManager.ConvertToBaseUnits("AUD/tonne")
      lowerThresholdPrice = 1900*AUDperTonne
      upperThresholdPrice = 4400*AUDperTonne
      # Where zinc is processed in Queensland to a metal content of at least 95%, the royalty payable is discounted by 35%.
      inStateAdjustment=0.65
    
    if(commodityPrice < lowerThresholdPrice):
      royaltyRate = 0.025
    elif(commodityPrice > upperThresholdPrice):
      royaltyRate = 0.05
    else:
      royaltyRate = 0.025+0.025*(commodityPrice-lowerThresholdPrice)/(upperThresholdPrice-lowerThresholdPrice)
      
      
    rv = value*royaltyRate
    if (processedState == "QLD"):
        rv*=inStateAdjustment


  elif(state == "TAS"):
    theUnitManager =  UnitManager()
    rv = value*0.019  # 1.9% of net sales 

    oneHundredThousandAUD = 100000* theUnitManager.ConvertToBaseUnits("AUD")
    if(np.isscalar(value)):
      if (value > oneHundredThousandAUD) and (profit > 0.0):
        rv += 0.4 * profit**2/value
    else:
      indx = np.logical_and(value > oneHundredThousandAUD, profit > 0.0)
      rv[indx] += 0.4 * profit[indx]**2/value[indx]
      
    #if rv > 0.0535 * value:  # cap
    rv = np.minimum(rv, 0.0535 * value  )
      
      
  elif(state == "SA"):
    if (commodity == "Au" or commodity =="Cu" or commodity =="Zn" or commodity =="Pb"):
      rv = 0.035 * value   # assumes refined value (concentrate is 5% of value)
    elif (commodity == "Ni"):
      rv = 0.035 * value
  
  
    
  return rv