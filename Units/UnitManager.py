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

import pint  
# Uses pint for unit conversion
# https://pint.readthedocs.io/en/latest/
# pip install pint

import os 



class UnitManager(object):
    """Singleton object to convert values into a consistent set of base units
    
    Usage:
      from Units.UnitManager import UnitManager
      theUnitManager = UnitManager()
      theUnitManager.ConvertToBaseUnits("50 m")
      fiftyFeetInBaseUnits = theUnitManager.ConvertToBaseUnits(50,"feet")
      valueInFeet = valueInBaseUnits*theUnitManager.ConvertTo("feet")
    """
    
    class __UnitManager:
      """ The unit manager is a singleton object."""
    
      def __init__(self):
            """Initialize the singleton object."""
            
            # would be better to define inline but seems pint does not allow yet
            dir_path = os.path.dirname(os.path.realpath(__file__)) 
            bluecapUnitsPath =  os.path.join(dir_path, "BlueCap_units.txt")  
            # uses m, kg, s, AUD by default
            # currency is defined in "currency.txt"
            self.ureg = pint.UnitRegistry(filename = bluecapUnitsPath)
  
      def ConvertToBaseUnits(self,value,units=""):
        """Convert from arbitrary units into the base units."""
        if(units):
          return value * self.ureg(units).to_base_units().magnitude
        else:
          try: 
            rv = self.ureg(value).to_base_units().magnitude
          except AttributeError:  # ints and floats don't convert
            rv = self.ureg(value)
          return rv
  
      def ConvertTo(self,units):
        """Convert from base units into the stated units."""
        return 1.0/self.ureg(units).to_base_units().magnitude

          
    instance = None    # this class variable points to the singleton object
    def __new__(cls): # __new__ always a classmethod
        if not UnitManager.instance:
            UnitManager.instance = UnitManager.__UnitManager()
        return UnitManager.instance
    def __getattr__(self, name):
        return getattr(self.instance, name)
    def __setattr__(self, name):
        return setattr(self.instance, name)
        