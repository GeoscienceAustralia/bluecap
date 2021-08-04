from bluecap_tests.Common import *



class Test_Regional_Solar(XmlRunNpyDiff_Base):  
  """Test regional calculation of hydrogen production from solar power."""
  _xmlFile = 'Solar_Regional.xml'
  _outputFile = 'Solar_Regional.out.npy'
  
class Test_Regional_Wind(XmlRunNpyDiff_Base):  
  """Test regional calculation of hydrogen production from wind power."""
  _xmlFile = 'Wind_Regional.xml'
  _outputFile = 'Wind_Regional.out.npy'

class Test_Regional_CSP(XmlRunNpyDiff_Base):  
  """Test regional calculation of hydrogen production from Concentrated Solar Power."""
  _xmlFile = 'CSP_Regional.xml'
  _outputFile = 'CSP_Regional.out.npy'

class Test_Regional_Gas(XmlRunNpyDiff_Base): 
  """Test regional calculation of hydrogen production from steam methane reformation.""" 
  _xmlFile = 'Gas_Regional.xml'
  _outputFile = 'Gas_Regional.out.npy'
  
class Test_Regional_BlackCoal(XmlRunNpyDiff_Base): 
  """Test regional calculation of hydrogen production from black coal gasification.""" 
  _xmlFile = 'BlackCoal_Regional.xml'
  _outputFile = 'BlackCoal_Regional.out.npy'
  
class Test_Regional_BrownCoal(XmlRunNpyDiff_Base): 
  """Test regional calculation of hydrogen production from brown coal gasification."""  
  _xmlFile = 'BrownCoal_Regional.xml'
  _outputFile = 'BrownCoal_Regional.out.npy'