from bluecap_tests.Common import *



class Test_Regional_CuAg(XmlRunNpyDiff_Base):  
  """Test regional economic fairways calculation for Copper/Silver developments."""
  _xmlFile = 'CuAg_Regional.xml'
  _outputFile = 'CuAg_Regional.out.npy'
  
class Test_Regional_PbZn(XmlRunNpyDiff_Base):  
  """Test regional economic fairways calculation for Lead/Zinc developments."""
  _xmlFile = 'PbZn_Regional.xml'
  _outputFile = 'PbZn_Regional.out.npy'

class Test_Regional_REO(XmlRunNpyDiff_Base):  
  """Test regional economic fairways calculation for coproduced Rare Earth Oxide."""
  _xmlFile = 'REO_Regional.xml'
  _outputFile = 'REO_Regional.out.npy'


class Test_Regional_Potash(XmlRunNpyDiff_Base):  
  """Test regional economic fairways calculation for KSO4 from Brine."""
  _xmlFile = 'Potash_Regional.xml'
  _outputFile = 'Potash_Regional.out.npy'
  
class Test_Regional_Sensitivity(XmlRunNpyDiff_Base):  
  """Test regional economic fairways sensitivity calculation."""
  _xmlFile = 'Sensitivity_Regional.xml'
  _outputFile = 'Sensitivity_Regional.out_snstvty_discount_rate.npy'
  
  
class Test_Regional_ComparativeSensitivity(XmlRunNpyDiff_Base):  
  """Test regional economic fairways comparative sensitivity calculation."""
  _xmlFile = 'ComparativeSensitivity_Regional.xml'
  _outputFile = 'ComparativeSensitivity_Regional.out_snstvty_RGB.npy'