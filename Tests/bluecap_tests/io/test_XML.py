from bluecap_tests.Common import *


class Test_XMLSchemaTest(XmlSchema_Base):  
  """Simple test of the XML Schema."""
  _xmlFile = 'XMLSchemaTest.xml'

class Test_SaveXMLTest(XmlRunDiff_Base):  
  """Simple test of save function.""" 
  _xmlFile = 'SaveXMLTest.xml'
  _outputFile = 'SaveXMLTest.A.out.xml'
  
