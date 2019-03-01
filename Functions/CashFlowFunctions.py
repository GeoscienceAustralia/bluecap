

def PriceIndexation(priceO,indexO,indexI):
  rv = priceO*indexI/indexO
  return rv
  

def PriceInflation(price,rate,years):
  rv = price*(1.0+rate)**years
  return rv
  


