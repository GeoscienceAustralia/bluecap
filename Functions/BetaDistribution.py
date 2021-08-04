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

from scipy import interpolate
from scipy.special import gamma, betainc

class ScaledBetaDistribution:

  def __init__(self,a,b,min=0.0,max=1.0):
    self.min = min
    self.max = max
    self.a = a
    self.b = b

  @staticmethod
  def FromMeanAndStd(mean,std=None,min=0.0,max=1.0):
    if (std is None):
      std = (max-min)/6.0
    # generate a beta distribution from the mean and std deviation
    dx = max-min + 1e-64
    scaledVar = (std/dx)**2
    scaledAv = (mean -min)/dx
    
    TT =  scaledVar + scaledAv**2

    cc,dd = ScaledBetaDistribution.FittedBetaParams(scaledAv,TT)
    if(np.isscalar(cc)):
      rv = ScaledBetaDistribution(cc,dd,min,max)
    elif(np.ndim(cc) == 1): #elif( isinstance(cc,np.ndarray) ):
      rv = np.array([ScaledBetaDistribution(cc[i],dd[i],min[i],max[i]) for i in range(len(cc)) ] )
    else: #elif( isinstance(cc,np.ndarray) ):
      cc=  cc.ravel()
      dd = dd.ravel()
      minr = min.ravel()
      maxr = max.ravel()
      rv = np.array([ScaledBetaDistribution(cc[i],dd[i],minr[i],maxr[i]) for i in range(len(cc)) ] ).reshape(mean.shape)
    return rv
    
  @staticmethod
  def FromMeanAndVariance(mean,var,min=0.0,max=1.0):
    # generate a beta distribution from the mean and variance
    dx = max-min + 1e-64
    scaledVar = var/dx**2
    scaledAv = (mean -min)/dx
    
    TT =  scaledVar + scaledAv**2

    cc,dd = ScaledBetaDistribution.FittedBetaParams(scaledAv,TT)
    if(np.isscalar(cc)):
      rv = ScaledBetaDistribution(cc,dd,min,max)
    elif(np.ndim(cc) == 1): #elif( isinstance(cc,np.ndarray) ):
      rv = np.array([ScaledBetaDistribution(cc[i],dd[i],min[i],max[i]) for i in range(len(cc)) ] )
    else: #elif( isinstance(cc,np.ndarray) ):
      cc=  cc.ravel()
      dd = dd.ravel()
      minr = min.ravel()
      maxr = max.ravel()
      rv = np.array([ScaledBetaDistribution(cc[i],dd[i],minr[i],maxr[i]) for i in range(len(cc)) ] ).reshape(mean.shape)
    return rv

  @staticmethod
  def FittedBetaParams(SS,TT):
    # generate beta parameters with first and second moments: SS and TT
    
    #a =  (SS-TT)*SS/(TT-SS**2+ 1e-64)
    #b = (SS-TT)*(1.-SS)/(TT-SS**2+ 1e-64)
    b = (SS-TT)/(TT-SS**2+ 1e-64)
    a = b*SS
    b -= a
  
    return a,b

  def __rmul__(self, other):
        return self.__mul__(other)

  def __mul__(self,Y):
    
    if( isinstance(Y,ScaledBetaDistribution) ):
        Zmin = self.min * Y.min
        Zmax = self.max * Y.max
    
        avA = self.mean()
        avB = Y.mean()
        av = avA*avB
    
        scaledAv = (av - Zmin)/(Zmax - Zmin + 1e-64)
    
        varA = self.var()
        varB = Y.var()
        var = varA*varB + varA*avB**2 + varB*avA**2
        scaledVar = var/(Zmax - Zmin+ 1e-64)**2
        
        TT =  scaledVar + scaledAv**2   # 2nd moment = Var + [E(X)]^2

        cc,dd = ScaledBetaDistribution.FittedBetaParams(scaledAv,TT)
    
        Z = ScaledBetaDistribution(cc,dd,Zmin,Zmax) 
    elif( isinstance(Y,np.ndarray) ):   
       # here we want to return a numpy array if we multiply a distribution by an array
       
       rv = np.full(Y.shape,self)*Y 
       return rv
    else:
        
        Z = ScaledBetaDistribution(self.a,self.b,Y*self.min,Y*self.max) 
    return Z
    
  def __lt__(self, other):
    # less than (based on expected value)
    if( isinstance(other,ScaledBetaDistribution) ):
      rv = self.mean() < other.mean()
    else:
      rv = self.mean() < other
    return rv
    
  def __gt__(self, other):
    # greater than (based on expected value)
    if( isinstance(other,ScaledBetaDistribution) ):
      rv = self.mean() > other.mean()
    else:
      rv = self.mean() > other
    return rv    
  
  def __ge__(self, other):
    # greater than or equal to (based on expected value)
    rv = not self.__lt__(other)
    return rv    
  
  def __le__(self, other):
    rv = not self.__gt__(other)
    return rv  
        
  def __radd__(self, other):
        return self.__add__(other)

  def __add__(self,Y):
    
    if( isinstance(Y,ScaledBetaDistribution) ):
    
        Zmin = self.min + Y.min
        Zmax = self.max + Y.max
    
        scaledMean = (self.mean() + Y.mean() -Zmin)/(Zmax-Zmin+ 1e-64)
        scaledVar = (self.var() + Y.var())/(Zmax-Zmin+ 1e-64)**2
        
        scaledSecond = scaledVar  + scaledMean**2
        cc,dd = ScaledBetaDistribution.FittedBetaParams(scaledMean,scaledSecond)
    
        Z = ScaledBetaDistribution(cc,dd,Zmin,Zmax) 
    elif( isinstance(Y,np.ndarray) ):   
        # here we want to return a numpy array if we add an array to a distribution
        Z = np.full(Y.shape,self) + Y 
    else:
        # Y is a scalar
        Z = ScaledBetaDistribution(self.a,self.b, self.min + Y, self.max + Y) 
      
    return Z

  def __neg__(self):
        return ScaledBetaDistribution(self.b, self.a,-self.max,-self.min)

  def __sub__(self, other):
      if( isinstance(other,ScaledBetaDistribution) ):
        rv = self + (-other)
      else:
        rv = ScaledBetaDistribution(self.a,self.b, self.min - other, self.max - other)
      return rv 

  def __rsub__(self, other):
      rv = other + (-self)
      return rv 
      
  
  def mean(self):
    # expected value
    unscaledAv = self.a/(self.a + self.b+ 1e-64)
    av = unscaledAv*(self.max-self.min) + self.min
    return av

  def median(self):
    # median value
    # if not (self.a > 1) & (self.b > 1):
      # print(self.a, self.b)
    # unscaledMd = (self.a - 1./3.) / (self.a + self.b - 2./3.)
    # md = unscaledMd*(self.max-self.min) + self.min
    # return md
    pass
    
  def var(self):
    # variance 
    denom = ((self.a+self.b)**2) * (self.a+self.b+1.0)+ 1e-64
    v = (self.a*self.b)/denom
    v *= (self.max-self.min)**2
    return v
    
  def pdf(self,x):
    # probability density for x
    p = 0.0
    if (x > self.min) and (x < self.max):
      scale = (self.max-self.min)+1e-64
      xx = (x-self.min)/scale
      p = xx**(self.a-1.0) * (1.0-xx)**(self.b-1.0) * gamma(self.a+self.b)/(gamma(self.a)*gamma(self.b) )/scale
    return p


# Functions for acting on numpy arrays

def GetExpectations(x):
  # test = np.array([ xi for xi in x.ravel()])
  # from matplotlib import pyplot as plt
  # plt.hist(test[-1])
  # print(test[-1].mean(), np.median(test[-1]))
  # print(np.sqrt(test[-1].var()))
  # plt.show()
  return np.array([ xi.mean() for xi in x.ravel()]).reshape(x.shape)

def GetVariances(x):
  return np.array([ xi.var() for xi in x.ravel()]).reshape(x.shape)

def GetMaximums(x):
  return np.array([ xi.max for xi in x.ravel()]).reshape(x.shape)

def GetMinimums(x):
  return np.array([ xi.min for xi in x.ravel()]).reshape(x.shape)

def GetPDFs(dists,xs):
  return np.array([ disti.pdf(xi) for disti,xi in zip(dists.ravel(),xs.ravel())]).reshape(xs.shape)


# 1d interpolation function
def UncertainInterp1d(xi,yi):
  
  mui = GetExpectations(yi)
  vari = GetVariances(yi)
  mini =  GetMinimums(yi)
  maxi =  GetMaximums(yi)
  
  mu_f = interpolate.interp1d(xi,mui)
  var_f = interpolate.interp1d(xi,vari)
  min_f = interpolate.interp1d(xi,mini)
  max_f = interpolate.interp1d(xi,maxi)
  
  def theInterpFun(x):
   
    doUncertainInterp = False
    if(np.ndim(x) > 0 and len(x) > 0):
      doUncertainInterp = isinstance(x.ravel()[0],ScaledBetaDistribution)
    
    if( doUncertainInterp ): 
      #when x is an array of distributions we integrate over the range of each distibution 
      # using Simpson's rule to evaluate the new distribution field
      quadWeights = np.array([1./3.,4./3.,2./3.,4./3.,2./3.,4./3.,1./3.]) 
      # Simpson's Rule NB quad-weights later gets scaled by dx
      nInterp = len(quadWeights)
      xmins = GetMinimums(x)
      xmaxs = GetMaximums(x)
      dx = (xmaxs-xmins)/(nInterp-1.0)
      rvInterp = 0
      for i in range(1,nInterp-1):
        # nb we assume that PDF(xmin) and PDF(xmax) = 0
        xis = xmins + i*dx   # quadrature points for distributions at each point
        pdfs = GetPDFs(x,xis) 
        rvInterp += (quadWeights[i]*dx*pdfs)*theInterpFun(xis)
    else:
      muInterp = mu_f(x)
      varInterp = var_f(x)
      stdInterp = varInterp**0.5
      minInterp = min_f(x)
      maxInterp = max_f(x)
    
      rvInterp = ScaledBetaDistribution.FromMeanAndStd(muInterp,stdInterp,min=minInterp, max=maxInterp)
    
    return rvInterp
  
  return theInterpFun
