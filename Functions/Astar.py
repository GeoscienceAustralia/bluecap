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


def FindNewUnvisitedNeighbours(x,y,visited,distances,heuristicDistances,parent, heuristicFn, costFn):
  nbrs = []
  nbrDistances = []
  nx,ny = visited.shape
  ddo = distances[x,y]
  for dx in [-1,0,1]:
    nbrx = x+dx 
    if (nbrx >= 0) and (nbrx < nx):
      for dy in [-1,0,1]:
        nbry = y+dy 
        if (nbry >= 0) and (nbry < ny) and (visited[nbrx,nbry] == 0) and ( dy is not 0 or dx is not 0):
          ddb = ddo + costFn(x,y,nbrx,nbry)
          if(distances[nbrx,nbry] > ddb):
             distances[nbrx,nbry] = ddb
             hd = ddb + heuristicFn(nbrx,nbry)
             heuristicDistances[nbrx,nbry] = hd
             nbrs.append((nbrx,nbry))
             nbrDistances.append([ddb,hd])
             parent[nbrx,nbry] =  x*ny + y
             
  return nbrs,nbrDistances



def AStar(nX,nY,initialX,initialY,targetX,targetY, heuristicFn, costFn):
    
    targetX,targetY = int(targetX),int(targetY)
    distances = 1e64*np.ones([nX,nY],dtype=np.float)
    visited = np.zeros([nX,nY],dtype=np.int)
    parent = -np.ones([nX,nY],dtype=np.int)
    heuristicDistances = np.zeros([nX,nY],dtype=np.float)

    fringe = {}
    

    activeX = int(initialX)
    activeY = int(initialY)
    
    # if initial node is target node then break
    if(activeX == targetX) and (activeY == targetY):
      path = []
      pathCost = 0.0
      return path, pathCost
    
    
    distances[activeX,activeY] = 0.0
    visited[activeX,activeY] = 1

    maxIters = 50000

    for iters in range(maxIters):
      updatedNbrs,nbrDistances = FindNewUnvisitedNeighbours(activeX,activeY,visited,distances,heuristicDistances,parent, heuristicFn, costFn)

      for i in range(len(updatedNbrs)):
        fringe[updatedNbrs[i]] = nbrDistances[i][1]
  
      fringe.pop((activeX,activeY),0)  # 0 prevents key error

      visited[activeX,activeY] = 1
      
      #break if no elements in fringe
      if(len(fringe) == 0):
        pathCost = 1e64
        break
  
      # get the heuristic distances associated with the fringe and return the minimum fringe
  
      # make the smallest the next active element
      activeX,activeY = min(fringe, key=fringe.get)
  
      # break if the target node is the smallest/active node
      if [activeX,activeY] == [targetX,targetY]:
        pathCost = fringe[ (activeX,activeY) ]
        break
    
      if iters is maxIters:
        print("Warning: Maximum number of iterations exceeded!")

    # recreate path
    path = []

    activeX,activeY = int(targetX),int(targetY)
    iX,iY = int(initialX),int(initialY)
    
    path.append( (activeX, activeY) )
    while([activeX,activeY] != [iX,iY]):
      pp  = parent[activeX,activeY]
      activeX  = pp//nY
      activeY  = pp%nY
      
      path.append( (activeX, activeY) )

    return path,pathCost
