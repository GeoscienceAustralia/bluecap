"""
Copyright (C) 2019-2024, Monash University, Geoscience Australia
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

# Calculate the pareto rank of every point in a set of points

def FindParetoRanks(pts):

    numDims = pts.shape[1]

    #####

    # sort the sum - later entries are not dominated by earlier ones
    
    sortedIds = pts.sum(1).argsort() 

    rankingsOrderedById = np.zeros(pts.shape[0],dtype=int)  # maps rank in sorted order back to id
    rankingsOrderedById[sortedIds] = np.array(range(pts.shape[0]),dtype=int)

    ##
    parentIds = np.array(range(pts.shape[0]),dtype=int)   # self parents
    paretoRank = pts.shape[0]*np.ones(pts.shape[0],dtype=int)


    # nb this is a slow way to find pareto optimal solutions:
    # we are looping from the least likely to dominate to the most likely and checking (then removing) earlier dominated nodes
    # going the other way (from most likely to least likely) is faster but only gives the pareto front
    # here we first build up a set of likely candidates to speed up the second step (assigning the nodes to particular pareto ranks)


    toSearch = np.array([],dtype=int)

    for id in sortedIds:


      toRemove = np.ones_like(toSearch,dtype=bool)
      for i in range(numDims):
        toRemove[toRemove] = pts[toSearch[toRemove],i] <= pts[id,i]
  
      parentIds[ toSearch[toRemove] ] = id
      
      #toSearch = toSearch[np.logical_not(toRemove)] 
      #toSearch = np.append(toSearch,id)
      toSearch = np.concatenate([ toSearch[np.logical_not(toRemove)]  ,[id]])


    # plot pareto front
    #if doPlots:
    #    pl.figure()
    #    for i in range(len(parentIds)):
    #          di = parentIds[i]
    #          pl.plot([pp[i,0],pp[di,0]],[pp[i,1],pp[di,1]],"k-")
    #          if(di == i):
    #            pl.plot([pp[i,0]],[pp[i,1]],"ro")
        


    # pareto front
    currentFront = np.array(toSearch) # np.where(parentIds == np.array(range(pts.shape[0])) )  # setting current front to pareto front. 

    activeRank = 0
    paretoRank[currentFront] = activeRank
    
    # next we want to iteratively remove the active front, 
    # then search those points connected to it to see if they are either in the new front, 
    # or if they are dominated by another node that is forward of the node they were dominated by - that has a rank less than the extracted ranks 
    # (they should not be dominated by any nodes behind the parent node by virtue of the previous search method)

    potentialNextRankNodes = np.where( np.logical_and(paretoRank[parentIds[sortedIds]] == activeRank, paretoRank[sortedIds] != 0  ) )[0]

    while(len(potentialNextRankNodes > 0)):
        
        potentialNodeIds = sortedIds[potentialNextRankNodes]
    
        for ii in range(len(potentialNodeIds)):
        
          nd = potentialNodeIds[ii]
          
          forwardNodes = np.array(potentialNodeIds[ii+1:])
          
          
          # Find dominating neighbours (i.e. all coords > than that of current node)
          
          i = 0
          while(forwardNodes.size > 0 and i < numDims):
            isDominatingNbr = pts[forwardNodes,i] >= pts[nd,i] 
            forwardNodes = forwardNodes[isDominatingNbr]
            
            i += 1
          
          # Check if the node is part of the new pareto rank (i.e. no dominating neighbours)
          if(forwardNodes.size > 0): 
            # i.e. dominated by one or more neighbouring nodes
        
            # get id of the closest dominating node (along mean axis)
            # set parent of this node to the (likely) closest (first) dominating neighbour
            parentIds[nd] =  forwardNodes[0] 
          else:
            # no other nodes dominate - this node is part of the new front
            paretoRank[nd] = activeRank +1

        activeRank += 1
        potentialNextRankNodes = np.where( paretoRank[parentIds[sortedIds]] == activeRank  )[0]
    
    return paretoRank,parentIds

