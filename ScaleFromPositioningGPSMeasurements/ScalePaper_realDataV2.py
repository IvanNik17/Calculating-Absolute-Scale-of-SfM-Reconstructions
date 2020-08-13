# -*- coding: utf-8 -*-
"""
Created on Fri Jun  1 12:43:04 2018

@author: ivan
"""

import numpy
import math


    
"""
## References
- [Umeyama's paper](Least-squares estimation of transformation parameters between two point patterns)
- [CarloNicolini's python implementation](https://gist.github.com/CarloNicolini/7118015)
"""    
    
def similarity_transform(from_points, to_points):
    
    assert len(from_points.shape) == 2, \
        "from_points must be a m x n array"
    assert from_points.shape == to_points.shape, \
        "from_points and to_points must have the same shape"
    
    N, m = from_points.shape
    
    mean_from = from_points.mean(axis = 0)
    mean_to = to_points.mean(axis = 0)
    
    delta_from = from_points - mean_from # N x m
    delta_to = to_points - mean_to       # N x m
    
    sigma_from = (delta_from * delta_from).sum(axis = 1).mean()
    sigma_to = (delta_to * delta_to).sum(axis = 1).mean()
    
    cov_matrix = delta_to.T.dot(delta_from) / N
    
    U, d, V_t = numpy.linalg.svd(cov_matrix, full_matrices = True)
    cov_rank = numpy.linalg.matrix_rank(cov_matrix)
    S = numpy.eye(m)
    
    if cov_rank >= m - 1 and numpy.linalg.det(cov_matrix) < 0:
        S[m-1, m-1] = -1
    elif cov_rank < m-1:
        raise ValueError("colinearility detected in covariance matrix:\n{}".format(cov_matrix))
    
    R = U.dot(S).dot(V_t)
    c = (d * S.diagonal()).sum() / sigma_from
    t = mean_to - c*R.dot(mean_from)
    
    return c*R,c,R, t


def scaleAndUncertanty(groundPoints, recPoints, stdDev):
    
    M_init,c_init,R, t_init = similarity_transform(recPoints,  groundPoints)
    stdDevXY = stdDev[0]
    stdDevZ = stdDev[1]

    
    allScales_newMethod_plus = numpy.zeros([len(groundPoints),3])
    allScales_newMethod_minus = numpy.zeros([len(groundPoints),3])



    for k in range(0,len(groundPoints)):
        for n in range(0,3):
            for l in range(0,2):
                
                groundPoints_test = []
                groundPoints_test = groundPoints.copy()
                
                if n == 0 or n == 1:
                    dimDev = stdDevXY
                else:
                    dimDev = stdDevZ
                
                if l == 0:
                    dimDev = dimDev
                else:
                    dimDev = -dimDev
                
                groundPoints_test[k,n] += dimDev

                M,c_new,R, t = similarity_transform(recPoints,  groundPoints_test)
                
                if l == 0:

                    allScales_newMethod_plus[k,n] = c_new
                else:
                    allScales_newMethod_minus[k,n] = c_new
                    
            

    allScales_newMethod_finalX = (allScales_newMethod_plus[:,0] - allScales_newMethod_minus[:,0])/(2*stdDevXY)
    allScales_newMethod_finalY = (allScales_newMethod_plus[:,1] - allScales_newMethod_minus[:,1])/(2*stdDevXY)
    allScales_newMethod_finalZ = (allScales_newMethod_plus[:,2] - allScales_newMethod_minus[:,2])/(2*stdDevZ)
    
    
    
    finalAll = numpy.zeros([1,3*len(allScales_newMethod_finalX)])
    
    for q in range(0,len(finalAll[0]),3):
        finalAll[0,q] = allScales_newMethod_finalX[int(q/3)]
        finalAll[0,q+1] = allScales_newMethod_finalY[int(q/3)]
        finalAll[0,q+2] = allScales_newMethod_finalZ[int(q/3)]
    
    firstSecondLine = numpy.ones((1,len(groundPoints)*3))*stdDevXY**2
    thirdLine = numpy.ones((1,len(groundPoints)*3))*stdDevZ**2
    combinedThree = numpy.concatenate([firstSecondLine,firstSecondLine,thirdLine])
    fullCombined = combinedThree.copy()
    for j in range(0,len(groundPoints)-1):
        fullCombined = numpy.concatenate([fullCombined,combinedThree])
    
    fullCombined[~numpy.eye(fullCombined.shape[0],dtype=bool) == True] = 0
    
    
    
    #scaleTest = numpy.linalg.multi_dot([finalAll, fullCombined,finalAll.T])
    scaleTest_temp = numpy.dot(fullCombined, numpy.transpose(finalAll))
    uncentainty = numpy.dot(finalAll, scaleTest_temp)
    
    
    return [c_init, math.sqrt(uncentainty)]
    


#if __name__ == '__main__': 
    
##    ax = initializeFigurePlot()
#    
##    axisSize = 500
#    
#    
##    groundPoints = numpy.loadtxt("coords.txt", delimiter=' ') 
#    
#    groundPoints = numpy.loadtxt("F:\FinalFiles\coords_noNames.txt", delimiter=' ')
#    
##    groundPoints = groundPoints[:,1:]
#    
#    recPoints = numpy.loadtxt(r"F:\FinalFiles\testBlade_cams.txt", delimiter=',')
##    recPoints = numpy.loadtxt("cameras_blade_new.txt", delimiter=',') 
##    recPoints = recPoints[:-1,10:13]
#
#
##    groundPoints = groundPoints[:18,:]
##    recPoints = recPoints[:18,:]
#    
##    groundPoints = groundPoints
#    
##    groundPoints = groundPoints*1000
## Find similarity matrix
#
#    
#
#    
#
#    M_init,c_init,R, t_init = similarity_transform(recPoints,  groundPoints)
##
##    
##
###                        
######   Fixed Points
##    
##    transformedPoints = M.dot(recPoints.T).T + t
##
##
##
##    maxError = numpy.abs(groundPoints - transformedPoints).max()
##    
##    meanError = numpy.abs(groundPoints - transformedPoints).mean()
##    
##    stdevError = numpy.abs(groundPoints - transformedPoints).std()
##    
##    
##   
##
###    plot3D_static(ax,groundPoints[:,0],groundPoints[:,1],groundPoints[:,2],[axisSize,axisSize,axisSize],'bo', "")
###    
###
###    
###    plot3D_static(ax,recPoints[:,0],recPoints[:,1],recPoints[:,2],[axisSize,axisSize,axisSize],'ro', "")
###
###    
###    plot3D_static(ax,transformedPoints[:,0],transformedPoints[:,1],transformedPoints[:,2],[axisSize,axisSize,axisSize],'r.', "")
##
###    
###    
#### Print stuff
##
##
##    print("Error in Transformed points - ", "Average Error = ", meanError, ", Max Error = ", maxError, ", Standard Dev = ", stdevError)
#    
#    
##  Adding noise and recompute 1000 times to get the standard deviation
#    
#    
#
##    numIter = 1000
##    iteratorN = 1
##    
##
##
##    allScales = []
##    for k in range(0, numIter, iteratorN):
##        
###        groundPoints_test = []
###        groundPoints_test = groundPoints.copy()
##        
##        groundPoints_test_X = groundPoints[:,0].copy()
##        groundPoints_test_Y = groundPoints[:,1].copy()
##        groundPoints_test_Z = groundPoints[:,2].copy()
##        
##
##        
##        addNoiseX_ground = [-0.02, 0.02]
##        addNoiseY_ground = [-0.02, 0.02]
##        addNoiseZ_ground = [-0.03, 0.03]
##
###    addNoiseX_ground = [-20, 20]
###    addNoiseY_ground = [-20, 20]
###    addNoiseZ_ground = [-30, 30]
##
##        groundPoints_test_X, ground_noiseArrX = addNoise(groundPoints_test_X,addNoiseX_ground)
##        groundPoints_test_Y, ground_noiseArrY = addNoise(groundPoints_test_Y,addNoiseY_ground)
##        groundPoints_test_Z, ground_noiseArrZ = addNoise(groundPoints_test_Z,addNoiseZ_ground) 
##        
##        
##        groundPoints_test = numpy.array([groundPoints_test_X,groundPoints_test_Y, groundPoints_test_Z]).T
##
##        M,c_t,R, t = similarity_transform(recPoints,  groundPoints_test)
##   
##        allScales.append(c_t)
##    allScales_ar = numpy.array(allScales)
###    Relative standard deviation
##    print(" RTK fake noise = " + str(allScales_ar.std()) + " " + str((allScales_ar.std()/allScales_ar.mean())*100) )
#
#
#    
#    
##    rtk_data = numpy.loadtxt("RTK_reading.txt")
##    rtk_data = rtk_data
##    allScales_real = []
##    
##    for j in range(0, numIter, iteratorN):
##        
##        groundPoints_test = []
##        groundPoints_test = groundPoints.copy()
##        
##        xNoise = numpy.random.choice(rtk_data[:,0], size = len(groundPoints_test[:,0]))
##        yNoise = numpy.random.choice(rtk_data[:,1], size = len(groundPoints_test[:,1]))
##        
##        groundPoints_test[:,0] += xNoise
##        groundPoints_test[:,1] += yNoise
##    
##        M,c_t,R, t = similarity_transform(recPoints,  groundPoints_test)
##        
##        allScales_real.append(c_t)
##        
##    allScales_real_ar = numpy.array(allScales_real)
##    
##    print(" RTK real = " + str(allScales_real_ar.std()) + " " + str((allScales_real_ar.std()/allScales_real_ar.mean())*100) )
#    
#    
#    
#    
#    stdDevXY = 0.017538494
#    stdDevZ = 0.02445698  
#    
#    allScales_newMethod_plus = numpy.zeros([len(groundPoints),3])
#    allScales_newMethod_minus = numpy.zeros([len(groundPoints),3])
#
#
#
#    for k in range(0,len(groundPoints)):
#        for n in range(0,3):
#            for l in range(0,2):
#                
#                groundPoints_test = []
#                groundPoints_test = groundPoints.copy()
#                
#                if n == 0 or n == 1:
#                    dimDev = stdDevXY
#                else:
#                    dimDev = stdDevZ
#                
#                if l == 0:
#                    dimDev = dimDev
#                else:
#                    dimDev = -dimDev
#                
#                groundPoints_test[k,n] += dimDev
#
#                M,c_new,R, t = similarity_transform(recPoints,  groundPoints_test)
#                
#                if l == 0:
#
#                    allScales_newMethod_plus[k,n] = c_new
#                else:
#                    allScales_newMethod_minus[k,n] = c_new
#                    
#            
#
#allScales_newMethod_finalX = (allScales_newMethod_plus[:,0] - allScales_newMethod_minus[:,0])/(2*stdDevXY)
#allScales_newMethod_finalY = (allScales_newMethod_plus[:,1] - allScales_newMethod_minus[:,1])/(2*stdDevXY)
#allScales_newMethod_finalZ = (allScales_newMethod_plus[:,2] - allScales_newMethod_minus[:,2])/(2*stdDevZ)
#
#
#
#finalAll = numpy.zeros([1,3*len(allScales_newMethod_finalX)])
#
#for q in range(0,len(finalAll[0]),3):
#    finalAll[0,q] = allScales_newMethod_finalX[int(q/3)]
#    finalAll[0,q+1] = allScales_newMethod_finalY[int(q/3)]
#    finalAll[0,q+2] = allScales_newMethod_finalZ[int(q/3)]
#
#firstSecondLine = numpy.ones((1,len(groundPoints)*3))*stdDevXY**2
#thirdLine = numpy.ones((1,len(groundPoints)*3))*stdDevZ**2
#combinedThree = numpy.concatenate([firstSecondLine,firstSecondLine,thirdLine])
#fullCombined = combinedThree.copy()
#for j in range(0,len(groundPoints)-1):
#    fullCombined = numpy.concatenate([fullCombined,combinedThree])
#
#fullCombined[~numpy.eye(fullCombined.shape[0],dtype=bool) == True] = 0
#
#
#
##scaleTest = numpy.linalg.multi_dot([finalAll, fullCombined,finalAll.T])
#scaleTest_temp = numpy.dot(fullCombined, numpy.transpose(finalAll))
#scaleTest = numpy.dot(finalAll, scaleTest_temp)
#print(" RTK analytical matrix = " + str(math.sqrt(scaleTest)))
#
##contr2NoiseX = 0
##contr2NoiseY = 0
##contr2NoiseZ = 0
##
##for i in range(len(allScales_newMethod_finalX)):
##    contr2NoiseX+=(allScales_newMethod_finalX[i]**2 * stdDevXY**2)
##    contr2NoiseY+=(allScales_newMethod_finalY[i]**2 * stdDevXY**2)
##    contr2NoiseZ+=(allScales_newMethod_finalZ[i]**2 * stdDevZ**2)
##    
##scaleNoise = contr2NoiseX + contr2NoiseY + contr2NoiseZ
##
##print(" RTK analytical loop = " + str(math.sqrt(scaleNoise)))
        