# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 10:53:12 2019

@author: ivan
"""

"""
Signature: mesh.ray.intersects_location(ray_origins,
                                        ray_directions,
                                        multiple_hits=True)
Docstring:
Return the location of where a ray hits a surface.
Parameters
----------
ray_origins:    (n,3) float, origins of rays
ray_directions: (n,3) float, direction (vector) of rays
Returns
---------
locations: (n) sequence of (m,3) intersection points
index_ray: (n,) int, list of ray index
index_tri: (n,) int, list of triangle (face) indexes
"""

import numpy as np

import pandas as pd

import trimesh

import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D

import math


def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])
    
    
# Raycast from each camera position towards the mesh - here you need camera position in 3D, camera forward direction, camera up direction, the LiDAR raycast angle of vision and the ammount of rays cast 
def getRayCastDist(inputMesh, cam_pos, cam_norm, cam_perp, rot_rayCast_max, rot_rayCast_delta):

    allHits_allCams = np.zeros([len(cam_pos), int((2*rot_rayCast_max)/rot_rayCast_delta)])
    
    allHitsAverage = []
    
    #  from the camera position cast rays at the forward direction
    for i in range(0, len(cam_pos)):
        
        ray_origins = np.array([cam_pos[i,:]])
    
        
        ray_perp = cam_perp[i,:]
        ray_norm = cam_norm[i,:]
        
    
        #  cast a number of rays equal to the field of view of the LiDAR and containing the angle density
        allHits_currCam = np.zeros([1, int((2*rot_rayCast_max)/rot_rayCast_delta)])
        count = 0
        for j in np.arange(-rot_rayCast_max, rot_rayCast_max, rot_rayCast_delta):
            rot_angle = np.deg2rad(j)
            
            # calculate the normal of the camera
            cam_norm_curr = np.dot(rotation_matrix(ray_perp, rot_angle), ray_norm)
            
            # calculate the ray direction
            ray_directions_curr = np.array([cam_norm_curr])
        
            # compute the intersection between the ray and the mesh and return the location in 3D, the index of the triangle
            locations, index_ray, index_tri = inputMesh.ray.intersects_location(
            ray_origins=ray_origins,
            ray_directions=ray_directions_curr)
            
            # check if there is a intersection, if there is calculate the distance and the angle between the camera and mesh
            if len(locations) != 0:
                allHitDist = np.sqrt(np.sum((locations-ray_origins)**2,axis=1)) 
                hitArg = np.argmin(allHitDist)
                
                allHits_currCam[0,count] = allHitDist[hitArg]
            
            
            count+=1
            

        
        
        
        # append and save all the hit distances and angles
        allHits_allCams[i] = allHits_currCam
        
        allHitsAverage.append(allHits_currCam[0,np.nonzero(allHits_currCam)[1]].mean() )
            
    # calcualte the average distance and angle
    allCams_meanRayDist = np.true_divide(allHits_allCams.sum(1),(allHits_allCams!=0).sum(1))  
    
    allHitsAverage = np.array(allHitsAverage)
    
    return allCams_meanRayDist    
    
    
folder = r"angel_rpLidar"



mesh = trimesh.load(folder + r'\currObj.obj')

#  initial required scalings to get the the unknown scale to a larger one, if the diffence between the real world one and it is too large the algorithm would fail.
initialScale = 1000
errorFixing = 1

additional_scale = 1

# load all necessary data
cam_pos = np.loadtxt(folder + r"\cameraPos.txt", delimiter=",")*initialScale*errorFixing * additional_scale
cam_norm = np.loadtxt(folder + r"\cameraNorm.txt", delimiter=",")

cam_perp = np.loadtxt(folder + r"\cameraPerp.txt", delimiter=",")

sparsePC = np.loadtxt(folder + r"\sparsePC.txt", delimiter=',')

sparsePC[:,1:] = sparsePC[:,1:] * initialScale*errorFixing * additional_scale

mesh.vertices *=errorFixing

mesh.vertices *=additional_scale

mesh.vertices *=initialScale


test = mesh.vertices





# Load real life LiDAR data and if there's nan values remove them

lidarDistances = np.loadtxt(folder + r"\LidarDist.txt",delimiter=",")
lidarDistances = lidarDistances.T
lidarDistances[lidarDistances == 0] = np.nan


# calculate the mean distances for each reading
lidarDistances_means = np.nanmean(lidarDistances, axis=0)
lidarDistances_means = lidarDistances_means


# calculated initial translation between the camera and lidar
dist_camLidar_real = np.array([0.2,150,-13])
#  initial rotation if between camera and lidar if needed
rot_camLidar_real = np.array([[1,0,0],[0,1,0],[0,0,1]]) 
scalingFactor_mean = 1

# number of iterations where the scale and uncertainty are recalculated
numIterations =10

scalingFactor_final = 1

rot_rayCast_max =10
rot_rayCast_delta = 1

# Lidar uncertainty in mm

distNoise =0.001773400863



#  number or lidar rays per capture
numLidarRays =24


numberOfCamsWithDist = len(lidarDistances_means)

allScales_plus = np.zeros([numberOfCamsWithDist,1])
allScales_minus = np.zeros([numberOfCamsWithDist,1])



init_cam_pos = cam_pos.copy()

allFakeDist = []
#  iterative approach for casting rays and adding and subtracting the distance noise for calculating the uncertainty
for k in range(0,len(lidarDistances_means)):
    for l in range(0,2):
        
        lidarDistances_means_test = lidarDistances_means.copy()
        
        lidarDistances_means_test = lidarDistances_means_test[:numberOfCamsWithDist]
        
        
        
        
        if l == 0:
            distNoise = distNoise
        else:
            distNoise = -distNoise

        
        mesh.vertices = test.copy()
        cam_pos = init_cam_pos
        
        lidarDistances_means_test[k] += distNoise

        cam_pos = cam_pos[:numberOfCamsWithDist]
        
        cam_norm = cam_norm[:numberOfCamsWithDist]
        cam_perp = cam_perp[:numberOfCamsWithDist]
        
        scalingFactor_final = 1

        scalingFactor_mean = 1
        
        scalingFactors = np.ones([len(cam_pos)])
        
        
        lidar_pos = cam_pos + cam_norm*dist_camLidar_real[2] + cam_perp*dist_camLidar_real[1] + np.cross(cam_norm,  cam_perp) *dist_camLidar_real[0]
        
        posNeg = 1
        
        prevRMSE = 9999
        i = 0
        #  do this for a number of iterations or when the RMSE becomes small enough
        for i in range(numIterations):
            

            
            lidar_pos = lidar_pos * scalingFactors[:, np.newaxis] 
            

            mesh.vertices*=scalingFactor_mean
            # if i == 0:
                
            #     lidar_pos = cam_pos + dist_camLidar_real
            # else:
            #     lidar_pos = cam_pos
            # allCams_meanRayDist = []
            
            
            
            allCams_meanRayDist = getRayCastDist(mesh,lidar_pos, cam_norm, cam_perp, rot_rayCast_max, rot_rayCast_delta)
                      
            scalingFactors = []
            scalingFactors = lidarDistances_means_test/allCams_meanRayDist

        
            scalingFactor_mean = scalingFactors.mean()
            currRMSE = np.sqrt(((allCams_meanRayDist - lidarDistances_means_test) ** 2).mean())
            print(currRMSE)
            
            if (currRMSE < prevRMSE):
                scalingFactor_final *= scalingFactor_mean
                
                prevRMSE = currRMSE
            else:
               
                break
            
        
    
        allFakeDist.append(allCams_meanRayDist.mean())    
        print("For lidar pos " + str(k) +  " final scaling factor " + str(scalingFactor_final))
        
        
        if l == 0:
            allScales_plus[k] = scalingFactor_final
        
        else:
            allScales_minus[k] = scalingFactor_final


# Use covariance propagation of distances uncertainty to the calculated scale
print("Average plus: " + str(allScales_plus.mean()) )

print("Average minus: " + str(allScales_minus.mean()))

allScales_finalDist = (allScales_plus - allScales_minus)/(2*np.abs(distNoise) )


matrixLine = np.ones((1,len(allScales_finalDist)))*((np.abs(distNoise)**2)/numLidarRays)

fullCombined = matrixLine.copy()
for j in range(0,len(allScales_finalDist)-1):
    fullCombined = np.concatenate([fullCombined,matrixLine])
#  covariance matrix
fullCombined[~np.eye(fullCombined.shape[0],dtype=bool) == True] = 0
             
scaleTest_temp = np.dot(fullCombined, allScales_finalDist)            
scaleTest = np.dot(np.transpose(allScales_finalDist), scaleTest_temp)    

scaleTest2 = np.transpose(allScales_finalDist).dot(fullCombined).dot(allScales_finalDist)
       

print(" RTK analytical matrix = " + str(scaleTest2))


contr2Noise = 0
for k in range(len(allScales_finalDist)):
    contr2Noise+=(allScales_finalDist[k]**2 * ((np.abs(distNoise)**2)/numLidarRays) )

    



# Visualize raycasts




#allHits_allCams = np.zeros([len(cam_pos), int((2*rot_rayCast_max)/rot_rayCast_delta)])
# 
#for i in range(0, len(cam_pos)):
#    
#    ray_origins = np.array([cam_pos[i,:]])
#
#    ray_directions = np.array([cam_norm[i,:]])
#    
#    ray_perp = cam_perp[i,:]
#    ray_norm = cam_norm[i,:]
#    
#
#    
#    allHits_currCam = np.zeros([1, int((2*rot_rayCast_max)/0.25)])
#    count = 0
#    for j in np.arange(-rot_rayCast_max, rot_rayCast_max, rot_rayCast_delta):
#        rot_angle = np.deg2rad(j)
#        
#        cam_norm_curr = np.dot(rotation_matrix(ray_perp, rot_angle), ray_norm)
#        
#        ray_directions_curr = np.array([cam_norm_curr])
#    
#    
#        locations, index_ray, index_tri = mesh.ray.intersects_location(
#        ray_origins=ray_origins,
#        ray_directions=ray_directions_curr)
#        
#        if len(locations) != 0:
#            allHitDist = np.sqrt(np.sum((locations-ray_origins)**2,axis=1)) 
#            hitArg = np.argmin(allHitDist)
#            
#            allHits_currCam[0,count] = allHitDist[hitArg]
#        
#        
#        ax.scatter(cam_pos[i,0], cam_pos[i,1], cam_pos[i,2], c='red', marker=".", edgecolors='none')
#        count+=1
#        
#        
#    #    ax.quiver(cam_pos[i,0], cam_pos[i,1], cam_pos[i,2], cam_norm[i,0], cam_norm[i,1], cam_norm[i,2],length=10, normalize = True)
#        
#        if len(locations) != 0:
#            ax.scatter(locations[hitArg,0], locations[hitArg,1], locations[hitArg,2], c='blue', marker="o", edgecolors='none')
#            ax.plot([cam_pos[i,0], locations[hitArg,0]], [cam_pos[i,1], locations[hitArg,1]], [cam_pos[i,2], locations[hitArg,2]])
#        else:
#            print("Camera " + str(i) + " does not hit")
#    
#    
#    
#    
#    allHits_allCams[i] = allHits_currCam
        
        
#allCams_meanRayDist = np.true_divide(allHits_allCams.sum(1),(allHits_allCams!=0).sum(1))        


