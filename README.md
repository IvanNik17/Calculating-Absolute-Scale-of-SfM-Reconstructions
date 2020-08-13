# Calculating Absolute Scale of SfM Reconstructions
 Repository for code featured in two papers for calculating the absolute scale of SfM reconstructions, as well as the uncertainty of said scale


# Scale From Distance Measurements

The python code used in the development of the paper **Nikolov, I., & Madsen, C. B. (2020). Calculating Absolute Scale and Scale Uncertainty for SfM Using Distance Sensor Measurements: A Lightweight and Flexible Approach. In Recent Advances in 3D Imaging, Modeling, and Reconstruction (pp. 168-192). IGI Global.**

**REQUREMENTS**
Numpy and [trimesh](https://github.com/mikedh/trimesh)

**Testing DATA**
An angel statue 3D reconstruction is given, as input data. The folder contains:
1. Images used to 3D reconstruct the statue
2. Camera normals, perpendiculars and positions for each image
3. Reconstructed mesh in .obj format
4. LiDAR distances captured from each camera position
5. Additional files containing separated the vertices, faces and normals, as well as a sparse point cloud

The transformation between the camera and LiDAR is given in the code, this is done once as a calibration step before capturing the data and the two should not be moved 
The uncertainty of the LiDAR is also given. The rpLiDAR A1 is used.


# Scale From Positioning GPS Measurements

The python code used in the development of the paper **Nikolov, I., & Madsen, C. B. (2019). Performance Characterization of Absolute Scale Computation for 3D Structure from Motion Reconstruction. In VISIGRAPP (5: VISAPP) (pp. 884-891).**

**REQUREMENTS**
Numpy

**Testing DATA**
For demonstrating the calculation of transformation matrix and the uncertainty, two files containing captured real world positions and calculated SfM positions of cameras are given. The uncertainty of the GPS used for the testing scenario is also given 
