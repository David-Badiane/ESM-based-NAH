# ESM-based-NAH
Equivalent Source Method for Near-field Acoustic Holography. Project for the Project Course of Polimi Music and Acoustic Engineering

[download acceleration and pressure measurements](https://polimi365-my.sharepoint.com/:f:/g/personal/10743504_polimi_it/ElZaXDe-48NIt7RkcwDgemUB_qcxJ5-4a2B1z9xCsb_G9A?e=9uiU2Q)

## Introduction

<p align="center"><img src="https://github.com/David-Badiane/ESM-based-NAH/blob/main/Figures/simulationSetupPicDiagram.png" width=70% height=70% centered>
  
NAH is a powerful inverse method for the estimation of the velocity and pressure field of a vibrating body from pressure measurements. The pressure field is captured in the close proximity of the object through a planar grid of microphones, called hologram. The pressure acquisitions are then used to recover the velocity field on the source by back-propagation (through Green's functions and Kirchhoff-Helmoltz (KH) integral theory). By taking near-field pressure measurements we capture the evanescent waves, which are essential to reconstruct the surface vibration velocity.

This projects adopts equivalent source method (ESM) based NAH to a violin back-plate. Since ESM based NAH avoids the discretization of the KH integral and can be applied to arbitrary geometries, it is particularly suited to our case study.
  
  The basic idea of ESM is that the acoustic field of the vibrating object can be approximated to the field produced by a set of virtual point sources, called equivalent sources, located on a surface right behind the vibrating object. Finding the complex weights of the equivalent sources from the pressure measurements means solving an undetermined problem. ESM is a two steps method: first, we obtain the weights through regularization techniques; then, the weights are back-propagated to the reconstruction surface by the means of the proper driving functions (based on 3D Green's functions) to estimate the vibration velocity of the object.
  
  <p align="center"><img src="https://github.com/David-Badiane/ESM-based-NAH/blob/main/Figures/ESM.png" width=70% height=60% centered>
  
  The free parameters of the method are the number of equivalent sources, their position and their distance with respect to the vibrating object. The choice of the free parameters has a strong impact on the method accuracy, so a careful and correct choice of those parameters is necessary. It is important to note that the parameters space of ESM is so wide that individuating the best solution is still an unresolved issue. Here a block diagram of the method, where Gp is the Green's function matrix hologram to virtual sources and Gv is the normal gradient of the Green's functions matrix virtual sources to reconstruction surface.

  In our work we propose the application of ESM on both synthetic and experimental data, taken by us at the anechoic chamber of Politecnico di Milano, Cremona, Italy and in the Musical Acoustics Lab of Politecnico di Milano, Museo del Violino - Cremona - Italy.  The code is MATLAB® based.

## Getting Started

### Dependencies

* Prerequisites: MATLAB® 2020 or later installed.
* Libraries : 
  * [Regularization Tools](https://it.mathworks.com/matlabcentral/fileexchange/52-regtools) by Per Christian Hansen - included in the repo.
  * [stlread](https://it.mathworks.com/matlabcentral/fileexchange/22409-stl-file-reader) function - included in the repo
  
### Executing programs

* **Main_experimental.m** - apply ESM on experimental data;
* **Main_synthetic.m** - apply ESM on synthetic data;
* **preprocessing_pressure.m** is used to compute the pressure data in the hologram point from measurement data;
* **velocityGroundTruth.m** computes the grid with the velocity ground truth value and coordinates, from the .csv file of measurement data, for every eigen frequencies;
* **violinMeshGenerator.m** generates violinMesh from large .stl file;
* **virtualPointsGenerator.m** generates the grid of virtual points of various geometries;
    
#### Main Synthetic
Applies ESM on simulated data.
Algorithm:
* function *importData* - Extracts the violin model geometry, velocity fieldsand hologram fields and points for each eigenfrequency;
* compute the *normal vectors* on the reconstruction surface, i.e.  violinMesh
* for each eigenfrequency, the user can choose to:
    * set retreat distance (RD) between virtual points and violin mesh and compute ESM estimation.  Save their metrics; 
    * optimize the estimation on the free parameters of the virtual points(RD - scaleX - scaleY), using the ground-truth data.  Save the optimization results;
    * see the optimization results, save their metrics;
ESM is carried on by the Matlab function *applyESM.m*.
    
#### Main Experimental
Applies ESM on experimental data.
The structure is the same as Mainsynthetic.m, but we have to treat different data, coming from our experimental measurements. The violin mesh is chosen from a set of grids located in the directory violinMeshes, which come from the laser scan of our test specimen. The hologram is filled with our experimental data. The velocity groundtruth is obtained from 27 measured points on the test sample through an hammer-accelerometer system.
 
#### preprocessing_pressure
Obtain the hologram entries for Main Experimental
The hologram geometry is modeled on the experimental setup and filled withexperimental data.  We recorded 9 channels (8 measurement mics + 1 reference mic) for 8 different heights. For each height, we recorded 7 independent takes of the plate excitation. 
In order, for each channel, the script:
* Imports the audio and force signals;
* Isolates each take synchronizing the signal with the reference mic acquisition;
* Applies exponential filters;
* Compute the H1 estimator;
* Perform  singular  values  decomposition  (SVD)  on  the  H1  estimator  toreduce noise;
* Peak analysis and pressure fields retrieval;
* save results in .csv and .mat files;
    
#### velocityGroundTruth
Obtain velocity ground−truth for Main Experimental.
The experimental velocity ground-truth is the H1 estimator of the mobility on 27 points of the violin back-plate. The  data  are  obtained  thourgh  an  hammer - accelerometer measurement system.
In order, for each point, the script:
* reads acceleration and force files;
* applies exponential filter;
* computes the mobility and the H1 estimator of the mobility; 
* applies SVD on the H1 estimator to reduce noise;
* peak analysis and velocity fields retrieval;
* save results in .csv and .mat files;    
    
#### violinMeshGenerator
Create and handle violin mesh. 
The violin meshes have been generated from a 4 million points .stl file obtained by laser scanning the violin back-plate.  
The code is divided in three sections that are meant to be opportunistically used:
* *first section* - read the .stl mesh and downsample it to a 128x128 gridwith the functiondownsamplingregular.m.  Save it back to a .csv file;
* *second section* - read a target grid and symmetrize it manually (user inputand figures) with respect to the y axis of the violin;
* *third section* - downsample target grid and save it
    
#### virtualPointsGenerator
Create and handle virtual sources grids.
This script generates grids of virtual points choosing between different geometries,  sparsity  levels  and  spatial  sampling  parameters.   The  meshes  are generated by user input, more in particular its possible to decide between the following geometries:
    - rectangular;                 - ellipsoidal;                    - circular + violin border;
    - ellipsoidal + violin border; - inner points + violin border;   - inner only
    - rectangular + violin border + in-ner points;
For each chosen geometry, the user can decide the level of sparsity of the grid,adjust its size and choose the downsampling parameters.  The functions used forthe virtual points generation are in the folder\functions\virtualPointsGenerators.
    
### Folders
* **Exp_Measurements**: folder containing the experimental data;
* **CSV**: synthetic data of pressure and velocity (.csv);
* **Data**: file .mat and .csv of sensible data of the code, useful to run it without reading all the experimental or synthetic data everytime the code is launched;
* **functions**: contains all the function used in the scripts;
* **VioliniMeshes**: different dimensions of grids for the violin mesh;
* **VPGrids**: virtual point grids;

### Important functions
* **applyESM.m** this function applies all ESM steps = (problem solution) + metrics + figures  
* **errorVelocity.m** this function computes the error through different metrics and save it into a struct.
* **getVelocityGroundtruth.m** this function calculates the array velocity into surface.
* **getVirtualPoints.m** this function calculates the virtual points grid.
* **Green_matrix.m** this funtion compute the Green's matrix.
* **importData.m** this function read the csv pressure and normal velocity extract their data and sotre the in a cell arrays
* **normalGradient.m** this functinon compute the normal gradient of the Green's matrix

## Help

Further infos can be found in the Documentation of our code, that you can find in this repository.

## Authors

David Giuseppe Badiane (davidgiuseppe.badiane@mail.polimi.it), Alessio Lampis (alessio.lampis@mail.polimi.it), Andrea Eugenio Losi(andreaeugenio.losi@mail.polimi.it)

