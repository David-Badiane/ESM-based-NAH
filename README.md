# ESM-based-NAH
Equivalent Source Method for Near-field Acoustic Holography. Project for the Project Course of Polimi Music and Acoustic Engineering

Raw measurements : https://polimi365-my.sharepoint.com/:f:/g/personal/10743504_polimi_it/ElZaXDe-48NIt7RkcwDgemUB_qcxJ5-4a2B1z9xCsb_G9A?e=9uiU2Q

## Introduction

<p align="center"><img src="https://github.com/David-Badiane/ESM-based-NAH/blob/main/Figures/simulationSetupPicDiagram.png" width=70% height=70% centered>
  
NAH is a powerful inverse method for the estimation of the velocity and pressure field of a vibrating body from pressure measurements. The pressure field is captured in the close proximity of the object through a planar grid of microphones, called hologram. The pressure acquisitions are then used to recover the velocity field on the source by back-propagation (through Green's functions and Kirchhoff-Helmoltz (KH) integral theory). By taking near-field pressure measurements we capture the evanescent waves, which are essential to reconstruct the surface vibration velocity.

This projects adopts equivalent source method (ESM) based NAH to a violin back-plate. Since ESM based NAH avoids the discretization of the KH integral and can be applied to arbitrary geometries, it is particularly suited to our case study.
  The basic idea of ESM is that the acoustic field of the vibrating object can be approximated to the field produced by a set of virtual point sources, called equivalent sources, located on a surface right behind the vibrating object. Finding the complex weights of the equivalent sources from the pressure measurements means solving an undetermined problem. ESM is a two steps method: first, we obtain the weights through regularization techniques; then, the weights are back-propagated to the reconstruction surface by the means of the proper driving functions (based on 3D Green's functions) to estimate the vibration velocity of the object.
  The free parameters of the method are the number of equivalent sources, their position and their distance with respect to the vibrating object. The choice of the free parameters has a strong impact on the method accuracy, so a careful and correct choice of those parameters is necessary. It is important to note that the parameters space of ESM is so wide that individuating the best solution is still an unresolved issue. 

    In our work we propose the application of ESM on both synthetic and experimental data, taken by us at the anechoic chamber of Politecnico di Milano, Cremona, Italy and in the Musical Acoustics Lab of Politecnico di Milano, Museo del Violino - Cremona - Italy. 
  The code is MATLAB® based.

## Getting Started

### Dependencies

* Prerequisites: MATLAB® 2020 or later installed.
* Libraries: Regularization Tools by Hansen.

### Executing program

* *Main_experimental.m* - apply ESM on experimental data.
* *Main_synthetic.m* - apply ESM on synthetic data.
* *platePoints.m* imports from the .stl scanned mesh of the violin.
* *preprocessing_pressure.m* is used to compute the pressure data in the hologram point from measurement data. Loads the audio files, cuts the different acquisition for every channel and apply to them a temporal exponential filter. From this cleaned files calculates the H1 estimator and extract the pressure field for every eigen freqencies.
* *symmetrizeViolinMesh.m* symmetrizes the point of the violin scanned grid.
* *velocityGroundTruth.m* computes the grid with the velocity ground truth value and coordinates, from the .csv file of measurement data, for every eigen frequencies.
* *violinMeshGenerator.m* creates the violin mesh from a grid of points.
* *virtualPointsGenerator.m* generates the grid of virtual points through the user input arguments.

### Folders
* **acq_11_06**: velocity acquisitions of the accelerometer (.csv).
* **CSV**: synthetic data of pressure and velocity (.csv).
* **Data**: file .mat and .csv of sensible data of the code, useful to run it without reading all the experimental or synthetic data everytime the code is launched.
* **functions**: contains all the function used in the scripts.
* **VioliniMeshes**: different dimensions of grids for the violin mesh.
* **VPGrids**: virtual point grids.

### Important functions
* *applyESM.m* this function applies all ESM steps = (problem solution) + metrics + figures  
* *errorVelocity.m* this function computes the error through different metrics and save it into a struct.
* *getVelocityGroundtruth.m* this function calculates the array velocity into surface.
* *getVirtualPoints.m* this function calculates the virtual points grid.
* *Green_matrix.m* this funtion compute the Green's matrix.
* *importData.m* this function read the csv pressure and normal velocity extract their data and sotre the in a cell arrays
* *normalGradient.m* this functinon compute the normal gradient of the Green's matrix

## Help

Any advise for common problems or issues.
```
command to run if program contains helper info
```

## Authors

David Giuseppe Badiane (davidgiuseppe.badiane@mail.polimi.it), Alessio Lampis (alessio.lampis@mail.polimi.it), Andrea Eugenio Losi(andreaeugenio.losi@mail.polimi.it)


## License

This project is licensed under the [NAME HERE] License - see the LICENSE.md file for details
