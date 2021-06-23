# ESM-based-NAH
Equivalent Source Method for Near-field Acoustic Holografy. Project for the Project Course of Polimi Music and Acoustic Engineering

Raw measurements : https://polimi365-my.sharepoint.com/:f:/g/personal/10743504_polimi_it/ElZaXDe-48NIt7RkcwDgemUB_qcxJ5-4a2B1z9xCsb_G9A?e=9uiU2Q

## Introduction

NAH is a powerful inverse method that allows the estimation of the velocity and pressure field of a vibrating body from pressure measurements. Typically, the pressure field is captured in the close proximity of the object through a planar grid of microphones, called hologram. The pressure acquisitions are successively used to recover the velocity field on the source by back-propagation. By taking near-field pressure measurements we capture the evanescent waves, which are essential to reconstruct the surface vibration velocity.

We adopt a computationally efficient approach to face the NAH which is represented by the equivalent source method (ESM), also known as wave superposition method. The basic idea is that the acoustic field of the vibrating object can be approximated to the field produced by a set of virtual point sources, called equivalent sources, located in a surface in the close proximity of the vibrating object. It's usually necessary to resort to regularization techiques to solve the inverse problem of finding the equivalent sources weights from the pressure measurements. Then, the weights are back-propagated to the reconstruction surface by the means of the proper driving functions to estimate the vibration velocity of the object. ESM is attractive because it avoids the Kirchhoff-Helmholtz integral discretization and does not loose the capability of reconstructing over non-separable geometries. The free parameters of the method are the number of equivalent sources, their position and their distance with respect to the vibrating object. Free parameters choice has a strong impact on the method accuracy so ESM involves a careful and correct choice of those parameters. It is important to note that the parameters space of ESM is so wide that individuating the best solution is still an unresolved issue. 

In our work we propose the application of ESM based NAH on MATLAB®  to a violin back plate in order to retrieve its normal velocity. The method has scripts applied on both synthetic and experimental data.

## Getting Started

### Dependencies

* Prerequisites: MATLAB® 2020 or later installed.
* Libraries: Regularization Tools by Hansen.

### Executing program

* *Main_experimental.m* applies the ESM method to experimental data.
* *Main_synthetic.m* apllies the ESM method to synthetic data.
* *platePoints.m* imports from the .csv file the plate points.
* *preprocessing_pressure.m* is used to compute the pressure data in the hologram point from measurement data. Loads the audio files, cuts the different acquisition for every channel and apply to them a temporal exponential filter. From this cleaned files calculates the H1 estimator and extract the pressure field for every eigen freqencies.
* *symmetrizeViolinMesh.m* symmetrizes the point of the violin scanned grid.
* *velocityGroundTruth.m* computes the grid with the velocity ground truth value and coordinates, from the .csv file of measurement data, for every eigen frequencies.
* *violinMeshGenerator.m* creates the violin mesh from a grid of points.
* *virtualPointsGenerator.m* generates the grid of virtual points through the user input arguments.

## Help

Any advise for common problems or issues.
```
command to run if program contains helper info
```

## Authors

Contributors names and contact info

ex. Dominique Pizzie  
ex. [@DomPizzie](https://twitter.com/dompizzie)

## Version History

* 0.2
    * Various bug fixes and optimizations
    * See [commit change]() or See [release history]()
* 0.1
    * Initial Release

## License

This project is licensed under the [NAME HERE] License - see the LICENSE.md file for details

## Acknowledgments

Inspiration, code snippets, etc.
* [awesome-readme](https://github.com/matiassingers/awesome-readme)
* [PurpleBooth](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2)
* [dbader](https://github.com/dbader/readme-template)
* [zenorocha](https://gist.github.com/zenorocha/4526327)
* [fvcproductions](https://gist.github.com/fvcproductions/1bfc2d4aecb01a834b46)
