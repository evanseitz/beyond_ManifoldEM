# Beyond ManifoldEM
This repository contains code for generating results shown in our paper **Beyond ManifoldEM: Geometric relationships between manifold embeddings of a continuum of 3D molecular structures and their 2D projections** (Seitz, Frank* and Schwander*): https://www.TBD.com. This work was developed in the Frank research group at Columbia University in collaboration with Peter Schwander at the University of Wisconsin-Milwaukee (UWM).

The algorithms presented here in their current form are developed for analyzing synthetic data. Custom synthetic datasets can be generated using code in our [previous repository](https://github.com/evanseitz/cryoEM_synthetic_continua) with additional information provided in the corresponding [previous paper](https://www.biorxiv.org/content/10.1101/864116v1): **Simulation of Cryo-EM Ensembles from Atomic Models of Molecules Exhibiting Continuous Conformations** (Seitz, Acosta-Reyes, Schwander, Frank*). We have also deposited our complete synthetic dataset in the [IEEE DataPort](https://ieee-dataport.org/documents/manifoldem-esper-data-and-code-repository), which is a mirror of the 2019 repository with the addition of several refinements for processing the provided data.

Users interested in applying these methods to experimental datasets should consult the [ManifoldEM Matlab suite](https://github.com/GMashayekhi/ManifoldEM_Matlab) or [ManifoldEM Python suite](https://github.com/evanseitz/ManifoldEM_Python), with the latter featuring a comprehensive GUI, user manual and robust automation strategies. As well, our [ManifoldEM ESPER method](https://github.com/evanseitz/ESPER) can be used to further refine ManifoldEM outputs given that certain prerequisites are met, with detailed information provided in our [ESPER paper](https://ieeexplore.ieee.org/document/9773954): **Recovery of Conformational Continuum From Single-Particle Cryo-EM Images: Optimization of ManifoldEM Informed by Ground Truth** (Seitz, Acosta-Reyes, Maji, Schwander, Frank*).

### Environment:
First, install [Anaconda](https://docs.anaconda.com/anaconda/install), and with Anaconda sourced, create a new Anaconda environment:

`conda create -n beyondMEM python=3`

Next, activate this environment via `conda activate beyondMEM`, and install the following packages:

- `pip install numpy`
- TBD

### Usage:
Detailed instructions and comments for all procedures are provided in the code. Within this repository, we have also supplied a sample of pristine data that can be used for reproducing outputs and figures.

## Attribution:
If this code is useful in your work, please cite:

TBD ZENODO LINK

### License:
Copyright (C) 2018-2023 Evan Seitz

The software, code sample and their documentation made available on this website could include technical or other mistakes, inaccuracies or typographical errors. We may make changes to the software or documentation made available on its web site at any time without prior notice. We assume no responsibility for errors or omissions in the software or documentation available from its web site. For further details, please see the LICENSE file.
