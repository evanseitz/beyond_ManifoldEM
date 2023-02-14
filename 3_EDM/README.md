# README
## Figures 7, 8

To generate Figure 7 and assets for Figure 8, run the following commands with the proper [Anaconda](https://docs.anaconda.com/anaconda/install) environment sourced:

- `python DifussionMaps.py`

Figure assets will render to the folder `_figure_assets_`.

All scripts and data are modifications of those found in the [ESPER repository](https://github.com/evanseitz/ManifoldEM_ESPER). In lieu of calculating distances manually, we have provided the distance matrix normally produced using `Distances.py`. However, this missing step can be reproduced using the EDM `.mrc` files located on our [IEEE repository](https://ieee-dataport.org/documents/manifoldem-esper-data-and-code-repository), which are too large (64 MB each) to store via GitHub. These exact EDM files can also be generated from scratch following directions in our [synthetic generation repository](https://github.com/evanseitz/cryoEM_synthetic_generation).