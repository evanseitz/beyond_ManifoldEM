# README
## Atomic-coordinate structures (ACS): Figures 6,8

To generate Figure 6 and assets for Figure 8, run the following commands with the proper [Anaconda](https://docs.anaconda.com/anaconda/install) environment sourced:

- `python DifussionMaps_PDB.py`

Figure assets will render to the folder `_figure_assets_`.

All scripts and data are modifications of those found in the [ESPER repository](https://github.com/evanseitz/ManifoldEM_ESPER). In lieu of calculating distances manually, we have provided the distance matrix normally produced using `GetDistances_PDB.py`. However, this missing step can be reproduced using the ACS `.pdb` files located in the `SS2/ACS/Pristine` folder. For generating PDB files from scratch, follow the instructions in our [synthetic generation repository](https://github.com/evanseitz/cryoEM_synthetic_generation).