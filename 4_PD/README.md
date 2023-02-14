# README
## Figure 11 and Figure 14

To generate Figures 11 and 14, first download the five .mrcs files in the `ManifoldEM_ESPER/0_Data_Inputs/2_PDs_2D/` folder and place them in the `SS2_PDs_Pristine` folder within the current directory.

Then run the following commands with the proper Anaconda environment sourced:

- `bash 1_Dist_Batch.sh`
- `bash 2_DM_Batch.sh`

Figure assets will render to the folders `Data_Distances` and `Data_Manifolds`.

All scripts and data are selected from the ESPER repository cited in our main text.