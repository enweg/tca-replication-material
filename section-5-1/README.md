# Replication of Section 5.1: Forward Guidance of Monetary Policy
## Data
- All data is contained in the 'data' folder. The data was obtained from the 
  replication material of McKay & Wolf (2023). 
    - Link to material: https://doi.org/10.3982/ECTA21045
    - Last Accessed: 17/06/2023 (dd/mm/yyyy)

## Replication

- Before replication, make sure that Julia is installed and that the 'DrWatson.jl' package is installed. The following instructions might be helpful. 
    - Julia installation: https://julialang.org/downloads/
    - DrWatson installation: https://juliadynamics.github.io/DrWatson.jl/stable

- The results of Section 5.1 can now be replicated by running the following scripts: 

1. Run the '01-results.jl' script to estimate the model and compute the transmission effects. The transmission effects are stored as JLD2 files in the 'output' directory. 
2. Run the '02-plots.jl' script to create the plots of Section 5.1. The plots are saved in the 'plots' folder as pdf files. 

## Other Directories and Files
- Helper functions are stored in the 'scripts' folder. These helper functions are used for the estimation of the VAR, the identification via internal instruments, the computation of transmission effects, and for plotting. 
## System Settings 
 - The scripts were run on a MacBook Air M1 running Julia 1.11.2.
 - The specific package versions can be found in the 'Project.toml' file. 
