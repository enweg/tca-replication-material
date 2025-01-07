# Replication of Section 5.2: Anticipation Effects of Government Spending
## Data

- All data is contained in the 'data' folder. The 'raw' folder within the 'data' folder contains all data obtained from other sources, while the data in the root 'data' folder contains compiled data. The 'README.md' in the 'raw' subfolder includes information about sources and access dates.

## Replication

- Before replication, make sure that Julia is installed and that the 'DrWatson.jl' package is installed. The following instructions might be helpful. 
   - Julia installation: https://julialang.org/downloads/
   - DrWatson installation: https://juliadynamics.github.io/DrWatson.jl/stable

- The results of Section 5.2 can now be replicated by running the following scripts: 

1. Run the '00-data-transforms.jl' script to create the compiled datasets 'spending_quarter.csv' and 'data_lp.csv' in the 'data' folder. The former file contains the constructed quarter military spending series, while the latter file contains all data needed for the local projections. The script also creates three plots in the 'plots' directory. These plots compare the constructed defense spending data to official FRED data (they are not included in the original paper). 
2. Run the '01-estimation.jl' script to estimate the local projections and calculate the transmission effects. The results are saved in the 'output' directory as JLD2 files.
3. Run the '02-plot.jl' file to create the plots of Section 5.2. Plots are saved in the 'plots' directory as PDF files.

## Other Directories and Files
- Helper functions are stored in the 'scripts' folder. These helper functions are used for the estimation of the LP, and for plotting. 
## System Settings 
 - The scripts were run on a MacBook Air M1 running Julia 1.11.2. 
 - Specific package versions can be found in the 'Project.toml' file. 
 - All packages in development (those loaded from GitHub) are fixed to a specific commit. To see the specific commit, type `]` followed by `status` in the Julia REPL.
