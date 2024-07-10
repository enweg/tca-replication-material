# Replication of Section 5.3: The Role of Inflation Expectations in DSGEs
## How to replicate? 
- To replicate the results, first open the 
  'run01_smets_wouters_tca.m' file and adjust the path to the Dynare installation 
  on line 10. 
- After adjusting the installation path of Dynare, running the entire file
  will replicate the results. Results are saved in the 'output' folder as 
  CSV files. 
- Figures can be obtained by running 'run02_smets_wouters_figures.m'
  The script saves four figures in PDF format in the 'plots' folder. 
  Each figure corresponds to the same analysis, with the only difference 
  being the maximum horizon for which the transmission effects are plotted.
- Figures can also be obtained by using the `run99_figures.jl` script, which
  uses Julia instead of Matlab to create figures. To run this file, DrWatson.jl
  and Julia must be installed. For guidance on this, see 
  - Julia installation: https://julialang.org/downloads/
  - DrWatson installation: https://juliadynamics.github.io/DrWatson.jl/stable

## Folder Structure
- 'Smets_Wouters_2007' contains the replication material for Smets & Wouters (2007)
- 'functions' contains all utility functions for transmission channel 
  analysis in DSGE models fitted using Dynare.
- 'output' contains the output of 'run01_smets_wouters_tca.m'.
- 'plots' contains all plots created using 'run02_smets_wouters_figures.m'
  and 'run99_figures.jl'.

## System Settings 
- The scripts were run on a MacBook Air M1 running Matlab R2023b 
  Update 4 (23.2.0.2428915) and Dynare 5.5 (arm64 version).
## Sources
- Dynare replication files for Smets & Wouters (2007) were obtained from 
  Johannes Pfeifer's replication repository (https://github.com/JohannesPfeifer/DSGE_mod). 

## References
- Smets, Frank and Wouters, Rafael (2007): "Shocks and Frictions in 
  US Business Cycles: A Bayesian DSGE Approach", American Economic Review, 
  97(3), 586-606
