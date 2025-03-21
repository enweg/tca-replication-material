# Replication of Section 5.3: The Role of Inflation Expectations in DSGEs
## How to replicate? 
- To replicate the results, first open the 'run001_sw2007.m' and 'run003_sw2007.m' files and adjust the path to the Dynare installation on line 12. 
- After adjusting the installation path of Dynare, running the entire files
  will replicate the results. Results are saved in the 'output' folder as 
  CSV files. 
- Figures can be obtained by running 'run002_sw2007_figures.m' and 'run004_sw2007_figures.m' which produce figures for the two decompositions separately. 
- Figures can also be obtained by using the `run99_figures.jl` script, which uses Julia instead of Matlab. This file also creates the figure of the paper in which both decompositions are depicted in a single figure. To run this file, Julia must be installed. For guidance on this, see 
  - Julia installation: https://julialang.org/downloads/

## Folder Structure
- 'SW2007' contains the replication material for Smets & Wouters (2007)
    - See the comments in 'SW2007.mod' for the original source and information 
      on calibration.
- 'functions' contains all utility functions for transmission channel 
  analysis in DSGE models fitted using Dynare.
- 'output' contains the output of 'run001_sw2007.m' and 'run003_sw2007.m'.
- 'plots' contains all plots created using 'run002_sw2007_figures.m', 'run004_sw2007_figures.m' and 'run99_figures.jl'.

## System Settings 
- The scripts were run on a MacBook Air M1 running Matlab R2024b 
  Update 5 (24.2.0.2863752) and Dynare 6.3 (arm64 version).
- The Julia version is 1.11.4; all package versions are fixed and can be found in the 'Project.toml' file. 
## Sources
- Dynare replication files for Smets & Wouters (2007) were obtained from Johannes Pfeifer's replication repository (https://github.com/JohannesPfeifer/DSGE_mod). 

## References
- Smets, Frank and Wouters, Rafael (2007): "Shocks and Frictions in US Business Cycles: A Bayesian DSGE Approach", American Economic Review, 97(3), 586-606
