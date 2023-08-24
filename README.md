# Adaptive sampling method to monitor “low-risk” pathways with limited surveillance resources

This repository contains the R implementation for the sampling method introduced in the paper *Adaptive sampling method to monitor “low-risk” pathways with limited surveillance resources*, by Thao P. Le, Thomas K. Waring, Howard Bondell, Andrew Robinson, and Christopher M. Baker.

**Note**: some figures require the **allocate.r** function which is taken from [*Allocating surveillance resources to reduce ecological invasions: maximizing detections and information about the threat* by Andrew Robinson, Mark A. Burgman, Rob Cannon](https://doi.org/10.1890/10-0195.1)

Aside from allocate.r, the code in this repository is written by Thao P. Le, Thomas K. Waring, and Christopher M. Baker.

The code is in R, and requires packages: tidyverse,  tidyr, memoise, and Rmpfr.

| Code file | Description |
| --------- | ---------- |
|**allocate.r** | Robinson et al's sampling method from [*Allocating surveillance resources to reduce ecological invasions: maximizing detections and information about the threat*](https://doi.org/10.1890/10-0195.1) |
|**core_functions.R** | Contains functions needed to calculate the minimum recommended sample size for the low risk sampling method |
|**generate_df.R** | Calls core_functions.R to create a dataframe (and csv file) containing combinations of parameters and subsequent recommended sample size |
| **plots_method.R** | Plots figures 3 and 8 in the paper |
| **plot_probability_status.R** | Plots figure 7 in the paper |
| **plots_scenario.R** | Runs very low leakage, low leakage, and high leakage scenarios, sampled with our method (Figures 4(a), 5(a), 6(a)) |
|**Robinson_2011_scenarios.R** | Runs very low leakage, low leakage, and high leakage scenarios, sampled with Robinson et al's method (Figures 4(b), 5(b), 6(b)) | 
| **fixed_600_scenario.R** | Runs very low leakage, low leakage, and high leakage scenarios, sampled using a fixed 600 sample volume (Figures 4(c), 5(c), 6(c)) | 
