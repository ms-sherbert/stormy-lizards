# stormy-lizards

Code and data repository for the manuscript entitled:

"Deacon MR, Herbert SM, Nelson NJ. Storm surge impacts on coastal lizard populations: a case study on how climate change may affect endemic species." Submitted to the New Zealand Journal of Ecology in July 2024. Please note that this manuscript and the pre-release version of this repository (0.9.1-alpha) has not yet been through peer review. 

All code is written for R version 4.3.1. We cannot guarantee that these codes will work with other versions of R. 

This code can be cited with the following doi: https://doi.org/10.5281/zenodo.12770298

## Repository contents

### Code

Code for all analyses, figures, and tabulated summary data reported in the manuscript or its Supplementary Material. This code pulls the raw data on lizard captures and survey covariates (in the `Data` subfolder), formats it correctly for further analysis, and writes tabulated data and analysis outputs to the `Outputs` subfolder. We have provided the tabulated data in `Outputs` so that these codes can be run independently. However, the code files were initially run in the following order:

1. `CPUE_Mt_calculations.R` Reads in raw data from the `Data` subfolder and calculates Mt+1 and CPUE for further analysis.
2. `Capture-time-inundation-analysis.R` Performs linear model analysis of the tabulated Mt+1 and CPUE data. Calls in helper functions from `Model-selection-functions.R`. Makes Figures 2 and 3. 
3. `Body-size-time-inundation-analysis.R` Performs linear model analysis of body sizes in the raw captures data. Calls in helper functions from `Model-selection-functions.R`. Makes Figure 4. 

### Data

`CC-grid-captures_S1-S8.csv` Lizard data collected in November 2017 - March 2021 from six monitoring grids. This file contains data from marked and observed lizards. 

`CMR_check_schedule.csv` Dates that each of the six lizard monitoring grids included in this study were checked (1 = checked, 0 = not checked), and the corresponding daily minimum and maximum temperatures recorded at Wellington airport weather station in degrees Celsius (temperature data sourced from https://cliflo.niwa.co.nz/). Breaks in consecutive dates that are otherwise not explained in the notes column were due to time spent monitoring lizards in an additional six (seasons 1-5) to ten (seasons 1-4) grids at Miramar and Baring Head for a wider lizard monitoring programme (Herbert 2020).   

`N-sessions.csv` are the number of capture sessions (aka secondary periods) for seasons 1 - 8 (2017-2021) used to calculate Mt+1 and CPUE. Summarizes the raw data in `CMR_check_schedule.csv`. 

`Site_inundation_and_trap_replacement.csv` A list of all lizard monitoring stations in the six grids, indicating which stations had the pitfall and/or Onduline ACOs replaced prior to the November 2020 survey. 1 = ACO or pitfall damaged or missing due to the storm surge and replaced, 0 = ACO or pitfall undamaged and not replaced.  

### Outputs

Contains tabulated data summaries exported from R codes. 
