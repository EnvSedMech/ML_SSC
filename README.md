Part 3A: Model Optimization using Landsat-derived SSC Data (See "Description_Part3A.docx" for more details)

Overview
This repository contains the code and data for the project titled "Development and Validation of an Integrated Modeling Approach to Reconstruct the Propagation of Fluvial Sediment Pulses After Dam Removals." The objective is to optimize a one-dimensional hydro-morphodynamic numerical model using Landsat-derived suspended sediment concentration (SSC) data.

Model Input:
1) Discharge (Q) extracted from "Q_h_list.rds"
2) Flow Depth (h) extracted from "Q_h_list.rds"

Model Calibration (Landsat-based SSC):
1) Landsat-based SSC data extraction, dependent on incident lighting conditions, cloud cover, and water pixel extraction efficiency (thresholding using NDWI).
2) Reflectance values used to derive SSC for water pixels in the study area.
3) Files "SSC_list.rds" and "dgo_list.rds" contain spatially explicit SSC values matched with the correct Finite volume location for the Numerical model.

Numerical Model:
1) Well-commented version shared in "Demo_w_desc.R."
2) Model Formulation:
	a. Governing equations: 1D unsteady St. Venant shallow water equations and Exner equations for conservation of sediment mass.
	b. Ability to simulate non-uniform bed load and suspended sediment transport rates.
	c. Ability to simulate changes in bed gradation (Castro Bolinaga, 2016).

Channel geometry:
1) Extracted using publicly available NAIP imagery.
2) Simplified using changepoint analysis.
3) Calibrated for Water Surface Elevation (m) and SSC (m³/m³).
4) Validated using remote sensing-derived Erosion and Deposition amounts (m).

Inverse Modeling Framework:
1) Relies on the conceptual framework validated for the Elwha Dam removals by Sharma, (2023) (also shown in the "Description_Part3A.docx" file)
2) Well-commented code available in "Demo_w_desc.R."
3) Parameter to be optimized: Upstream sediment supply concentration (Cs in m³/m³).
4) Objective function: RMSE.
5) MCMC used with adaptive scaling for efficient exploration of Cs parameter space.
6) Results saved at each time step.
7) Multiple scripts are run to spatiotemporally reconstruct sediment pulses using remote sensing, numerical modeling, and machine learning. "Parallel.R" demonstrates a simple method for running multiple scripts in parallel.

License
This project is licensed under the Creative Commons Zero V1.0 Universal License.

References
Castro Bolinaga, C. F. (2016). Numerical Simulation of the Propagation of Fine-Grained Sediment Pulses in Alluvial Rivers. Virginia Tech.
Sharma, A. (2023). Development and Validation of an Integrated Modeling Approach to Reconstruct the Propagation of Fluvial Sediment Pulses After Dam Removals.
