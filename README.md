# ACT-Baseflow-Separation

A generic hydrological analysis workflow for separating baseflow from total streamflow using multiple separation methods. The workflow compares "Without Development" (WOD) scenario data against empirical flow measurements for any gauge site.

## Overview

This repository contains a Jupyter notebook workflow that performs baseflow separation analysis to distinguish groundwater contributions (baseflow) from total streamflow. The analysis compares two datasets:

1. **WOD (Without Development) Data**: Flow observations representing a "without development" scenario (pre-development or naturalized conditions)
2. **Empirical Flow Data**: Measured flow data representing actual observed conditions

The primary objective is to understand how baseflow characteristics differ between these two scenarios, providing insights into the hydrological behavior and the impact of development on baseflow patterns.

**Note**: This workflow is designed to work with any gauge site. Example data is provided for gauges 410734 and 410761, but the workflow can be adapted for any gauge by updating the data file paths.

## Purpose

The fundamental goal of this baseflow separation workflow is to improve the determination of baseflow environmental water requirements for the SWIOD review and future hydrological modeling in the ACT. By comparing "Without Development" (WOD) scenarios with empirical flow data, this analysis provides critical insights into:

- Baseflow characteristics under naturalized conditions vs. current observed conditions
- Seasonal patterns and variability in baseflow contributions
- The impact of development on baseflow regimes
- Data-driven support for environmental water requirement assessments

The outputs from this workflow inform hydrological modeling decisions and environmental water management strategies for the Australian Capital Territory.

## What is Baseflow Separation?

Baseflow separation is a hydrological technique used to partition total streamflow into two components:
- **Baseflow**: The portion of streamflow that comes from groundwater discharge and other delayed sources (sustained flow during dry periods)
- **Quickflow/Runoff**: The portion of streamflow that comes from direct surface runoff and interflow (rapid response to precipitation)

This separation is essential for:
- Understanding groundwater contributions to streamflow
- Water resource management and planning
- Environmental flow assessments
- Hydrological modeling and analysis

## Baseflow Separation Methods

The workflow uses the `baseflow` Python library, which implements multiple separation algorithms. While the analysis focuses primarily on the **UKIH method**, the following methods are available:

- **UKIH**: United Kingdom Institute of Hydrology method
- **Local Minimum**: Local minimum method
- **Fixed Interval**: Fixed interval method
- **Sliding Interval**: Sliding interval method
- **LH**: Lyne-Hollick method
- **Chapman**: Chapman method
- **CM**: CM method
- **Boughton**: Boughton method
- **Furey**: Furey method
- **Eckhardt**: Eckhardt filter method
- **EWMA**: Exponentially Weighted Moving Average
- **Willems**: Willems method

## Workflow Components

### 1. Data Loading
- Loads WOD (Without Development) scenario data from CSV files
- Loads empirical flow data from CSV files
- Processes datetime indices and handles missing values
- **Note**: File paths should be updated in the notebook to match your gauge site data

### 2. Baseflow Separation
- Applies `baseflow.single()` to both datasets
- Calculates baseflow using multiple methods (with focus on UKIH)
- Generates time series of separated baseflow components

### 3. Data Export
- Merges original flow data with baseflow results
- Exports results to CSV files:
  - `merged_baseflow_results_wod.csv`
  - `merged_baseflow_results_empirical.csv`

### 4. Comparative Analysis

#### Monthly Statistical Analysis
- **Monthly KDE (Kernel Density Estimation) plots**: Compare the distribution of baseflow values by month between WOD and empirical data
- **Monthly histograms**: Visualize baseflow frequency distributions
- **Statistical summaries**: Median, 10th percentile, and 90th percentile values for each month

#### Baseflow Proportion Analysis
- Calculates the proportion of baseflow to total flow (UKIH baseflow / total flow)
- Compares monthly median proportions between datasets
- Generates KDE plots showing the distribution of baseflow proportions by month

#### Flow Comparison
- Compares median total flows and baseflow by month
- Can include reference thresholds for context (e.g., 150 ML/day)
- Visualizes both total flow and baseflow components

### 5. Visualizations

The workflow generates several types of visualizations:

- **Interactive Plotly charts**: Time series of flow and baseflow
- **Matplotlib plots**: Monthly comparisons, statistical distributions, and KDE plots
- **Comparative plots**: Side-by-side analysis of WOD vs. empirical data

## Key Outputs

1. **Time Series Data**: Complete time series of baseflow separated from total flow
2. **Monthly Statistics**: Median, percentiles, and distribution characteristics by month
3. **Proportion Analysis**: Baseflow as a percentage of total flow
4. **Visual Comparisons**: Graphical representations of differences between WOD and empirical datasets

## Dependencies

The workflow requires the following Python packages:

- `pandas`: Data manipulation and analysis
- `matplotlib`: Static plotting
- `plotly`: Interactive visualizations
- `seaborn`: Statistical visualizations
- `numpy`: Numerical computations
- `scipy`: Statistical functions (KDE)
- `baseflow`: Baseflow separation algorithms

## Data Requirements

The workflow expects two data files in the `Data/` directory for each gauge site:

- **WOD data file**: Flow data representing "without development" scenario (e.g., `WOD_<gauge_id>.csv`)
- **Empirical data file**: Measured flow data representing observed conditions (e.g., `<gauge_id>_empirical_flow_<period>.csv`)

**Example files provided**:
- Gauge 410734: `WOD_410734.csv` and `410734_empirical_flow_<period>.csv`
- Gauge 410761: `WOD_410761.csv` and `410761_empirical_flow_1992_2022.csv`

Both files should have:
- Date/time as the index (first column)
- Flow values in ML/day or ML/month units
- Consistent date formatting

To use the workflow with a different gauge site, simply update the file paths in the notebook to point to your gauge's data files.

## Usage

1. Ensure all dependencies are installed
2. Place your gauge site data files in the `Data/` directory
3. Open `Baseflow_Separation.ipynb` and update the file paths to match your gauge site data
4. Run the notebook cells sequentially
5. Review outputs and exported CSV files

**For new gauge sites**: Update the CSV file paths in the data loading cells to reference your gauge's WOD and empirical flow data files.

## Repository Structure

```
ACT-Baseflow-Separation/
├── Baseflow_Separation.ipynb    # Main analysis notebook
├── README.md                     # This file
├── .gitignore                    # Git ignore rules
└── Data/                         # Data directory (ignored by git)
    ├── WOD_410734.csv            # Example: WOD data for gauge 410734
    ├── 410734_empirical_flow_*.csv
    ├── WOD_410761.csv            # Example: WOD data for gauge 410761
    └── 410761_empirical_flow_1992_2022.csv
```

## Notes

- The analysis focuses on the **UKIH method** for baseflow separation, though other methods are calculated
- Results are exported to CSV files for further analysis or reporting
- The workflow includes extensive visualization for exploratory data analysis
- Monthly analysis helps identify seasonal patterns in baseflow behavior
- **WOD** stands for "Without Development" scenario, representing pre-development or naturalized flow conditions
- This workflow is generic and can be applied to any gauge site by updating the data file paths

## License

This repository is part of the ACT Government hydrological modeling work.
