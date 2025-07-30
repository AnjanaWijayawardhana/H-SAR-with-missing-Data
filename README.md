# H-SAR-with-Missing-Data

This repository contains R code for the simulation studies and real-world examples presented in the manuscript.

---

##  Folder Structure

Each simulation or real-world example is located in a separate folder, containing the following key files:

- **`Source.R`** – Contains the core R code for running the algorithms.  
- **`Implement.R`** – Implements and runs all algorithms for the selected scenario.

---

## Running the Scripts

To reproduce results from the simulation section of the manuscript (e.g., *Simulation: Fit SEM under MAR with n = 5,041 and 90% missing data*):

1. **Download** the relevant folder `Simulations`.
2. **Set the working directory** to that folder in R or RStudio.
3. **Open** the `Implement.R` script.
4. **Specify** the model type, number of observations (`n`), number of simulations (`N`), and the missing data percentage.
5. **Run** the script to execute the simulation.

>  For other simulations, simply repeat the steps above using the appropriate folder and parameters.

To reproduce results from the Real-world example section of the manuscript (e.g., *Fit SEM for Lucas County house price under MAR with 90% missing data*):
1. **Download** the relevant folder `Simulations`.
2. **Set the working directory** to that folder in R or RStudio.
3. **Open** the `Implement.R` script.
4. **Specify** the model type, number of observations (`n`), number of simulations (`N`), and the missing data percentage.
5. **Run** the script to execute the simulation.

---

## Environment Requirements

###  R Version

- Tested on **R 4.4.2**  
- Also compatible with **xxxxx**

### Required Packages

Install the necessary packages by running:

```r
install.packages(c("MASS", "mvtnorm", "truncnorm", "BayesLogit", "Matrix", "tidyverse", "ggplot2", "patchwork"))
