---
output:
  html_document: default
  pdf_document: default
---
## 📚 Table of Contents

------------------------------------------------------------------------

### [1. `extraction_from_big_datasets`]{style="color:#1F77B4;"}

> Implements the big data extraction simulations described in **Section 5.1**.

#### 📄 `sec51_extraction_sim.R`

**Main simulation script**

- Conducts the simulation study.
- Generates the boxplots for:
  - **nRMSPE**
  - **Computation time**

#### ▶️ How to run

1. **Install the required packages**

   Before running the simulation, please install all required R packages.
   The package list is given at the top of `sec51_extraction_sim.R`.
   In particular, the script uses packages such as:

   - `pacman`
   - `Metrics`
   - `MultiRNG`
   - `hetGP`
   - `mvtnorm`
   - `SPlit`
   - `class`
   - `mixAK`
   - `stringr`
   - `mined`
   - `twinning`
   - `supercompress`
   - `svMisc`
   - `tidyverse`
   - `gss`
   - `FNN`
   - `tmvtnorm`
   - `np`
   - `twingp`
   - `GpGp`
   - `fields`
   - `doParallel`
   - `foreach`

   The script also requires the `GEMSS` package. If it is not available locally, install it from GitHub first.

2. **Load the required packages and source helper functions**

   After installing the packages, run the package-loading commands at the beginning of
   `sec51_extraction_sim.R`, and then source

   - `all_functions.R`

   which contains all helper functions used in the simulation.

3. **Set the simulation scenario**

   The script uses four indices to control the simulation setting:

   - `index1`: true response function
   - `index2`: distribution of **X**
   - `index3`: covariance (kernel) function
   - `index4`: signal-to-noise ratio (**SNR**)

   Specifically,

   - `fun_type = c("brat", "robot", "g_fun", "welch")[index1]`
   - `X_type = c("Uniform", "Mvt", "MVNormal")[index2]`
   - `Cov_Fun = c("Matern3_2", "Matern5_2")[index3]`
   - `SNR = c(2, 10)[index4]`

   **Example:**  
   `(index1, index2, index3, index4) = (1, 1, 1, 1)` corresponds to

   - Brat function
   - Uniform \(X\)
   - Matérn 3/2 kernel
   - SNR = 2
   
   **Reproducing Figure 2 in the paper.**  
To reproduce Figure 2, run `sec51_extraction_sim.R` under the following four settings:
`(1,1,1,1)`, `(1,1,2,1)`, `(1,1,1,2)`, and `(1,1,2,2)`, shown from left to right in the figure.

4. **Set the data sizes and replication number**

   The simulation script allows the user to specify:

   - `n_train`: full training data size
   - `n_test`: testing data size
   - `ns_set`: subdata sizes
   - `dup`: number of replicates

   For example, the setting in the script uses

   - `n_train = 10000`
   - `n_test = 10000`
   - `ns_set = c(100, 200, 300, 400)`
   - `dup = 10`

   while the paper setting (commented out in the script) uses larger values.

5. **Run the simulation**

   After setting the indices and sample sizes, run `sec51_extraction_sim.R`.
   The script will generate the simulation results and produce the corresponding boxplots
   for **nRMSPE** and **computation time**.

#### 🛠️ `all_functions.R`

**Supporting functions**

- Contains all helper functions required for the simulation.
- Includes:
  - response functions
  - baseline methods
    - Chang (2023)
    - ASMECr
  - evaluation metrics

------------------------------------------------------------------------

### [2. `elimination_of_redundant_data_points`]{style="color:#D62728;"}

> Implements the redundant data removal simulations described in **Section 5.2**.

#### 📄 `sec52_removal_sim.R`

**Backward elimination script**

- Runs the backward elimination procedure.
- Generates:
  - **Figure 3a** through **Figure 4c**

#### ▶️ How to run

1. **Install the required packages**

   Before running the script, install the required R packages.  
   The script uses the following packages:

   - `tidyverse`
   - `hetGP`
   - `twinning`
   - `Metrics`
   - `GEMSS`
   - `ContourFunctions`

   If the `GEMSS` package is not available locally, install it first (e.g., from GitHub).

2. **Load the required packages**

   At the beginning of `sec52_removal_sim.R`, load the required packages:

   - `library(tidyverse)`
   - `library(hetGP)`
   - `library(twinning)`
   - `library(Metrics)`
   - `library(GEMSS)`
   - `library(ContourFunctions)`

3. **Set the working directory**

   Set the working directory to

   - `elimination_of_redundant_data_points`

   so that the script can correctly read the provided CSV files.

4. **Prepare the data**

   For exact reproducibility of the paper figures, the script reads the following files directly:

   - `data_train.csv`
   - `data_test.csv`
   - `validation_index.csv`

   These files contain the training data, testing data, and validation split used in the paper.

   The script also includes commented code showing how the original Dropwave dataset was generated and how the validation split was constructed.  
   For reproducing the paper figures exactly, use the provided CSV files rather than regenerating the data.

5. **Set the simulation configuration**

   The main simulation setting in the script is controlled by:

   - `Cov_Fun <- "Matern5_2"`: covariance/kernel function
   - `ngrid <- 100`: grid resolution for contour plots

   The plotting grid is then constructed by

   - `grid.x <- seq(0, 1, len = ngrid)`
   - `grid <- expand.grid(x1 = grid.x, x2 = grid.x) %>% as.matrix()`

6. **Generate Figures 3a through 4c**

   Running `sec52_removal_sim.R` produces the figures in the following order:

   - **Figure 3a**: contour plot of the Dropwave function
   - **Figure 3b**: GP prediction performance using training + validation data
   - **Figure 4a**: GP prediction performance using training data only
   - **Removing process**: run
     - `result <- gemss_remove(...)`
   - **Figure 4c**: predictive \(R^2\) plot during the removal process
   - **Figure 4b**: GP prediction after eliminating the selected redundant points

7. **Reproducing the paper figures**

   To reproduce the figures in Section 5.2 exactly, run `sec52_removal_sim.R` with the provided data files and the default settings in the script, in particular:

   - `Cov_Fun <- "Matern5_2"`
   - `ngrid <- 100`

#### 📂 Data files

- `data_test.csv`
- `data_train.csv`
- `validation_index.csv`

**Purpose:**

- These three CSV files contain the generated data used in the paper.
- They are provided to ensure the **exact reproducibility** of the plots in Section 5.2.
- The original code used to generate these data is also included as **comments within the script**.

------------------------------------------------------------------------

### [3. `sacros_dataset`]{style="color:#2CA02C;"}

> Implements the real-world data extraction application described in **Section 6.1**.

#### 📄 `sec61_sarcos.R`

**SARCOS application script**

- Conducts data extraction on the **SARCOS dataset**.
- Generates the boxplots for:
  - **nRMSPE** (see Figure 5)
  - **Computation time**

#### ▶️ How to run

1. **Install the required packages**

   Before running the script, install the required R packages.

   The script uses packages such as:

   - `Metrics`
   - `MultiRNG`
   - `hetGP`
   - `mvtnorm`
   - `SPlit`
   - `class`
   - `mixAK`
   - `stringr`
   - `mined`
   - `twinning`
   - `supercompress`
   - `svMisc`
   - `tidyverse`
   - `gss`
   - `FNN`
   - `tmvtnorm`
   - `np`
   - `twingp`
   - `GpGp`
   - `fields`
   - `doParallel`
   - `foreach`
   - `GEMSS`

   If the `GEMSS` package is not available locally, install it first (e.g., from GitHub).

2. **Load the required packages**

   The script uses `pacman::p_load(...)` to load all required packages automatically.

3. **Set the working directory**

   Set the working directory to

   - `sacros_dataset`

   so that the script can correctly read the SARCOS training and testing files.

4. **Prepare the data**

   The script reads

   - `sarcos_inv.csv`
   - `sarcos_inv_test.csv`

   and then defines

   - `X`, `Y`: training predictors and responses
   - `X_test`, `Y_test`: testing predictors and responses

5. **Set the experimental configuration**

   The main settings in the script are:

   - `ns_set`: subdata sizes  
     - example in the script: `c(100, 200, 300, 400)`
     - larger paper setting (commented out): `c(250, 500, 750, 1000)`
   - `dup`: number of replicates  
     - example in the script: `dup = 5`
     - paper setting (commented out): `dup = 50`
   - `index1`: covariance-function index
     - `index1 = 1` corresponds to `Cov_Fun = "Matern3_2"`
     - `index1 = 2` corresponds to `Cov_Fun = "Matern5_2"`
   - `cl`: number of CPU cores for parallel computing

6. **Run the simulation**

   After setting `index1`, `ns_set`, `dup`, and the number of cores, run `sec61_sarcos.R`.

   In each replicate, the script:

   - uses Twinning to select a manageable subset for GP parameter estimation,
   - fixes the estimated GP parameters,
   - compares the following methods:
     - `GEMSS`
     - `Chang (2023)`
     - `ASMECr`
     - `Supercompress`
     - `SRS`
     - `twinGP`
     - `GpGp`

7. **Generate the plots**

   At the end of the script, two boxplots are produced:

   - **nRMSPE boxplot**
   - **log10(Time in minutes) boxplot**

#### 🔁 Reproducing Figure 5 in the paper

To reproduce the main SARCOS prediction figure in the paper, run the script twice:

- once with `index1 = 1` (`Matern3_2`)
- once with `index1 = 2` (`Matern5_2`)

These two settings correspond to the two kernel choices reported for the SARCOS application.

#### 📂 Data files

- `sarcos_inv.csv`
- `sarcos_inv_test.csv`

**Purpose:**

- These are the training and testing datasets for the **SARCOS kinematics application**.
- They are used directly by `sec61_sarcos.R` to reproduce the real-data extraction results.

------------------------------------------------------------------------

### [4. `real_estate_valuation_dataset`]{style="color:#9467BD;"}

> Implements the real-world data removal application described in **Section 6.2**.

#### 📄 `sec62_real_estate.R`

**Real estate valuation script**

- Runs the redundant data removal process.
- Generates:
  - **Figure 6a** through **Figure 6d**

#### ▶️ How to run

1. **Install the required packages**

   Before running the script, install the required R packages.

   The script uses:

   - `tidyverse`
   - `hetGP`
   - `twinning`
   - `Metrics`
   - `GEMSS`
   - `ContourFunctions`

   If the `GEMSS` package is not available locally, install it first (e.g., from GitHub).

2. **Load the required packages**

   At the beginning of `sec62_real_estate.R`, load the required packages using:

   - `library(tidyverse)`
   - `library(hetGP)`
   - `library(twinning)`
   - `library(Metrics)`
   - `library(GEMSS)`
   - `library(ContourFunctions)`

3. **Set the working directory**

   Set the working directory to

   - `real_estate_valuation_dataset`

   so that the script can correctly read the dataset and reproducibility index file.

4. **Prepare the data**

   The script reads

   - `dataset.csv`
   - `index_for_reproduce.csv`

   Then it:

   - extracts the relevant columns,
   - log-transforms the response (`unit.price`),
   - reconstructs the exact training / validation / testing split used in the paper.

   The script also includes commented code showing the original data-splitting procedure.
   For exact reproduction of the paper figures, use the provided `index_for_reproduce.csv`.

5. **Set the analysis configuration**

   The main settings in the script are:

   - `Cov_Fun <- "Gaussian"`: covariance/kernel function
   - `ngrid <- 80`: grid resolution used for contour plots

   The plotting grid is then constructed from the latitude and longitude ranges of the data.

6. **Generate Figures 6a through 6d**

   Running `sec62_real_estate.R` produces the figures corresponding to the paper in the following order:

   - **Figure 6a**: GP prediction performance using training + validation data
   - **Figure 6b**: GP prediction performance using training data only
   - **Figure 6c**: GP prediction after removing the selected redundant points
   - **Figure 6d**: predictive \(R^2\) plot during the removal process

7. **Run the removal procedure**

   The main removal step is

   - `result <- gemss_remove(...)`

   where the training set is reduced using the validation set to assess predictive performance during elimination.

#### 🔁 Reproducing the paper figures

To reproduce the figures in Section 6.2 exactly, run `sec62_real_estate.R` with the default settings in the script, in particular:

- `Cov_Fun <- "Gaussian"`
- `ngrid <- 80`

together with the provided reproducibility file:

- `index_for_reproduce.csv`

#### 📂 Data files

- `dataset.csv`

**Purpose:**

- Contains the **Real Estate Valuation dataset**.

#### 📄 `index_for_reproduce.csv`

**Reproducibility file**

- Contains the exact indices used to partition the data into:
  - training set
  - validation set
  - testing set
- Included to guarantee **full reproducibility** of the paper figures.

------------------------------------------------------------------------
