### Table of contents



##### **extraction\_from\_big\_datasets:**

Implements the big data extraction simulations described in Section 5.1.

###### **sec51\_extraction\_sim.R**:

* The main script to conduct the simulation and generate the boxplots for nRMSPE and computation time. To handle the wide variety of simulation settings conveniently, the script uses a set of indices (index1 through index4) to assign the response function type, X distribution, kernel type, and signal-to-noise ratio (SNR). For example, an index set of (1, 1, 1, 1) corresponds to the Brat function, Uniform X, Matern 3/2 kernel, and SNR = 2.

###### **all\_functions.R**:

* Contains all helper functions required for the simulation, including the response functions, baseline methods (e.g., Chang (2023), ASMECr), and evaluation metrics.



##### **elimination\_of\_redundant\_data\_points:**

Implements the redundant data removal simulations described in Section 5.2.

###### **sec52\_removal\_sim.R:**

* Runs the backward elimination process and generates Figures 3a through 4c.

###### **data\_test.csv, data\_train.csv, validation\_index.csv:**

* These three CSV files contain the generated data used in the paper. They are provided to ensure the exact reproducibility of the plots in Section 5.2. The original code used to generate this data is also included as comments within the script.



##### **sacros\_dataset:**

Implements the real-world data extraction application described in Section 6.1.

###### **sec61\_sarcos.R**:

* Conducts the data extraction on the SARCOS dataset and generates the boxplots for nRMSPE and computation time.

###### **sarcos\_inv.csv, sarcos\_inv\_test.csv:**

* The training and testing datasets for the SARCOS kinematics application.



##### **real\_estate\_valuation\_dataset:**

Implements the real-world data removal application described in Section 6.2.

###### **sec62\_real\_estate.R**:

* Runs the redundant data removal process and generates Figures 6a through 6d.

###### **dataset.csv**:

* The Real Estate Valuation dataset.

###### **index\_for\_reproduce.csv**:

* Contains the exact indices used to partition the data into training, validation, and testing sets to guarantee reproducibility.

