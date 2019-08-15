# Quantis
Helps to visualise the effluent data of experimental setups.
Creates input files for the CHEMKINParallellizer and COILSIM1DParallellizer for the performed experiments.
Enables comparison between simulation and experiment.

## Requirements
Install python 3+ via anaconda: https://www.continuum.io/downloads
Install the Quantis conda environment from the "environment.yaml" file located in the main folder of the Quantis:
```
$ conda env create -f environment.yaml
```

Add python to your user path (if you didn't check to add it via Anaconda).
- **Windows 7**  
 Right-click on your C:\ Drive properties -> Advanced system settings -> Add your administrator password -> Environment Variables -> Add python to the 'Path' system variables Example
``` 
 C:\Users\\*username*\AppData\Local\Continuum\Anaconda3
 C:\Users\\*username*\AppData\Local\Continuum\Anaconda3\Scripts
``` 
 

### Necessary input
The input for the `GC_Processor_HDF5` of **Quantis** are the integrated RGA, GCxGC, (LOA) data (resulting excel files).
These are stored in the appropriate folders named RGA, GCxGC, and (optionally) LOA.
The program gets his input from two '.csv' files: conditions and gen_configs (Examples of both can be found in the working example folder.)
If the older '*name*.xml' file is given as input, the two '.csv' files are partially created.

The 'conditions_*name*.csv' has the following outline:
Each row depicts a single datapoint of the experiment and each column represents the following information about the datapoint:
 - Name
 - Time
 - FeedRate (g/h) of the hydrocarbon in: 'Feed_*Hydrocarbon*_FeedRate'
 - FeedRate (g/h) of H2O in: 'Feed_H2O_FeedRate'
 - FileName of the corresponding GCxGC data 
 - FeedRate (g/h) of the Primary Internal Standard in: 'Feed_PrimaryInternalStandard_FeedRate'
 - Name of the Primary Internal Standard
 - FileName of the corresponding RGA-FID (detector 1)
 - FileName of the corresponding RGA-TCD-L (detector 2)
 - FileName of the corresponding RGA-TCD-R (detector 3)
 - Info of the Secondary Internal Standard
 - Name of the Secondary Internal Standard
 - Temperature
 - Info of the Tertiary Internal Standard
 - Name of the Tertiary Internal Standard
 - The GCxGC correction factor
 
 The 'gen_configs_*name*.csv' has the following outline:
 - Name of the experiment
 - The FileName of the experiment: *YearMonthDay*-*Name of the experiment*
 - FileLocation of the RGA
 - FileLocation of the GCxGC FID
 - FileLocation of the GCxGC SCD
 - FileLocation of the GCxGC NCD
 - FileLocation of the logging files
 - FileLocation of the Calibration factors
 - FileName for RGA FID Calibration factors
 - FileName for RGA TCD Calibration factors
 - FileName for GCxGC Calibration factors
 - Name of the reactor
 - Usage of SCD: TRUE or FALSE
 - Usage of LOA: TRUE or FALSE
 - Usage of NCD: TRUE or FALSE
 - Name of the hydrocarbon feed
 - Name of the diluent
 - Start time of the experiment: *DD/MM/YYYY HH:MM*
 - End time of the experiment: *DD/MM/YYYY HH:MM*
 - Name of the Researcher who performed the experiments
 - Name of the Researcher who performed the processing
 - Confidential data: yes or no
 - Published data: yes or no
 - Reliability factor: 1-5
 - Remarks
 - Absolute pressure (bara)
 - Availability of a mass spectrum: yes or no
 - Purpose/campaign of the experiments
 

**deprecated:** The '*name*.xml' file has the following outline.
Every injection is denoted as a *datapoint*.
And every *datapoint* has a
- Name
- Temperature
- RGA-TCD-L filelocation
- RGA-TCD-R filelocation
- RGA-FID filelocation
- GCxGC filelocation
- Date
- Time
- Feed (can be multiple)
  - Name
  - FeedRate in g/h
- PrimaryInternalStandard
  - Name
  - FeedRate in g/h
- SecondaryInternalStandard
  - Name
  - Info
- TertiaryInternalStandard
  - Name
  - Info
- QuaterneryInternalStandard (optional)
  - Name
  - Info
The '*.xml*' file further has a switch to enable LOA data (*process_LOA*),
and the locations of the response libraries for GCxGC (*RF-GCxGC*), RGA-TCD (*RF-RGA-TCD*), and RGA-FID (*RF-RGA-FID*). 
An example of this '.xml' file can be found in the GC_processor GitHub.

### Usage
**Quantis** works by double-clicking the run.bat file and small GUI appears that asks for what part of the program should run, a command prompt window will appear showing the progress and/or warnings. An HDF5 file container is created that holds all data: raw, processed, processing conditions and (if ran) simulation results. Attached to the file container are the conditions of the experiment as metadata which enables easier filtering of the database.

4 modes are available for the program:
- Process the results to create the HDF5 file
- Visualize the experimental results
- Simulate the experimental conditions in either CHEMKIN or COILSIM1D with the CHEMKINParallellizer or COILSIMParallellizer respectively.
- Compare the simulations with the experiments
The first input field is to insert the already created HDF5 file or the input.csv and gen_configs.csv. The latter input is to include the CHEMKIN reaction network names. Mind that the reaction networks should be stored in the folder 'CHEMKIN networks' and an InChI converter should be stored in the folder 'InChI_converters'. Note that performing CHEMKIN and/or COILSIM1D simulations requires additional programs (not public).

#### Create HDF5 file
The library response files are necessary to create the HDF5 file. These response files are stored centrally under 'N:\Project\918-Experimental set-ups\GC_library' or under a local folder named 'GC_library'. The GCxGC library is a shared Excel document where every researcher can add their missing components. Remind that this document is used by other researchers as well, so make sure that the data that is added is correct!

The processed GC data is stored on the N: drive in the mirrorData and is read from there or from a local folder called 'GCxGC'. The folders of all the GC's should have the same folder name. (e.g. the folder in which the data is stored for the RGA should have the same name as the folder that includes the GCxGC-FID data)

By applying the date and time of injection, the program searches for the logging files on the N: drive in the setup folder (mirrored from the lab pc every night) or on a local folder 'logging'.

#### Visualize results
Several graphs are being generated to conclude what injections are trustworthy.
- Mass balance difference with 100%
- Relative difference between TCD and FID
- Elemental balance differences
- H2 versus C2H4 (wt.%)
- H2 versus C2H4+C3H6 (wt.%)
- C2H4, C3H6 and C4- yield
- Principal component analysis (PCA) yields
- PCA temperature profiles
- Centered temperature profiles
- Correlation matrix
- Scatter matrix between major components
- Major products vs temperature (wt.%)
- Minor products vs temperature (wt.%)

#### Simulate experiment
The CHEMKINParallellizer (not public) is being run. The folder 'N:\Project\918-Experimental set-ups\ChemkinParallellizer' is copied to the users desktop and from there the parallellizer is run.

The COILSIM1DParallellizer (not public) is also being run.
#### Compare simulations and experiment
Creates a comparison file in the folder 'Comparison' and comparison figures (conversion vs T, compound vs conversion, parity plot per compound) are created. Simulation summaries need to be stored in the 'simulationSummaries' folder.

### The HDF5 file
See for more information on HDF5, the website: https://www.hdfgroup.org/.

The processed yields file includes several sheets
- Results (raw results from RGA and GCxGC, sums FID > TCD > GCxGC) 'raw_yields'
- Normalised results wet/dry (mathematically normalized to 100) 'norm_yields(_wet)'
- C/H/O/N/S balance (based on raw results) '*C*_balance'
- mass bal from C norm (mass balance based on the closed C balance) 'C_norm_yields'
- norm C/H/O/N/S balance (based on normalised results) '*C*_balance_norm'
- Mol% dry/wet (molar yields based on normalized results) 'mol_perc(_wet)'
