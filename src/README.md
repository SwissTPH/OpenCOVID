# OpenCOVID

An individual-based model of SARS-CoV-2 transmission and COVID-19 disease dynamics developed at the Swiss Tropical and Public Health institute (Swiss TPH). A stochastic, discrete-time, setting-agnostic model with age, risk-group, contact network, waning immunity, and viral variant structure. Numerous interventions can be simulated, including vaccination (plus booster doses), pre-exposure prophylaxis, testing and diagnosis, treatment, isolation, and non-pharmaceutical interventions (such as physical distancing policies of varying intensities).

#### Version 4.0 (beta) (stable with R version 4.1)

## The repository

If you will be regularly contributing to model development, `git clone` this repository to your desired location. If you wish to edit, commit, and push any model code, please branch off from master first. If you would just like to experiment then please `git fork` instead.

## Analysis files

For applications, most of the interaction you will need to have with the model will be via an 'analysis file'. Analysis files are written in the YAML language. An analysis file details model parameters and metaparameters to be simulated, along with any alternative scenarios. The general principle is that one analysis file is sufficient to define all scenarios needed to answer a research question. 

For the sake of consistency and readability, a default file sets a value for each required model parameter, with the purpose of an analysis file to be to overwrite any such defaults. This overwriting can occur for all scenarios in the analysis or only certain scenarios.

Two analysis files exist in the repository as standard:
1. `/config/default.yaml`: A full analysis file with all default model parameters. The values in this file should not be altered. 
2. `/input/demo.yaml`: A demo analysis file showing some basic functionality. This file shows an example of an 'array scenario' (discussed below).

To create a new application, create a new YAML file and save this file in the `/input/` directory with a (machine readable) file name of your choice (ensuring the file has a `.yaml` extension). The contents of this directory are git-ignored (aside from `demo.yaml`), so it is perfectly acceptable to create and store any number of analysis files. To instruct the model to use a particular analysis file, set the value of `o$analysis_name` in `options.R` to the name of the analysis file (without the `.yaml` extension).

It is good practice to familiarise yourself with the `/config/default.yaml` file to 1) have an understanding of which model parameters, metaparameters, and options are used by default, and 2) check syntax and formatting. If certain default values are sufficient, it is unnecessary (but not an error) to define these parameters in your analysis file(s).

Note that it is generally not necessary to preserve ordering in analysis files.

## The pipeline

### Getting started
- All model analyses are launched from `launch.R`
  - Preferred usage is to 'source' this file (without 'echo' is ideal)
- Alternative command line usage is to `cd` to the code directory then call `sh bash_launch.sh`
- Any required packages not currently installed are automatically downloaded and installed the very first time you run the model
  - If several large packages need to be installed this process may take a while (potentially 30-60 minutes)
  - See `dependencies.R` for a complete list of packages used by the model
  
### Directory structure
- This is all taken care of when `launch.R` (or `bash_launch.sh`) is called
- The model will automatically re-point the current working directory
- The model will also automatically create the necessary output directory structure
  
### Step 0: Running a test simulation
- Set `do_step = 0` in `launch.R` and 'source' to run a test simulation with baseline model parameters
- Once finished (this may take several minutes depending on settings), see `~\output\3_figures\demo\` for default output plots

### Step 1: Generate a model calibration
- Set `do_step = 1` in `launch.R` and 'source' to calibrate the model
- There are three broad options for model calibration available (by setting `calibration_type` in your analysis file):
  1. Calibrate the model to a user-defined effective reproduction number at the start of the simulation period (`calibration_type: "r_user"`, the default)
  2. Calibrate the model to an effective reproduction number that is calculated from setting-specific epidemiological data (`calibration_type: "r_data"`)
  3. Calibrate the model to a set of setting-specific epidemiological metrics (such as confirmed cases, hospitalisations, and deaths) over time (`calibration_type: "epi_data"`)

### Step 2: Run scenarios
- Set `do_step = 2` in `launch.R` and 'source' to simulate the baseline and any user-defined alternative scenarios detailed in the analysis file
  
### Step 3: Plot all your results
- Set `do_step = 3` in `launch.R` and 'source' to generate a series of plots that examine the analyses you have run
- See `results.R` to see which results are to be generated
- All plotting functionality itself resides in `plotting.R`
- Plotting toggles can be turned on and off, see `plot_xxx` in `options.R`

## Calibration methodology

For calibration type 1 (`"r_user"`) and type 2 (`"r_data"`), it is generally sufficient to calibrate only one model parameter: the average number of daily contacts (the default setting). For calibration type 3 (`"epi_data"`), it may be necessary to calibrate multiple model parameters in order to achieve a good fit between model output and the epidemiological data.

The general calibration process is the same for all three approaches: 
<ol type="i">
<li>A number of parameter sets are sampled from the input parameter space using Latin Hypercube sampling</li>
<li>These parameters sets are simulated for a sufficient period</li>
<li>The quality of fit between the model output and the target value(s) is calculated for each parameter set</li>
<li>A surrogate model emulator is trained to learn the relationship between parameter values and quality of fit</li>
<li>Adaptive sampling is performed using the trained emulator and an expected improvement acquisition function</li>
<li>Steps ii-iv are then repeated for the newly sampled parameter sets attained from step v</li>
<li>Finally, the optimal parameter set that achieves the best quality of fit is determined via a stochastic decent algorithm that operates on the (re-)trained model emulator</li>
</ol>

The best-fitting parameter set determined from this process is then used when simulating scenarios (pipeline step 2).

Note that adaptive sampling can be performed multiple times if desired (set `adaptive_sampling::rounds` in your analysis file). In such a case, steps ii-v are repeated.

## Parallel computing

Due to the potentially large number of individual-based model simulations required in 'steps' 1 and 2, it is necessary for OpenCOVID to be run using parallel computing. Ideally this would be done on a cluster computing framework. 

Out of the box, OpenCOVID uses the UniBas (University of Basel) cluster to run simulations in parallel. To run the full model pipeline outside of the University of Basel network, you will first need to configure parallel simulations, ideally on a cluster. 

Note that step 0 simulates the model locally (as this is just one single simulation), and therefore does not require parallel computing.

## Final notes
- It is totally acceptable to run the full pipeline from end to end by setting `do_step = 1 : 3` in `launch.R`
- These notes are just to get you started, of course
- For a full introduction to the model and it's capabilities contact the development team

## Release details

#### V4.0 release date: 07-10-2022

#### Development and maintenance:
* Andrew J. Shattock (andrewjames.shattock@swisstph.ch)

#### Contributors:
* Epke A. Le Rutte
* Max Richter
* Cassandra Alvarado
* Sherrie L. Kelly
* Nakul Chitnis
* Melissa A. Penny

#### Previous contributors:
* Robert P. Duenner
* Swapnoleena Sen

