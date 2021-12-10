# OpenCOVID
An individual-based model of SARS-CoV-2 transmission and COVID-19 disease dynamics developed at Swiss TPH. A stochastic, discrete-time model with age, risk-group, network, and viral variant structure. From version 2.0 onwards, the model is setting-agnostic.

##### Version 2.0 (stable with R versions 3.6-4.1)

## The repository
If you will be regularly contributing to model development, `git clone` this repository to your desired location. If you wish to edit, commit, and push any model code, please branch off from master first. If you would just like to experiment then please `git fork` instead.

## Analysis files
For applications, most of the interaction you will need to have with the model will be via an 'analysis file'. Analysis files are written in the YAML language. An analysis file details model parameters and metaparameters to be simulated, along with any alternative scenarios. The general principle is that one analysis file is sufficient to define all scenarios needed to answer a research question. 

For the sake of consistency and readability, a default file sets a value for each required model parameter, with the purpose of an analysis file to be to overwrite any such defaults. This overwriting can occur for all scenarios in the analysis or only certain scenarios.

Two analysis files exist in the repository as standard:
1. `/config/default.yaml`: A full analysis file with all default model parameters. The values in this file should not be altered. 
2. `/input/demo.yaml`: A demo analysis file showing some basic functionality. This file shows an example of an 'array scenario' (discussed below).

To create a new application, create a new YAML file and save this file in the `/input/` directory with a (machine readable) file name of your choice (ensuring the file has a `.yaml` extension). The contents of this directory are git-ignored (aside from `demo.yaml`), so it is perfectly acceptable to create and store any number of analysis files. To instruct the model to use a particular analysis file, set the value of `o$analysis_name` in `options.R` to the name of the analysis file (without the `.yaml` extension).

It is good practice to familiarise yourself with the `/config/default.yaml` file to 1) have an understanding of which model parameters, metaparameters, and options are used by default, and 2) check synatx and formatting. If certain default values are sufficient, it is unecessary (but not an error) to define these parameters in your analysis file(s).

Note that it is generally not necessary to preserve ordering in analysis files.

### Time parameters

```
n_days: <integer>
```

Parmeter `n_days` defines the number of days for which to forward project model outcomes.

### Demographic properties

```
population_size: 100000
risk_groups:
  comorbidities:
    age_lower: <integer>
    age_upper: <integer>
    probability: <float between 0 and 1>
  healthcare_worker:
    age_lower: <integer>
    age_upper: <integer>
    probability: <float between 0 and 1>
```
We recommend using a `population_size` of 10,000 to 100,000 for applications. Note that smaller populations enable the master to run faster, while larger populations produce higher quality results.

### Network properties

```
population_size: 100000
network_structure: <"random", "age", "layers>
network_layers:
- <"household">
- <"school">
- <"workplace">
contact_matrix_countries: <country_name(s)>
```

### Model metrics

```
model_metrics:
- <metric_name>: <yes/no>
  by: <none/age/variant/vaccine_priority>
```

OpenCOVID can provide the following output metrics, by default these are presented as values per 100,000 people per day:
- all_new_infections
- new_local_infections
- new_importations
- confirmed
- deaths
- hospital_beds
- hospital_admissions
- icu_beds
- icu_admissions
- currently_infected
- currently_infectious
- currently_symptomatic
- currently_isolated
- recovered
- n_vaccinated
- total_vaccinated
- n_doses
- n_infections*
- variant_prevalence*
- R_effective*
- seroprevalence*
- pop_susceptibility*
- pop_prevalence*
- seasonality*

Unless denoted with * all metrics can be dissagregated by age, variant infected with, and vaccine priority group.

Metrics denoted with * cannot be further dissagregated. The 'by' key is not required for these metrics.

### Vaccination

### Viral variants

### Scenarios

#### Single scenario

#### Array scenario

## The pipeline

### Getting started
- All model analyses are launched from `launch.R`
  - Preferred usage is to 'source' this file (without echo is ideal)
  - Alternative command line usage is to `cd` to the code directory then call `sh bash_launch.sh`
- Any required packages not currently installed are automatically downloaded and installed the very first time you run the model
  - If several large packages need to be installed this process may take a while (potentially 30-60 minutes)
  - See `dependencies.R` for a complete list of packages used by the model
  
### Directory structure
- This is all taken care of when `launch.R` (or `bash_launch.sh`) is called
  - The model will automatically re-point the current working directory
  - The model will also automatically create the necessary output directory structure
  
### Running a test simulation
- Set `do_step = 0` in `launch.R` and 'source' to run a test simulation with default model parameters
  - Once finished (this may take several minutes depending on settings), see `~\output\3_figures\demo\` for default output plots

### Generate a model calibration
- Set `do_step = 1` in `launch.R` and 'source' to calibrate the model
- Step 1 simulates the model for several thousand sampled parameter sets then trains a model emulator to interpolate the entire solution space
  - Due to the large number of model simulations required, this process uses the UniBas cluster and cannot be done locally
  - Even with the cluster, this process can take several hours - under default settings allow it to run overnight
  - Once finished, see `~\output\3_figures\demo\` for emulator performance plots
- Step 2 applies an MCMC algorithm to identify the global minima of the emulated space and produce parameter posteriors
  - Once finished, see `~\output\3_figures\demo\` for 'quick' plots to gauge calibration quality

### Run scenarios
- Set `do_step = 2` in `launch.R` and 'source' to simulate the calibrated baseline and any user-defined alternative scenarios
  
### Plot all your results
- Set `do_step = 3` in `launch.R` and 'source' to generate a series of plots that examine the analyses you have run
  - See `results.R` to see which results are to be generated
  - All plotting functionality itself resides in `plotting.R`
  - Plotting toggles can be turned on and off, see `plot_xxx` in `options.R`

### Final notes
- It is totally acceptable to run the full pipeline from end to end by setting `do_step = 1 : 3` in `launch.R`
- These notes are just to get you started, of course
  - For a full introduction to the model and it's capabilities contact the development team

## Release details

#### V2.0 release date: 09-12-2021

#### Development and maintenance:
* Andrew J. Shattock (andrewjames.shattock@swisstph.ch)

#### Contributors:
* Epke A. Le Rutte
* Robert P. Duenner
* Max Richter
* Sherrie L. Kelly
* Swapnoleena Sen
* Nakul Chitnis
* Melissa A. Penny

