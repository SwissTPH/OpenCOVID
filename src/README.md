---
output:
  html_document: default
  word_document: default
---
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
- Step 1 determines the average number of daily contacts required to achieve a pre-defined effective reproduction number at the start of the simulation period
- This is achieved using a stochastic decent algorithm that evalutes the effective reproduction number over the first few time steps of the model

### Run scenarios
- Set `do_step = 2` in `launch.R` and 'source' to simulate the baseline and any user-defined alternative scenarios
- Due to the potentially large number of model simulations required, this process uses the UniBas cluster to run simulations in parallel and cannot be done locally
- To run step 2 outside of the University of Basel network, you will first need to configure parallel simualtions, ideally on a cluster
  
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

#### V2.0 release date: 10-12-2021

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

