# OpenCOVID

An individual-based model of SARS-CoV-2 transmission and COVID-19 disease dynamics developed at Swiss TPH. A stochastic, discrete-time model with age, risk-group, network, viral variant, and geo-spatial structure.

##### Version 1.0 (stable with R version 3.6.0)

### Before you begin
- If you will be regularly contributing to model development, `git clone` this repository to your desired location
  - If you wish to edit, commit, and push any model code, please branch off from master first
- If you would just like to experiment then please `git fork` instead

### Getting started
- All model analyses can be launched from `launch.R`
  - Preferred usage is to 'source' this file (without echo is ideal)
  - Alternative command line usage is to `cd` to the code directory then call `sh bash_launch.sh`
- Any required packages not currently installed are automatically downloaded and installed the very first time you run the model
  - If several large packages need to be installed this process may take a while (potentially 30-60 minutes)
  - See `dependencies.R` for a complete list of packages used by the model
  
### Directory structure
- This is all taken care of when `launch.R` (or `bash_launch.sh`) is called
  - The model will automatically re-point the current working directory
  - The model will also automatically create the necessary output directory structure
  
### Exploring the data
- Set `do_step = 0` in `launch.R` and 'source' to explore the data sources available
  - Once finished (~1 minute), see `~\output\4_figures\my_analysis\` for data visualisation plots
  - All the data used by the model is publicly available so access rights *shouldn't* be an issue
  
### Running a test simulation
- Set `do_step = NA` in `launch.R` and 'source' to run a test simulation with default model parameters
  - Once finished (5-30 minutes depending on settings), see `~\output\4_figures\my_analysis\` for default output plots
  
### Getting acquainted with key model options
- For applications, most of the interaction you will need to have with the model will be through `options.R` and several key input csv and xlsx files
- The file `options.R` contains numerous settings that you may want to explore, key settings include:
  - `analysis_name` := Point to a directory to store model output
  - `calibration_name` := Point to a new or previously created model calibration
  - `cantons` := Which canton(s) you wish to model
  - `scenarios` and `strategies` := Alternative user-defined past and future scenarios to model
  - `max_population_size` := Number of individuals to simulate in the model (< for speed, > for quality)
  
### Generate a model calibration
- Set `do_step = 1 : 2` in `launch.R` and 'source' to calibrate the model
- Step 1 simulates the model for several thousand sampled parameter sets then trains a model emulator to interpolate the entire solution space
  - Due to the large number of model simulations required, this process uses the UniBas cluster and cannot be done locally
  - Even with the cluster, this process can take several hours - under default settings allow it to run overnight
  - Once finished, see `~\output\4_figures\my_analysis\` for emulator performance plots
- Step 2 applies an MCMC algorithm to identify the global minima of the emulated space and produce parameter posteriors
  - Once finished, see `~\output\4_figures\my_analysis\` for 'quick' plots to gauge calibration quality
- See `~\config\model_parameters.xlsx` for which parameters are subject to calibration
- See `~\config\calibration_weights` for which metrics (and respective weightings) are considered in the likelihood function

### Run some alternative scenarios
- Set `do_step = 3` in `launch.R` and 'source' to simulate the calibrated baseline and any user-defined alternative scenarios
  - See `scenarios.R` for an example set of scenarios that can be simulated
  - See `strategies.R` for example strategies that can be simulated - strategies are essentially glorified, full-factorial scenarios
  - Either of these files can be adapted to represent scenarios that help to answer your research question(s)
  
### Plot all your results
- Set `do_step = 4` in `launch.R` and 'source' to generate a series of plots that examine the analyses you have run
  - See `results.R` to see which results are to be generated
  - All plotting functionality itself resides in `plotting.R`
  - Plotting toggles can be turned on and off, see `plot_xxx` in `options.R`

### Final notes
- It is totally acceptable to run the full pipeline from end to end by setting `do_step = 1 : 4` in `launch.R`
- These notes are just to get you started, of course
  - For a full introduction to the model and it's capabilities contact the development team

#### V1.0 release date: 2021-02-19

#### Authors:
* Andrew J. Shattock
* Epke A. Le Rutte
* Robert P. Duenner
* Nakul Chitnis
* Melissa A. Penny

#### Contributors:
* Swapnoleena Sen
* Sherrie L. Kelly

