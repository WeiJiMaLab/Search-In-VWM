# Color localization data and code

This repository recreates the figures in the paper: "A Computational Approach to Search in Visual Working Memory"

## Data

Run the function `All_Data.m` in the data-analysis folder. This script loads in `Data_Wen-Chuang.mat` to reformat the raw data into `alldata.mat`, which is a 1 * Nsubj struct with the  field `data`, which each in turn have the following fields, all with size 1 * Ntrials (here 640), with the exception of stim_order and col_val, with size 8*640. The data contains the following fields:

- `set_size`: N in [2, 4, 6, 8]
- `stim_order`: stimulus indices organized from top center, going clockwise; 1 is always the place holder for the target stimuli, and 2 for the non-target/distractor
- `spatial_dist`: between target and nontarget, in [1,2,3,4]
- `col_val`: color values represented in radians from [-180°,180°]
- `response`: 1 (correct) or 0 (incorrect)
- `reaction_time` (seconds)
- `col_dist`: absolute circular_distance between target and non-target colors within domain of [0°,180°] (uniformly distributed)

For creation of simulated data, use the script `sim_data_all_models.m` that creates a simulated data structure like the real data. All the model fitting scripts described below should work for the simulated data. To check the parameter recovery, use the script `param_recovery.m`. `TODO for ADITI: ADD THESE SCRIPTS IN`

## Model fitting

#### Resource Modeling
We explore 4 models in this approach:

- Equal Precision (EP) resource model
- Equal Precision model with fixed item limit (EPF) resource+slots model  
- Variable Precision (VP) resource model
- Variable Precision model with fixed item limit (VPF) resource+slots model

The code for the modeling code for these models can be found in the `resource_parametric` folder. Note that for all the scripts here, an identifying tag variable `mi` is used to determine which model will be run (EP, EPF, VP, or VPF).

The function `Master_fitting_real.m` in the `resource_parametric` folder calls the function `fit_pred.m` and finds for each subject the parameters that maximize the loglikelihood `LL.m` and calculates the model predictions with the function `Pred_change.m`.

Both `LL.m` and  `Pred_change.m` call the function `calc_prob_corr.m`. Additional section inside the modeling folder `model_fits_cluster2` stores parameter values for all 4 models, 14 subjects, and all runs (currently set to 20).

Note that tau and lambda_alpha are optimized in log-space, while Jbar is optimized in regular exponentiated form. Both optimization are done with BADS.

##### Resource Modeling (non-parametric)

We also test a non-parametric instance of the best performing paramatric resource model from above (in this case, the VP model), and determine whether having no assumption that encoding precision decreases with set size according to a power law function. To run modeling for this, run the `Master_Fitting.m` script in the `resource_non-parametric` folder.

#### Resource-rational Modeling

In this approach, we use the encoding and decision-making rule according to the most optimal model of the resource models (VP) -- meaning, we assume variability in precision as represented by the parameter 'tau' and do not set an item limit.

We consider two variations of the resource-rational model---one with a linear neural cost, and the other with a neural cost defined by a power law function. The model-fitting scripts are in the `resource-rational_linear-cost` and `resource-rational_power-law` folders respectively.

In both folders, `optim_normative.m` holds the code for the first round of optimization (tau and lambda_alpha). This script then calls `LL_cost.m` which holds the observer's optimization of precision. `LL_cost.m` calls `calc_prob_corr.m` in its script to calculate the total cost value of the observer. The only difference in these two approaches is the cost function the observer is ultimately trying to minimize (which is defined in the `Fun_Jbars_optimizeNEW_pow.m` script)

## Analysis and visualization scripts

All the plotting scripts are within one folder called `plotting`. Below is a list of the plotting script name and the figures it maps to:

- `analysis_all_summary_statistics.m` provides an analysis of proportion correct behavior across trials, with each experimental condition: set size, color distance, and spatial distance. 
- `analysis_all_parametric.m` plots the parametric resource models. Note the script asks for an `mi` variable to be set which specifies the model index that the script will produce a figure for. This model index can be EP, EPF, VP, or VPF. 
- `analysis_all_model_non_parametric.m` produces the modeling fits for the non-parametric VP resource model.
- `analysis_all_model_comparison.m` uses the summary paramater mat files for each model fit in the `model_params_summary` folder (which is automatically saved if you run the modeling fits using any of the above scripts), and creates a model comparision plot with the files that are loaded in the script. Note that by default, all the comparisions are relative to the third model---but this can be changed through the `compare_index` variable in this script. `TODO FOR ADITI -- CLEAN THIS UP.`

