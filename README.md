# Multi-pathogen situational assessment and forecasting of respiratory disease in New Zealand -- 2025 results

Code used in the ACEFA NZ winter situational assessment program 2025, as described in the paper `Multi-pathogen situational assessment and forecasting of respiratory disease in New Zealand`.

## Epidemic trend analysis model
 R code for the Bayesian P-spline model, used to analyse past and current trends in the case and hospitalisation time series, is in the folder '

## Forecasting model

Matlab code for the forecasting model is in the folder `forecasting_code`.

To run the weekly forecasts there are two steps:
* Step 1. Run the script `runRetrospective` to create retrospective forecasts from historic data files.
This will run forecasts for a series of weekly origin dates, using the closest available prior data file for each origin date (checking first if a forecast already exsists and only running those that don't).
Forecast results will be saved in a series of .mat files (one for each origin date and pathogen) in the folder `retrospective/`.
* Step 2. Run the script `plotRetrospective` to plot summary graphs for the retrospective forecasts, and calculate and plot forecast scores (CRPS on log(1+x) transoformed data). This reads in the .mat files produced in Step 1 and saves figures in the folder `retrospective_figures`.

Note: this public repo does not contain the input data for confidentiality reasons, but does contain the .mat files produced by running Step 1. Hence, running step 2 will produce graphs of the forecast model outputs like those that appear in the paper, but without data to compare to. 

To rerun the forecasting model, the appropriate input data will need to be stored in the folder `processed-data/2025/`.

