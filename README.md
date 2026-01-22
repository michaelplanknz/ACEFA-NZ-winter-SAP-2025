# ACEFA-NZ-winter-SAP-2025

Code used in the ACEFA NZ winter situational assessment program 2025, as described in the paper `Multi-pathogen situational assessment and forecasting of respiratory disease in New Zealand`.


# Forecasting model

Matlab code for the forecasting model is in the folder `forecasting_code`.

The run the weekly forecasts, data will need to be saved in the folder `processed-data/2025/`.

Run the script `runRetrospective` to create retrospective forecasts from historic data files.
This will run forecasts for a series of weekly origin dates, using the closest available prior data file for each origin date (checking first if a forecast already exsists and only running those that don't).
Forecast results will be saved in a series of .mat files (one for each origin date and pathogen) in the folder `retrospective/`.

Run the script `plotRetrospective` to plot summary graphs for the retrospective forecasts, and calculate and plot forecast scores (CRPS on log(1+x) transoformed data). Figures will be saved in the folder `retrospective_figures`.
