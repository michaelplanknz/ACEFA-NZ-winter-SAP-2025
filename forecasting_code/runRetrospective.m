clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For reproducibility
rng(50522);

% Range of file dates (datestamp on folder/file name) to use
% 2025-06-19 is the earliest data in a consistent format
fileDates = datetime(2025, 6, 19):7:datetime(2025, 10, 9);        

% Season start - earliest origin date to run
firstOrigin = datetime(2025, 5, 4);

% Use contemporary data 
% - set to true to use the closest relevant file
% date for each origin date
% - set to false to use the latest file date for all origins
useContempData = true;

% Set to true to rerun all forecasts on the latest file, or false to read
% in previously saved forecasts where available (and rerun if not)
rerunFlag = false;

% Input/output filenames and locations
fileNames.dataFolder = "../processed-data/2025/";
fileNames.outputFolder = "../retrospective/";    

% Specify location and pathogen names
location_name = "NZ";
pathogen_name = ["SARSCOV2", "flu", "RSV"];


% Input filename identifiers (if no hospital data available, set fNameData_hosp to
% "")
fileNames.date_info = "date-information-";

% quantiles for plotting (choose 5 levels, with 0.5 as the middle one)
qt = [0.05, 0.25, 0.5, 0.75, 0.95]; 

    
% Get date information for the most recent file date
date_info_latest = getDateInfo(fileNames, fileDates(end));



nPathogens = length(pathogen_name);
for iPathogen = 1:nPathogens
    % Get pathogen specific data input settings
    [fileNames.cases, fileNames.hosp, test_types, useHospAsCases, pathogen_name_full] = getPathogenInputSettings(location_name, pathogen_name(iPathogen));

    
    % Get last origin date for this pathogen
    lastOrigin = date_info_latest.origin_date(date_info_latest.location == location_name & date_info_latest.pathogen == pathogen_name_full);

    % Set the last origin to be the Sunday (denoted by weekday=1) on or before the last specified orgin date
    lastOrigin = lastOrigin - (weekday(lastOrigin)-1);

    % Setup vector of origin dates to work with
    origins = fliplr(lastOrigin:-7:firstOrigin);
    nOrigins = length(origins);

    % Get model parameters
    par = getPar(location_name, pathogen_name(iPathogen));
    
    % Loop through origin dates
    for iOrigin = 1:nOrigins
        doOneOrigin(fileNames, useHospAsCases, fileDates, location_name, pathogen_name(iPathogen), pathogen_name_full, test_types, origins(iOrigin), par, qt, useContempData, rerunFlag);
    end
end
