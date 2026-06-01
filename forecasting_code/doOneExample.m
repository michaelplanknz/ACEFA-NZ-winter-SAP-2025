clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For reproducibility
rng(50522);

% File date (datestamp on folder/file name) 
fileDate = datetime(2025, 7, 31);        

% Input/output filenames and locations
fileNames.dataFolder = "../processed-data/2025/";
fileNames.outputFolder = "../outputs/"; 

% Specify location and pathogen names
location_name = "NZ";
pathogen_name = ["SARSCOV2", "flu", "RSV"];

% Input filename identifiers (if no hospital data available, set fNameData_hosp to "")
fileNames.date_info = "date-information-";

    
% Get date information
date_info = getDateInfo(fileNames, fileDate);

nPathogens = length(pathogen_name);
for iPathogen = 1:nPathogens
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set origin date
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    % Get pathogen specific data input settings
    [fileNames.cases, fileNames.hosp, test_types, useHospAsCases, pathogen_name_full] = getPathogenInputSettings(location_name, pathogen_name(iPathogen));
    
    originDate(iPathogen) = date_info.origin_date(date_info.location == location_name & date_info.pathogen == pathogen_name_full);
    
    % Get model parameters
    par = getPar(location_name, pathogen_name(iPathogen));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Run model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    fprintf('%s/%s, origin %s, file %s\n', location_name, pathogen_name(iPathogen), originDate(iPathogen), fileDate)
    particles = runOneForecast(fileNames, useHospAsCases, fileDate, location_name, pathogen_name(iPathogen), test_types, originDate(iPathogen), par);
    
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Format and save output 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
    writeOutput(particles, par, location_name, pathogen_name(iPathogen), fileNames, originDate(iPathogen), fileDate);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot graphs 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Note: it is possible to call this function with a later value of fileDateData to plot and score forecast
    % against subsequent data
    fileDateData = fileDate;
    plotGraphs(location_name, pathogen_name(iPathogen), test_types, fileNames, useHospAsCases, originDate(iPathogen), fileDate, fileDateData);
end

