clear
close all

% Script to copy the required NZ data from the ACEFA repo into the
% current directory
% NB Make sure these files are not uploaded to the public repo by including
% *count*.csv in the .gitignore 

% Range of file dates (datestamp on folder/file name) needed
fileDates = [datetime(2025, 6, 5):7:datetime(2025, 10, 9), datetime(2025, 10, 23)];    

% Name of local directory in which to put data files
local_dir = "2025\";

% Path to ACEFA directory in which data is stored
ACEFA_dir = "V:\mpl31\ACEFA_forecasting_repos\forecast-hub\processed-data\2025\";

nFileDates = length(fileDates);

% Loop through the set of file dates
for iFileDate = 1:nFileDates
    dateString = string(datetime(fileDates(iFileDate), "Format", "yyyy-MM-dd"));

    % Set source and destination directories
    source_dir = ACEFA_dir + dateString + "\";
    dest_dir = local_dir + dateString + "\";

    % Specify file names for hospitalisation and case data (SARSCOV2 only)
    dates_fName = "date-information-" + dateString + ".csv";
    hosp_fName = "NZ-hospitalisations-count-" + dateString + ".csv";
    hosp_fName_old = "SARSCOV2-NZ-hospitalisations-count-" + dateString + ".csv";
    case_fName = "SARSCOV2-all-tests-case-count-" + dateString + ".csv";

    % Copy the date information across 
    copyfile(source_dir+dates_fName, dest_dir);

    % Copy the hopsitalisation data across (as this only contains NZ data)
    % Try the standard file name first
    if exist(source_dir+hosp_fName, 'file' )
        copyfile(source_dir+hosp_fName, dest_dir);
    % If the standard file name doesn't exist, try the old file name which
    % was used in some of the early rounds
    elseif exist(source_dir+hosp_fName_old, 'file' )
        copyfile(source_dir+hosp_fName_old, dest_dir);
    else
        fprintf('Warning: no hopsitalisation data found for file date %s\n', dateString)
    end

    % Read in the case data
    tbl = readtable(source_dir+case_fName);

    % Filter to NZ data only
    NZflag = tbl.location == "NZ";
    tbl_NZ = tbl(NZflag, :);

    % Write new file with NZ data in destination directory
    writetable(tbl_NZ, dest_dir+case_fName);
end







