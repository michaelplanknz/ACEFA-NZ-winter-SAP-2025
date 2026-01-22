function processed = readData(fileNames, useHospAsCases, fileDate, location_name, pathogen_name, test_types )

% Set full input filenames
fileDateString = datestr(fileDate, 'YYYY-mm-DD');

fName_cases = fileNames.dataFolder + fileDateString + "/" + fileNames.cases + fileDateString + ".csv";
fName_hosp = fileNames.dataFolder + fileDateString + "/" + fileNames.hosp + fileDateString + ".csv";

% Check if input files exist, otherwise return empty output variable
if exist(fName_cases, "file") || exist(fName_hosp, "file")
    
    if ~useHospAsCases
        
        fprintf('   Reading data from %s\n', fName_cases);
        
        % Read ACEFA data on cases
        opts = detectImportOptions(fName_cases);
        opts = setvartype(opts, {'location', 'pathogen'}, 'categorical');
        if ismember('test_type', opts.VariableNames)
            testTypeDataFlag = 1;
            opts = setvartype(opts, {'test_type'}, 'categorical');
        else
            testTypeDataFlag = 0;
        end
        data_cases = readtable(fName_cases, opts);
        
        
        
        % Filter to specified test types if test type data is available
        if testTypeDataFlag
            % Unstack case data by test type
            data_unstacked = unstack(data_cases, "cases", "test_type");
            if test_types == "PCR"
                data_unstacked.cases = data_unstacked.PCR;
            elseif test_types == "PCR_RAT"
                data_unstacked.cases = data_unstacked.PCR + data_unstacked.RAT;
            elseif test_types == "all"
                data_unstacked.cases = data_unstacked.PCR + data_unstacked.RAT + data_unstacked.Unknown;
            else
                error('test_types must be either "PCR", "PCR_RAT", or "all"\n')
            end
        else
            % No test type data available - just use all cases
            data_unstacked = data_cases;
            if test_types ~= "all"
                fprintf('   WARNING: test_types is set to %s, but no test_type field found in data - using all tests \n', test_types)
            end
        end
        
        % filter to specified location and pathogen and sort by notification date
        casesSorted = sortrows(data_unstacked(data_unstacked.location == location_name & data_unstacked.pathogen == pathogen_name, :), 'notification_date') ;
    else
        casesSorted = [];
    end
    
    
    % Read ACEFA data on hospitalisations if available
    if fileNames.hosp ~= ""
        fprintf('   Reading data from %s\n', fName_hosp);
        opts = detectImportOptions(fName_hosp);
        opts = setvartype(opts, {'location', 'pathogen'}, 'categorical');
        data_hosp = readtable(fName_hosp, opts);
    
        % filter to specified location and pathogen and sort by notification date
        hospSorted = sortrows(data_hosp(data_hosp.location == location_name & data_hosp.pathogen == pathogen_name, :), 'admission_date');
    else 
        hospSorted = [];
    end
    
    
    
    if ~isempty(casesSorted) & ~isempty(hospSorted)
        % Join case and hospitalisation tables
        processed = outerjoin(casesSorted, hospSorted, 'LeftKeys', {'notification_date', 'location', 'pathogen'}, 'RightKeys', {'admission_date', 'location', 'pathogen'}, 'MergeKeys', true, 'LeftVariables', {'notification_date', 'location', 'pathogen', 'cases'}, 'RightVariables', {'admission_date', 'location', 'pathogen', 'hospitalisations'});
        processed = renamevars(processed, {'notification_date_admission_date'}, {'date'});
    elseif ~isempty(casesSorted)
        processed = renamevars(casesSorted, {'notification_date'}, {'date'});
    elseif ~isempty(hospSorted)
        % No case data - use hospitalisation data in place of cases
        processed = renamevars(hospSorted, {'admission_date', 'hospitalisations'}, {'date', 'cases'});
    else
        error('no case or hospitalisation data found')
    end
    
    
    % Remove all data earlier than the last "NA" for cases
    % This deals with the fact that NZ SARI data has long runs of "NA" over
    % summer
    % "NA" for hospitalisations is OK - the model should deal with these (and
    % there will sometimes be some "NA"s at the end for hospitalisations
    % due to data truncation
    iLastNA = find(isnan(processed.cases), 1, 'last');
    if ~isempty(iLastNA)
        processed = processed(iLastNA+1:end, :);
    end

else
    processed = [];
end

