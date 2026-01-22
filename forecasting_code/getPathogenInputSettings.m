function  [fNameData_cases, fNameData_hosp, test_types, useHospAsCases, pathogen_name_full] = getPathogenInputSettings(location_name, pathogen_name)

% Return data input settings for a specified location_name and pathogen_name 
% OUTPUTS: fNameData_case - file name identifier for case data
%          fNameData_hosp - file name identifier for hospitalisations data
%          (NZ only, otherwise "")
%          test_types - string set to "PCR" to target PCRs, "PCR_RAT" to target PCR+RATs, or "all" to target all test types
%          useHospAsCases - flag variable that is true if hospitalisation
%          are to be used as the target in place of cases (for NZ SARI
%          data)
%          pathogen_name_full - pathogen name to be queries in the date
%          information file (suffixed with "_NZhosp" to target SARI data)


if location_name == "NZ"
    fNameData_hosp = "NZ-hospitalisations-count-";
else
    fNameData_hosp = "";
end

if pathogen_name == "SARSCOV2"
    fNameData_cases = pathogen_name + "-all-tests-case-count-";
    test_types = "PCR";
    % Targeting case data
    useHospAsCases = false;
    pathogen_name_full = pathogen_name;
else
    fNameData_cases = pathogen_name + "-case-count-";
    test_types = "all";
    if location_name == "NZ"
        % Targeting NZ SARI data
        useHospAsCases = true;
        pathogen_name_full = pathogen_name + "_NZhosp";
    else
        % Targeting case data
        useHospAsCases = false;
        pathogen_name_full = pathogen_name;
    end
end



