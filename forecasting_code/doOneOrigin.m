function doOneOrigin(fileNames, useHospAsCases, fileDates, location_name, pathogen_name, pathogen_name_full, test_types, originDate, par, qt, useContempData, rerunFlag)

% number of samples to save (for scoring)
nSav = 2000;           

if useContempData
    % Select the file date relevant for this origin date
    fileDateCurrent = getRightFileDate(originDate, fileDates, fileNames, location_name, pathogen_name_full);
else
    % Use most recent file date
    fileDateCurrent = fileDates(end);
end


% Create a strting from the file data and origin date for putting in file names
fileDateString = datestr(fileDateCurrent, 'YYYY-mm-DD');
originDateString = datestr(originDate, 'YYYY-mm-DD');

% Check if previously saved forecast available
fName = fileNames.outputFolder + "origin-" + originDateString + "-file-" + fileDateString + "-" + location_name + "-" + pathogen_name + "-quantiles.mat";

fprintf('%s/%s, origin %s, file %s - ', location_name, pathogen_name, originDate, fileDateCurrent)
if exist(fName, 'file') & ~rerunFlag
     fprintf('found saved forecast at %s\n', fName)
else
    % No matching forcast file found - run new forecast
    fprintf('running model... \n')

    % Run model
    particles = runOneForecast(fileNames, useHospAsCases, fileDateCurrent, location_name, pathogen_name, test_types, originDate, par);
    
    % Randomly select trajectories for saving
    if nSav < par.nParticles
        iSav = randsample(par.nParticles, nSav);
    else
        iSav = 1:par.nParticles;
        nSav = par.nParticles;
    end

    % Extract quantiles and samples for saving
    [t, qts, samples] = getParticleQuantiles(particles, originDate, qt, iSav );

    % Save results as a .mat file
    fprintf('...writing to %s\n', fName)
    save(fName, 't', 'qts', 'samples', 'par');
end   

