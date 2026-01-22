function writeOutput(particles, par, useHospAsCases, location_name, pathogen_name, fileNames, originDate, fileDate )

% quantiles for plotting (choose 5 levels, with 0.5 as the middle one)
qt = [0.05, 0.25, 0.5, 0.75, 0.95]; 

% Parameters for formatting ACEFA outputs
nSav = 2000;            % number of particles to save
backHorizon = 14;       % number of days to backcast relative to origin date
            
% Create a strting from the file date for putting in file names
fileDateString = datestr(fileDate, 'YYYY-mm-DD');
originDateString = datestr(originDate, 'YYYY-mm-DD');



% Randomly select trajectories for saving
if nSav < par.nParticles
    iSav = randsample(par.nParticles, nSav);
else
    iSav = 1:par.nParticles;
    nSav = par.nParticles;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save quantles of hidden states etc. for plotting and subsequent
% comparison with future data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('   Formatting particle results...  ')

% Extract quantiles for plotting
[t, qts, samples] = getParticleQuantiles(particles, [], qt, iSav);


% Set output .mat file name
fOutMat = fileNames.outputFolder + "origin-" + originDateString + "-file-" + fileDateString + "-" + location_name + "-" + pathogen_name + "-quantiles.mat";

% Save results as a .mat file
fprintf('writing to %s\n ', fOutMat)
save(fOutMat, 't', 'qts', 'samples', 'par');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save ACEFA formatted results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('   Formatting ACEFA results...  ')

tInd = find( particles.t >= originDate-backHorizon & particles.t <= originDate+par.timeHorizon );
nDays = length(tInd);
nRows = nSav*nDays;

% First part of output table: daily case incidence including DOW effect
% (case_incidence)
outTab1.round_id = repmat( datetime(fileDate, 'Format', 'yyyy-MM-dd') , nRows, 1 );
outTab1.origin_date = repmat( datetime(originDate, 'Format', 'yyyy-MM-dd') , nRows, 1 );
if useHospAsCases
    target_name = "hosp";
else
    target_name = "case_incidence";
end
outTab1.target = repmat(target_name, nRows, 1);
outTab1.horizon = repmat(reshape(int64(days(particles.t(tInd)-originDate)), nDays, 1), nSav, 1);
outTab1.location = repmat(location_name, nRows, 1);
outTab1.pathogen = repmat(pathogen_name, nRows, 1);
outTab1.output_type = repmat("sample", nRows, 1);
outTab1.output_type_id = repelem( int64(1:nSav)', nDays, 1 );
outTab1.value = reshape( int64(particles.Ct(iSav, tInd ))', nSav*nDays, 1);

% Second part of output table: daily case incidence with DOW effect removed
% (case_incidence_smoothed)
outTab2.round_id = repmat( datetime(fileDate, 'Format', 'yyyy-MM-dd') , nRows, 1 );
outTab2.origin_date = repmat( datetime(originDate, 'Format', 'yyyy-MM-dd') , nRows, 1 );
if useHospAsCases
    target_name = "hosp_smoothed";
else
    target_name = "case_incidence_smoothed";
end
outTab2.target = repmat(target_name, nRows, 1);
outTab2.horizon = repmat(reshape(int64(days(particles.t(tInd)-originDate)), nDays, 1), nSav, 1);
outTab2.location = repmat(location_name, nRows, 1);
outTab2.pathogen = repmat(pathogen_name, nRows, 1);
outTab2.output_type = repmat("sample", nRows, 1);
outTab2.output_type_id = repelem( int64(1:nSav)', nDays, 1 );
outTab2.value = reshape( int64(particles.Cndw(iSav, tInd ))', nSav*nDays, 1);

% Combine into a single table
outTab = [struct2table(outTab1); struct2table(outTab2)];


% Write additional output data for hospitalisation forecast if available
if ~useHospAsCases & any(~isnan(particles.At))
    % Output with DOW included
    outTab1.round_id = repmat( datetime(fileDate, 'Format', 'yyyy-MM-dd') , nRows, 1 );
    outTab1.origin_date = repmat( datetime(originDate, 'Format', 'yyyy-MM-dd') , nRows, 1 );
    outTab1.target = repmat("hosp", nRows, 1);
    outTab1.horizon = repmat(reshape(int64(days(particles.t(tInd)-originDate)), nDays, 1), nSav, 1);
    outTab1.location = repmat(location_name, nRows, 1);
    outTab1.pathogen = repmat(pathogen_name, nRows, 1);
    outTab1.output_type = repmat("sample", nRows, 1);
    outTab1.output_type_id = repelem( int64(1:nSav)', nDays, 1 );
    outTab1.value = reshape( int64(particles.At(iSav, tInd ))', nSav*nDays, 1);

    % Output with DOW removed
    outTab2.round_id = repmat( datetime(fileDate, 'Format', 'yyyy-MM-dd') , nRows, 1 );
    outTab2.origin_date = repmat( datetime(originDate, 'Format', 'yyyy-MM-dd') , nRows, 1 );
    outTab2.target = repmat("hosp_smoothed", nRows, 1);
    outTab2.horizon = repmat(reshape(int64(days(particles.t(tInd)-originDate)), nDays, 1), nSav, 1);
    outTab2.location = repmat(location_name, nRows, 1);
    outTab2.pathogen = repmat(pathogen_name, nRows, 1);
    outTab2.output_type = repmat("sample", nRows, 1);
    outTab2.output_type_id = repelem( int64(1:nSav)', nDays, 1 );
    outTab2.value = reshape( int64(particles.Andw(iSav, tInd ))', nSav*nDays, 1);
    
    % Append to existing table
    outTab = [outTab; struct2table(outTab1); struct2table(outTab2)];
end



% Set output file name
fOut = fileNames.outputFolder + fileDateString + "-" + location_name + "-" + pathogen_name + "-forecast.parquet";

% Write output parquet file
fprintf('writing to %s\n', fOut)
parquetwrite(fOut, outTab);


