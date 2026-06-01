function writeOutput(particles, par, location_name, pathogen_name, fileNames, originDate, fileDate )

% Number of particles to save
nSav = 2000;           

% quantiles for plotting (choose 5 levels, with 0.5 as the middle one)
qt = [0.025, 0.25, 0.5, 0.75, 0.975]; 

            
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


