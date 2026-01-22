function plotGraphs(location_name, pathogen_name, test_types, fileNames, useHospAsCases, originDate, fileDate, fileDateData)

% Read back in a previously saved forecast and plot alongside data (which can be
% the same data or subsequent)

% INPUTS:
%    location_name - string specifying jursidction
%    pathogen_name - string specifying pathogen "SARSCOV2", "flu" or "RSV"
%    test_types - string specify which test types are targeted by the
%    forecast "PCR", "PCR_RAT" or "all"
%    fileNames - structure with folder and file name info
%    useHospAsCases - flag variable indicating whether to use hospital
%    admissions in place of cases as forecast target
%    originDate - datetime specifying the originDate of the forecast (used
%    to identify the name of hte file containing the forecast results)
%    fileDate - datetime specifying the fileDate used to generate the forecast (used
%    to identify the name of hte file containing the forecast results)
%    fileDateData - datetime specifying the fileDate of the file containing
%    the data to be plotted


% colour settings
darkBlue = [0.1 0.2 0.6];

% Days to plot back from origin date
plotBackDays = 90; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load forecast results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a strting from the file date for putting in file names
originDateString = datestr(originDate, 'YYYY-mm-DD');
fileDateString = datestr(fileDate, 'YYYY-mm-DD');

% .mat file name containing quantiles of relevant variables for plotting
fIn = fileNames.outputFolder + "origin-" + originDateString + "-file-" + fileDateString + "-" + location_name + "-" + pathogen_name + "-quantiles.mat";

% Load file with specified variables
load(fIn, 't', 'qts', 'samples', 'par');


% Use default values of parameters that did not exist when forecast was made
if ~ismember('caseIgnoreDays', fields(par))
    par.caseIgnoreDays = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read epi data
processed = readData(fileNames, useHospAsCases, fileDateData, location_name, pathogen_name, test_types);

% Call getTrainingData with orignDate = max(processed.date) to return all
% data as 'training' data (for plotting)
[~, tData, nCases, nHosp] = getTrainingData(processed, max(processed.date), par.maxTimeBack, par.caseIgnoreDays, par.hospIgnoreDays);


% Calculate growth rate via Wallinga-Lipsitch relation
r_qt = calcGrowth_WL(qts.Rt, par.GTD);

% Calculate 7 day moving average of cases
nCasesSmoothed = smoothdata(nCases, 'movmean', 7);

% Discard last part of smoothed data due to incomplete data for
% smoothing (and truncation according to par.caseIgnoreDays)
nCasesSmoothed(end-(3+par.caseIgnoreDays-1):end) = nan;







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tPlot = [max(tData(1), originDate - plotBackDays), t(end)];

h = figure;
if location_name == "NZ" & pathogen_name == "SARSCOV2"
    h.Position = [      212         274        1359         674];
    tiledlayout(2, 3, "TileSpacing", "compact");
else
    h.Position = [      212         274        920         674];
    tiledlayout(2, 2, "TileSpacing", "compact");
end

nexttile
myModelPlot(t, qts.Rt, darkBlue);
xline(originDate, 'k--')
yline(1, 'k--')
xlim(tPlot)
ylabel('R(t)')
grid on

nexttile
myModelPlot(t, qts.Ct, darkBlue);
plot(tData, nCases, 'r.')
xline(originDate, 'k--')
xlim(tPlot)
if location_name == "NZ" & pathogen_name ~= "SARSCOV2"
    ylabel('daily admissions')
else
    ylabel('daily PCR notifications')
end
grid on



if location_name == "NZ" & pathogen_name == "SARSCOV2"
    nexttile
    myModelPlot(t, qts.Pt, darkBlue);
    xline(originDate, 'k--')
    xlim(tPlot)
    ylim([0 inf])
    ylabel('case-hospitalisation ratio')
    grid on
end




nexttile
myModelPlot(t, r_qt, darkBlue);
xline(originDate, 'k--')
yline(0, 'k--')
xlim(tPlot)
ylabel('r(t) (days^{-1})')
grid on


nexttile
myModelPlot(t, qts.Cndw, darkBlue);
plot(tData, nCasesSmoothed, 'r.')
xline(originDate, 'k--')
xlim(tPlot)
if useHospAsCases
    ylabel('daily admissions (smoothed)')
else
    ylabel('daily PCR notifications (smoothed)')
end
grid on



if location_name == "NZ" & pathogen_name == "SARSCOV2"
    nexttile
    myModelPlot(t, qts.At, darkBlue);
    plot(tData, nHosp, 'r.')
    xline(originDate, 'k--')
    xlim(tPlot)
    ylabel('daily admissions')
    grid on
end

sgtitle(sprintf('%s, %s, origin date %s', location_name, pathogen_name, originDate))

% Save output figure (unless working in scratch mode)
if ~contains(fileNames.outputFolder, "scratch")
    fName = sprintf('../figures/origin-%s-file-%s-%s-%s.png', datetime(originDate, 'Format', 'yyyy-MM-dd'), datetime(fileDate, 'Format', 'yyyy-MM-dd'), location_name, pathogen_name);
    saveas(h, fName);
end




