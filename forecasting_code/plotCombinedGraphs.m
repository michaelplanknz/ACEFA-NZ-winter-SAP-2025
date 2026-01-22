function plotCombinedGraphs(location_name, pathogen_name, fileNames, originDate, fileDate, fileDateData)

% Read back in a previously saved forecasts and plot alongside data (which can be
% the same data or subsequent)

% INPUTS:
%    location_name - string specifying jursidction
%    pathogen_name - string or vector of string specifying pathogen names
%    to plot
%    fileNames - structure with folder and file name info
%    originDate - vector of datetimes specifying the originDate of the forecast (used
%    to identify the name of hte file containing the forecast results) for
%    each pathogen
%    fileDate - datetime specifying the fileDate used to generate the forecast (used
%    to identify the name of hte file containing the forecast results)
%    fileDateData - datetime specifying the fileDate of the file containing
%    the data to be plotted



% Choose whether to plot data as raw or smoothed (weekly moving avg)
smoothDataFlag = false;

% Choose whether to plot model outputs with/without DOW 
removeModelDOW = false;

% colour settings
darkBlue = [0.1 0.2 0.6];

% Days to plot back from origin date
plotBackDays = 90; 

nPathogens = length(pathogen_name);

% Latest origin date (for setting plot ranges)
latestOrigin = max(originDate); 

% Plot figure
h = figure;
h.Position = [ 19    19   1100   320*nPathogens];
tiledlayout(nPathogens, 2, "TileSpacing", "compact");
iTile = 1;
for iPathogen = 1:nPathogens
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load forecast results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Get pathogen specific data input settings
    [fileNames.cases,  fileNames.hosp, test_types, useHospAsCases, ~] = getPathogenInputSettings(location_name, pathogen_name(iPathogen));


    % Create a strting from the file date for putting in file names
    originDateString = datestr(originDate(iPathogen), 'YYYY-mm-DD');
    fileDateString = datestr(fileDate, 'YYYY-mm-DD');
    
    % .mat file name containing quantiles of relevant variables for plotting
    fIn = fileNames.outputFolder + "origin-" + originDateString + "-file-" + fileDateString + "-" + location_name + "-" + pathogen_name(iPathogen) + "-quantiles.mat";
    
    % Load file with specified variables
    load(fIn, 't', 'qts', 'samples', 'par');
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Read epi data
    processed = readData(fileNames, useHospAsCases, fileDateData, location_name, pathogen_name(iPathogen), test_types);
    
    % Call getTrainingData with orignDate = max(processed.date) to return all
    % data as 'training' data (for plotting)
    [~, tData, nCases, nHosp] = getTrainingData(processed, max(processed.date), par.maxTimeBack, par.caseIgnoreDays, par.hospIgnoreDays);
    
    
    % Calculate 7 day moving average of cases
    nCasesSmoothed = smoothdata(nCases, 'movmean', 7);
    nHospSmoothed = smoothdata(nHosp, 'movmean', 7);  

    % Discard last part of smoothed data due to incomplete data for
    % smoothing (and truncation according to par.caseIgnoreDays)
    nCasesSmoothed(end-(3+par.caseIgnoreDays-1):end) = nan;
    nHospSmoothed(end-(3+par.hospIgnoreDays-1):end) = nan;
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plotting
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tPlot = latestOrigin + [-plotBackDays, par.timeHorizon];
    
    % Skip left hand plot if use modelling hospitalisations only
    if useHospAsCases
        iTile = iTile+1;
    end
    
    nexttile(iTile);
    iTile = iTile+1;
    
    if removeModelDOW
        y = qts.Cndw;
    else
        y = qts.Ct;
    end
    myModelPlot(t, y, darkBlue);

    if smoothDataFlag
        plot(tData, nCasesSmoothed, 'r.')
    else
        plot(tData, nCases, 'r.')
    end
    xline(originDate(iPathogen), 'k--')
    xlim(tPlot)
    grid on
    if pathogen_name(iPathogen) == "SARSCOV2"
        pathogen_title = "SARS-CoV-2";
    elseif pathogen_name(iPathogen) == "flu"
        pathogen_title = "Influenza";
    elseif pathogen_name(iPathogen) == "RSV"
        pathogen_title = "RSV";
    end
    if useHospAsCases
        ylabel('daily admissions (smoothed)')
        title(pathogen_title + " hospitalisations")
    else
        ylabel('daily PCR notifications (smoothed)')
        title(pathogen_title + " cases")
    end


    if ~useHospAsCases
        nexttile(iTile);
        iTile = iTile+1;
 
        if removeModelDOW
            y = qts.Andw;
        else
            y = qts.At;
        end
        myModelPlot(t, y, darkBlue);

        if smoothDataFlag
            plot(tData, nHospSmoothed, 'r.')
        else
            plot(tData, nHosp, 'r.')
        end
        xline(originDate(iPathogen), 'k--')
        xlim(tPlot)
        grid on
        ylabel('daily admissions (smoothed)')
        title(pathogen_title + " hospitalisations")
    end

end

% Save output figure (unless working in scratch mode)
if ~contains(fileNames.outputFolder, "scratch")      
    fName = sprintf('../figures/file-%s-%s.png', datetime(fileDate, 'Format', 'yyyy-MM-dd'), location_name);
    saveas(h, fName);
end


    