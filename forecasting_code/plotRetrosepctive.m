clear 
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose whether to plot data as raw or smoothed (weekly moving avg)
smoothDataFlag = false;

% Choose whether to plot model outputs with/without DOW
removeModelDOW = false;


% Range of file dates (datestamp on folder/file name) to use
% 2025-06-19 is the earliest data in a consistent format
fileDates = datetime(2025, 6, 19):7:datetime(2025, 9, 25); 

% If later data is available for evaluation data purposes, specify the file
% date here (otherwise set to the latest data in fileDates):
%evalFileDate = fileDates(end);
evalFileDate = datetime(2025, 10, 23);

% Season start - earliest origin date to run
firstOrigin = datetime(2025, 5, 4);

% Use contemporary data 
% - set to true to use the closest relevant file
% date for each origin date
% - set to false to use the latest file date for all origins
useContempData = true;

% Input/output filenames and locations
fileNames.dataFolder = "../processed-data/2025/";
fileNames.outputFolder = "../retrospective/";    
fileNames.figureFolder = "../retrospective_figures/";    


% Specify location and pathogen names
location_name = "NZ";
pathogen_name = ["SARSCOV2", "flu", "RSV"];

% Origin-file min delay
% For Covid - the origin date is typically 4 days before the file date, so
% take the next available file date after the origin, i.e. delay = 0
% For flu & RSV, the origin date is typically 11 days before the file date,
% so take the next available file date at least 7 days after the origin,
% i.e. delay = 7
originFileDelay = [0, 7, 7];

% Input filename identifiers (if no hospital data available, set fNameData_hosp to
% "")
fileNames.date_info = "date-information-";



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in model outputs and plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% colour settings
nCols = 5;
clrs = colororder;
letters = ["(a)", "(b)", "(c)", "(d)"];


% Set up mutli-pathogen figure
h = figure(100);
h.Position = [ 680   467   913   511];
tiledlayout(2, 2, "TileSpacing", "compact");

% Get date information for the most recent file date
date_info = getDateInfo(fileNames, fileDates(end));

nPathogens = length(pathogen_name);
for iPathogen = 1:nPathogens
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get last origin date
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    % Get pathogen specific data input settings
    [fileNames.cases, fileNames.hosp, test_types, useHospAsCases, pathogen_name_full] = getPathogenInputSettings(location_name, pathogen_name(iPathogen));

    % Get last origin date for this pathogen
    lastOrigin = date_info.origin_date(date_info.location == location_name & date_info.pathogen == pathogen_name_full);

    % Set the last origin to be the Sunday (denoted by weekday=1) on or before the last specified orgin date
    lastOrigin = lastOrigin - (weekday(lastOrigin)-1);

    origins = fliplr(lastOrigin:-7:firstOrigin);
    nOrigins = length(origins);

    % Get model parameters
    par = getPar(location_name, pathogen_name(iPathogen));

    % Initialise array for storing forecast scores
    scoreCases = nan(nOrigins, par.timeHorizon);
    scoreHosp = nan(nOrigins, par.timeHorizon);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read data (from specified file date for evaluation data) for plotting
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % Read epi data
    processed = readData(fileNames, useHospAsCases, evalFileDate, location_name, pathogen_name(iPathogen), test_types);

    if ~isempty(processed)
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

        % Date range for plotting
        tPlot = [max(tData(1), origins(1) - 28), origins(end)+par.timeHorizon];

    else
        % If no input data exists, make empty arrays so model can still be
        % plotted without data
        tData = [];
        nCases = [];
        nHosp = [];
        nCasesSmoothed = [];
        nHospSmoothed = [];
        tPlot = [origins(1) - 28, origins(end)+par.timeHorizon];
    end


    % Set up 2x2 set of axes to plot into one origin date at a time
    h = figure;
    h.Position = [ 50    50   1600   800];
    tiledlayout(2, 2, "Tilespacing", "compact");
    iTile = 1;

    for iOrigin = 1:nOrigins

        if useContempData
            % Select the file date relevant for this origin date
            fileDateCurrent = getRightFileDate(origins(iOrigin), fileDates, fileNames, location_name, pathogen_name_full);
        else
            % Use most recent file date
            fileDateCurrent = evalFileDate;
        end
        
        % Create a strting from the file data and origin date for putting in file names
        fileDateString = datestr(fileDateCurrent, 'YYYY-mm-DD');
        originDateString = datestr(origins(iOrigin), 'YYYY-mm-DD');

        % Load relevant forecast results
        fName = fileNames.outputFolder + "origin-" + originDateString + "-file-" + fileDateString + "-" + location_name + "-" + pathogen_name(iPathogen) + "-quantiles.mat";
        fileRead = load(fName, 't', 'qts', 'samples', 'par');

        %  Store the relevant variables in a structure array
        results(iOrigin).t = fileRead.t;
        results(iOrigin).qts = fileRead.qts;
        results(iOrigin).samples = fileRead.samples;

        % Calculate forecast score (if data is available)
        if ~isempty(tData)
            tDataEnd = max(tData(~isnan(nCases)));
            if tDataEnd >= origins(iOrigin) + par.timeHorizon
                ind1 = results(iOrigin).t > origins(iOrigin) & results(iOrigin).t <= origins(iOrigin)+par.timeHorizon;
                ind2 = tData > origins(iOrigin) & tData <= origins(iOrigin)+par.timeHorizon;
                assert(sum(ind1) == par.timeHorizon);
                assert(sum(ind2) == par.timeHorizon);
                % Use log transformed data (log(1+x)) - see Bosse et al, https://doi.org/10.1371/journal.pcbi.1011393
                logTransFlag = true;
                scoreCases(iOrigin, :) = calcCRPS(results(iOrigin).samples.Ct(:, ind1)', nCases(ind2)', logTransFlag);
                scoreHosp(iOrigin, :) = calcCRPS(results(iOrigin).samples.At(:, ind1)', nHosp(ind2)', logTransFlag);
            end
        end

        % Plotting
        nexttile(iTile);
        iTile = mod(iTile, 4)+1;

        % Colour index
        iCol = mod(floor((iOrigin-1)/4), nCols)+1;

        % Select variable to plot
        if removeModelDOW
             y = results(iOrigin).qts.Cndw;
        else
            y = results(iOrigin).qts.Ct;
        end

        % Plot forecast
        myModelPlot(results(iOrigin).t, y, clrs(iCol, :));
    end

    % Plot a copy of the same data into each axis
    yUpper = zeros(4, 1);
    for iTile = 1:4
        nexttile(iTile);
        if smoothDataFlag
            plot(tData, nCasesSmoothed, 'k.')
            ylbl = " (smoothed)";
        else
            plot(tData, nCases, 'k.')
            ylbl = "";
        end
        xlim(tPlot)
        if useHospAsCases
            ylabel("daily admissions" + ylbl)
        else
            ylabel("daily PCR notifications" + ylbl')
        end
        % Record the upper limit of the y-axis in each tile to make them
        % consistent
        h = gca;
        yUpper(iTile) = h.YLim(2);
    end
    % Set y-axis limits
    yMax = max(yUpper);
    for iTile = 1:4
        nexttile(iTile);
        ylim([0 yMax])
        grid on
    end
    if pathogen_name(iPathogen) == "SARSCOV2"
        pathogen_title = "SARS-CoV-2";
    elseif pathogen_name(iPathogen) == "flu"
        pathogen_title = "Influenza";
    elseif pathogen_name(iPathogen) == "RSV"
        pathogen_title = "RSV";
    end
    if useHospAsCases
       sgtitle(pathogen_title + " hospitalisations");
       fName = pathogen_name(iPathogen) + "_hosp.png";
    else
       sgtitle(pathogen_title + " cases");
       fName = pathogen_name(iPathogen) + "_cases.png";
    end
    saveas(h, fileNames.figureFolder+fName);



    % Make hospitalisation plot as well 
    if ~useHospAsCases
        h = figure;
        h.Position = [ 100    50   1600   800];
        tiledlayout(2, 2, "Tilespacing", "compact");
        iTile = 1;
    
        for iOrigin = 1:nOrigins
            nexttile(iTile);
            iTile = mod(iTile, 4)+1;

            % Colour index
            iCol = mod(floor((iOrigin-1)/4), nCols)+1;

            % Select variable to plot
            if removeModelDOW
                 y = results(iOrigin).qts.Andw;
            else
                y = results(iOrigin).qts.At;
            end

            % Plot forecast
            myModelPlot(results(iOrigin).t, y, clrs(iCol, :));
        end
    
        yUpper = zeros(4, 1);
        for iTile = 1:4
            nexttile(iTile);
            % Plot data
            if smoothDataFlag
                plot(tData, nHospSmoothed, 'k.')
                ylbl = " (smoothed)";
            else
                plot(tData, nHosp, 'k.')
                ylbl = "";
            end
            xlim(tPlot)
            ylabel("daily admissions" + ylbl)
            % Record the upper limit of the y-axis in each tile to make them consistent
            h = gca;
            yUpper(iTile) = h.YLim(2);
        end
        % Set y axis limits
        yMax = max(yUpper);
        for iTile = 1:4
            nexttile(iTile);
            ylim([0 yMax]);
            grid on
        end
        if pathogen_name(iPathogen) == "SARSCOV2"
            pathogen_title = "SARS-CoV-2";
        elseif pathogen_name(iPathogen) == "flu"
            pathogen_title = "Influenza";
        elseif pathogen_name(iPathogen) == "RSV"
            pathogen_title = "RSV";
        end
        sgtitle(pathogen_title + " hospitalisations");
        fName = pathogen_name(iPathogen) + "_hosp.png";
        saveas(h, fileNames.figureFolder + fName);


    end


    % Plot scores
    tScore = 1:par.timeHorizon;

    h = figure(100);
    nexttile(iPathogen);
    plot(tScore, nanmean(scoreCases, 1), 'LineWidth', 2)
    hold on
    plot(tScore, nanmean(scoreHosp, 1), 'LineWidth', 2)
    ylim([0 0.65])
    grid on
    xlabel('time horizon (days)')
    ylabel('mean CRPS')
    if ~all(all(isnan(scoreHosp)))
        legend('cases', 'hospitalisations', 'location', 'southeast')
    end
    title(letters(iPathogen) + " " + pathogen_title);
end

fName = "scores.png";
saveas(h, fileNames.figureFolder+fName);

