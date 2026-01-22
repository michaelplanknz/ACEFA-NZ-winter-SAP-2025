function particles = runOneForecast(fileNames, useHospAsCases, fileDate, location_name, pathogen_name, test_types, originDate, par)

% Code to read data for a specified location and pathogen and specified
% file date, and run the forecast model for a specified origin date

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Read epi data
processed = readData(fileNames, useHospAsCases, fileDate, location_name, pathogen_name, test_types);

if originDate > max(processed.date)
    fprintf('   WARNING: origin date is after the most recent data for this location-pathogen which is on %s\n', max(processed.date(processed.location == location_name & processed.pathogen == pathogen_name)));
end

% Extract training data (filtered to relevant dates) to feed to model
[~, tData, nCases, nHosp] = getTrainingData(processed, originDate, par.maxTimeBack, par.caseIgnoreDays, par.hospIgnoreDays);


% Uncomment this code to produce a temporary plot of the data that is being
% fed to the particle filter
% h = figure(100);
% h.Position = [     38         250        1103         695];
% tiledlayout(2, 1, "TileSpacing", "compact");
% nexttile;
% plot(processed.date, processed.cases, 'k.', tData, nCases, 'r.-' )
% xline(originDate, 'k--')
% xline(fileDate, 'g--')
% h = gca;
% xl = h.XLim;
% xl(1) = min(tData)-14;
% h.XLim = xl;
% ylabel('cases')
% grid on
% nexttile;
% if ismember('hospitalisations', fieldnames(processed))
%     plot(processed.date, processed.hospitalisations, 'k.', tData, nHosp, 'r.-' )
%     xline(originDate, 'k--')
%     xline(fileDate, 'g--')
%     h = gca;
%     xl = h.XLim;
%     xl(1) = min(tData)-14;
%     h.XLim = xl;
%     ylabel('admissions')
%     grid on
% end
% sgtitle(sprintf('%s/%s data being used for forecast, file dates %s', location_name, pathogen_name, fileDate))
% drawnow


% Run particle filter
particles = runPF_hosp(tData, nCases, nHosp, par);


%close(100);

