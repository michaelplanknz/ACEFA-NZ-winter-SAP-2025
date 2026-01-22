function [filteredData, tData, nCases, nHosp] = getTrainingData(processed, originDate, maxTimeBack, caseIgnoreDays, hospIgnoreDays)

% Filter data to relevant dates for training
filteredData = processed(processed.date >= originDate-maxTimeBack, :);
nRows = height(filteredData);


% Extract relevant data for nCases and nHosp in the training period
indTraining = (filteredData.date <= originDate);
tData = filteredData.date(indTraining)';
nCases = (filteredData.cases(indTraining))';
if ismember('hospitalisations', fieldnames(filteredData))
    nHosp = filteredData.hospitalisations(indTraining)';
else
    filteredData.hospitalisations = nan(nRows, 1);
    nHosp = nan(size(nCases));
end

% Pad with NaNs if origin date is after most recent data
nPad = days(originDate - max(tData));
tData = [tData, max(tData)+1:max(tData)+nPad];
nCases = [nCases, nan(1, nPad)];
nHosp = [nHosp, nan(1, nPad)];

% Replaced "ignored" data with nan
nCases(tData > originDate-caseIgnoreDays) = nan;
nHosp(tData > originDate-hospIgnoreDays) = nan;

% Check dates are unique and consecutive
if ~isequal(tData, tData(1):tData(end) )
    error('data needs to have unique and consecutive dates')
end
