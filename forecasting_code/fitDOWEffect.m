function DOWEffect = fitDOWEffect(t, nCases)

maxWeeks = 16;

% Index of last available data
iLast = find(~isnan(nCases), 1, 'last');

% Keep an integer number of weeks (up to specified maximum)
nToKeep = min(7*floor(iLast/7), 7*maxWeeks);

indKeep = iLast-nToKeep+1:iLast;
t = t(indKeep);
nCases = nCases(indKeep);

DOW = weekday(t);

sumCases = zeros(1, 7);
for iDay = 1:7
    sumCases(iDay) = sum(nCases(DOW == iDay));
end


DOWEffect = 7*sumCases/sum(nCases);
