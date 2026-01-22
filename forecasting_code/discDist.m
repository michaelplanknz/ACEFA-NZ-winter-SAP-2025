function f = discDist(pdfFunc, a, b )

% Calculate a discrete probability mass function from a continuous PDF on
% non-negative values using the Method of Cori et al 2013 - see also
% Thompson et al 2024
% The discrete PMF will be truncated (at the maximum value b) and
% normalised to sum to 1
% If the starting value a is 1 instead of 0, the first element (probability of 1) will
% include all mass between 0 and 1, plus a share of mass between 1 and 2

nn = a:b;

nValues = length(nn);
f = zeros(1, nValues);


if ~ismember(a, [0 1])
    error('Lower limit for PMF discretisation (a) needs to be either 0 or 1')
end

for iValue = 1:nValues
    % If the starting value a is 1 instead of 0, set the first element
    % (probability of 1) to include all mass between 0 and 1, plus a share of mass between 1 and 2
    if iValue == 1 & a == 1
        ig1 = @(x)(pdfFunc(x) );
    else
        ig1 = @(x)(pdfFunc(x).*(1-nn(iValue)+x ) );
    end
    ig2 = @(x)(pdfFunc(x).*(1+nn(iValue)-x ) );
    f(iValue) = quad(ig1, nn(iValue)-1, nn(iValue)) + quad(ig2, nn(iValue), nn(iValue)+1);
end
f = f/sum(f);




