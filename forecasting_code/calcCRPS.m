function score = calcCRPS(xs, y, logTransFlag)

% Compute the CRPS from a set of samples from a distribution (e.g. forecast) for a given data point xData
%
% USAGE: score = calcCRPS(xSample, xData)
%
% INPUTS: xs - n x m matrix of samples of x, with each row representing a
% different distribution (e.g. different forecast output or time point) and
% each column representing a different sample
%         y - n x 1 array of data points against which the CRPS will be
%         computed for each set of samples
%
% OUTPUTS: score - n x 1 array of scores corresponding to the n data points and forecast outputs

% Remove any rows where xData or a value of xSample is nan (values of score for
% these rows will be returned as nan)
[nRows, nPoints] = size(xs);
nanFlag = sum(isnan(xs), 2) > 0 | isnan(y);
xs(nanFlag, :) = [];
y(nanFlag) = [];

% Number of rows after removal of rows containing NaN
[ny, ~] = size(xs);

% Apply log transform if applicable
if logTransFlag
    if any(xs(:) < 0) || any(y < 0)
        error('Cannot apply log transform to negative forecast samples or negative data')
    end
    xs = log(1+xs);
    y = log(1+y);
end

% Ensure each row of samples is sorted in ascending order
xs = sort(xs, 2);                             

% index counter corrresponding to the column number of xs (dividing this by nPoints gives the corresponding values of the ECDF in the intevral immediately to the right of each x value)
iCol = repmat(1:nPoints, ny, 1 );          

% Flag indicating whether y is great than each vlaue of x in the array
y_above_x_flag = y > xs;

% kCrit is the index of the largest value of x for which y>x
% Either kCrit = 0 if y is less than all x values, kCrit = n if y is
% greater than all x values , or otherwise kCrit is such that x_k < y <
% x_k+1
kCrit = max(iCol.*y_above_x_flag, [], 2);

% Value of the two summands in the equation for CRPS
summand1 = (1-2*iCol).*xs;
summand2 = (2*(nPoints-iCol)+1).*xs;

% Create a combined summand selecting summand1 for values of x < y
% and summand2 for values of x > y
summand_combined = y_above_x_flag.*summand1 + (1-y_above_x_flag).*summand2;

% Any scores where xData or a value of xs is nan are returned as nan
score = nan(nRows, 1);

% Evaluate the equation for CRPS for the valid rows
score(~nanFlag) = 1/nPoints^2 * sum(summand_combined, 2) + (2*kCrit/nPoints-1).*y;



