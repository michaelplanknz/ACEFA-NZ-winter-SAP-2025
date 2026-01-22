function [xOut, yOut] = getFillArgs(x, yLow, yHi)

% Return arguments for a fill graph to plot a band between yLow and yHi against x
% The outputs can be used in a 'fill' call as follows:
% >> fill(xOut, yOut, ...)

xOut = [x, fliplr(x)];
yOut = [yLow, fliplr(yHi)];


