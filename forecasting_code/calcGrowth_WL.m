function r = calcGrowth_WL(R, GTD)

nPoints = 1000;
pad = 1e-5;

x = 1:length(GTD);

% Normalise GTD and calculate its mean
GTD = GTD/sum(GTD);
GTD_mean = sum(x.*GTD);


% Calculate minimum and maximum valuues of r by inverting the
% Wallinga-Lipsitch relation
opts = optimset('Display', 'off');
myFn = @(y, R)(1/R - sum( exp(-y*x).*GTD  ) );

Rmin = min(min(R));
y0 = (Rmin-1)/GTD_mean;     % initial estimate for r based on an exponential GTD
rMin = fsolve(@(y)myFn(y, Rmin), y0, opts);

Rmax = max(max(R));
y0 = (Rmax-1)/GTD_mean;       % initial estimate for r based on an exponential GTD
rMax = fsolve(@(y)myFn(y, Rmax), y0, opts);

% Create an array of values for r and apply the WL relation to get
% corresponding values of R
ra = linspace(rMin-pad, rMax+pad, nPoints);
Ra = calcR_WL(ra, GTD);

% Invert via linear interpolation to find r 
r = interp1(Ra, ra, R);

