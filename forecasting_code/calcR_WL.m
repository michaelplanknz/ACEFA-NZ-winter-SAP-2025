function R = calcR_WL(r, GTD)


% Normalise GTD
GTD = GTD/sum(GTD);

% Render GTD and x as vectors along the 3rd dimension, so that summing
% along the 3rd dimensions will return a 2D array
nx = length(GTD);
GTD = reshape(GTD, 1, 1, nx);
x = reshape( 1:length(GTD), 1, 1, nx);

% Calculate R using Wallinga-Lipsitch relation
R = 1./sum( exp(-r.*x).*GTD, 3 );


