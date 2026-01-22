function [t, qts, samples] = getParticleQuantiles(particles, tMin, qt, iSav)

% Get quantiles and samples from particles for saving 
% iSav is the index of the particles to be saved

% Truncate to dates after tMin only (if specified)
if ~isempty(tMin)
    ind = particles.t >= tMin;
else
    ind = true(1, length(particles.t));
end

t = particles.t(ind);

% Calculate quantiles of each variable at levels specified by qt and store
% in qts
qts.Rt = quantile(particles.Rt(:, ind), qt);
qts.Ct = quantile(particles.Ct(:, ind), qt);
qts.Pt = quantile(particles.Pt(:, ind), qt);
qts.At = quantile(particles.At(:, ind), qt);
qts.Ct_smoothed = quantile(smoothdata(particles.Ct(:, ind), 2, 'movmean', 7), qt);
qts.Cndw = quantile(particles.Cndw(:, ind), qt);
qts.Andw = quantile(particles.Andw(:, ind), qt);



% Save outputs from the first N particles to enable scoring
samples.Ct = particles.Ct(iSav, ind);
samples.Cndw = particles.Ct(iSav, ind);
samples.At = particles.Ct(iSav, ind);
samples.Andw = particles.Ct(iSav, ind);

