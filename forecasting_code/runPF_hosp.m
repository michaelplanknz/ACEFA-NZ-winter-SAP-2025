function particles = runPF_hosp(t, nCases, nHosp, par)

% Function to run the particle filter on given input data
%
% USAGE: particles = runPF(t, nCasesLoc, nInfImp, par)
%
% INPUTS: t - vector of daily dates defining the model simulation period
%         nCases - corresponding vector of daily cases
%         par - structure of model parameters with fields
%         - par.nParticles - number of particles to use
%         - par.GTD - vector of probability masses for generation time t = 1, 2, .. days
%         - par.RTFD - vector of probability masses for infection to reporting time t = 0, 1, 2, .. days
%         - par.R_shape, par.R_scale - shape and scale parameters for the
%         prior distribution for initial Rt
%         - par.deltat, par.sigmat - vectors (of same length as t) of means and s.d. for the daily change in Rt 
%         - par.obsModel - noise model for observed daily cases data, either
%         "bin" or "negbin" 
%         - par.kObs (if obsModel is "negbin") - dispersion parameter for
%         observed daily cases (inf for Poisson dist)
%         - par.pReport - probability of infecitons being reported as cases
%         - par.resampleLag - fixed lag (number of time steps) for bootstrap particle filter resampling
%
% OUTPUTS: particles - a structure with the following fields
%          t - vector of dates corresponding to the columns of the
%          following variables
%          Rt - matrix of reproduction numbers - (i,j) element corresponds to particle i on day j
%          It - matrix of daily infections  - (i,j) element corresponds to particle i on day j
%          Zt - matrix of infections by assigned date of report (independent of whether they actually reported as a case or not)  - (i,j) element corresponds to particle i on day j
%          Ct - matrix of simulated daily cases  - (i,j) element corresponds to particle i on day j
%          LL - particle marginal log likelihood - can be used in PMMH to fit fixed parameters

% Fit day of the week multipliers
DOWEffect_cases = fitDOWEffect(t, nCases);
DOWEffect_hosp = fitDOWEffect(t, nHosp);

% Extend time vector by specified time horizon
tExt = t(1):t(end) + par.timeHorizon;  
nSteps = length(tExt);
iLastData = length(nCases);

% Get day of the week indices for date in tExt
DOW = weekday(tExt);

% Estimate hospitalisation probability as a crude ratio for initial
% condition
pHospEst = sum(nHosp(1:par.pHospInitDays))/sum(nCases(1:par.pHospInitDays));

tInit = max([length(par.GTD), length(par.RTD), length(par.HTD)]);        % initialisation period (when not running renewal equation)
                              

% Initialise arrays for renewal equation particle filter
Rt = zeros(par.nParticles, nSteps);
It = zeros(par.nParticles, nSteps);
Zt = zeros(par.nParticles, nSteps);
Pt = zeros(par.nParticles, nSteps);
Ht = zeros(par.nParticles, nSteps);
if par.Rt_method == "GP_diff"
    dR = zeros(par.nParticles, nSteps);
end

% Initialise array for log mean weights 
lmw = zeros(1, nSteps);


% Draw inital value of Rt from prior (lognormal with same mean and var as
% GP if using the "GP" method, or with specified mu and sigma otherwise)
if par.Rt_method == "GP"
    Rt(:, tInit) = lognrnd(0, par.GP_sd, par.nParticles, 1);
else
    ln_Mu = log(par.R0_mean^2 / sqrt(par.R0_sd^2 + par.R0_mean^2));
    ln_Sigma = sqrt(log(par.R0_sd^2/par.R0_mean^2 + 1));
    Rt(:, tInit) = lognrnd(ln_Mu, ln_Sigma, par.nParticles, 1);
end

% Assumne the same value for t=tInit-1 (required for the one-step
% autogregression model)
Rt(:, tInit-1) = Rt(:, tInit);




% Draw initial value of Pt from prior
[sh, sc] = gamShapeScale(pHospEst, pHospEst*par.pHospCV);
Pt(:, tInit) = gamrnd(sh, sc, par.nParticles, 1);
    
% Initialise infections in the burn-in period
nCasesSmoothed = smoothdata(nCases, 'movmean', 7);
RTmean = round(sum( par.RTD .* (0:length(par.RTD)-1) ));
It(:, 1:tInit) = poissrnd( repmat(nCasesSmoothed(1+RTmean:tInit+RTmean)/par.pReport, par.nParticles, 1) ); % offset by RT mean and inflate by factor of pReport to account approximately for reporting lag and under-reporting



% Loop through time steps
for iStep = tInit+1:nSteps

   % Draw current day's value for reproduction number Rt 
   if par.Rt_method == "RW"
       Rt(:, iStep) = Rt(:, iStep-1) .* (Rt(:, iStep-1)./Rt(:, iStep-2)).^par.Beta .* lognrnd(0, par.sigmaR, par.nParticles, 1);
   elseif par.Rt_method == "GP"
       Rt(:, iStep) = exp( GPstep(tExt(tInit:iStep-1), log(Rt(:, tInit:iStep-1)), tExt(iStep), par.GP_covFn, par.GP_sd, par.GP_autoCorrTime, par.sNoise   )  );
   elseif par.Rt_method == "GP_diff"
       dR(:, iStep) = GPstep(tExt(tInit:iStep-1), dR(:, tInit:iStep-1), tExt(iStep), par.GP_covFn, par.GP_diff_sd, par.GP_autoCorrTime, par.sNoise  );
       Rt(:, iStep) = Rt(:, iStep-1).*exp(dR(:, iStep));
   else
        error("par.Rt_method must be eiter RW (random walk) or GP (Gaussian process)");
   end

   % Draw current day's value for hospitalisation probability Pt
   Pt(:, iStep) = Pt(:, iStep-1) .* lognrnd(0, par.sigmaP, par.nParticles, 1);

   % Compute It according to renewal equation
   ind = iStep-1:-1:iStep-length(par.GTD);
   It(:, iStep) = poissrnd(  Rt(:, iStep) .* sum(par.GTD.*It(:, ind), 2 ) );     

    % Deterministic convolution for number of infections by notification date (Zt) for efficiency
   ind = iStep:-1:iStep-length(par.RTD)+1;
   Zt(:, iStep) = sum(par.RTD.*It(:, ind), 2);

   % Deterministic convolution for number of infections by admission date (Ht) for efficiency
   ind = iStep:-1:iStep-length(par.HTD)+1;
   Ht(:, iStep) = sum(par.HTD.*It(:, ind), 2);

   % Calculate weights for particle resampling (only during period for which data is available)
   if iStep <= iLastData
       if ~isnan(nCases(iStep))
           meanCases = par.pReport * DOWEffect_cases(DOW(iStep)) * Zt(:, iStep);
            if isfinite(par.kObs)   
               weights_cases = nbinpdf(nCases(iStep), par.kObs, par.kObs./(meanCases+par.kObs));
            else
               weights_cases = poisspdf(nCases(iStep), meanCases );
            end
       else
          weights_cases = ones(par.nParticles, 1);
       end   
        if ~isnan(nHosp(iStep))
            meanHosp = par.pReport*DOWEffect_hosp(DOW(iStep))*Pt(:, iStep).*Ht(:, iStep);
            if isfinite(par.kHosp)
               weights_hosp = nbinpdf(nHosp(iStep), par.kObs, par.kObs./(meanHosp+par.kObs));
            else
               weights_hosp = poisspdf(nHosp(iStep), meanHosp);
            end
        else
            weights_hosp = ones(par.nParticles, 1);
        end
        weights = weights_cases .* weights_hosp;
        lmw(iStep) = log(mean(weights));

        % Resample particles according to weights
        if any(isnan(weights))
            pause
        end
        resampInd = randsample(par.nParticles, par.nParticles, true, weights);

        iResample = max(1, iStep-par.resampleLag);
        Rt(:, iResample:end) = Rt(resampInd, iResample:end);
        It(:, iResample:end) = It(resampInd, iResample:end);
        Zt(:, iResample:end) = Zt(resampInd, iResample:end);
        Pt(:, iResample:end) = Pt(resampInd, iResample:end);
        Ht(:, iResample:end) = Ht(resampInd, iResample:end);
   end

end


 % Generate samples from the reported case and admissions distributions to construct
 % prediction intervals.
 % Separately output values with and without the DOW effect
meanCases_noDOW = par.pReport*Zt;
meanCases_withDOW = DOWEffect_cases(DOW).*meanCases_noDOW;
if isfinite(par.kObs)  
    Ct = nbinrnd(par.kObs, par.kObs./(meanCases_withDOW+par.kObs)); 
    Cndw = nbinrnd(par.kObs, par.kObs./(meanCases_noDOW+par.kObs)); 
else  
    Ct = poissrnd(meanCases_withDOW);
    Cndw = poissrnd(meanCases_noDOW);
end
meanHosp_noDOW = par.pReport*Pt.*Zt;
meanHosp_withDOW = DOWEffect_hosp(DOW).*meanHosp_noDOW;
if isfinite(par.kObs)  
    At = nbinrnd(par.kObs, par.kObs./(meanHosp_withDOW+par.kObs)); 
    Andw = nbinrnd(par.kObs, par.kObs./(meanHosp_noDOW+par.kObs)); 
else  
    At = poissrnd(meanHosp_withDOW);
    Andw = poissrnd(meanHosp_noDOW);
end



% Store outputs in a structure called particles, discsrding initialisation period plus one additional week of burn in:
particles.t = tExt(tInit+1:end);
particles.Rt = Rt(:, tInit+1:end);
particles.It = It(:, tInit+1:end);
particles.Zt = Zt(:, tInit+1:end);
particles.Ct = Ct(:, tInit+1:end);
particles.Cndw = Cndw(:, tInit+1:end);
particles.Pt = Pt(:, tInit+1:end);
particles.Ht = Ht(:, tInit+1:end);
particles.At = At(:, tInit+1:end);
particles.Andw = Andw(:, tInit+1:end);

particles.LL = sum(lmw);


