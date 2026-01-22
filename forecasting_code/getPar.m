function par = getPar(location_name, pathogen_name)

% number of particles to simulate - usually 1e5 for full forecast
par.nParticles = 1e5;      

% fixed resampling lag (days)
par.resampleLag = 42;       



par.maxTimeBack = 250;              % time to start model runs (days before origin date)
par.timeHorizon = 28;               % time to run forward (days after origin date)

% Specify model type for Rt
par.Rt_method = "GP";               % Set to "RW" for random walk or "GP" for Gaussian process on Rt or "GP_diff" for Gaussian process on daily change in Rt 

% RW parameters for Rt
par.sigmaR = 0.015;                % 0.015 S.D. in random walk step for log Rt
par.Beta = 0;                       % autoregression coefficient for Rt (if 

% GP parameters for Rt
par.GP_sd = 0.1;                        % GP s.d. for "GP" method (s.d. in Rt)
par.GP_diff_sd = 0.01;                  % GP s.d. for "GP_diff" method (s.d. in daily change in Rt)
par.GP_autoCorrTime = 30;
par.GP_covFn = @(t, SD, L)covFnMatern(days(t), SD, L, 5/2);
par.sNoise = 0.001;                                    % "observation" noise level
                                                                            

par.pHospCV = 0.025;                % CV in initial condition for pHosp
par.pHospInitDays = 21;             % Initialise pHosp with a mean equal to the ratio of admissions to cases over this number of days at the start of the run
par.sigmaP = 0.01;                  % S.D. in random walk step for log(pHosp)  
    
% Get shape and scape parameters for the intial value of Rt (when NOT using the "GP" method)
%[par.R_shape, par.R_scale] = gamShapeScale(1, 1 );      % Inputs to gamShapeScale are the mean and s.d.
par.R0_mean = 1;
par.R0_sd = 0.2;

par.pReport = 1;          % case reporting probability
par.kObs = 25;             % 25 overdispersion parameter for daily observed cases (set to inf for a Poisson distribution)
par.kHosp = 25;             % 25 overdispersion parameter for daily hospital admissions (set to inf for a Poisson distribution)

% Pathogen-specific generation time distribution parameters

if pathogen_name == "SARSCOV2"
    % Number of days to ignore recent case/hospitalisation data 
    % NB This is relative to origin date, so if the data is truncated before being downloaded, this will not truncate it further 
    par.caseIgnoreDays = 0;
    par.hospIgnoreDays = 5;

    GTmax = 15; 
    % Mean 3.3, s.d. 3.5
    GTD_shape = 0.89;
    GTD_scale = 1/0.27;
    incMean = 4.4;
    incSD = 2.5;
elseif pathogen_name == "flu"
    % Number of days to ignore recent case/hospitalisation data 
    % NB in NZ "case data" is actually hospitalisations for flu/RSV 
    par.caseIgnoreDays = 0;
    par.hospIgnoreDays = 0;

    GTmax = 15;        
    % Mean 2.6, s.d. 1.3
    GTD_shape = 4;
    GTD_scale = 1/1.54;
    incMean = 1.7;
    incSD = 1.0;
elseif pathogen_name == "RSV"
    % Number of days to ignore recent case/hospitalisation data 
    % NB in NZ "case data" is actually hospitalisations for flu/RSV 
    par.caseIgnoreDays = 0;
    par.hospIgnoreDays = 0;
    
    GTmax = 15; 
    % Mean 7.5, s.d. 2.1
    GTD_shape = 12.76;
    GTD_scale = 1/1.7;
    incMean = 5.0;
    incSD = 2.0;
end


% Infection to report time assumed to be a gamma distribution whose mean and
% variance are the sum of the incubation period and onset to report mean and
% variance respectively 
RTmax = 25;
onsetToReportMean = 1.9;
onsetToReportSD = 1.8;
RTmean = incMean + onsetToReportMean;       
RTsd = sqrt(incSD^2 + onsetToReportSD^2);  

% Infection to admission time assumed to be a gamma distribution defined similarly by incubation period, onset to report and report to admission times
HTmax = 25;
reportToAdmitMean = 0;
reportToAdmitSD = 2.1;
HTmean = RTmean + reportToAdmitMean;        
HTsd = sqrt(RTsd^2 + reportToAdmitSD^2);


% Form generation time PMF vector
pdfFnGTD = @(x)(gampdf(x, GTD_shape, GTD_scale ));
par.GTD = discDist( pdfFnGTD, 1, GTmax ); % Probability mass of discretised GT distribution on integers 1, 2, ...


% Form report time PMF vector
[RTD_shape, RTD_scale] = gamShapeScale(RTmean, RTsd);       
pdfFnRTD = @(x)(gampdf(x, RTD_shape, RTD_scale));
par.RTD = discDist( pdfFnRTD, 0, RTmax);

% Form hospitalisation time PMF vector
[HTD_shape, HTD_scale] = gamShapeScale(HTmean, HTsd);       
pdfFnHTD = @(x)(gampdf(x, HTD_shape, HTD_scale));
par.HTD = discDist( pdfFnHTD, 0, HTmax);


