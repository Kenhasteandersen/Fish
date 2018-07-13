% ==================================================
% Set basic parameters
% ==================================================
function param = baseparameters(wInf)

A = 4.39;
fc = 0.2;

nSpecies = length(wInf);
% -------------------------------------------------
% Numerical parameters:
% -------------------------------------------------
param.nSpecies = nSpecies;   % no. of species
param.fGridexp = 0.2;        % Expansion of grid
param.nGridR = 20;          % No. of grid points for the resource
param.dt = 0.25;             % Time step
param.tEnd = 100;             % End time
param.iSave = 1;             % how often to save results (in iterations)
param.wMax = 2*max(wInf);    % Maximum size of grid
param.bVerbose = true;       % Whether to write progress to the screen while running

% -------------------------------------------------
% General parameters:
% -------------------------------------------------
param.f0 = 0.6;              % Equilibrium feeding level
kappa = 5;                % The factor determining the total biomass
                             % in the system. This factor is essentially
                             % irrelevant, except for the determining
                             % absolute concentrations or total biomass
% -------------------------------------------------
% Species-specific parameters
% -------------------------------------------------
param.wInf = shiftdim(wInf); % Theoretical asymptotic size
param.mu0prefactor = .5;      % Natural mortality is Z0prefactor*Winf^Z0exponent
param.mu0exponent = -1/4;
param.muS0 = 0.1;            % Prefactor for starvation mortality

% Growth parameters:
param.alpha = 0.6;           % Assimilation efficiency
param.n = 3/4;               % Scaling of intake
param.h = A/(param.alpha*(param.f0-fc));                % Prefactor intake
param.k = 0;                 % Activity
param.p = 0.75;              % Exponent of std. metabolism
param.ks = fc*param.alpha*param.h; % Factor of std. metabolism set 
                                    % to give a critical feeding level of 0.2
% Calculate the critical feeding level and write it on the screen:
param.fc = param.ks(1)/(param.alpha*param.h(1));
fprintf('Critical feeding level: %f\n',param.fc);
% -------------------------------------------------
% Predator-prey encounter:
% -------------------------------------------------
param.q = 0.8;               % Scaling of search volume
param.beta = 700;            % Predation/prey mass ratio
param.sigma = 1.3;           % Width of size preference
%
% Here gamma is calculated such that the initial feeding level is close to
% the feeding level specificed in the parameter f0est (Hartvig et al 2011, eq. 16):
%
lambda = 2+param.q-param.n;
alphae = sqrt(2*pi)* ...
		    param.sigma*param.beta^(lambda-2)* ...
		    exp((lambda-2)^2*param.sigma^2/2);
param.gamma = param.f0*param.h / (alphae*kappa*(1-param.f0));
param.v = (param.wInf/param.wInf(end)) .^ 0.0;  % Vulnerabilities
% -------------------------------------------------
% Reproduction:
% -------------------------------------------------
param.w0 = 0.001;           % Weigth to recruit at (egg size)
param.alphaMature = 0.25;    % Fraction of wInf to mature at
param.eRepro = .02;           % Efficiency of gonad production
% -------------------------------------------------
% Resource spectrum
% -------------------------------------------------
param.typeResource = 1;      % Type of dynamics: semi-chemostat
param.rR = 4;                % Factor for growth rate
param.lR = 0.25;            % Exponent for growth rate
param.Rmin = 0.001;         % Minimum abundance
param.kR = -2-param.q+param.n; % Exponent of carrying capacity spectrum
param.kappaR = kappa;       % Factor for carrying capacity spectrum
param.wRcut = 1;            % Cut off of recource spectrum
% -------------------------------------------------
% Recruitment:
% -------------------------------------------------
param.nRecruitmentType = 2; % 0: Physiological
                            % 1: Fixed N(w0) given in fN0
                            % 2: Beverton-Holt with limit given in Rmax
                            % 3: Fixed flux given in fFixedR
param.fRandomRecruitment = 0; % Amount of noise on recruitment
%
% Calculate the maximum recruitment for the Beverton-Holt function:
%
% 1. Calculate dW:
tmpA = param.wInf(1);
tmpB = (log10(param.wInf(end))-log10(param.wInf(1)))/(param.nSpecies-1);
dW = tmpB*tmpA*10.^(tmpB*((1:param.nSpecies)-1))';
% 2. Calculate the physiological mortality parameter a:
alphap = param.f0 * param.h(1)*param.beta^(2*param.n-param.q-1) * ...
	 exp((2*param.n*(param.q-1)-param.q^2+1)*param.sigma^2/2);
a = alphap / (param.alpha(1)*param.h(1)*param.f0-param.ks)
% 3. Calculate the expected Rmax offset by a factor (the "kappa"):
param.Rmax = 1000*50*param.kappaR*(param.alpha*param.f0*param.h*param.w0.^param.n-param.ks*...
    param.w0.^param.p) * param.wInf.^(param.n*2-param.q-3+a).*dW;
% -------------------------------------------------
% Initial spectrum is just a powerlaw N propto w^(-n-a): 
% -------------------------------------------------
w = makegrid(param);
param.Ninit = repmat(w.^(-param.n-a), param.nSpecies,1);
for i = 1:nSpecies
  param.Ninit(i,w>param.alphaMature*param.wInf(i)) = 0;
  if param.nSpecies==1
      param.Ninit = param.Ninit*param.kappaR*w(1).^(-2-param.q+param.n)/param.Ninit(1);
  else
      param.Ninit(i,:) = param.Ninit(i,:)/param.Ninit(i,1)*param.Rmax(i)/...
          (param.alpha*param.f0*param.h*param.w0.^param.n-param.ks*...
          param.w0.^param.p);
  end
end
% -------------------------------------------------
% Fishing is  specified as a "trawl" selectivity pattern through
% F0 and wFstart. An alternive is to define a function which returns the
% fishing mortality as a function of (param,iSpecies,w).
% -------------------------------------------------
param.bFishingEcosystem = false;
param.F0 = 0 * ones(nSpecies,1);   % Overall fishing pressure on each species
param.wFstart = 0.05.*param.wInf;  % Start of fishing on each species.
% Example of a function which specifices a trawl fishing pattern:
% param.funcFishing = @(param,iSpecies,w) 0*w(w>0.05*param.wInf(iSpecies))