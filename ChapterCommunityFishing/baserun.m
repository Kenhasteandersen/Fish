addpath('~/Source/sizespectrum/svn/devel/')

param = baseparameters(logspace(1,6,10));
param.F0 = 0.1 * ones(1,param.nSpecies);
S = IterateSpectrum(param);
PlotResult(param,S)

%%
% Check Fmsy:
%
Fmsy = zeros(1,param.nSpecies);
param.tEnd = 100;
parfor i = 1:param.nSpecies
  Fmsy(i) = calcFmsy(param,i)
end

clf
semilogx(param.wInf, Fmsy, param.wInf, 2*param.wInf.^(-0.25))
%%
% Check cascade:
%
param.F0 = 0.15 * ones(1,param.nSpecies);
param.wFstart = 5000 * ones(1,param.nSpecies);
SF = IterateSpectrum(param,S);
PlotTwo(param,S,SF)
%%
% Kappa sensitivity
%
Kappa = logspace(-2,2,10);

clf
w = S.w;
W = param.wInf;
param.tEnd = 100;
param.bVerbose = false;

parfor i = 1:length(Kappa)
  par = param;
  par.Rmax = Kappa(i) * param.Rmax;
  SS(i) = IterateSpectrum(par, S);
end

for i = 1:length(Kappa)
  ix = SS(i).nSave*0.5:SS(i).nSave;
  subplot(3,1,1)
  loglog(w, mean(SS(i).Ntot(ix,:)).*w.^2)
  hold on
  
  subplot(3,1,2)
  loglog(W, mean(SS(i).Rp(ix,:)./SS(i).R(ix,:)))
  hold on
  
  subplot(3,1,3)
  semilogx(w, squeeze(mean(SS(i).f(ix,1,:))))
  hold on
end
