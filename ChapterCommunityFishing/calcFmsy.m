function Fmsy = calcFmsy(param, iSpecies)
param.bVerbose = false;
S0 = IterateSpectrum(param);

Fmsy = fminbnd(@fun, 0,0.6, optimset('MaxFunEvals',10));

  function obj = fun(F)
    param.F0(iSpecies) = F;
    S = IterateSpectrum(param,S0);
    Y = calcYield2(param,S);
    obj = -Y(iSpecies);
  end
end
