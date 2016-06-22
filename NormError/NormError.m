function [e] = NormError(N) 

[fHat, fx] = GetFourierCoefficients('piecewise', N);
[S_Nf, x] = ComputeFourierReconstruction(fHat);
error = (fx(x) - S_Nf);

h = x(2) - x(1);

e = sqrt(h) *norm(error);

return



