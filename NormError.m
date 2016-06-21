function [e] = NormError(N)

%% returns the norm error 


[fHat, fx] = GetFourierCoefficients('piecewise', N)
[S_Nf, x] = ComputeFourierReconstruction(fHat)

error = (fx(x) - S_Nf)
h = (2*pi)/1024
e = h*norm(error)

return



