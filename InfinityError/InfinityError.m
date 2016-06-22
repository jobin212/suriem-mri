function [e] = InfinityError(N) 

[fHat, fx] = GetFourierCoefficients('piecewise', N);
[S_Nf, x] = ComputeFourierReconstruction(fHat);
error = abs(fx(x) - S_Nf);

e = max(error);

%figure;
%plot(x, S_Nf);
figure;
plot(x, S_Nf);
return
