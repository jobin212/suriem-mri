function [e] = NormErrorPronyEstJumps(N) 

%%estimated jmp info
[est_jump_vals,est_jump_locs] = PronyMethod(N,9,20);

%%same method now as NormErrorTrueJumps
[fHat, fx] = GetFourierCoefficients('piecewise', N);
[~, x] = ComputeFourierReconstruction(fHat);
S_Nf_edge = EdgeEnhancedReconstruction(fHat, real(est_jump_vals), est_jump_locs);

%error = abs(fx(x) - S_Nf_edge);
error = (fx(x) - S_Nf_edge);
h = x(2) - x(1);

%e = sqrt(h)*norm(max(error));
e = sqrt(h)*norm(error);


return
