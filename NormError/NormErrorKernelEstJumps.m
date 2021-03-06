function [e] = NormErrorKernelEstJumps(N) 

[fHat, fx] = GetFourierCoefficients('piecewise', N);
[~, x] = ComputeFourierReconstruction(fHat);

[kernel_est_jump_heights, kernel_est_jump_locs] = KernelEstPiecewise(N);

S_Nf_edge = EdgeEnhancedReconstruction(fHat, kernel_est_jump_heights, x(kernel_est_jump_locs));

figure;
plot(x, fx(x), x, S_Nf_edge);


error = (fx(x) - S_Nf_edge);

h = x(2) - x(1);

e = sqrt(h) *norm(error);

return