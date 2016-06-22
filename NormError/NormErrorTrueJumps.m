function [e] = NormErrorTrueJumps(N) 

true_jump_locs = [-3*pi/4, -pi/2, -pi/4, pi/8, 3*pi/8, 3*pi/4];
true_jump_heights = [3/2, -3/2, (7/4) - ((-pi/4)/2) + sin((-pi/4)- (1/4)),-((7/4) - ((pi/8)/2) + sin((pi/8)- (1/4))), ((33/32)*pi -5), -((33/16)*pi-5)];

[fHat, fx] = GetFourierCoefficients('piecewise', N);
[~, x] = ComputeFourierReconstruction(fHat);
S_Nf_edge = EdgeEnhancedReconstruction(fHat, true_jump_heights, true_jump_locs);


error = (fx(x) - S_Nf_edge);

h = x(2) - x(1);

e = sqrt(h) *norm(error);

return

