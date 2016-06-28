%get coefficients
FncType = 'box';
N = 20;
sN = num2str(N);
M = 20;
sM = num2str(M);
leg = strcat('N=M=',sN); 

ErrType = 'absolute';
ReconstructionType = lower('standard');

[fHat, fxy] = Get2DFourierCoefficients(FncType, N, M);

[S_NMf, x, y] = Compute2DFourierReconstruction(fHat, ReconstructionType);

[error, xx, yy] = Get2DError(FncType, ErrType, S_NMf, x, y);

mesh(xx, yy, S_NMf);
legend(leg);

figure;
mesh(xx, yy, error);
legend(leg);




