%%%the amount of exponenents we want to take 
exp = 7;
ErrType = 'infinity';
FncType = 'box';


%fourier coefficients
k = 50*2.^(0:exp);

[fHat, fxy] = Get2DFourierCoefficients(FncType, N, M);

[S_NMf, x, y] = Compute2DFourierReconstruction(fHat, ReconstructionType);

error_vector = zeros(size(k));

for i = 1:length(k)
    error_vector(i) = Get2DError(k(i));
    
end
