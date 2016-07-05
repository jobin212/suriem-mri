clear;

%%%the amount of exponenents we want to take 
runs = 4;
ErrType = '2norm';
FncType = 'circle';
ReconstructionType = 'circle-true-jumps';


%fourier coefficients
k = 50*2.^(0:runs-1);

error_vector = zeros(size(k));

for i = 1:length(k)
    N = k(i);
    M = N;
    
    [fHat, fxy] = Get2DFourierCoefficients(FncType, N, M);

    [S_NMf, x, y] = Compute2DFourierReconstruction(fHat, ReconstructionType);
    
    [xx, yy] = meshgrid(x , y);
    
    abs_error = abs(fxy(xx,yy) - S_NMf);
    
    
    
    error_vector(i) = GetError(ErrType, abs_error(:), x);
    
end

loglog(k, error_vector)
hold on;
loglog(k, k.^(-1))
