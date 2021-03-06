clear;

%%%the amount of exponenents we want to take 
exp = 3;
ErrType = '2norm';
FncType = 'circle';
ReconstructionType = 'true-jumps';


%fourier coefficients
k = 50*2.^(0:exp);

error_vector = zeros(size(k));

for i = 1:length(k)
    N = k(i);
    M = N;
    
    [fHat, fxy] = Get2DFourierCoefficients(FncType, N, M);

    [S_NMf, x, y] = Compute2DFourierReconstruction(fHat, ReconstructionType);
    
    [xx, yy] = meshgrid(x , y);
    
    
    
    
    abs_error = abs(fxy(xx,yy) - S_NMf);
    
    cross_section = abs(abs_error(N+1, :));
    
    error_vector(i) = GetError(ErrType, cross_section, x);
    
end



loglog(k, error_vector);


