%%%the amount of exponenents we want to take 
exp = 7;
ErrType = '2norm';
FncType = 'box';
ReconstructionType = 'standard';


%fourier coefficients
k = 50*2.^(0:exp);

error_vector = zeros(size(k));

for i = 1:length(k)
    [fHat, fxy] = Get2DFourierCoefficients(FncType, k(i), k(i));

    [S_NMf, x, y] = Compute2DFourierReconstruction(fHat, ReconstructionType);
    
    [xx, yy] = meshgrid(x , y);
    
    abs_error = abs(fxy(xx,yy) - S_NMf);
    
    error_vector(i) = Get2DError(ErrType, abs_error, x, y);
    
end

loglog(k ,error_vector);
