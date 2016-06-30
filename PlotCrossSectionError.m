clear;

%get coefficients%get coefficients
FncType = 'circle';
N = 100;
sN = num2str(N);
M = N;
sM = num2str(M);
leg = strcat('N=M=',sN); 

ReconstructionType = lower('standard');

[fHat, fxy] = Get2DFourierCoefficients(FncType, N, M);

%overlay cross sections to compare error

[S_NMf, x, y] = Compute2DFourierReconstruction(fHat, ReconstructionType);

[xx, yy] = meshgrid(x,y);

abs_error = abs(S_NMf - fxy(xx,yy));



%{

%%%the amount of exponenents we want to take 
exp = 5;
ErrType = '2norm';
FncType = 'box';
ReconstructionType = 'standard';


%fourier coefficients
k = 50*2.^(0:exp);

error_vector = zeros(size(k));

for i = 1:length(k)
    N = k(i);
    M = N;
    
    [fHat, fxy] = Get2DFourierCoefficients(FncType, N, M);

    [S_NMf, x, y] = Compute2DFourierReconstruction(fHat, ReconstructionType, 2*N+1, 2*M+1);
    
    [xx, yy] = meshgrid(x , y);
    
    abs_error = abs(fxy(xx,yy) - S_NMf);
    
    
    
    error_vector(i) = GetError(ErrType, abs_error(:), x);
    
end

plot(log2(k),log2(error_vector));

%}
