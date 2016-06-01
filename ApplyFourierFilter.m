 function fHat_sigma = ApplyFourierFilter(fHat, p)                  

%% Apply Fourier (Exponential, Low-Pass) Filter
% Applies a Fourier-space exponential low-pass filter
%
% Usage:    fHat_sigma = ApplyFourierFilter(fHat, sigma)
%
%   Inputs:
%       fHat        - the Fourier coefficients (complex vector)
%       p           - filter order (even integer)
%
%   Outputs:
%       fHat_sigma  - filtered Fourier coefficients (complex vector)
%

% Initialization
% Convert to a column vector
fHat = fHat(:);

% No. of Fourier modes
M = length(fHat);           % note: we assume there are 2N+1 total 
                            % coefficients f^_{-N},.....,f^_N
N = (M-1)/2;
% Fourier modes
k = (-N:N).';

% Filter function
sig = exp( -36*( k/N ).^p );

% Apply Fourier filter
fHat_sigma = fHat .* sig;


return