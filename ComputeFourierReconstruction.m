function [S_Nf, x] = ComputeFourierReconstruction(fHat)                  

%% Compute Fourier Partial Sum Approximation
% Script to compute the Fourier partial sum approximation of a function
%
% Usage:    S_Nf = ComputeFourierReconstruction(fHat)
%
%   Inputs:
%       fHat        - the Fourier coefficients (complex vector)
%
%   Outputs:
%       S_Nf        - the Fourier partial sum approximation (real vector)
%       x           - the grid on which the approximation is computed
%                     (real vector)
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

% Generate equispaced grid in [-pi, pi)
ngrid = 2^10;
x = -pi + (2*pi/ngrid)*( 0:ngrid-1 ).';

% Compute Fourier partial sum
% Note: this is not very efficient - later, you will use FFTs to compute
% this.
S_Nf = exp( 1i*x*k.' )*fHat;

% We will mainly consider only real functions. The following ensures tiny 
% complex values due to roundoff errors and such are not included
S_Nf = real(S_Nf);


return