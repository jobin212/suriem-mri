function [S_Nf, x] = ComputeFourierReconstruction(fHat, varargin)                  

%% Compute Fourier Partial Sum Approximation
% Script to compute the Fourier partial sum approximation of a function
%
% Usage:    S_Nf = ComputeFourierReconstruction(fHat)
%
%   Inputs:
%       fHat        - the Fourier coefficients (complex vector)
%
%   Optional Inputs:
%       ngrid       - number of grid points to compute reconstruction on
%                     (integer)
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

% Process optional arguments
if( nargin == 2 )           % reconstruction grid length 
    ngrid = varargin{1};
else
    ngrid = 3*(2*N+1);      % reconstruct on 3x oversampled grid
end

% Generate equispaced grid in [-pi, pi)
x = -pi + (2*pi/ngrid)*( 0:ngrid-1 ).';

% Compute Fourier partial sum
% % Note: this is not very efficient - later, you will use FFTs to compute
% % this.
S_Nf = exp( 1i*x*k.' )*fHat;

%{
% To account for Matlab/FFTW's fftshift behavior
fHat = fHat .* exp( -1i*pi*k );

% Oversampled reconstruction
fk              = zeros(ngrid, 1);
fk(1:N+1)       = fHat(N+1:end);
fk(end-N+1:end) = fHat(1:N);

% use FFT
S_Nf = ngrid*ifft(fk);
%}

% We will mainly consider only real functions. The following ensures tiny 
% complex values due to roundoff errors and such are not included
S_Nf = real(S_Nf);


return