function [S_NMf, x, y] = Compute2DFourierReconstruction(fHat, varargin)                  

%% Compute Fourier Partial Sum Approximation
% Script to compute the Fourier partial sum approximation of a function
%
% Usage:    S_NMf = ComputeFourierReconstruction(fHat)
%
%   Inputs:
%       fHat        - the Fourier coefficients (complex vector)
%
%   Optional Inputs:
%       ngrid       - number of grid points to compute reconstruction on
%                     (integer)
%       mgrid       - number of grid points to compute reconstruction on 
%
%   Outputs:
%       S_NMf        - the Fourier partial sum approximation (real vector)
%       x           - the grid on which the approximation is computed
%                     (real vector)
%



% No. of Fourier modes
NM = size(fHat);           % note: we assume there are 2N+1 total 
                            % coefficients f^_{-N},.....,f^_N
N = (NM(1)-1)/2;
M = (NM(2)-1)/2;

% Fourier modes
kn = (-N:N).';
km = (-M:M).';

% Process optional arguments
if( nargin == 3 )           % reconstruction grid length 
    ngrid = varargin{1};
    mgrid = varargin{2};
else
    ngrid = 3*(2*N+1);      % reconstruct on 3x oversampled grid
    mgrid = 3*(2*M+1);
end


% Generate equispaced grid in [-pi, pi)
x = -pi + (2*pi/ngrid)*( 0:ngrid-1 ).';
y = -pi + (2*pi/mgrid)*( 0:mgrid-1 ).';


S_NMf = zeros(length(x), length(y));

%outer product
F1 = exp(1i * x * kn.');
F2 = exp(1i * y * km.');

S_NMf = F1 * fHat * F2.';

% We will mainly consider only real functions. The following ensures tiny 
% complex values due to roundoff errors and such are not included
S_NMf = real(S_NMf);


return