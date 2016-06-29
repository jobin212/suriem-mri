function [fHat, fxy ] = Get2DFourierCoefficients(FncType, N, M)
%% Generate Fourier Series Coefficients of 2D Test Function
% Script to generate the Fourier coefficients of a 2D test function 
% TODO: Circle function
% find paper
% 2D error with one jump
% construct along rays or in different ways!
%
% Usage:    [fHat, fxy] = GetFourierCoefficients(FncType, N, M)
%
%   Inputs:
%       FncType     - the type of test function to generate (string)
%                     
%       N           - Number of Fourier coefficients to compute (integer)
%                     (this generates f^_{-N},.....,f^_N)
%
%       M           - Number of Fourier coefficients to compute (integer)
%
%   Outputs:
%       fHat        - the Fourier coefficients (complex vector)
%       fxy          - function handle to construct (physical-space)
%                     function on a grid (usage: fxy(z), where z is a
%                     matrix
%                     of grid points)
%

% Initialization
FncType = lower(FncType);

% Fourier coefficients to compute
fHat = zeros(2*N+1,2*M+1);

% Compute Fourier coefficients
% Note: We will assume that the function is 2pi periodic
switch FncType

    % Smooth function
    % f = exp( sin(x) )
    % Will compute the Fourier coefficients numerically
    case ('box')
        
        kN = -N:N;
        kM = -M:M;
        
        [kx, lx] = meshgrid(kN, kM);
        fHat = (1/(pi^2)) .* sin(kx).*sin(lx) ./ (kx.*lx);
        
        fHat(kx==0) = sin(lx(kx==0)) ./ (pi^2 .* lx(kx==0));
        fHat(lx==0) = sin(kx(lx==0)) ./ (pi^2  .*kx(lx==0));
        
        fHat(M+1, N+1) = 1 / (pi^2);
               
        % Function handle
        fxy = @(x,y)    0 + ...
                        1 * ( (x >= -1) & (x <= 1) & (y >= -1 ) & (y <= 1));
    %
    case ('circle')
        kN = -N:N;
        kM = -M:M;
       
        [kx, lx] = meshgrid(kN, kM);
        
        rx = sqrt(kx.^2 + lx.^2);
        
        %fix equation!
        fHat = besselj(1, rx) ./ (2*pi * rx); 
        
        fHat(rx == 0) = 1/(2*pi);
        
        % Function handle
        fxy = @(x,y)    0 + ...
                        1 * (sqrt(x.^2 + y.^2) <= 1);


end

return


