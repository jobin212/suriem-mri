function [fHat, fxy ] = Get2DFourierCoefficients(FncType, N, M)
%% Generate Fourier Series Coefficients of 2D Test Function
% Script to generate the Fourier coefficients of a 2D test function 
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
fHat = zeros(N,M);

% Compute Fourier coefficients
% Note: We will assume that the function is 2pi periodic
switch FncType

    % Smooth function
    % f = exp( sin(x) )
    % Will compute the Fourier coefficients numerically
    case ('box')
        %{
        % equispaced grid in [-pi, pi)        
        ngrid = 2^12;
        x = -pi + (2*pi/ngrid)*(0:ngrid-1).';
        % function definition
        fx = exp( sin(x) );
        % Fourier coefficients
        fHat = ( exp(-1i*k*x.') * fx )/ngrid; 
         %}
        
        %do we need ngrid stuff?        
       %
        
        %calculate fHat values
        for kx = -N:N
            for lx = -M:M
                fHat(kx+N+1,lx+M+1) = (1/(4*pi^2)) * sin(kx)*sin(lx) / (kx*lx);
                
                if(kx == 0 && lx == 0)
                    fHat(N+1, M+1) = 1 / (pi^2);
                elseif(kx == 0)
                    fHat(kx+N+1, lx+M+1) = sin(lx) / (2 * pi^2 * lx);                             
                elseif (lx == 0)
                    fHat(kx+N+1, lx+M+1) = sin(kx) / (2 * pi^2  * kx); 
                end
                
            end
        end
        
        

        
        
        % Function handle
        fxy = @(x,y)    0 + ...
                        1 * (-1 <= x & 1 >= x & -1 <= y & y >= 1);


end

return


