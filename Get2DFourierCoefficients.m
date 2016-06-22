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
fHat = zeros(2*N+1,2*M+1);

% Compute Fourier coefficients
% Note: We will assume that the function is 2pi periodic
switch FncType

    % Smooth function
    % f = exp( sin(x) )
    % Will compute the Fourier coefficients numerically
    case ('box')
        
        %calculate fHat values
        %{
        for kx = -N:N
            for lx = -M:M
                fHat(kx+N+1,lx+M+1) = (1/(pi^2)) * sin(kx)*sin(lx) / (kx*lx);
                
                if(kx == 0 && lx == 0)
                    fHat(N+1, M+1) = 1 / (pi^2);
                elseif(kx == 0)
                    fHat(kx+N+1, lx+M+1) = sin(lx) / (pi^2 * lx);                             
                elseif (lx == 0)
                    fHat(kx+N+1, lx+M+1) = sin(kx) / (pi^2  * kx); 
                end
                
            end
        end
        %}
        
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


end

return


