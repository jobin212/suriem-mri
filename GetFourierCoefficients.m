function [fHat, fx] = GetFourierCoefficients(FncType, N)                                       

%% Generate Fourier Series Coefficients of 1D Test Function
% Script to generate the Fourier coefficients of a 1D test function 
%
% Usage:    fHat = GetFourierCoefficients(FncType, N)
%
%   Inputs:
%       FncType     - the type of test function to generate (string)
%                     ( 'smooth', 'piecewise' )
%       N           - Number of Fourier coefficients to compute (integer)
%                     (this generates f^_{-N},.....,f^_N)
%
%   Outputs:
%       fHat        - the Fourier coefficients (complex vector)
%       fx          - function handle to construct (physical-space)
%                     function on a grid (usage: fx(y), where y is a vector
%                     of grid points)
%

% Initialization
FncType = lower(FncType);

% Fourier coefficients to compute
k = (-N:N).';

% Compute Fourier coefficients
% Note: We will assume that the function is 2pi periodic
switch FncType

    % Smooth function
    % f = exp( sin(x) )
    % Will compute the Fourier coefficients numerically
    case ('smooth')
        % equispaced grid in [-pi, pi)        
        ngrid = 2^12;
        x = -pi + (2*pi/ngrid)*(0:ngrid-1).';
        % function definition
        fx = exp( sin(x) );
        % Fourier coefficients
        fHat = ( exp(-1i*k*x.') * fx )/ngrid;    
        
        % Function handle
        fx = @(y) exp( sin(y) );
        
        
    % Piecewise-smooth function (see definition of g in Exercise B)
    case ( 'piecewise' )
        % Fourier modes computed analytically using Mathematica
        fHat = exp(3*pi*1i*k/4).*(3./(4*pi*1i*k)) + ...
            exp(1i*pi*k/2).*(-3./(4*pi*1i*k)) + ...
            exp(-3*1i*pi*k/4).*(5./(1i*k) - ...
            (33*pi)./(16*1i*k) -11./(4*1i*1i*k.*...
            k))/(2*pi) + exp(-3*1i*pi*k/8).*...
            (-5./(1i*k) + (33*pi)./(32*1i*k) + ...
            11./(4*1i*1i*k.*k))/(2*pi) + ...
            exp(-1i*pi*k/8).*(-7./(4*1i*k) + ...
            pi./(16*1i*k) + 1./(2*1i*1i*k.*k) ...
            - sin(pi/8 - 1/4)*(k./(1i*(k.*k-1)))...
            + cos(pi/8 - 1/4)*(1./(k.*k-1)) )/(2*pi) + ...
            exp(1i*pi*k/4).*(7./(4*1i*k) + ...
            pi./(8*1i*k) - 1./(2*1i*1i*k.*k) +...
            sin(-pi/4 - 1/4)*(k./(1i*(k.*k-1))) ...
            - cos(-pi/4 -1/4)*(1./(k.*k-1)) )/(2*pi);

        % The DC component
        if( sum(k==0) == 1 )
            fHat(k==0) = ( (-pi/2+3*pi/4)*1.5 - ...
                5*(3*pi/4-3*pi/8) + 11*(9*pi*pi/16 - 9*pi*pi/64)/8 ...
                + 7*(pi/8+pi/4)/4 - .25*(pi*pi/64 - pi*pi/16) - ...
                cos(pi/8 - .25) + cos(-pi/4 -.25) )/(2*pi);
        end
        
        % Components at k=+/-1
        if( sum(k==1) == 1 )
            fHat(k==1) = -( 48 - 160*(-1)^(1/8) - (56+16i)*...
                (-1)^(3/8) -88*(-1)^(5/8) + (100-60i)*sqrt(2) - 8*...
                sqrt( 28-45i ) + (4-4i)*((1-1i)+sqrt(2))*exp(1i/4) + ...
                (33*(-1)^(1/8)+2*(-1)^(3/8)-(35-35i)*sqrt(2))*pi + ...
                6i*pi*exp(-1i/4) )/(64*pi);
        end
        if( sum(k==-1) == 1 )
            fHat(k==-1) = exp(-1i/4)*( (-4-4i)*((1+1i)+sqrt(2)) ...
                + 6i*exp(1i/2)*pi + exp(1i/4)*(-8*(6+(2+7i)*(-1)^(1/8) +...
                11*(-1)^(3/8) + 20*(-1)^(7/8) + (8+5i)*sqrt(2)) + ( 2*...
                (-1)^(5/8) +33*(-1)^(7/8) +(35+35i)*sqrt(2))*pi) )/(64*pi);
        end        
                    
  
        % Function handle
        fx = @(y)   1.5*(y>=-3*pi/4).*(y<-pi/2) + ...
                    ( 7/4 - y/2 + sin(y-1/4) ).*(y>=-pi/4).*(y<pi/8) + ...
                    ( y*11/4 - 5 ).*(y>=3*pi/8).*(y<3*pi/4);
                
end    

return