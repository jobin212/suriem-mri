function sig = confac( fourModes, cfacType )

%% Generate concentration factors
% Script to generate concentration factors on a given set of Fourier modes
%
% Usage:    sig = confac(fourModes)
%
%   Inputs:
%       fourModes - the modes on which to compute the Fourier coefficients
%                   (real vector)
%       cfacType - type of concentration factor
%                   'Trig' - Trigonometric
%                   'Poly' - First order polynomial
%                   'Exp' - Sixth order exponential
%
%   Outputs:
%       sig - the concentration factors (real vector)
%

switch ( cfacType )

    case 'Trig'
        % Gibbs/Trigonometric concentration factor
        Sipi = 1.85193705198247;            % normalizing constant
        sig = pi*sin( pi*abs(fourModes)/max(fourModes) )/Sipi;

    case 'Poly'
        % Polynomial factor
        sig = ( pi*abs(fourModes)/max(fourModes) );

    case 'Exp'
        % Exponential factor (alpha=6)
        alpha = 6;
        tau = linspace(1/max(fourModes), 1-1/max(fourModes), 1000);
        res = tau(2)-tau(1);
        const = pi/( res*sum(exp(1./(alpha*tau.*(tau-1)))) );
        sig = const*( abs(fourModes)/max(fourModes) ).* ...
                exp(1./(alpha*(abs(fourModes)/max(fourModes)).*...
                                ((abs(fourModes)/max(fourModes))-1)));
        sig(fourModes==0) = 0; sig(1) = 0; sig(end) = 0;
        
end

return


