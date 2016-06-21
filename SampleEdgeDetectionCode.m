% Sample edge detection code

clear; close all; clc

% for pretty pictures
set(0,'defaulttextinterpreter','latex');
set(0,'defaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultLegendInterpreter','latex')

% Generate Fourier coefficients
N = 50;
[fk, f] = GetFourierCoefficients( 'piecewise', N );

% Concentration factors
k = (-N:N).';             % We have 2N+1 coefficients
sig_p = confac(k, 'Poly');
sig_t = confac(k, 'Trig');
sig_e = confac(k, 'Exp');

% Compute jump approximation
jmp_cfs = (1i*sign(k).*sig_t) .* fk;
[jmp_apprx, x] = ComputeFourierReconstruction(jmp_cfs);

% Extract jump locations and values
[jmp_ht, jmp_locs] = findpeaks( abs(jmp_apprx), 'MinPeakHeight', 0.5, ...
                            'MinPeakDistance', 10/(2*N+1));
% since we may have negative jumps...                        
jmp_ht = jmp_ht .* sign( jmp_apprx(jmp_locs) );

%first part of sum
coeff_diff = fk - jmp_cfs;

% Plot
plot(x, f(x), 'k', 'linewidth', 2); hold on
plot(x, jmp_apprx, 'r', 'linewidth', 2); 
stem(x(jmp_locs), jmp_ht, 'b'); grid; xlim([-pi pi])
xlabel '$x$'; ylabel '$\tilde S_N^\sigma[f](x)$'
title 'Approximating Jumps from Fourier Data'
legend( '$f$', '$\tilde S_N^\sigma[f]$', 'Jumps')


%%%%true values for piecewise jump locations, heights
%true_jump_locs = [-3*pi/4, -pi/2, -pi/4, pi/8, 3*pi/8, 3*pi/4]
%true_jump_heights = [3/2, -3/2, (7/4) - ((-pi/4)/2) + sin((-pi/4)- (1/4)),
%-((7/4) - ((pi/8)/2) + sin((pi/8)- (1/4))), ((33/32)*pi -5), -((33/16)*pi
%-5)];
