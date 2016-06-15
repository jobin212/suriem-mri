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

% Plot
plot(x, f(x), 'k', 'linewidth', 2); hold on
plot(x, jmp_apprx, 'r', 'linewidth', 2); 
stem(x(jmp_locs), jmp_ht, 'b'); grid; xlim([-pi pi])
xlabel '$x$'; ylabel '$\tilde S_N^\sigma[f](x)$'
title 'Approximating Jumps from Fourier Data'
legend( '$f$', '$\tilde S_N^\sigma[f]$', 'Jumps')


x_jump = zeros(6,1);

%define location of jumps: jump(0) = x_0, jump(1) = x_1, ... and jump
%values 
x_jump(1) = (-3*pi/4);
x_jump(2) = (-pi/2);
x_jump(3) = (-pi/4);
x_jump(4) = (pi/8);
x_jump(5) = (3*pi/8);
x_jump(6) = (3*pi/4);


%define jump function, define values: r_jump(0) = [r](x_0), r_jump(1) = [r](x_1)
r_jump = zeros(size(x_jump));
r_jump(1) = (3/2);
r_jump(2) = (-3/2);
r_jump(3) = (7/4) - ((-pi/4)/2) + sin((-pi/4)- (1/4));
r_jump(4) = -((7/4) - ((pi/8)/2) + sin((pi/8)- (1/4)));
r_jump(5) = ((33/32)*pi -5);
r_jump(6) = -((33/16)*pi -5);


