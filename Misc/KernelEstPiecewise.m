function [kernel_est_jump_heights, kernel_est_jump_locs] = KernelEstPiecewise(N)

%clear; close all; clc

% Generate Fourier coefficients
%N = 50;
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
[kernel_est_jump_heights, kernel_est_jump_locs] = findpeaks( abs(jmp_apprx), 'MinPeakHeight', 0.5, ...
                            'MinPeakDistance', 10/(2*N+1));
% since we may have negative jumps...                        
kernel_est_jump_heights = kernel_est_jump_heights .* sign( jmp_apprx(kernel_est_jump_locs) );

return


