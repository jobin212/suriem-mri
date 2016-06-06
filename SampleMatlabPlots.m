% Model solutions from Matlab problem set (Exercise B)

clear; close all; clc

%% Problem 1.1 (1)

% Generate Fourier coefficients
N = 32;
fhat_k = GetFourierCoefficients( 'smooth',    N );
ghat_k = GetFourierCoefficients( 'piecewise', N );

% Plot coefficients (in log scale)
k = (-N:N).';             % We have 2N+1 coefficients

figure(1); 
semilogy( k, abs(fhat_k), 'r-+', 'linewidth', 2); 
            % 'r-+' means plot in red color using a solid line with a +
            % marker
            % the linewidth argument makes the line thicker
hold on     % this allows us draw multiple plots in the same graph
semilogy( k, abs(ghat_k), 'b-o', 'linewidth', 2);
legend( {'$|\hat f_k|$', '$|\hat g_k|$'}, 'interpreter', 'latex', ...
                'fontsize', 14 );   
            % you can typeset the legend in LaTeX at a given fontsize
xlabel( '$k$', 'interpreter', 'latex', 'fontsize', 14 );
            % as before, you can label the axes in LaTeX at a desired
            % fontsize
ylabel( '$|\hat f_k| \, \rm{and} \, |\hat g_k|$', ...
                'interpreter', 'latex', 'fontsize', 14 );
title( 'Fourier Coefficients of $f$ and $g$', ...
                'interpreter', 'latex', 'fontsize', 14 )
grid on;    % this superimposes a set of gridlines
xlim([-N N])% this limits the plot to the interval [-N,N] 


%% Problem 1.1 (3)

% Fourier partial sum approcimation of f

% First, get Fourier coefficients for N = 4, 8, 16
[fhat_4, fx]  = GetFourierCoefficients( 'smooth',  4 );
fhat_8        = GetFourierCoefficients( 'smooth',  8 );
fhat_16       = GetFourierCoefficients( 'smooth', 16 );
% Note: fx is a function handle which contains a description of f (we will
% use it soon...)

% Compute and plot approximation for N = 8
[S_Nf_8, x] = ComputeFourierReconstruction( fhat_8 );
% Note: Here, x is the equispaced grid in [-pi, pi) where we plot our
% reconstruction.

figure; plot( x, fx(x), 'k', 'linewidth', 2 ); hold on
plot( x, S_Nf_8, 'r--', 'linewidth', 2);
xlabel x; ylabel S_Nf(x); title 'Partial Fourier Sum Reconstruction'
legend( 'f', 'S_Nf' ); grid

% Let's now make an error plot
S_Nf_4 = ComputeFourierReconstruction( fhat_4 );
S_Nf_16 = ComputeFourierReconstruction( fhat_16 );

figure; semilogy( x, abs( fx(x) - S_Nf_4 ), 'r', 'linewidth', 2 ); hold on
semilogy( x, abs( fx(x) - S_Nf_8 ), 'b', 'linewidth', 2 );
semilogy( x, abs( fx(x) - S_Nf_16 ), 'k', 'linewidth', 2 );
xlabel x; grid on
ylabel( '$|f(x) - S_Nf(x)|$', 'interpreter', 'latex', 'fontsize', 14 );
legend( {'$S_4f$', '$S_8f$', '$S_{16}f$'}, 'interpreter', 'latex' );
axis([-pi pi 1e-16 1e0])
title 'Absolute Error Plot'
