%This is the same as inputting adjusted coefficient arrays into
%reconstruciton function

%%not plotting correctly

%%TODO generalize for J jumps

N = 100;
M = N+50;

%define intervals
k0 = -N:N;
k1 = -M:-N-1;
k2 = N+1:M;

%define jumps
%jump discontinuities
J = 1;
x_0 = 0;

%define jump function
r_jump = 1;


%calculate coefficient
coefficients_k0 =  1 ./(pi*2i*k0);
%don't forget to adjust for null values!
coefficients_k0(N+1) = 0;

%calculate coefficient estimates 
coefficients_est_k1 = (r_jump .* exp(-1*k1*x_0*1i)) ./ (2i*pi*k1);
coefficients_est_k2 = (r_jump .* exp(-1*k2*x_0*1i)) ./ (2i*pi*k2);

%combine coefficients
coefficients = [coefficients_est_k1, coefficients_k0, coefficients_est_k2];

%reconstruction
[reconstruction_mn, domain_mn] = ComputeFourierReconstruction(coefficients);
[reconstruction_n, domain_n] = ComputeFourierReconstruction(coefficients_k0);

plot(domain_mn, reconstruction_mn, domain_n, reconstruction_n), title('Reconstruction');
legend('mn', 'n');



