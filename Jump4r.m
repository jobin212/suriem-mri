%the higher k values appear to give a less accurate picture

N = 50;
M = 100;

k1 = -M:-N-1;
k2 = N+1:M;

k = [k1, k2];

%calculate coefficients according to equation
coefficients =  1 ./(pi*2i*k);

%jump discontinuities
J = 1;
x_0 = 0;

%define jump function
r_jump = 1;

%calculate coefficient estimates 
coefficients_est = (r_jump .* exp(-1*k*x_0*1i)) ./ (2i*pi*k)

%function call to plot reconstruction
PlotReconstruction(coefficients)
hold on
PlotReconstruction(coefficients_est)


