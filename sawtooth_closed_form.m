%%%%%%the real fourier coefficients up to N
N = 50;
k = -N:N;

coefficients =  1 ./(pi*2i*k);

%replace k = 0 with integral value
coefficients(N+1) = 0;

%Reconstruction
[reconstruction, domain] = ComputeFourierReconstruction(coefficients);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%jump discontinuities
J = 1;
x_jump = zeros(J);

%define values: x(0) = x_0, x(1) = x_1, ...
x(1) = 1;

%define jump function
r_jump = zeros(size(x_jump));

%define values: r_jump(0) = [r](x_0), r_jump(1) = [r](x_1)
r_jump(1) = 1;

%define sawtooth function g_j
%J + 1 piecewise smooth for J jumps
region1 = (domain < x(1))
region2 =  (domain > x(1))

r = zeros(size(domain));

r(region1) = ((- pi - domain(region1))/ (2*pi))
r(region2) = ((pi - domain(region2)) / (2*pi))

%plot
plot(domain, r);
