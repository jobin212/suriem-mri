


%jump discontinuities
J = 1;
x_jump = zeros(J);

%define values: x(0) = x_0, x(1) = x_1, ...
x(0) = 0;

%define jump function
r_jump = zeros(size(x_jump));

%define values: r_jump(0) = [r](x_0), r_jump(1) = [r](x_1)
r_jump(0) = 1;

%define sawtooth function g_j
%
region1 = zeros(J +1)

