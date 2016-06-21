%%%edge enhacned reconstruction has better code!

%calculate the coefficients for j jumps
N = 50;
k = -N:N;

%jump discontinuities
J = 1;
x = zeros(J);
%x(0) = x_0, x(1) = x_1, ...

%define jump function
r_jump = zeros(size(x));
%r_jump(0) = [r](x_0), r_jump(1) = [r](x_1)


%calculate coefficient estimates
sum = 0;
coefficients_est_k = zeros(size(k));
for y = 1:size(k)
    for z = 1:size(r_jump)
        coefficients_est_k(y) = coefficients_est_k(y) + (r_jump(z) * exp(-1i* x(z)*y))/(2i*pi*y);      
    end
end
        