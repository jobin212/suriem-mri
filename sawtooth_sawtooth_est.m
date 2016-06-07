%%%%%%the real fourier coefficients up to N
N = 50;
K = -N:N;

coefficients_est =  1 ./(pi*2i*K);

%replace k = 0 with integral value
coefficients_est(N+1) = 0;

%%%%%%%%%%%%%%%%%% initialize jump function

%jump discontinuities
JUMPS = 1;
x_jump = zeros(JUMPS);

%define location of jumps: jump(0) = x_0, jump(1) = x_1, ...
jump(1) = 0;

%define jump function, define values: r_jump(0) = [r](x_0), r_jump(1) = [r](x_1)
r_jump = zeros(size(x_jump));
r_jump(1) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%% calculate f_k^est
%calculate coefficient estimates
sum = 0;
coefficients_est_k = zeros(size(K));


%%!! cannot use size(K) for index?
for k = 1:(2*N+1)
    for j = 1:JUMPS
        coefficients_est_k(k) = coefficients_est_k(k) + (r_jump(j) * exp(-1i* jump(j)*k))/(2i*pi*k);  
    end
end

coefficients_est_k(N + 1) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% reconstruct new difference of coefficients

coefficients_diff = coefficients_est - coefficients_est_k;

[reconstruction_combined, domain] = ComputeFourierReconstruction(coefficients_diff);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculating second half of the sum

jump_sum = zeros(size(domain));

%define sawtooth function g_j
g = zeros(JUMPS, 1024);

%J + 1 piecewise smooth for J jumps
region1 = (domain < jump(1));
region2 =  (domain > jump(1));

g(1,:) = zeros(size(domain));


g(1,region1) = ((- pi - domain(region1))/ (2*pi));
g(1,region2) = ((pi - domain(region2)) / (2*pi));


%%caculate sum
for x = 1:size(jump_sum)
    for j = 1:JUMPS
        jump_sum(x) = jump_sum(x) + (r_jump(j) * g(1,x));
    end
end

reconstruction_final = reconstruction_combined + jump_sum;

%plot
plot(domain, reconstruction_final);

