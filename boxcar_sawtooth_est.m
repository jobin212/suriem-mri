%%%%%%the real fourier coefficients up to N
N = 50;
K = -N:N;

%same formula
coefficients_est =  sin(K*pi/2) ./(pi^2 *K);

%replace k = 0 with integral value
coefficients_est(N+1) = 0;
[reconstruction_orig, ~] = ComputeFourierReconstruction(coefficients_est);

%%%%%%%%%%%%%%%%%% initialize jump function

%jump discontinuities
JUMPS = 2;
x_jump = zeros(JUMPS,1);

%define location of jumps: jump(0) = x_0, jump(1) = x_1, ...
x_jump(1) = -pi/2;
x_jump(2) = pi/2;

%define jump function, define values: r_jump(0) = [r](x_0), r_jump(1) = [r](x_1)
r_jump = zeros(size(x_jump));
r_jump(1) = 1/pi;
r_jump(2) = -1/pi;


%calculate coefficient estimate
coefficients_est_k = zeros(size(K));

for k = -N:N
    for j = 1:JUMPS
        coefficients_est_k(k+N+1) = coefficients_est_k(k+N+1) + (r_jump(j) * exp(-1i* x_jump(j)*k))/(2i*pi*k);  
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

%%anonymous function
g = @(x) ((-pi-x)/(2*pi) .* (x<0)) + ((pi -x)/(2*pi) .* (x>0));

%%caculate sum
for x = 1:size(jump_sum)
    for j = 1:JUMPS
        jump_sum(x) = jump_sum(x) + (r_jump(j) * g(domain(x) - x_jump(j)));
    end
end

reconstruction_final = reconstruction_combined + jump_sum;

%plot
plot(domain, reconstruction_final, 'k', domain, reconstruction_orig, 'r');
legend('New Reconstruction', 'Old Reconstruction');


