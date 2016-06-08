%%%%%%the real fourier coefficients up to N
N = 50;
K = -N:N;

%same formula
coefficients_est =  GetFourierCoefficients('piecewise', N);

%replace k = 0 with integral value
[reconstruction_orig, domain] = ComputeFourierReconstruction(coefficients_est);

edge_coefficients = coefficients_est .* (1i*sign(K) .* pi .*abs(K) ./ N).';

[edge_approx, domain] = ComputeFourierReconstruction(edge_coefficients);
plot(domain, reconstruction_orig, 'k');
hold on;
stem(domain, edge_approx .* (abs(edge_approx) > 1));
legend('Function', 'Conjugate Sum');
title('Conjugate Sum with Polynomial Kernel');


figure;

%%%%%%%%%%%%%%%%%% initialize jump function

%jump discontinuities
JUMPS = 6;
x_jump = zeros(JUMPS,1);

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


%calculate coefficient estimate
coefficients_est_k = zeros(size(K));

for k = -N:N
    for j = 1:JUMPS
        coefficients_est_k(k+N+1) = coefficients_est_k(k+N+1) + (r_jump(j) * exp(-1i* x_jump(j)*k))/(2i*pi*k);  
    end
end

coefficients_est_k(N + 1) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% reconstruct new difference of coefficients

coefficients_diff = coefficients_est.' - coefficients_est_k;

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
g = @(x) ((-pi-x)/(2*pi) .* (x<0)) + ((pi -x)/(2*pi) .* (x>=0));
%h = @(y) (pi -y) /(2*pi) + floor(y/(2*pi));

%plot(domain, g(domain), 'k-*', domain, h(domain), 'r-o');



%%caculate sum
for x = 1:size(jump_sum)
    for j = 1:JUMPS
        jump_sum(x) = jump_sum(x) + (r_jump(j) * g(domain(x) - x_jump(j)));
    end
end

reconstruction_final = reconstruction_combined + jump_sum;






%%error
gg = zeros(size(domain))

%region1 = ((domain < (-3*pi /4)) | (domain >= (3 * pi /4)) | (domain >= -pi /2 & domain < -pi / 4) | (domain >= pi /8 & domain< 3* pi / 8));
%gg(region1) = 0;

region2 = (domain >= (-3*pi/4) & domain< -pi /2);
gg(region2) = 3/2;

region3 = (domain >= -pi /4 & domain < pi /8);
gg(region3) = 7/4 - domain(region3)/2 + sin(domain(region3) - (1/4));

region4 = (domain >= 3*pi /8 & domain < 3*pi /4);
gg(region4) = (11/4) * domain(region4) - 5;


%plot
plot(domain, reconstruction_final, 'k', domain, reconstruction_orig, 'r', domain, gg, 'm');
legend('New Reconstruction', 'Old Reconstruction', 'Original Function');
title('Reconstruction');

figure;

err_final = abs(reconstruction_final - gg);
err_orig = abs(reconstruction_orig - gg);

semilogy(domain, err_final, 'k', domain, err_orig, 'r');
title('Error plot');
legend('New Reconstruction', 'Old Reconstruction');








