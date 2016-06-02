N = 50;
M = 100;

k1 = -M:-N-1;
k2 = N+1:M;

k = [k1, k2]

%jump discontinuities
J = 1;
x_0 = 0;

%define jump function
r_jump = 1

%calculate coefficients
coefficients = (r_jump .* exp(k*x_0*i)) ./ (2i*pi*k)

%plot (code copied from Jump3r.m)
%compute r(x)
region1 = (domain < 0);
region2 = (domain > 0);

r = zeros(size(domain));

r(region1) = ((- pi - domain(region1))/ (2*pi))
r(region2) = ((pi - domain(region2)) / (2*pi))


%calculate error
error = abs(r - reconstruction);


plot(domain, reconstruction), title('Reconstruction');
figure;
plot(domain, error), title('Error Function');

