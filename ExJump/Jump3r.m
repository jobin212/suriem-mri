%%%%%%what are the approximation properties?
N = 50;
k = -N:N;

coefficients =  1 ./(pi*2i*k);

%replace k = 0 with integral value
coefficients(N+1) = 0;

%Reconstruction
[reconstruction, domain] = ComputeFourierReconstruction(coefficients);


%compute r(x)
region1 = (domain < 0);
region2 = (domain > 0);

r = zeros(size(domain));

r(region1) = ((- pi - domain(region1))/ (2*pi))`
r(region2) = ((pi - domain(region2)) / (2*pi))


%calculate error
error = abs(r - reconstruction);

%plot
plot(domain, r, 'k', domain, reconstruction, 'r','LineWidth',2 ), title('r(x) and its reconstruction');
legend('r(x)', 'reconstruction');

%figure;
%plot(domain, error), title('Error Function');
