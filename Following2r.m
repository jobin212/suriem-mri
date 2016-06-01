N = 50;
k = -N:N;

coefficients =  1 ./(pi*2i*k);

%replace k = 0 with integral value
coefficients(N+1) = 0;


[reconstruction, domain] = ComputeFourierReconstruction(coefficients);


%compute pi(x)
region1 = (domain < 0);
region2 = (domain > 0);

p = zeros(size(domain));

p(region1) = ((- pi - domain(region1))/ (2*pi))
p(region2) = ((pi - domain(region2)) / (2*pi))


%calculate error
error = abs(p - reconstruction);


plot(domain, reconstruction), title('Reconstruction');
figure;
plot(domain, error), title('Error Function');
