N = 50;

%coefficients = 1:2*N+1;


k = -N:N;
coefficients = (sin(k*pi/2) ./ (pi^(2)*k));

%replace k = 0 with 0 (evalute original integral with k = 0)
coefficients(N+1) = 0;


[reconstruction, domain] = ComputeFourierReconstruction(coefficients);


%compute pi(x)
region1 = (domain <= -pi/2 | domain >= pi/2);
region2 = (domain > -pi/2 & domain < pi/2);

p = zeros(size(domain));
p(region2) = 1/ pi;


%calculate error
error = abs(p - reconstruction);


plot(domain, reconstruction);
figure;
plot(domain, error);
