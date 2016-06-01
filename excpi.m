N = 50;

coefficients = 1:2*N+1

for i = 1:(2*N+1)
    k = coefficients(i) - N -1;
    
    %%insert formula here
    coefficients(i) = ((-1)^(k)) / (pi^(2)*k);
end
  
%replace k = 0 with 0 (evalute original integral with k = 0)
coefficients(N+1) = 0

[reconstruction, domain] = ComputeFourierReconstruction(coefficients);


%%

plot(domain, reconstruction);

