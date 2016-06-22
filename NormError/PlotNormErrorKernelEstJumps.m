%%%the amount of exponenents we want to take 
exp = 7;

%fourier coefficients
k = 200*2.^(0:exp);

error_vector = zeros(size(k));

for i = 1:length(k)
    error_vector(i) = NormErrorKernelEstJumps(k(i));
end

figure;
plot(log2(k), log2(error_vector));
