%%%the amount of exponenents we want to take 
exp = 4;

%fourier coefficients
k = 50*2.^(0:exp);

error_vector = zeros(size(k));

for i = 1:length(k)
    error_vector(i) = NormErrorTrueJumps(k(i));
end

figure;
plot(log2(k), log2(error_vector));
