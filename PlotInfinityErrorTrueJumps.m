%%%the amount of exponenents we want to take 
exp = 7;

%fourier coefficients
k = 50*2.^(0:exp);

error_vector = zeros(size(k));

for i = 1:length(k)
    error_vector(i) = InfinityErrorTrueJumps(k(i));
end

figure;
loglog(k, error_vector, k, k.^(-1))
