N = 100
K = [-N :N]

f = zeros(size(K))

for i = K
    f(101 + i) = 2*(sin(i)) / (pi * i);
end

f

jobini = ComputeFourierReconstruction(f)


semilogy(jobini)