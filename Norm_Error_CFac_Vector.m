n = [0:1:5]
k = 50*2.^n

ev = zeros(size(k))

for i = [1:6]
    ev(i) = Norm_Error_CFac(k(i));
end

figure;

lev = log2(ev)
plot(lev)