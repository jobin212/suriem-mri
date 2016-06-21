n = [0:1:5]
k = 50*2.^n

ev = zeros(size(k))

for i = [1:6]
    ev(i) = NormError(k(i));
end

lk = log(k) 
lev = log(ev)
plot(lk, lev);

figure
loglog(k , ev, k, k.^-1, '--', k, 1e-1*k.^-0.5, '-.')
legend('error', 'first order slope', 'N^-0.5')

%lev = log2(ev)
%plot(lev)