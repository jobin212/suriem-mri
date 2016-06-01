s16 = GetFourierCoefficients('smooth', 16);
s64 = GetFourierCoefficients('smooth', 64);
s256 = GetFourierCoefficients('smooth', 256);

xyz = size(s16);
yy = size(s64);
zz = size(s256);


s16r = ComputeFourierReconstruction(s16);
s64r = ComputeFourierReconstruction(s64);
s256r = ComputeFourierReconstruction(s256);


%size of array
ss16 = size(s16r);
ss64 = size(s64r);
ss256 = size(s256r);

%adjust x values so that we can subtract f and reconstruction
x16r = linspace(-pi, pi, ss16(1));
x64r = linspace(-pi, pi, ss64(1));
x256r = linspace(-pi, pi, ss256(1));
x = x16r

%calculate f values
%we are lucky because size(x16r) = size(x64r) = size(x256r)
gg = zeros(ss16)

region1 = ((x < (-3*pi /4)) | (x > (3 * pi /4)) | (x >= -pi /2 & x < -pi / 4) | (x >= pi /8 & x < 3* pi / 8));
gg(region1) = 0;

region2 = (x > (-3*pi/4) & x <= -pi /2);
gg(region2) = 3/2;

region3 = (x >= -pi /4 & x < pi /8);
gg(region3) = 7/4 - x(region3)/2 + sin(x(region3) - (1/4));

region4 = (x >= 3*pi /8 & x < 3*pi /4)
gg(region4) = (11/4) * x(region4) - 5
 


size(s16r)
size(gg)



%calculate error funciton
err16 = abs(gg - s16r)
err64 = abs(gg - s64r)
err256 = abs(gg - s256r)



plot(x, err16-err64)


%it appears that adding points does not noticeably improve the graph
%plot
semilogy(x16r, err16, x16r, err64, x16r, err256)

title('Error plot for N = 16');
legend('N = 16', 'N = 64', 'N = 256');

xyz
yy
zz