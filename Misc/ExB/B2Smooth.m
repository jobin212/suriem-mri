s4 = GetFourierCoefficients('smooth', 4);
s8 = GetFourierCoefficients('smooth', 8);
s16 = GetFourierCoefficients('smooth', 16);

s4r = ComputeFourierReconstruction(s4);
s8r = ComputeFourierReconstruction(s8);
s16r = ComputeFourierReconstruction(s16);

%{
semilogy(abs(s4r))
hold on
semilogy(abs(s8r))
hold on
semilogy(abs(s16r))
hold on
%}

%{
semilogy(abs(s4)), title('Smooth with N = 4, 8, 16')
hold on 
semilogy(abs(s8))
hold on
semilogy(abs(s16))
%}

%size of array
ss4 = size(s4r);
ss8 = size(s8r);
ss16 = size(s16r);

%adjust x values so that we can subtract f and reconstruction
x4r = linspace(-pi, pi, ss4(1));
x8r = linspace(-pi, pi, ss8(1));
x16r = linspace(-pi, pi, ss16(1));

%calculate f values
%we are lucky because size(x4r) = size(x8r) = size(x16r)
ff = exp(sin(2*x4r));

%calculate error funciton
err4 = abs(ff' - s4r)
err8 = abs(ff' - s8r)
err16 = abs(ff' - s16r)


%plot
semilogy(x4r, err4, x4r, err8, x4r, err16)

title('Error plot for N = 4');
legend('N = 4', 'N = 16', 'N = 32');

