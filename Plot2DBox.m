%get coefficients
N = 100;
M = 100;
[fHat, fxy] =Get2DFourierCoefficients('box', N, M);

%compute reconstruction
[S_NMf, x, y] = Compute2DFourierReconstruction(fHat);

%use meshgrid for plotting
[xx, yy] = meshgrid(x,y);

%original function
mesh(xx, yy, fxy(xx,yy));
legend('Original Function fxy');


figure;
mesh(xx,yy,S_NMf.');
legend('Reconstruction N=M=100');

figure;
error = abs(S_NMf.' - fxy(xx,yy));
mesh(xx,yy,error);
legend('Error plot');

