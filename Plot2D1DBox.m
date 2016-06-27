%get coefficients
N = 35;
M = 35;
[fHat, fxy] =Get2DFourierCoefficients('box', N, M);

%compute reconstruction
[S_NMf, x, y] = Compute2D1DFourierReconstruction(fHat);

%use meshgrid for plotting
[xx, yy] = meshgrid(x,y);

figure;
mesh(xx,yy,S_NMf.');
%imagesc(x, y, S_NMf);