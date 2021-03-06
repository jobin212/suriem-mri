%get coefficients
FncType = 'box';
N = 100;
sN = num2str(N);
M = N;
sM = num2str(M);
leg = strcat('N=M=',sN); 

ReconstructionType = 'standard';

[fHat, fxy] = Get2DFourierCoefficients(FncType, N, M);

%overlay cross sections to compare error

[S_NMf, x, y] = Compute2DFourierReconstruction(fHat, ReconstructionType);

[xx, yy] = meshgrid(x,y);

abs_error = abs(S_NMf - fxy(xx,yy));


mesh(xx, yy, S_NMf);
%legend(leg);
xlabel('x')
ylabel('y')

figure;
mesh(xx, yy, abs_error);
zlim([0 0.6]);
view(45, 65);
legend(leg);





