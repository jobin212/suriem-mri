function [est_jmpvals,est_jmplocs] = PronyMethod(N,nJumps,ndata)

%% Prony's Method - First Implementation!

%clear; close all; clc

%%Getting Fourier Coefficients
%N         = N;
k         = (-N:N).';
[fk, f]   = GetFourierCoefficients('piecewise', N);

% True jump information
% true_jmpvals = [        3/2 ...
%                        -3/2 ...
%                         f(-pi/4) ...
%                         -(7/4-pi/16+sin(pi/8-1/4)) ...
%                         f(3*pi/8) ...
%                         -((11/4)*(3*pi/4)-5) ...
%                                 ].';
%                             
% true_jmplocs = [-3*pi/4 -pi/2 -pi/4 pi/8 3*pi/8 3*pi/4].';


% Multiply out the (2i pi k) factor
d_fk      = (2i*pi*k) .* fk;

% Jump parameters
%nJumps = 06;         % Max. expected no. of jumps
%ndata  = 20;         % No. of data/coefficient measurements to use
M      = 01;         % We will use Fourier modes -(N+M-1:N+M+ndata-2)


% Set up linear prediction model
LinPredictionMat = zeros(ndata, nJumps);
for ix = 1:ndata
    LinPredictionMat(ix,:) = d_fk( M+(nJumps-1)+(ix-1) : -1 : M+(ix-1) );
%     disp( M+(nJumps-1)+(ix-1) : -1 : M+(ix-1) )   % Indexing
end

model_measurements = d_fk( M+nJumps : M+nJumps+ndata-1 );
% disp( (M+nJumps : M+nJumps+ndata-1).' )           % Indexing

% Solve linear prediction model
poly_cfs = LinPredictionMat\model_measurements;

% Getting the roots of the associated polynomial
poly_roots = roots([1; -poly_cfs]);

% Get Jump Locations
est_jmplocs = unwrap( -angle(poly_roots) );
est_jmplocs = sort( est_jmplocs )              % Sort

% Print these out
% fprintf( '\n Here are the true jump locations \n' );
% disp( true_jmplocs.' )
% 
% fprintf( ' Here are the estimated jump locations \n' );
% disp( est_jmplocs.' )

% % Error
% fprintf( ' Max absolute error in jump locations is %3.3e\n\n', ...
%                norm(true_jmplocs-est_jmplocs, inf) );

% Jump Heights
% Set up and solve linear system
fit_cfs = (abs(k)>0.75*N);       % Which coefficients to use for data fit?
sysMat      = exp( -1i*k(fit_cfs)*est_jmplocs.' );
est_jmpvals = sysMat \ d_fk(fit_cfs)

% Print these out
% fprintf( '\n Here are the true jump heights \n' );
% disp( true_jmpvals.' )
% 
% fprintf( ' Here are the estimated jump heights \n' );
% disp( real(est_jmpvals.') )

% % Error
% fprintf( ' Max absolute error in jump values is %3.3e\n\n', ...
%                norm(true_jmpvals-est_jmpvals, inf) );
 
%%%Plotting Original Function
 %z = linspace(-pi,pi,4096);
  
 %figure;

% %%Plotting Error in Jump Locs and Vals
% err_jmplocs = abs(est_jmplocs-true_jmplocs);
% err_jmpvals = abs(est_jmpvals-true_jmpvals);
% semilogy(z, err_jmplocs, 'r', z, err_jmpvals,'b');

%%Plotting Function with Edge 
%%see joe_piecewise
%plot(z,EdgeEnhancedReconstruction(fk,real(est_jmpvals),est_jmplocs),'r',z,f(z),'b')
%legend('Reconstruction', 'Function');
%title('Partial Sum Using Prony Approximations');

%figure;

%semilogy(abs(f(z).'-EdgeEnhancedReconstruction(fk,real(est_jmpvals),est_jmplocs)),'k')
%title('Partial Sum Error (N=50)');

return
 
 
 