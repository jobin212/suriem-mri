clear; close all; clc

% for pretty pictures
set(0,'defaulttextinterpreter','latex');
set(0,'defaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',14);
% set(0,'TicklabelInterpreter','latex');
set(0,'DefaultLegendInterpreter','latex')

% Make convolution movie/animation

% Grid
ngrid = 1000;
x = -pi + (2*pi/ngrid)*(0:ngrid-1).';

% Function definitions
f1 = @(x) double(abs(x)<=pi/2);
% General convolution
% f2 = @(x) (x>=-pi/2).*(x<=pi/2) + (x>=3*pi/2).*(x<=5*pi/2) + ...
%           (x>=-5*pi/2).*(x<=-3*pi/2);
% Edge detection
N = 50;
k = (-N:N).';
Sipi = 1.85193705198247;            % normalizing constant
sig = pi*sin( pi*abs(k)/N )/Sipi;
f2 = @(x) real( exp(1i*x*k.')*( 1i*sign(k).*sig/(2*pi) ) );

% Plot f_1
figure; subplot(3, 1, 1); plot( x, f1(x), 'b' ); grid on; xlabel '$x$'
ylabel '$f_1(u)$'; axis([-4 4 0 1.25])

% Store result here
g = zeros(ngrid, 1);
F(ngrid) = struct('cdata',[],'colormap',[]);

for ix = 1:ngrid
    y = x(ix);
    
    % Flip and Shift second function
    subplot(3,1,2); plot(x,f2(y-x),'r'); grid; 
    xlabel '$x$'; ylabel '$f_2(x-u)$'; 
    axis([-4 4 1.25*min(f2(y-x)) 1.25*max(f2(y-x))])

    % Evaluate Convolution
    g(ix) = sum( f1(x) .* f2(y-x) )*(2*pi/ngrid);
    subplot(3,1,3); plot(x(1:ix), g(1:ix),'k'); grid; 
%     axis([-4 4 0 1.1*pi])
    axis([-4 4 -1.25 1.25])
    xlabel '$x$'; ylabel '$(f_1*f_2)(x)$'; drawnow;
    
    % Make movie
    F(ix) = getframe(gcf);
end

% Export movie
movie2avi(F, 'edge_detection_animation.avi');