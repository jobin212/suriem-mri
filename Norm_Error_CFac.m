function [e] = Norm_Error_CFac(N)


clear; close all; clc

SampleEdgeDetectionCode;
f;
x_jump = zeros(6);
x_jump(1) = (-3*pi/4);
x_jump(2) = (-pi/2);
x_jump(3) = (-pi/4);
x_jump(4) = (pi/8);
x_jump(5) = (3*pi/8);
x_jump(6) = (3*pi/4);

r_jump = zeros(6)
 r_jump(1) = (3/2);
 r_jump(2) = (-3/2);
 r_jump(3) = (7/4) - ((-pi/4)/2) + sin((-pi/4)- (1/4));
 r_jump(4) = -((7/4) - ((pi/8)/2) + sin((pi/8)- (1/4)));
 r_jump(5) = ((33/32)*pi -5);
 r_jump(6) = -((33/16)*pi -5);
%x;
%S_Nf_edge = EdgeEnhancedReconstruction(fk, jmp_ht, x(jmp_locs));
S_Nf_edge = EdgeEnhancedReconstruction(fk, r_jump, x_jump);
%figure;
%plot(x, S_Nf_edge)

error = (f(x) - S_Nf_edge);
h = (2*pi)/1024;
e = h*norm(error);

return