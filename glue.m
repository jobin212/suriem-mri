clear; close all; clc

SampleEdgeDetectionCode;

r = EdgeEnhancedReconstruction(fk, jmp_ht, x(jmp_locs));
figure;
plot(x, r)

error = abs(f(x) - r);
semilogy(x, error);
title('Error using Trigonometric Concentration Kernel for N=50')

