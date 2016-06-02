function done = PlotReconstruction(Coefficients) 


[reconstruction, domain] = ComputeFourierReconstruction(Coefficients);

plot(domain, reconstruction), title('Reconstruction');

done = true;


return;