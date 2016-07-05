function resdl = FitJumpValues( jmp_hts, k, fk, jmp_locs )
%FitJumps Estimate Fourier Coefficients from Jump Information

% Compute estimated Fourier coefficients
fk_est = zeros(length(k), 1);
% This is a sum over the total number of jumps
for ix = 1:length(jmp_locs)
    fk_est = fk_est + jmp_hts(ix)*exp(-1i*k*jmp_locs(ix))./(2i*pi*k);
end

% Compute residual
resdl = (fk_est-fk);

end

