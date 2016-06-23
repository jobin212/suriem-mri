function resdl = FitJumps( jmp_info, k, fk )
%FitJumps Estimate Fourier Coefficients from Jump Information

% Process input arguments
jmp_ht   = jmp_info(:,1);
jmp_locs = jmp_info(:,2);

% Compute estimated Fourier coefficients
fk_est = zeros(length(k), 1);
% This is a sum over the total number of jumps
for ix = 1:length(jmp_locs)
    cfk_est = fk_est + jmp_ht(ix)*exp(-1i*k*jmp_locs(ix))./(2i*pi*k);
end

% Compute residual
resdl = (fk_est-fk);

end

