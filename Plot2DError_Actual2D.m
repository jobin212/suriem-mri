clear; close all; clc

% for pretty pictures
set(0,'defaulttextinterpreter','latex');
set(0,'defaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultLegendInterpreter','latex')

%%%the amount of exponenents we want to take 
ErrType = '2norm';
FncType = 'circle';
% ReconstructionType = 'standard';

% ReconstructionType = 'box-true-jumps';
% ReconstructionType = 'box-prony-jumps';
% ReconstructionType = 'box-conc-jumps';

ReconstructionType = 'circle-true-jumps';
% ReconstructionType = 'circle-prony-jumps';
% ReconstructionType = 'circle-conc-jumps';
%order of errors decreased

%fourier coefficients
k = 2.^(5:8);

error_vector = zeros(size(k));

for i = 1:length(k)
    fprintf( '\n No computing error for N=%d', k(i) );
    N = k(i);
    M = N;
%     x0 = 0.10;

    kn = (-N:N).';
    km = (-M:M).';
    
    ngrid = 2*(2*N+1);
    mgrid = 2*(2*M+1);
%     mgrid = 2000;
    
    x = -pi + (2*pi/ngrid)*( 0:ngrid-1 ).';
    y = -pi + (2*pi/mgrid)*( 0:mgrid-1 ).';

%     x0 = x(3*N+4);
    
    [fHat, fxy] = Get2DFourierCoefficients(FncType, N, M);
    
    %%take from  Compute2DFourierReconstruction)  
    if(strcmp(ReconstructionType, 'standard'))
        %outer product
        F1 = exp(1i * x * kn.');
        F2 = exp(1i * y * km.');

        S_NMf = F1 * fHat * F2.';
        S_NMf = real( S_NMf );
        
        
    else
        %%calculate    
        
        S_NMf = zeros(ngrid, mgrid);
        
        for igrid = 1:ngrid
            comp_exp = exp(1i * x(igrid) * kn.');
        
            CFR = comp_exp * fHat;
        
            jmp_heights = [];
            jmp_locs = [];


            switch(ReconstructionType)
                case('box-true-jumps')
                    jmp_heights = [0];
                    jmp_locs = [0];

                    if(abs(x(igrid)) <= 1) 
                        jmp_heights = [-1 1].';
                        jmp_locs = [1 -1].';
                    end;
                case('circle-true-jumps')
                    jmp_heights = [0];
                    jmp_locs = [0];

                    if( abs(x(igrid)) <= 1) 
                        jmp_locs = [-1*sqrt(1- x(igrid)^2) 1*sqrt(1-x(igrid)^2)].';
                        jmp_heights = [1 -1].';
                    end

                case('box-prony-jumps')

                    jmps = 0;

                    if( x(igrid) <= 1 && x(igrid) >= -1)
                        jmps = 2;
                    end


                    [jmp_locs, jmp_heights] = FindJumps(CFR, 'prony', false, [], jmps);

                case('circle-prony-jumps')
                    jmps = 0;

                    if( x(igrid) <= 1 && x(igrid) >= -1)
                        jmps = 2;
                    end

                    [jmp_locs, jmp_heights] = FindJumps(CFR, 'prony', false, [], jmps);

                case('box-conc-jumps')
                    [jmp_locs, jmp_heights] = FindJumps(CFR, 'conc', true);

                case('circle-conc-jumps')
                    [jmp_locs, jmp_heights] = FindJumps(CFR, 'conc', true);
            end;  




            S_NMf(igrid,:) =  EdgeEnhancedReconstruction(CFR, jmp_heights, jmp_locs, mgrid);

        end 
    end
    
    [xgrid, ygrid] = ndgrid(x, y);
    f = fxy(xgrid, ygrid);
    
    switch ErrType
        case '2norm'
            h = y(2)-y(1);
            error_vector(i) = h*norm(f(:)-S_NMf(:));
        case 'infinity'
            error_vector(i) = norm(f(:)-S_NMf(:), inf);
    end

end

fprintf('\n');

figure(2); loglog(k, error_vector, '-+'); hold on;
normFac1 = 0.75*error_vector(1)/( k(1)^-0.5 );
normFac2 = 0.75*error_vector(1)/( k(1)^-1 );
normFac3 = 0.75*error_vector(1)/( k(1)^-2 );
loglog(k, normFac1*k.^-0.5, '--', k, normFac2*k.^-1, ':', ...
                                  k, normFac3*k.^-2, '-.')
legend('Error', '$\mathcal O(N^{-0.5})$', '$\mathcal O(N^{-1})$', ...
                '$\mathcal O(N^{-2})$')
axis([1e1 2e3 1e-6 1e-1]); grid