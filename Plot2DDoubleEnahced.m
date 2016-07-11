clear;

%%%the amount of exponenents we want to take 
runs = 4;
ErrType = '2norm';
FncType = 'box';
ReconstructionType = 'box-prony-jumps';
ErrTitle = strcat(ErrType, ' ', FncType, ' ', ReconstructionType, ' 3*N+2');
%order of errors decreased

%fourier coefficients
k = 50*2.^(0:runs-1);

error_vector = zeros(size(k));

for i = 1:length(k)
    N = k(i);
    M = N;
    ix = 3*N + 2;
    kn = (-N:N).';
    km = (-M:M).';
    
    ngrid = 3*(2*N+1);
    mgrid = 3*(2*M+1);
    
    x = -pi + (2*pi/ngrid)*( 0:ngrid-1 ).';
    y = -pi + (2*pi/mgrid)*( 0:mgrid-1 ).';
    
    [fHat, fxy] = Get2DFourierCoefficients(FncType, N, M);
    
    %%take from  Compute2DFourierReconstruction)  
    if(strcmp(ReconstructionType, 'standard'))
        %outer product
        F1 = exp(1i * x * kn.');
        F2 = exp(1i * y * km.');

        temp = F1 * fHat * F2.';
        
        S_NMf = temp(ix, :).';

    else
        %%calculate    
        
        comp_exp = exp(1i * x(ix) * (-N:N));
        
        Edge_Enhanced_fHat = zeros(size(fHat));
        
        jmp_heights = [];
        jmp_locs = [];
        
        for jx = 1:length(fHat)
            jmp_heights = [0];
            jmp_locs = [0];

            if(abs(x(ix)) <= 1) 
                jmp_heights = [-1 1].';
                jmp_locs = [1 -1].';
            end;
            Edge_Enhanced_fHat = EdgeEnhancedReconstruction(fHat(jx, :), [jmp_heights], [jmp_locs]);
        end;
            
        
         CFR = comp_exp * Edge_Enhanced_fHat;
        
        

        
        switch(ReconstructionType)
            case('box-true-jumps')
                jmp_heights = [0];
                jmp_locs = [0];

                if(abs(x(ix)) <= 1) 
                    jmp_heights = [-1 1].';
                    jmp_locs = [1 -1].';
                end;
            case('circle-true-jumps')
                jmp_heights = [0];
                jmp_locs = [0];

                if( abs(x(ix)) <= 1) 
                    jmp_locs = [-1*sqrt(1- x(ix)^2) 1*sqrt(1-x(ix)^2)].';
                    jmp_heights = [1 -1].';
                end

            case('box-prony-jumps')

                jmps = 0;

                if( x(ix) <= 1 && x(ix) >= -1)
                    jmps = 2;
                end


                [jmp_locs, jmp_heights] = FindJumps(CFR, 'prony', false, [], jmps);

            case('circle-prony-jumps')
                jmps = 0;

                if( x(ix) <= 1 && x(ix) >= -1)
                    jmps = 2;
                end

                [jmp_locs, jmp_heights] = FindJumps(CFR, 'prony', false, [], jmps);                

            case('box-conc-jumps')
                [jmp_locs, jmp_heights] = FindJumps(CFR, 'conc', true);

            case('circle-conc-jumps')
                [jmp_locs, jmp_heights] = FindJumps(CFR, 'conc', true);
        end;  




        S_NMf =  EdgeEnhancedReconstruction(CFR, jmp_heights, jmp_locs);


    end 
    
    f = fxy(x, x(ix));
    abs_error = abs(f - S_NMf);
    
    
    %{
    figure;
    plot(x, S_NMf);
    hold on;
    plot(x, f);
    legend('reconstruction' ,'f');
    
    figure;
    plot(x, abs_error);
    title('abs_error');
    %}
    
        
    error_vector(i) = GetError(ErrType, abs_error, x);
    
end

figure;
loglog(k, error_vector)
hold on;
loglog(k, k.^(-1))
hold on;
loglog(k, k.^(-2))
legend(ErrTitle, 'k^{-1}', 'k^{-2}')
ylim([1e-4, 1e0]);