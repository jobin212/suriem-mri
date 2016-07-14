%get coefficients
%discretize
%cross section near edges and center
FncType = 'box';
N = 100;
sN = num2str(N);
M = N;
sM = num2str(M);
leg = strcat('N=M=',sN); 

ReconstructionType = 'box-true-jumps';

[fHat, fxy] = Get2DFourierCoefficients(FncType, N, M);


kn = (-N:N).';
km = (-M:M).';

ngrid = 3*(2*N+1);
mgrid = 3*(2*M+1);

x = -pi + (2*pi/ngrid)*( 0:ngrid-1 ).';
y = -pi + (2*pi/mgrid)*( 0:mgrid-1 ).';
z = -pi + (2*pi/(2*N+1))*( 0:(2*N)).';

jmp_heights = [];
jmp_locs = [];

CFR = zeros(2*N+1, 1);

S_NMf = zeros(length(x), length(y));

for ix = 1:length(x)
    for jx = 1:length(fHat)
        jmp_heights = [0];
        jmp_locs = [0];

        if(abs(z(jx)) < 1) 
            %jmp_heights = [-1 1].';
            p = jx - M - 1;
            jmp_heights = (sin(p) / (pi * p)) * ([-1 1]).';
            if(p == 0 ) 
                jmp_heights = [-1/pi 1/pi];
            end;
            jmp_locs = [1 -1].';
        end;

        Edge_Enhanced_fHat = EdgeEnhancedReconstruction(fHat(:, jx), jmp_heights, jmp_locs);
        CFR(jx) = Edge_Enhanced_fHat(ix);
    end;

    comp_exp = exp(1i * x(ix) * (-N:N));

    %CFR = comp_exp * fHat;
    
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
    end;

    S_NMf(:, ix) =  EdgeEnhancedReconstruction(CFR, jmp_heights, jmp_locs);

    %F1 = exp(1i * x * kn.');
    %S_NMf = F1 * CFR.';
end;

[xx, yy] = meshgrid(x,y);

abs_error = abs(S_NMf - fxy(xx,yy));


mesh(xx, yy, real(S_NMf));
legend(leg);
xlabel('x')
ylabel('y')

figure;
mesh(xx, yy, abs_error);
legend(leg);





