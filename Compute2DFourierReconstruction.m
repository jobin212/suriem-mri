function [S_NMf, x, y] = Compute2DFourierReconstruction(fHat, ReconstructionType, varargin)                  

%% Compute Fourier Partial Sum Approximation
% Script to compute the Fourier partial sum approximation of a function
%
% Usage:    S_NMf = ComputeFourierReconstruction(fHat)
%
%   Inputs:
%       fHat        - the Fourier coefficients (complex vector)
%
%   Optional Inputs:
%       ngrid       - number of grid points to compute reconstruction on
%                     (integer)
%       mgrid       - number of grid points to compute reconstruction on 
%
%   Outputs:
%       S_NMf        - the Fourier partial sum approximation (real vector)
%       x           - the grid on which the approximation is computed
%                     (real vector)
%



% No. of Fourier modes
NM = size(fHat);           % note: we assume there are 2N+1 total 
                            % coefficients f^_{-N},.....,f^_N
N = (NM(1)-1)/2;
M = (NM(2)-1)/2;

% Fourier modes
kn = (-N:N).';
km = (-M:M).';

% Process optional arguments
if( nargin == 4 )           % reconstruction grid length 
    ngrid = varargin{1};
    mgrid = varargin{2};
else
    ngrid = 3*(2*N+1);      % reconstruct on 3x oversampled grid
    mgrid = 3*(2*M+1);
end


% Generate equispaced grid in [-pi, pi)
x = -pi + (2*pi/ngrid)*( 0:ngrid-1 ).';
y = -pi + (2*pi/mgrid)*( 0:mgrid-1 ).';

S_NMf = zeros(length(x), length(y));


if(strcmp(ReconstructionType, 'standard'))
    %outer product
    F1 = exp(1i * x * kn.');
    F2 = exp(1i * y * km.');

    S_NMf = F1 * fHat * F2.';

else
    %%calculate
    for ix = 1:length(x)
        %%use edge enhanced 
        comp_exp = exp(1i * x(ix) * (-N:N));
        CFR = comp_exp * fHat;
        
        
        %S_NMf(:, ix) = ComputeFourierReconstruction(CFR);

        jmp_heights = [];
        jmp_locs = [];

        %for box!
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
                [jmp_locs, jmp_heights] = FindJumps(CFR, 'conc', false);
                
            case('circle-conc-jumps')
                [jmp_locs, jmp_heights] = FindJumps(CFR, 'conc', false);
        end;  
        
        
        

        S_NMf(:, ix) = EdgeEnhancedReconstruction(CFR, jmp_heights, jmp_locs);       

    end
end


% We will mainly consider only real functions. The following ensures tiny 
% complex values due to roundoff errors and such are not included
S_NMf = real(S_NMf);


return
