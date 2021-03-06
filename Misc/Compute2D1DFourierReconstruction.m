function [S_NMf, x, y] = Compute2D1DFourierReconstruction(fHat, ReconstructionType, varargin)                  

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
if( nargin == 3 )           % reconstruction grid length 
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

%%calculate
for ix = 1:length(x)
    comp_exp = exp(1i * x(ix) * (-N:N));
    CFR = comp_exp * fHat;
    %S_NMf(:, ix) = ComputeFourierReconstruction(CFR);
    
    jmp_heights;
    jmp_locs;
    
    %for box!
    switch(ReconstructionType)
        case('true-jumps')
            jmp_heights = [0];
            jmp_locs = [0];

            if(abs(x(ix)) <= 1) 
                jmp_heights = [-1 1].';
                jmp_locs = [1 -1].';
            end;
        case('prony-jumps')
            jmps = 0;
            if( x(ix) <= 1 && x(ix) >= -1)
                jmps = 2;
            end
            

            [jmp_heights, jmp_locs] = FindJumps(CFR, 'prony', false, [], jmps);
            
        case('conc-jumps')
            [jmp_heights, jmp_locs] = FindJumps(CFR, 'conc', false);
    end;  
        
    S_NMf(:, ix) = EdgeEnhancedReconstruction(CFR, jmp_heights, jmp_locs);       
    
end

% We will mainly consider only real functions. The following ensures tiny 
% complex values due to roundoff errors and such are not included
S_NMf = real(S_NMf);


return






%{
%%calculate
for ix = 1:length(x)
    for iy = 1:length(y)
        
        
        for ik = -N:N
            for il = -M:M
                
                S_NMf(ix , iy) = S_NMf(ix, iy) + fHat(ik + N + 1, il + M + 1) * exp(1i * ik * x(ix) + 1i * il * y(iy));
            end
        end
    end
end
%}
