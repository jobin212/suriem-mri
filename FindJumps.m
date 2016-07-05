function [jmp_loc, jmp_ht] = FindJumps(fk, varargin)

%% Find Jumps given Fourier Coefficients
% Script to find jump locations and heights in a piecewise-smooth function
% given Fourier coefficients
%
% Usage:    [jmp_loc, jmp_ht] = FindJumps(fk, varargin)
%
%   Inputs:
%       fk          - the Fourier coefficients (complex vector)
%
%   Optional Inputs:
%       alg         - algorithm to use (string)
%                       - prony  (Prony-based edge detection)
%                       - conc   (using concentration kernels)
%       refine      - refine jumps using Matlab's fsolve? (boolean)
%       ngrid       - no. of grid points to use when evaluating Fourier 
%                     partial sums (integer)
%       njumps      - for Prony-based method - estimated no. of jumps
%
%   Outputs:
%       jmp_loc     - estimated jump locations (real vector)
%       jmp_ht      - estimated jump heights (real vector)
%



%% Initialization

% Convert to a column vector
fk = fk(:);

% No. of Fourier modes
M = length(fk);           % note: we assume there are 2N+1 total 
                            % coefficients f^_{-N},.....,f^_N
N = (M-1)/2;
% Fourier modes
k = (-N:N).';

% Process optional arguments
if( nargin>=2 )
    alg = varargin{1};
else
    alg = 'prony';              % by default, Prony-based edge detection
end

if( nargin>=3 )
    refine = varargin{2};
else
    refine = false;             % by default, no refinement of jumps
end

if( nargin>=4 )
    ngrid = varargin{3};
else
    ngrid = 3*(2*N+1);          % by default, use 3x oversampled grid
end

if( nargin>=5 )
    njumps = varargin{4};
else
    njumps = 10;                % default estimate of no. of jumps
end


%% Jump Detection

switch lower(alg)
    
    % Prony-based jump detection
    case 'prony'
        
        % Multiply out the (2i pi k) factor
        d_fk      = (2i*pi*k) .* fk;

        % Jump parameters
        ndata  = 20;         % No. of data/coefficient measurements to use
        M      = 01;         % We will use Fourier modes -(N+M-1:N+M+ndata-2)


        % Set up linear prediction model
        LinPredictionMat = zeros(ndata, njumps);
        for ix = 1:ndata
            LinPredictionMat(ix,:) = ...
                            d_fk( M+(njumps-1)+(ix-1) : -1 : M+(ix-1) );
        end

        model_measurements = d_fk( M+njumps : M+njumps+ndata-1 );

        % Solve linear prediction model
        poly_cfs = LinPredictionMat\model_measurements;

        % Getting the roots of the associated polynomial
        poly_roots = roots([1; -poly_cfs]);

        % Get Jump Locations
        jmp_loc = unwrap( -angle(poly_roots) );
        jmp_loc = sort( jmp_loc );              % Sort

        
        % Jump Heights
        % Set up and solve linear system
        fit_cfs = (abs(k)>0.75*N);       % coefficients to use for data fit
        sysMat  = exp( -1i*k(fit_cfs)*jmp_loc.' );
        jmp_ht  = real( sysMat \ d_fk(fit_cfs) );
        
        
    % Concentration kernel based jump detection
    case 'conc'
        
        % Concentration factor
        sig = confac(k, 'Trig');
        
        % Compute jump approximation
        jmp_cfs = (1i*sign(k).*sig) .* fk;
        [jmp_apprx, x] = ComputeFourierReconstruction(jmp_cfs, ngrid);

        % Extract jump locations and values
        [~, jmp_loc] = findpeaks( abs(jmp_apprx), ...
                'MinPeakHeight', 0.5, 'MinPeakDistance', 05/(2*N+1));
        jmp_ht = jmp_apprx(jmp_loc);
        jmp_loc = x(jmp_loc);         % convert from index to absolute 
                                      % location

end


% Refine jump information?
if( refine )
    % We will use the high freq. coefs (k>=M) to fit data
%         M = 50;
    M = round(3*N/4);

    options = optimoptions('fsolve', ...
                            'Algorithm', 'levenberg-marquardt', ...
                            'StepTolerance', 1e-10, ...
                            'FunctionTolerance', 1e-10, ...
                            'MaxIterations', 1e3, ...
                            'MaxFunctionEvaluations', 1e3, ...
                            'Display','off' );
    [refined_jmps, rsdl, flag, dtls]  = fsolve( @(jmp_info) ...
                    FitJumps( jmp_info, k(abs(k)>=M), fk(abs(k)>=M) ), ...
                    [jmp_ht(:) jmp_loc(:)], options );
%         [refined_jmps, rsdl, flag, dtls]  = fsolve( @(jmp_info) ...
%                 FitJumpValues( jmp_info, k(abs(k)>=M), fk(abs(k)>=M), ...
%                     x(jmp_locs) ), jmp_ht(:), options );

    jmp_ht   = real( refined_jmps(:,1) );
    jmp_loc  = real( refined_jmps(:,2) );

end


return