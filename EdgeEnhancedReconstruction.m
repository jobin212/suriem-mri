function S_Nf_edge = EdgeEnhancedReconstruction(fhat, ...
                            jmpvals, jmplocs, varargin)

%% Compute Edge Enhanced Fourier Reconstructions
% Use jump locations and jump values to estimate Fourier coefficients. Use 
% these to compute enhanced Fourier reconstructions.
%
% Usage:    S_Nf_edge = EdgeEnhancedReconstruction(fhat, jmpvals, jmplocs)
%
%   Inputs:
%       fhat        - the Fourier coefficients (complex vector)
%       jmpvals     - size of jump in function (real vector)
%       jmplocs     - locations of jumps in function (real vector)
%
%   Outputs:
%       S_Nf_edge   - edge enhanced Fourier reconstruction (real vector)
%

% Initialization
% Convert to a column vector
fhat = fhat(:);

% No. of Fourier modes
M = length(fhat);           % note: we assume there are 2N+1 total 
                            % coefficients f^_{-N},.....,f^_N
N = (M-1)/2;
% Fourier modes
k = (-N:N).';

% Optional arguments
if( nargin>3 )
    ngrid = varargin{1};
end


% Compute estimated Fourier coefficients
fhat_est = zeros(2*N+1, 1);
% This is a sum over the total number of jumps
for ix = 1:length(jmplocs)
    fhat_est = fhat_est + jmpvals(ix)*exp(-1i*k*jmplocs(ix))./(2i*pi*k);
end
fhat_est(k==0) = 0;             % for k==0

% Define a mean-shifted ramp function
ramp = @(y) (pi-y)/(2*pi) + floor( y/(2*pi) );

% Fourier reconstruction
if( nargin>3 )
    [S_Nf_edge, x] = ComputeFourierReconstruction(fhat-fhat_est, ngrid);
else
    [S_Nf_edge, x] = ComputeFourierReconstruction(fhat-fhat_est);
end


% Edge-based correction term
edge_correction = zeros(length(x), 1);
for ix = 1:length(jmplocs)
    edge_correction = edge_correction + jmpvals(ix)*ramp(x-jmplocs(ix));
end
S_Nf_edge = S_Nf_edge + edge_correction;

return