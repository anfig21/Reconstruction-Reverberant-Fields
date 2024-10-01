function [H,uk] = dictionaryPW(c,f,r,N)
%[H,uk] = dictionaryPW(c,f,r0,N) Obtains the H matrix for a sound field at
% a position r0 comprised of N plane waves with propagation vector sampled
% along the surface of a unit sphere.
%   Input:
%       - c         : speed of sound (m/s). Scalar
%       - f         : frequency (Hz). Nf x 1
%       - r         : measurement point. 3 x M

%       - N         : number of plane waves. Scalar
%   Output:
%       - H         : dictionary. M x N x Nf
%       - uk        : propagation unit vector. 3 x N
%
% Author: Antonio Figueroa Durán
% Date: February 2022

%% ERROR HANDLING
if nargin < 4, error('dictionaryPW Error: Not enough input parameters.'), end

%% MAIN CODE
Nf = length(f);
[~,M] = size(r);

% Spatial sampling - Fibonacci Algorithm
uk = -fibonacciSampling(N);      % uk point inwards the sphere

% Propagation vector
k = (2*pi*f)/c;
kv = nan(3,N,Nf);
for ii = 1:Nf, kv(:,:,ii) = k(ii)*uk; end

% Dictionary
H = nan(M,N,Nf);
% if Nf > 1, d = waitbar(0,'Loading...0\%','Name','Building dictionary...'); end
for ii = 1:Nf
    H(:,:,ii) = exp(-1i*r'*squeeze(kv(:,:,ii)));
    % if Nf > 1, waitbar(ii/Nf,d,strcat("Loading... ",string(round(100*ii/Nf,2)),"\,\%")); end
end
% if Nf > 1, delete(d), end

% Normalisation
H = H./vecnorm(H,2,2);

end