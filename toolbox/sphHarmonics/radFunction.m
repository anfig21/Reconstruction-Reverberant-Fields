function h = radFunction(k,r,Nmax,config)
%h = radFunction(k,r,Nmax,config) Obtain the spherical Hankel function
%depending on the spherical microphone array configuration.
%   Input:
%       - k         : wavenumber. nF x 1
%       - r         : distance in m. nR x 1
%       - Nmax      : maximum SH order. Scalar
%       - config    : microphone array configuration:
%                       'open' - Open configuration. Returns Hankel
%                           function of the first kind (Bessel)
%                       'closed' - Closed configuration. Returns Hankel
%                           function of the third kind (Hankel)
%   Output:
%       - h         : radial functions. nF x nR x Nmax
%
% Author: Antonio Figueroa Durán
% Date: July 2023

%% ERROR HANDLING
arguments
    k (:,1) double {mustBeNonempty,mustBeNonnegative}
    r (1,:) double {mustBePositive}
    Nmax (1,1) {mustBeInteger,mustBeNonnegative}
    config char {mustBeMember(config,{'open','closed'})}
end

%% MAIN CODE
kr = k*r;

nF = length(k);
nR = length(r);

h = nan(nF,nR,Nmax);
for nn = 0:Nmax
    % Cylindrical Bessel function (Hankel function of the first kind)
    h(:,:,nn+1) = besselj(nn+0.5,kr);
    
    if strcmp(config,'closed')
        % Closed array configuration:
        %   Neumann function (Hankel function of the second kind)
        h(:,:,nn+1) = h(:,:,nn+1) - 1i*bessely(nn+0.5,kr);
    end
end

% Spherical Hankel function
h = h.*sqrt(pi./(2*kr));
h(isnan(h)) = 0;
h(1,1) = 1;
h = squeeze(h);

end