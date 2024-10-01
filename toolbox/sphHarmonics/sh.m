function Y = sh(Omega,Nmax,typeBasis,normalisation)
%Y = sh(Omega,Nmax,typeBasis,normalisation) Obtain the spherical harmonics
%up to order Nmax.
%   Input:
%       - Omega         : angular directions (theta,phi) in rad. Theta is
%                           defined from zenith. Phi is defined from x-axis
%                           towards y-axis. nDirections x 2
%       - Nmax          : maximum SH order. Scalar
%       - typeBasis     : 'complex' or 'real' Spherical Harmonics
%       - normalisation : Legendre polynomials normalisation:
%                           'unnorm' - Associated Legendre Functions
%                                      (default)
%                           'sch' - Schmidt Seminormalisation
%                           'norm' - Normalisation according to L_2 space
%   Output:
%       - Y             : spherical harmonics. nDirections x (Nmax+1)^2
%
% Author: Antonio Figueroa Durán
% Date: July 2023

%% ERROR HANDLING
arguments
    Omega (:,2) double {mustBeNonempty}
    Nmax (1,1) {mustBeInteger,mustBeNonnegative}
    typeBasis char {mustBeMember(typeBasis,{'complex','real'})} = 'complex'
    normalisation char {mustBeMember(normalisation,{'unnorm','sch','norm'})} = 'unnorm'
end

%% MAIN CODE
% Dimensions
Nangles = size(Omega,1);
Nsignals = (Nmax+1)^2;

% Extract coordinates
theta = Omega(:,1);
phi = Omega(:,2);

% Memory allocation
Y = nan(Nangles,Nsignals);
for n = 0:Nmax  % Order
    m = 0:n;    % Degree

    % Legendre polynomials
    Leg = legendre(n,cos(theta),normalisation);
    
    % Spherical Harmonics
    Y(:,(n:2*n)+n^2+1) = sqrt((2*n+1).*factorial(n-m)./(4*pi*factorial(n+m))).*Leg.'.*exp(1i*m.*phi);   % m = 0:n
    if n ~= 0, Y(:,(1:n)+n^2) = flip((-1).^m(2:end).*conj(Y(:,(n+1:2*n)+n^2+1))); end                   % m = -n:-1

    if strcmp(typeBasis,'real') && numel(m) > 1
        Y(:,(n+1:2*n)+n^2+1) = sqrt(2)*(-1).^m(2:end).*real(Y(:,(n+1:2*n)+n^2+1));
        Y(:,(1:n)+n^2) = sqrt(2)*(-1).^m(end:-1:2).*imag(Y(:,(1:n)+n^2));
    end
end

end