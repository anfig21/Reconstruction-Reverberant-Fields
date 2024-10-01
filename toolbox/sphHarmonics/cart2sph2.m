function [R,THETA,PHI] = cart2sph2(X,Y,Z)
%[R,THETA,PHI] = cart2sph2(X,Y,Z) Transform cartesian to spherical
%coordinates. 
%   Input:
%       - X         : x dimension. N x 1
%       - Y         : y dimension. N x 1
%       - Z         : z dimension. N x 1
%   Outputs:
%       - R         : distance to origin. N x 1
%       - THETA     : elevation (from zenith). N x 1
%       - PHI       : azimuth (from X-axis). N x 1
%
% Author: Antonio Figueroa Durán
% Date: June 2023

%% ERROR HANDLING
if nargin < 3, error('cart2sph2 Error: Not enough input parameters.'), end
if ~isequal(size(X),size(Y),size(Z))
    error('cart2sph2 Error: Input dimensions not matching.')
end

%% MAIN CODE
X = X(:);
Y = Y(:);
Z = Z(:);

R = sqrt(X.^2 + Y.^2 + Z.^2);
THETA = acos(Z./R);
% PHI = sign(Y).*acos(X./sqrt(X.^2 + Y.^2));
PHI = (sign(Y)<0)*2*pi + atan2(Y,X);

end