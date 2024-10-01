function [X,Y,Z] = sph2cart2(R,THETA,PHI)
%[X,Y,Z] = sph2cart2(R,THETA,PHI) Transform spherical to cartesian
%coordinates. 
%   Input:
%       - R         : distance to origin. N x 1
%       - THETA     : elevation (from zenith). N x 1
%       - PHI       : azimuth (from X-axis). N x 1
%   Outputs:
%       - X         : x dimension. N x 1
%       - Y         : y dimension. N x 1
%       - Z         : z dimension. N x 1
%
% Author: Antonio Figueroa Durán
% Date: June 2023

%% ERROR HANDLING
if nargin < 3, error('sph2cart2 Error: Not enough input parameters.'), end
if ~isequal(size(R),size(THETA),size(PHI))
    error('sph2cart2 Error: Input dimensions not matching.')
end

%% MAIN CODE
R = R(:);
THETA = THETA(:);
PHI = PHI(:);

X = R.*sin(THETA).*cos(PHI);
Y = R.*sin(THETA).*sin(PHI);
Z = R.*cos(THETA);

end