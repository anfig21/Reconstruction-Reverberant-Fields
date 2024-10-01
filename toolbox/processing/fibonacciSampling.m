function uk = fibonacciSampling(N)
%uk = fibonacciSampling(N) Spherical uniform sampling along the surface of
%a sphere using the Fibonacci ratio.
%Ref:
%   [1] Gonzalez, Alvaro. “Measurement of Areas on a Sphere Using Fibonacci
%       and Latitude-Longitude Lattices”, International Association for
%       Mathematical Geosciences, 2009.
%   Input:
%       - N         : number of plane waves. Scalar
%   Output:
%       - uk        : propagation unit vector. 3 x N
%
% Author: Antonio Figueroa Durán
% Date: March 2022

%% ERROR HANDLING
if nargin < 1, error('fibonacciSampling Error: Not enough input parameters.'), end

%% MAIN CODE
% N must be odd
even = false;
if mod(N,2) == 0, N = N+1; even = true; end
P = (N-1)/2;

% Fibonacci golden ratio
phi = (1 + sqrt(5))/2;

% Uniform sphere sampling algorithm
lat = nan(1,N);
lon = nan(1,N);
i = -P:P;
for ii = 1:N
    lat(ii) = asin(2*i(ii)/N)+pi/2;
    lon(ii) = 2*pi*mod(i(ii),phi)/phi;
    if lon(ii) < -pi, lon(ii) = lon(ii) + 2*pi; end
    if lon(ii) > pi, lon(ii) = lon(ii) - 2*pi; end
end

% Cartesian coordinates
z = cos(lat);
x = cos(lon).*sin(lat);
y = sin(lon).*sin(lat);

% Propagation unit vector
uk = [x; y; z];
if even, uk = uk(:,1:end-1); end

end

