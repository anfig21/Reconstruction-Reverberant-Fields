function [idx,fNyq] = nearestNeighboursSph(micPos,M,c)
%[idx,fNyq] = nearestNeighboursSph(micPos,M,c) Resamples the spherical
%microphone array described by micPos taking M positions pseudo-uniformly
%using a nearest neighbour algorithm. Microphones inside the sphere are
%kept.
%   Input:
%       - micPos    : microphone positions of the spherical array. nMic x 3
%       - M         : number of samples. Scalar
%       - c         : speed of sound. Scalar
%
%   Output:
%       - idx       : indices of the optimised subsampling. 1 x M
%       - fNyq      : Nyquist frequency for uniform sampling. Scalar
%
% Author: Antonio Figueroa Durán
% Date: November 2022

%% ERROR HANDLING
N = size(micPos,1);
if M > N, M = N; warning('nearestNeighboursSph Warning: Reached maximum number of microphones.'), end
if M < 1, error('nearestNeighboursSph Error: M must be a positive integer.'), end

%% MAIN CODE
r0 = mean(micPos);
micPos = micPos-r0;
radius = mean(vecnorm(micPos,2,2));

extPos = micPos(vecnorm(micPos,2,2) > radius*0.9,:);
intPos = micPos(vecnorm(micPos,2,2) < radius*0.9,:);

radius = mean(vecnorm(extPos,2,2));

uk = radius*fibonacciSampling(M-size(intPos,1))';

idxExt = nearestNeighbor(delaunayTriangulation(extPos),uk);
idxInt = find(ismember(micPos,intPos,'rows'));

idx = [idxExt; idxInt];

% Nyquist frequency
Nambisonics = floor(sqrt(numel(idxExt))-1);
k = Nambisonics/radius;
fNyq = c*k/(2*pi);

end

