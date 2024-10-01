function [] = drawRoom(lx,ly,lz)
%[] = drawRoom(lx,ly,lz) Draws the room walls
%   Input:
%       - lx        : x dimension
%       - ly        : y dimension
%       - lz        : z dimension
%
% Author: Antonio Figueroa Durán
% Date: April 2022

%% ERROR HANDLING
if nargin < 3, error('drawRoom Error: Not enough input parameters.'), end
if lx < 0 || ly < 0 || lz < 0, error('drawRoom Error: Room dimensions need to be positive scalars.'), end

%% MAIN CODE
x = lx*[0 0 0 0 1 0; 0 0 1 1 1 0; 0 1 1 1 1 1; 0 1 0 0 1 1];
y = ly*[1 0 1 1 1 1; 0 0 1 1 0 1; 0 0 0 0 0 1; 1 0 0 0 1 1];
z = lz*[1 0 1 0 1 0; 1 1 1 0 1 1; 0 1 1 0 0 1; 0 0 1 0 0 0];

fill3(x,y,z,'k','FaceAlpha',5e-2)

end

