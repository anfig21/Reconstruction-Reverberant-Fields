function vertices = genTetrahedron()
%vertices = genTetrahedron() Generates coordinates for a regular
%tetrahedron of side length 1
%   Output:
%       - vertices   : coordinates of the vertices. 4 x 3
%
% Author: Antonio Figueroa Durán
% Date: June 2023

%% MAIN CODE
vertices = [1 0 -1/sqrt(2);
    -1 0 -1/sqrt(2);
    0 1 1/sqrt(2);
    0 -1 1/sqrt(2)]/2;

end