function DR = colormapDR(cData)
%DR = colormapDR(cData) Obtaines the dynamic range for a 2D color plot
%   Input:
%       - cData     : data. 2D matrix
%   Output:
%       - DR        : dynamic range. Scalar
%
% Author: Antonio Figueroa Durán
% Date: February 2023

%% ERROR HANDLING
if imag(cData) ~= 0, error('colormapDR Error: Input data must be real.'), end

%% MAIN CODE
minV = min(cData,[],'all');
maxV = max(cData,[],'all');

% Dinamic Range
DR = max(abs([minV maxV]));

end

