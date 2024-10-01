function applyColorbarProperties(c,label)
%applyColorbarProperties(c,label)) Applies label configuration to colorbar
%axis
%   Input:
%       - c         : colorbar handler
%       - label     : colorbar label. String
%   Output:
%       none
%
% Author: Antonio Figueroa Durán
% Date: February 2022

c.Label.String = label;
c.Label.Interpreter = 'Latex';

end

