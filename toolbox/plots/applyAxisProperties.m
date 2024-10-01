function applyAxisProperties(ax)
%applyAxisProperties(ax) Reads axis configuration from plotParams.mat and
%applies to axis in plot ax
%   Input:
%       - ax        : plot handler
%   Output:
%       none
%
% Author: Antonio Figueroa Durán
% Date: February 2022

% Check parameters existence
[filePath,~,~] = fileparts(mfilename('fullpath'));
if ~isfile([filePath,'/plotParams.mat'])
    warning(['Axis properties have not been defined. ' ...
        'Please, run loadPlotParams.m'])
    return
end

load([filePath,'/plotParams.mat'],'params');

% Check if params.plot.axis is a valid structure
if ~isfield(params.plot,'axis'), return, end
if ~isstruct(params.plot.axis), return, end

% Get the parameters from the params.plot.axis structure
fieldName = fieldnames(params.plot.axis);

for p_idx = 1:length(fieldName)
    s = fieldName{p_idx};
    if isprop(ax,s), set(ax, s, params.plot.axis.(s))
    else, warning('Property "%s" not found for axis configuration', s)
    end 
end

end