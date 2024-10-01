function applyLegendProperties(h)
%applyLegendProperties(leg) Reads axis configuration from plotParams.mat
%and applies to legend in plot leg
%   Input:
%       - h         : figure handler
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

% Check if params.plot.legend is a valid structure
if ~isfield(params.plot,'legend'), return, end
if ~isstruct(params.plot.legend), return, end

leg = findobj(h,'type','Legend');

% Get the parameters from the params.plot.legend structure
fieldName = fieldnames(params.plot.legend);

for p_idx = 1:length(fieldName)
    s = fieldName{p_idx};
    if isprop(leg,s), set(leg, s, params.plot.legend.(s))
    else, warning('Property "%s" not found for legend configuration', s)
    end 
end

end