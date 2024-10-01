%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                           FIGURE PARAMETERS
%
% Antonio Figueroa Durán
% anfig@elektro.dtu.dk
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INTERPRETER
set(0,'DefaultTextInterpreter', 'latex')
set(0,'DefaultLegendInterpreter', 'latex')
set(0,'DefaultAxesTickLabelInterpreter', 'latex')
set(0,'DefaultColorbarTickLabelInterpreter', 'latex')

%% PARAMETERS DEFINITION
% Line properties
params.plot.linewidth = 4;
params.plot.markersize = 5;

set(0,'DefaultLineLineWidth',params.plot.linewidth)
set(0,'DefaultLineMarkerSize',params.plot.markersize)

% Color
params.plot.nColors = 5;
params.plot.colormap = cbrewer2('qual', 'Set1', params.plot.nColors);
params.plot.nColors2 = 11;
params.plot.colormap2 = cbrewer2('qual', 'Paired', params.plot.nColors2);
params.plot.nColors3 = 11;
params.plot.colormap3 = cbrewer2('div', 'Spectral', params.plot.nColors3);
params.plot.markers = {'square','^', 'v', 'o', 'd'};

% Axis properties
params.plot.axis.fontsize = 20;%28;

% Legend properties
params.plot.legend.fontsize = 20;%28;
params.plot.legend.location = 'best';
params.plot.legend.edgecolor = 'w';

% Ticks properties
params.plot.tick.TickLabelInterpreter = 'latex';

% Default cm width of IEEE two-column papers
params.plot.size.onecolumnwidth = 8.5985;
params.plot.size.twocolumnwidth = 17.797;

%% OVERWRITE STRUCTURE IN CASE OF CHANGE
[filePath,~,~] = fileparts(mfilename('fullpath'));
save([filePath,'/plotParams'],'params');

clear params