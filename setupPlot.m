function setupPlot(Data,flagS,flagH,flagT)
%setupPlot(Data,flagS,flagH,flagT) Plots the measurement setup, the
%frequency response and the RIR at the reference microphones.
%   Input:
%       - Data      : data structure.
%       - flagS     : plots setup
%                       'true'
%                       'false' (Default value)
%       - flagH     : plots frequency response
%                       'true'
%                       'false' (Default value)
%       - flagT     : plots RIR at reference line
%                       'true'
%                       'false' (Default value)
% Author: Antonio Figueroa Durán
% Date: August 2023

%% ERROR HANDLING
% plotFlag default value
if nargin < 1, error('setupPlot Error: Not enough input parameters.'), end
if nargin < 4, flagT = false; end
if nargin < 3, flagH = false; end
if nargin < 2, flagS = false; end

%% MAIN CODE
% Time vector
T = [15 40]*1e-3;

%% PLOT: SETUP
if flagS
    figure
    scatter3(Data.Ref.pos(:,1),Data.Ref.pos(:,2),Data.Ref.pos(:,3)), hold on
    scatter3(reshape(Data.Mic.pos(:,1,:),[],1),reshape(Data.Mic.pos(:,2,:),[],1),reshape(Data.Mic.pos(:,3,:),[],1))
    scatter3(Data.Source.pos(1),Data.Source.pos(2),Data.Source.pos(3),200,'filled',"MarkerFaceColor",'#7E2F8E')
    NBIRoom(), axis tight
    xlabel('$x$ / m'), ylabel('$y$ / m'), zlabel('$z$ / m')
    legend('Reference Line','Eigenmike','Source')
    applyAxisProperties(gca)
    applyLegendProperties(gcf)
end

%% PLOT: FREQUENCY RESPONSE
if flagH, plotFreqResponse(Data), end

%% PLOT: REFERENCE RIR
if flagT
    idx = Data.t > T(1) & Data.t < T(2);

    RIRPlot = Data.Ref.h(idx,:)./max(abs(Data.Ref.h(idx,:)),[],'all');
    figure, hold on
    s = pcolor(Data.Ref.pos(:,1),Data.t(idx)*1e3,RIRPlot);
    set(s,'edgecolor','none')
    xlabel('$x$ / m'), ylabel('Time / ms')
    axis([0 Data.D(1) T(1)*1e3 T(2)*1e3])
    DR = colormapDR(RIRPlot);
    colorbarpwn(-DR,DR)
    c = colorbar;
    applyColorbarProperties(c,'Normalised $h(t)$')
    applyAxisProperties(gca)
end
end

