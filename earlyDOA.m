function DOA = earlyDOA(Data,Early,f,N,plotFlag,micIdx)
%Direct = earlyDOA(Data,Early,f,plotFlag) Estimate the Direction-of-Arrival
%(DOA) of an individual reflection using SRP-PHAT.
%   Input:
%       - Data      : raw data. Structure
%       - Early     : structure with RIR of dim. [nSamples x nArrays*32]
%       - f         : frequency span. 1 x Nf
%       - plotFlag  : 'true' to plot DOA estimation
%                     'false' (Default value)
%   Output:
%       - Early     : DOA Estimation. Structure
%
% Author: Antonio Figueroa-Duran
% Date: August 2023

%% ERROR HANDLING
if nargin < 3, error('directSoundDOA Error: Not enough input parameters.'), end
if nargin < 5, plotFlag = false; end
if nargin < 6, micIdx = []; end

%% MAIN CODE
if any(micIdx)
    SphH_ii = Early.Mic.H(ismember(Data.f,f),:);

    % DOA Estimation via SRP-PHAT
    DOA = dirDOA_SRP_PHAT(SphH_ii,Data.Mic.pos(:,:,micIdx),f,Data.Fs,Data.c,N);
else
    DOA = nan(3,Data.Mic.nArrays);
    for ii = 1:Data.Mic.nArrays
        SphH_ii = Early.Mic.H(ismember(Data.f,f),:,ii);

        % DOA Estimation via SRP-PHAT
        DOA(:,ii) = dirDOA_SRP_PHAT(SphH_ii,Data.Mic.pos(:,:,ii),f,Data.Fs,Data.c,N);
    end
end

% Include frequency vector in structure
Early.DOA.f = f;

%% PLOT
if plotFlag
    % 3-D Estimation
    figure
    scatter3(Data.Mic.R0(:,1),Data.Mic.R0(:,2),Data.Mic.R0(:,3)), hold on   % em32
    scatter3(Data.Source.pos(1),Data.Source.pos(2),Data.Source.pos(3),200,'filled')
    quiver3(Data.Mic.R0(:,1),Data.Mic.R0(:,2),Data.Mic.R0(:,3),...
        DOA(1,:)',DOA(2,:)',DOA(3,:)',1.8,'Linewidth',4)

    drawRoom(Data.D(1),Data.D(2),Data.D(3)), axis equal
    xlabel('$x$ / m'), ylabel('$y$ / m'), zlabel('$z$ / m')
    legend('Eigenmike','Source','DOA Estimation')
    applyAxisProperties(gca)
    applyLegendProperties(gcf)
end

end

