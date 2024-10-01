function range = earlyRange(Data,Early,plotFlag)
%range = earlyRange(Data,Early,plotFlag) Estimates the Range of the 
%apparent origin of the corresponding reflection based on DOA
%triangulation.
%   Input:
%       - Data          : raw data. Structure
%       - Early         : early part of the RIR. Structure
%       - plotFlag      : 'true' to plot Range estimation
%                         'false' (Default value)
%   Output:
%       - range         : range estimation. 3 x 1
%
% Author: Antonio Figueroa Durán
% Date: August 2023

%% ERROR HANDLING
if nargin < 2, error('earlyRange Error: Not enough input parameters.'), end
if ~isfield(Early,'DOA'), error('earlyRange Error: DOA has not been estimated.'), end
if nargin < 3, plotFlag = false; end

%% MAIN CODE
% Direction of line connecting the closest points
Idx = nchoosek(1:Data.Mic.nArrays,2);
Idx_i = Idx(:,1);
Idx_j = Idx(:,2);
n = cross(Early.DOA.Mode(:,Idx_i),Early.DOA.Mode(:,Idx_j));

range_ij = nan(3,numel(Idx_i));
for ij = 1:numel(Idx_i)
    n_ii = n(:,ij);

    % Parameter for closest points on each line
    k_i = dot(cross(Early.DOA.Mode(:,Idx_j(ij)),n_ii),(Data.Mic.R0(Idx_j(ij),:)-Data.Mic.R0(Idx_i(ij),:)))/dot(n_ii,n_ii);
    k_j = dot(cross(Early.DOA.Mode(:,Idx_i(ij)),n_ii),(Data.Mic.R0(Idx_j(ij),:)-Data.Mic.R0(Idx_i(ij),:)))/dot(n_ii,n_ii);

    % Closest points on each line
    y_iClosest = Early.DOA.Mode(:,Idx_i(ij))*k_i+Data.Mic.R0(Idx_i(ij),:)';
    y_jClosest = Early.DOA.Mode(:,Idx_j(ij))*k_j+Data.Mic.R0(Idx_j(ij),:)';

    % Source estimation
    range_ij(:,ij) = (y_iClosest+y_jClosest)/2;
end

% Remove outliers using Mean Absolute Deviation & obtain centroid
centroid = median(range_ij,2,'omitnan');
distMatrix = vecnorm(range_ij-centroid);
finalIdx = distMatrix < 1.5*mad(distMatrix);

range = median(range_ij(:,finalIdx),2,'omitnan');

%% PLOT
if plotFlag
    figure
    scatter3(Data.Mic.R0(:,1),Data.Mic.R0(:,2),Data.Mic.R0(:,3)), hold on   % em32
    scatter3(Data.Source.pos(1),Data.Source.pos(2),Data.Source.pos(3),200,'filled')
    scatter3(range(1),range(2),range(3),150,'filled')
    scatter3(range_ij(1,:),range_ij(2,:),range_ij(3,:),50,'filled')
    drawRoom(Data.D(1),Data.D(2),Data.D(3)), axis equal
    xlabel('$x$ / m'), ylabel('$y$ / m'), zlabel('$z$ / m')
    legend('Eigenmike','Source','Range Estimation','Estimation Candidates')
    applyAxisProperties(gca)
    applyLegendProperties(gcf)
end

end

