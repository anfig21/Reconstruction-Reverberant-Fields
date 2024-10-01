function Rec = clusterCoefficients(Data,Rec,R0,plotFlag)
%Rec = clusterCoefficients(Data,Rec,R0,plotFlag) Obtain the coefficients
%for a cluster of point sources.
%   Input:
%       - Data          : raw data. Structure
%       - Rec           : windowed RIR section to reconstruct. Structure
%       - R0            : estimated position of the source. 3 x 1
%       - plotFlag      : 'true' to plot coefficients
%                         'false' (Default value)
%   Output:
%       - Rec           : windowed RIR section with coefficients. Structure
%
% Author: Antonio Figueroa Durán
% Date: February 2024

%% ERROR HANDLING
if nargin < 4, plotFlag = false;    % plotFlag default value
elseif nargin < 3, error('clusterCoefficients Error: Not enough input parameters.'), end

%% MAIN CODE
% Combine data from arrays
pos = reshape(permute(Data.Mic.pos,[2 1 3]),[],Data.Mic.nM*Data.Mic.nArrays)';
P = Rec.Mic.H(ismember(Data.f,Rec.f),:,:);
P = reshape(P,[],Data.Mic.nM*Data.Mic.nArrays).';

% Tetrahedron of edge length 15 cm
L = 15e-2;      % Tetrahedron edge length (Loudspeaker dimensions: 21x34x32 cm)
Rec.rs = L*genTetrahedron+R0(:)';

% Plot cluster geometry
% figure, scatter3(Rec.rs(:,1),Rec.rs(:,2),Rec.rs(:,3)), axis equal

N = size(Rec.rs,1);             % Number of point sources
Nf = length(Rec.f);             % Number of frequency bins

k = (2*pi*Rec.f)/Data.c;        % Propagation vector

% Euclidean distance matrix
d = pdist2(pos,Rec.rs);         % mics x point sources

% Coefficient estimation
x = zeros(N,Nf);                % Allocate memory
Res = nan(size(P,1),Nf);
for ii = 1:Nf
    H = (1./d).*exp(-1i*d*k(ii));       % Dictionary (point sources)
    r = P(:,ii);                        % Initialise residual

    % Sequential LS
    for jj = 1:N
        x(jj,ii) = H(:,jj)\r;                   % LS for coefficient estimation
        r = P(:,ii) - H(:,1:jj)*x(1:jj,ii);     % Update residual
    end

    % Reconstruction at measurement positions
    Res(:,ii) = H*x(:,ii);
end

% Save coefficients
Rec.x = zeros(N,Data.Nsamples/2);
Rec.x(:,ismember(Data.f,Rec.f)) = x;

% Return reconstruction at measurement positions
Rec.Rec.H = zeros(size(Rec.Mic.H));
Rec.Rec.H(ismember(Data.f,Rec.f),:,:) = reshape(Res.',[],Data.Mic.nM,Data.Mic.nArrays);

%% PLOT
% Coefficients
if plotFlag
    figure
    plot(Data.f,abs(Rec.x)'), grid on
    xlabel('Frequency / Hz'), ylabel('$|\hat{x}|$')
    applyAxisProperties(gca)

    figure
    scatter3(Rec.rs(:,1),Rec.rs(:,2),Rec.rs(:,3),[],log10(sum(abs(Rec.x),2)),'filled'), hold on
    scatter3(R0(1),R0(2),R0(3),'filled')
    if size(Rec.rs,1) < 15, text(Rec.rs(:,1),Rec.rs(:,2),Rec.rs(:,3),num2str((1:size(Rec.rs,1))'),"FontSize",18), end
    xlabel('x / m'), ylabel('y / m'), zlabel('z / m')
    colorbar
    axis equal
    axis([min(Rec.rs(:,1))-0.5 max(Rec.rs(:,1))+0.5 min(Rec.rs(:,2))-0.5 ...
        max(Rec.rs(:,2))+0.5 min(Rec.rs(:,3))-0.5 max(Rec.rs(:,3))+0.5])
    applyAxisProperties(gca)
end

end