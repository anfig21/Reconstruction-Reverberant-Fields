function Rec = kernelReconstructionOverlap(Data,Rec,r,plotFlag)
%Rec = kernelReconstructionOverlap(Data,Rec,r,plotFlag) Reconstruct the RIR
%contained in Rec at points r using the Bessel function as kernel. The
%reconstruction is performed using the Representer Theorem, with aid of
%Anchor Points (synthesised sound field samples based on data). The
%time-frequency processing is done using overlapping windows.
%   Input:
%       - Data          : raw data. Structure
%       - Rec           : windowed RIR section to reconstruct. Structure
%       - r             : reconstruction points. R x 3
%       - plotFlag      : 'true' to plot reconstructed RIR
%                         'false' (Default value)
%   Output:
%       - Rec           : reconstructed RIR. Structure
%
% Author: Antonio Figueroa Durán
% Date: June 2023

%% ERROR HANDLING
if nargin < 4, plotFlag = false;    % plotFlag default value
elseif nargin < 3, error('kernelReconstructionOverlap Error: Not enough input parameters.'), end

%% MAIN CODE
% Dimensions
Nr = size(r,1);

% Cut-off frequency: beyond this freq correlation is not relevant
fCut = Data.c/min(pdist(r),[],'all');

%% OVERLAP PARAMETERS
overlap = 0.5;     % Overlap rate [0,1)
nLim = ceil(Rec.T*Data.Fs);
nLim(nLim == 0) = 1;

% Stack arrays data
h = reshape(Rec.Mic.h(nLim(1):nLim(2),:,:),[],Data.Mic.nM*Data.Mic.nArrays);
pos = reshape(permute(Data.Mic.pos,[2 1 3]),[],Data.Mic.nM*Data.Mic.nArrays)';

M = 2^8;                    % samples per window (256 samples ~ 5.3 ms)
L = ceil(M*overlap);
[Nx,Nmic] = size(h);

% Perfect recovery condition - kSTFT must be an integer
kSTFT = (Nx-L)/(M-L);
Ny = (ceil(kSTFT)+1)*(M-L)+L;   % Add 1 frame to avoid non-causal artifacts
h = vertcat(h,zeros(Ny-Nx,Nmic));

% Segmentation & windowing
[hSTFT,fSTFT] = stft(h,Data.Fs,'Window',hanning(M), ...
                            'OverlapLength',L, ...
                            'FFTLength',M, ...
                            'FrequencyRange','onesided');

[Nf,nFrames,~] = size(hSTFT);

% Reconstruction frequencies
fRec = fSTFT(fSTFT >= Rec.f(1) & fSTFT <= Rec.f(end));
k = 2*pi*fRec/Data.c;      % Wavenumber

%% SYNTHESISED PRESSURE POINTS
dMin = Data.c/min([fCut max(fRec)]);    % 1 wavelength

% Define Pressure Points coverage
sppSpace = [min(r)-dMin/2; max(r)+dMin/2];
sppDist = abs(sppSpace(2,:)-sppSpace(1,:));

%% RECONSTRUCTION
rng('default')
HRec = zeros(Nf,nFrames,Nr);
for mm = 1:nFrames
    disp(['Frame number ' num2str(mm) '/' num2str(nFrames)])

    % Frame m: initialise
    Hm = squeeze(hSTFT(:,mm,:));
    sppPos = [];
    sppN = 0;

    % Regularisation
    for ii = 1:length(fRec)
        % Input data
        Hmi = Hm(fSTFT == fRec(ii),:);

        if fRec(ii) <= fCut
            % Synthetised pressure points
            sppNf = ceil(sppDist/dMin);

            if sum(sppNf) > 0
                dMinSPP = sppDist./sppNf;
                xx = dMinSPP(1)*(0.5+(0:sppNf(1)-1)); if isempty(xx), xx = 0; end
                yy = dMinSPP(2)*(0.5+(0:sppNf(2)-1)); if isempty(yy), yy = 0; end
                zz = dMinSPP(3)*(0.5+(0:sppNf(3)-1)); if isempty(zz), zz = 0; end

                [XX,YY,ZZ] = meshgrid(xx,yy,zz);
                sppPos = [XX(:) YY(:) ZZ(:)]+sppSpace(1,:);
                
                sppN = size(sppPos,1);
            end

            % Gaussian distribution parameter estimation
            [meanReal,stdReal] = normfit(real(Hmi(:)));
            [meanImag,stdImag] = normfit(imag(Hmi(:)));

            meanReal = meanReal*ones(sppN,1);
            stdReal = 2*stdReal*ones(sppN,1);
            meanImag = meanImag*ones(sppN,1);
            stdImag = 2*stdImag*ones(sppN,1);

            % Anchor pressure (synthetised)
            hVirtual = normrnd(meanReal,stdReal)+1j*normrnd(meanImag,stdImag);
    
            % Update input data
            Hmi = vertcat(Hmi(:),hVirtual);
    
            % Euclidean distance matrices
            d_MNr = pdist2(r,[pos; sppPos])';
            d_MM = squareform(pdist([pos; sppPos]));

            % Bessel Kernel
            K_MM = var(Hmi)*sinc(k(ii)*d_MM/pi);
            K_MNr = var(Hmi)*sinc(k(ii)*d_MNr/pi);
    
            % Parameter estimation: Representer Theorem
            [U,s,V] = csvd(K_MM);
            [regParameter,~,~] = gcv(U,s,Hmi);
            [x,~,~] = tikhonov(U,s,V,Hmi,regParameter);
            
            x(isnan(x)) = 0;
            
            % Reconstruction
            HRec(fSTFT == fRec(ii),mm,:) = x.'*K_MNr;
        else
            pdReal = fitdist(real(Hmi(:)),'Normal');
            pdImag = fitdist(imag(Hmi(:)),'Normal');

            HRec(fSTFT == fRec(ii),mm,:) = sqrt(2)*(random(pdReal,Nr,1)+1j*random(pdImag,Nr,1));
        end
    end

    % GAIN ADJUSTMENT
    total_m = sum(abs(Hm),'all');
    total_r = sum(abs(squeeze(HRec(:,mm,:))),'all');
    G = (Nr/Nmic)*(total_m/total_r);
    HRec(:,mm,:) = G*HRec(:,mm,:);
end

% Time domain
hRec = istft(HRec,Data.Fs,'Window',hanning(M), ...
                            'OverlapLength',L, ...
                            'FFTLength',M, ...
                            'ConjugateSymmetric',true, ...
                            'FrequencyRange','onesided');
hRec = hRec(1:Nx,:);

% Total signal length
Rec.h = zeros(size(Rec.Mic.h,1),Nr);
Rec.h(nLim(1):nLim(2),:) = hRec;

% Frequency domain
[Rec.H,~] = fftUniBi(Rec.h);


%% PLOT: REFERENCE RIR
if plotFlag
    T = [5 50]*1e-3;
    t = Data.t(Data.t >= T(1) & Data.t <= T(2));

    figure, hold on
    s = pcolor(r(:,1),t*1e3,Rec.h(ismember(Data.t,t),:));
    set(s,'edgecolor','none')
    xlabel('$x$ / m'), ylabel('Time / ms')
    xlim([0 Data.D(1)])
    DR = colormapDR(Rec.h(ismember(Data.t,t),:));
    colorbarpwn(-DR,DR)
    c = colorbar;
    xline(min(Data.SphL.pos(:,1))), xline(max(Data.SphL.pos(:,1)))
    xline(min(Data.SphR.pos(:,1))), xline(max(Data.SphR.pos(:,1)))
    applyColorbarProperties(c,'Normalised $h(t)$')
    applyAxisProperties(gca)
end

end