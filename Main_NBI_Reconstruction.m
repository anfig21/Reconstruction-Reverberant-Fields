%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                  ROOM IMPULSE RESPONSE RECONSTRUCTION
%
% Dataset: Niels Bohr Institute (NBI) - classroom
%          Multiple Eigenmikes em32
%
% Antonio Figueroa-Duran
% anfig@dtu.dk
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear, close all

addpath(genpath('toolbox'))

% Download data from:
addpath(genpath('M:\Data\2023 NBI Data'))

% Matlab dependencies:
% - Statistics and Machine Learning Toolbox

loadPlotParams

%% INITIAL PARAMETERS
Data.T = 1;         % Pre-processing T = 1 s
Data.Fs = 48e3;
Data.t = 0:1/Data.Fs:Data.T-1/Data.Fs;
Data.D = [6.7920,9.4750,3.7050];       % Room dimensions (simplified)
Data.Nsamples = Data.Fs*Data.T;
Data.f = Data.Fs/Data.Nsamples*(0:Data.Nsamples/2-1);
Data.Temp = 24.3;
Data.Humi = 50.9;
Data.p0 = 100658;

[Data.rho,Data.c,~,Data.gamma,~,~,~,~,~,~] = amb2prop(Data.p0,Data.Temp,...
    Data.Humi,1);

% Select Eigenmikes
Data.Mic.IdxArrays = 2:3:8;
Data.Mic.nArrays = numel(Data.Mic.IdxArrays);

% Eigenmike specifications
Data.Mic.nM = 32;
Data.Mic.rArray = 0.042;         % em32 radius in m

% Window parameters (for Early Reflections extraction)
propagationTime = 2*Data.Mic.rArray/Data.c;    % Propagation through sphere in s
Peaks.wLength = 2^(ceil(log2(propagationTime*Data.Fs))+2);  % Propagation delay through sphere + margin -> 64 samples (default)

% Early reflections threshold
zeta = 10^(-30/20);

clear propagationTime

%% DATA ACQUISITION
Data = dataAcquisitionNBI(Data);

%% DATA HANDLING
Data = dataHandling(Data);

%% SETUP PLOT
% Flags
% - Setup
% - Frequency response
% - RIR reference line
setupPlot(Data,false,false,true);

%% RECONSTRUCTION
% Reference line: straight line to avoid artifact visualisation
Rec.recPosition = Data.Ref.pos;
Rec.recPosition(:,2) = mean(Rec.recPosition(:,2))*ones(size(Rec.recPosition,1),1);
Rec.recPosition(:,3) = mean(Rec.recPosition(:,3))*ones(size(Rec.recPosition,1),1);

%% ------------ LOCALISATION: DIRECT SOUND FIELD ------------ %%
% Peak detection: Settings
Peaks.Settings.SlopeThreshold = 0;
Peaks.Settings.AmpThreshold = max(Data.Mic.h.^2,[],'all')*zeta;
Peaks.Settings.smoothwidth = 0;
Peaks.Settings.peakgroup = 1;
Peaks.Settings.smoothtype = 1;

[~,Early.Peaks.direct] = peakDetection(Data.Mic,Data.t,Peaks.Settings);

% Windowing
Direct.T = [zeros(Data.Mic.nArrays,1) Early.Peaks.direct+1e-3];
Direct = windowRIR(Data,Direct.T(:,1),Direct.T(:,2));

% DOA Estimation
f = Data.f(2e3 <= Data.f & Data.f <= 5e3);
f = f(1:100:end);
N = 6e3;
Direct.DOA.Mode = earlyDOA(Data,Direct,f,N);

clear f N

%% Range Estimation (DOA-triangulation)
Direct.Range.pos = earlyRange(Data,Direct,true);

% Error
disp(['Position estimation off by ' num2str(norm(Direct.Range.pos-Data.Source.pos')) ' m'])

%% Adjust pre-delays (Mic & Ref)
% Based on the delay between onset and estimated Time-of-Flight
peak1_N = nan(Data.Mic.nArrays,1);
delay1_N = nan(Data.Mic.nArrays,1);

for ii = 1:Data.Mic.nArrays
    IR = squeeze(Data.Mic.h(:,:,ii));
    IR = IR./max(abs(IR),[],"all");
    ene = IR.^2;

    % Direct sound delay: maximum energy
    [~,idxPeak] = max(ene,[],'all');
    [peak1_N(ii),idxMic1] = ind2sub(size(ene),idxPeak);

    % Expected peak
    delay1_N(ii) = ceil(vecnorm(Direct.Range.pos'-squeeze(Data.Mic.pos(idxMic1,:,ii)))*Data.Fs/Data.c);
end

diff_N = peak1_N-delay1_N;

% Eigenmikes
for ii = 1:Data.Mic.nArrays
    if diff_N(ii) > 0
        Data.Mic.h(:,:,ii) = cat(1,Data.Mic.h(diff_N(ii):end,:,ii),zeros(diff_N(ii)-1,Data.Mic.nM,1));
        Direct.Mic.h(:,:,ii) = cat(1,Direct.Mic.h(diff_N(ii):end,:,ii),zeros(diff_N(ii)-1,Data.Mic.nM,1));
    else
        Data.Mic.h(:,:,ii) = cat(1,zeros(abs(diff_N(ii)),Data.Mic.nM,1),Data.Mic.h(1:end-abs(diff_N(ii)),:,ii));
        Direct.Mic.h(:,:,ii) = cat(1,zeros(abs(diff_N(ii)),Data.Mic.nM,1),Direct.Mic.h(1:end-abs(diff_N(ii)),:,ii));
    end
    [Data.Mic.H(:,:,ii),~] = fftUniBi(Data.Mic.h(:,:,ii));
    [Direct.Mic.H(:,:,ii),~] = fftUniBi(Direct.Mic.h(:,:,ii));
end

% Reference
delayRef_N = ceil(vecnorm(Direct.Range.pos-Data.Ref.pos.',2,1)*Data.Fs/Data.c);
[~,peakRef_N] = max(abs(Data.Ref.h));
diffRef_N = round(mode(peakRef_N-delayRef_N));

if diffRef_N > 0
    Data.Ref.h = cat(1,Data.Ref.h(diffRef_N:end,:),zeros(diffRef_N-1,size(Data.Ref.h,2)));
else
    Data.Ref.h = cat(1,zeros(abs(diffRef_N),size(Data.Ref.h,2)),Data.Ref.h(1:end-abs(diffRef_N),:));
end
[Data.Ref.H,~] = fftUniBi(Data.Ref.h);

clear peak1_N delay1_N diff_N ii delayRef_N peakRef_N diffRef_N ene IR idxPeak idxMic1

%% ------------ RECONSTRUCTION: DIRECT SOUND ------------ %%
% Window properties - tukey window
peak_Samp = ceil(Data.Fs*vecnorm(Direct.Range.pos.'-Data.Mic.R0,2,2)/Data.c);   % Peak wrt array center
delay = peak_Samp+[-Peaks.wLength Peaks.wLength]/2;
Peaks.w = tukeywin(Peaks.wLength);
Peaks.wInv = 1-Peaks.w;

% Frequency span of the reconstruction
Rec.Direct.f = Data.f(Data.Mic.BW(1) <= Data.f & Data.f <= Data.Mic.BW(2)*1.1);

% Memory allocation
nIdx = nan(Data.Mic.nArrays,Peaks.wLength);  % Index for windowed samples
Rec.Direct.Mic.h = zeros(Data.Nsamples,Data.Mic.nM,Data.Mic.nArrays);
Rec.Direct.Mic.H = nan(size(Rec.Direct.Mic.h,1)/2,size(Rec.Direct.Mic.h,2),size(Rec.Direct.Mic.h,3));

% Window direct part - Tukey window
for ii = 1:Data.Mic.nArrays
    nIdx(ii,:) = delay(ii,1):delay(ii,2)-1;
    Rec.Direct.Mic.h(nIdx(ii,:),:,ii) = Peaks.w.*Data.Mic.h(nIdx(ii,:),:,ii);

    % Frequency domain
    [Rec.Direct.Mic.H(:,:,ii),~] = fftUniBi(Rec.Direct.Mic.h(:,:,ii));
end

% Coefficient estimation
Rec.Direct = clusterCoefficients(Data,Rec.Direct,Direct.Range.pos);

% Residual in time domain
Hdiff = Rec.Direct.Mic.H - Rec.Direct.Rec.H;
hDiff = nan(size(Data.Mic.h));
for ii = 1:Data.Mic.nArrays
    P2 = [real(Hdiff(1,:,ii)); Hdiff(2:end,:,ii)/2];
    P2 = [P2; flip(conj(P2))];
    hDiff(:,:,ii) = ifft(P2*Data.Nsamples,[],1,'symmetric');

    % Post-window: ensures zero (low) values at edges
    hDiff(:,:,ii) = cat(1,zeros(nIdx(ii,1)-1,Data.Mic.nM), ...
        Peaks.w.*hDiff(nIdx(ii,:),:,ii), ...
        zeros(Data.Nsamples-nIdx(ii,end),Data.Mic.nM));
end

%% Synthetise residual @ mic positions
Rec.Res.Mic.h = Data.Mic.h;     % Initialise residual
for ii = 1:Data.Mic.nArrays
    % Inverse Tukey window
    Rec.Res.Mic.h(nIdx(ii,:),:,ii) = Peaks.wInv.*Rec.Res.Mic.h(nIdx(ii,:),:,ii);
    Rec.Res.Mic.h(nIdx(ii,:),:,ii) = Rec.Res.Mic.h(nIdx(ii,:),:,ii)+Peaks.w.*hDiff(nIdx(ii,:),:,ii);

    % Frequency domain
    [Rec.Res.Mic.H(:,:,ii),~] = fftUniBi(Rec.Res.Mic.h(:,:,ii));
end

%% Extrapolation
Rec.Direct = reconstructReflection(Data,Rec.Direct,Rec.recPosition);

% Clear workspace
clear P2 ii delay nIdx peak_Samp hDiff Hdiff

%% ------------ ITERATIVE RECONSTRUCTION: EARLY REFLECTIONS ------------ %%
% Directional Decomposition: DOA & Range
f = Data.f(2e3 <= Data.f & Data.f <= 4e3);
f = f(1:50:end);
N = 4e3;

% Peak detection algorithm with mic 1 (front-most)
[Early.Peaks.peaks,~] = peakDetection(Rec.Res.Mic,Data.t,Peaks.Settings);
Early.Peaks.peaks(Early.Peaks.peaks <= max(Early.Peaks.direct)+1e-3) = nan;
Early.totalPeaks = 0;
Early.Peaks.modelled = [];

Rec.Early.h = zeros(Data.Nsamples,size(Rec.recPosition,1));
while any(Early.Peaks.peaks>0,'all') %& min(Early.Peaks.peaks,[],'all') < 60e-3
    % ----- DOA estimation ----- %
    Early.totalPeaks = Early.totalPeaks + 1;
    [peak_ii,micIdx] = min(Early.Peaks.peaks,[],'all');
    Early.Peaks.modelled = [Early.Peaks.modelled peak_ii];  % Backlog
    [~,micIdx] = ind2sub(size(Early.Peaks.peaks),micIdx);

    disp(['Reflection number ' num2str(Early.totalPeaks) ' at ' num2str(peak_ii*1e3) ' ms from array ' num2str(micIdx)])

    % Leading microphone: windowing
    delay_ii = ceil(peak_ii*Data.Fs)+[-Peaks.wLength Peaks.wLength]/2;
    nIdx_ii = delay_ii(1):delay_ii(2)-1;
    rir_ii.Mic.h = zeros(Data.Nsamples,Data.Mic.nM);
    rir_ii.Mic.h(nIdx_ii,:) = Peaks.w.*Rec.Res.Mic.h(nIdx_ii,:,micIdx);
    [rir_ii.Mic.H,~] = fftUniBi(rir_ii.Mic.h);

    % DOA estimation
    rir_ii.DOA.Mode = earlyDOA(Data,rir_ii,f,N,[],micIdx);

    % ----- Range estimation ----- %
    delay = peak_ii - Early.Peaks.direct(micIdx);
    range = Data.c*delay + vecnorm(Direct.Range.pos.'-Data.Mic.R0(micIdx,:));
    rir_ii.pos = rir_ii.DOA.Mode.*range + Data.Mic.R0(micIdx,:);

    % ----- Reconstruction  ----- %
    % Frequency span of the reconstruction
    Early_ii.f = Data.f(Data.Mic.BW(1) <= Data.f & Data.f <= Data.Mic.BW(2)*1.1);

    % Define delay
    peak_Samp = ceil(Data.Fs*vecnorm(rir_ii.pos-Data.Mic.R0,2,2)/Data.c);
    delay = peak_Samp+[-Peaks.wLength Peaks.wLength]/2;
    nIdx = nan(Data.Mic.nArrays,Peaks.wLength);
    for ii = 1:Data.Mic.nArrays
        % Window RIR
        nIdx(ii,:) = delay(ii,1):delay(ii,2)-1;
        Early_ii.Mic.h(:,:,ii) = cat(1,zeros(nIdx(ii,1)-1,Data.Mic.nM), ...
            Peaks.w.*Rec.Res.Mic.h(nIdx(ii,:),:,ii), ...
            zeros(Data.Nsamples-nIdx(ii,end),Data.Mic.nM));

        % Frequency domain
        [Early_ii.Mic.H(:,:,ii),~] = fftUniBi(Early_ii.Mic.h(:,:,ii));
    end

    % Coefficient estimation
    Early_ii = clusterCoefficients(Data,Early_ii,rir_ii.pos);

    % Residual in time domain
    Hdiff = Early_ii.Mic.H - Early_ii.Rec.H;
    hDiff = nan(size(Data.Mic.h));
    for ii = 1:Data.Mic.nArrays
        P2 = [real(Hdiff(1,:,ii)); Hdiff(2:end,:,ii)/2];
        P2 = [P2; flip(conj(P2))];
        hDiff(:,:,ii) = ifft(P2*Data.Nsamples,[],1,'symmetric');

        % Post-window: ensures zero (low) values at edges
        hDiff(:,:,ii) = cat(1,zeros(nIdx(ii,1)-1,Data.Mic.nM), ...
            Peaks.w.*hDiff(nIdx(ii,:),:,ii), ...
            zeros(Data.Nsamples-nIdx(ii,end),Data.Mic.nM));
    end

    % Synthetise residual
    for ii = 1:Data.Mic.nArrays
        % Inverse Tukey window
        Rec.Res.Mic.h(nIdx(ii,:),:,ii) = Peaks.wInv.*Rec.Res.Mic.h(nIdx(ii,:),:,ii);
        Rec.Res.Mic.h(nIdx(ii,:),:,ii) = Rec.Res.Mic.h(nIdx(ii,:),:,ii)+Peaks.w.*hDiff(nIdx(ii,:),:,ii);
    
        % Frequency domain
        [Rec.Res.Mic.H(:,:,ii),~] = fftUniBi(Rec.Res.Mic.h(:,:,ii));
    end

    % Reconstruction
    Early_ii = reconstructReflection(Data,Early_ii,Rec.recPosition);

    Rec.Early.h = Rec.Early.h+Early_ii.h;

    % Peak detection algorithm with mic 1 (front-most)
    [Early.Peaks.peaks,~] = peakDetection(Rec.Res.Mic,Data.t,Peaks.Settings);

    % Discard already-modelled peaks
    Early.Peaks.peaks(Early.Peaks.peaks <= max(Early.Peaks.modelled)) = nan;
end

% Clear workspace
clear ii rir_ii delay f P2 Early_ii range nIdx peak_Samp micIdx peak_ii

%% ------------ KERNEL RIDGE REGRESSION ------------ %%
% Frequency span of the reconstruction
Rec.Kernel.f = Data.f(Data.Mic.BW(1) <= Data.f & Data.f <= Data.Mic.BW(2)*1.1);

% Generate window
Rec.Kernel.T = [0 1];
Rec.Kernel.w = createWindow(Data,Rec.Kernel.T(1),Rec.Kernel.T(2));

for ii = 1:Data.Mic.nArrays
    % Window residual
    Rec.Kernel.Mic.h(:,:,ii) = Rec.Kernel.w.*Rec.Res.Mic.h(:,:,ii);

    % Frequency domain
    [Rec.Kernel.Mic.H(:,:,ii),~] = fftUniBi(Rec.Kernel.Mic.h(:,:,ii));
end

% Reconstruction via Representer Theorem
Rec.Kernel = kernelReconstructionOverlap(Data,Rec.Kernel,Rec.recPosition);

%% KERNEL: FADE-IN
% In accordance to causality property
Rec.Kernel.hFI = Rec.Kernel.h;

[~,peakIdx] = max(abs(Rec.Direct.h));
Rec.Kernel.hFI(repmat(Data.t',1,numel(peakIdx)) < Data.t(peakIdx-1)) = 0;

%% ------------ MERGED ROOM IMPULSE RESPONSE ------------ %%
% Direct sound + Kernel
Rec.Total.h = Rec.Direct.h+Rec.Early.h+Rec.Kernel.hFI;

% Frequency domain
[Rec.Direct.H,~] = fftUniBi(Rec.Direct.h);
[Rec.Early.H,~] = fftUniBi(Rec.Early.h);
[Rec.Kernel.H,~] = fftUniBi(Rec.Kernel.hFI);
[Rec.Total.H,~] = fftUniBi(Rec.Total.h);

%% SAVE RIR AS WAV
% Rec.Total.h = Rec.Total.h/max(abs(Rec.Total.h),[],'all');
% audiowrite('rirReconstructed.wav',Rec.Total.h,Data.Fs);

%% ------------ PLOT: MERGED ROOM IMPULSE RESPONSE ------------ %%
% REFERENCE LINE
T = [15 40]*1e-3;

dataPlot = Rec.Total.h(Data.t>T(1) & Data.t<T(2),:);
timePlot = Data.t(Data.t>T(1) & Data.t<T(2));

figure, hold on
s = pcolor(Data.Ref.pos(:,1),timePlot*1e3,dataPlot);
set(s,'edgecolor','none')
xlabel('$x$ / m'), ylabel('Time / ms')
axis([0 Data.D(1) T(1)*1e3 T(2)*1e3])
DR = colormapDR(dataPlot);
colorbarpwn(-DR,DR)
c = colorbar;
for ii = 1:Data.Mic.nArrays
    xline(min(Data.Mic.pos(:,1,ii))), xline(max(Data.Mic.pos(:,1,ii)))
end
clim([-1 1])
applyColorbarProperties(c,'Normalised $h(t)$')
applyAxisProperties(gca)

% Clear workspace
clear s

%% DIFFERENT CONTRIBUTIONS
figure
subplot(131)
dataPlot = Rec.Direct.h(Data.t>T(1) & Data.t<T(2),:);
s = pcolor(Data.Ref.pos(:,1),timePlot*1e3,dataPlot);
set(s,'edgecolor','none')
xlabel('$x$ / m'), ylabel('Time / ms'), title('Direct sound')
axis([0 Data.D(1) T(1)*1e3 T(2)*1e3])
colorbarpwn(-DR,DR)
colorbar
for ii = 1:Data.Mic.nArrays
    xline(min(Data.Mic.pos(:,1,ii))), xline(max(Data.Mic.pos(:,1,ii)))
end
applyAxisProperties(gca)
clim([-1 1])

subplot(132)
dataPlot = Rec.Early.h(Data.t>T(1) & Data.t<T(2),:);
s = pcolor(Data.Ref.pos(:,1),timePlot*1e3,dataPlot);
set(s,'edgecolor','none')
xlabel('$x$ / m'), title('Early reflections')
axis([0 Data.D(1) T(1)*1e3 T(2)*1e3])
colorbarpwn(-DR,DR)
colorbar
for ii = 1:Data.Mic.nArrays
    xline(min(Data.Mic.pos(:,1,ii))), xline(max(Data.Mic.pos(:,1,ii)))
end
applyAxisProperties(gca)
clim([-1 1])

subplot(133)
dataPlot = Rec.Kernel.hFI(Data.t>T(1) & Data.t<T(2),:);
s = pcolor(Data.Ref.pos(:,1),timePlot*1e3,dataPlot);
set(s,'edgecolor','none')
xlabel('$x$ / m'), title('Kernel processing')
axis([0 Data.D(1) T(1)*1e3 T(2)*1e3])
colorbarpwn(-DR,DR)
c = colorbar;
for ii = 1:Data.Mic.nArrays
    xline(min(Data.Mic.pos(:,1,ii))), xline(max(Data.Mic.pos(:,1,ii)))
end
applyColorbarProperties(c,'Normalised $h(t)$')
applyAxisProperties(gca)
clim([-1 1])
