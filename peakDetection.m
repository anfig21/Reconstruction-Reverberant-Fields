function [peaks,direct] = peakDetection(Res,t,Settings)
%[peaks,direct] = peakDetection(Data,t,Settings) Detect the main
%peaks for each individual spherical array using the frontmost microphone.
%The peak detection is performed over the squared RIR.
%   Input:
%       - Res       : data structure with RIR of dim [nSamples x nArray*32]
%       - t         : time vector. nSamples x 1
%       - Settings  : peak detection algorithm settings. Structure
%   Output:
%       - peaks     : detected peak delays. [nPeaks x nArray]
%       - direct    : direct sound TOA. [nArray x 1]
%
% Author: Antonio Figueroa-Duran
% Date: August 2023

%% ERROR HANDLING
if nargin < 3, error('peakDetection Error: Not enough input parameters.'), end

%% MAIN CODE
% Peak detection using the closest microphone to the center of the array
idxMic = 1;
nArrays = size(Res.h,3);

IR = squeeze(Res.h(:,idxMic,:));
ene = IR.^2;                % Energy

% Peak detection
direct = nan(nArrays,1);
peaksCell = cell(nArrays,1);
for ii = 1:nArrays
    peaks_ii = findpeaksx(t,ene(:,ii),Settings.SlopeThreshold, ...
        Settings.AmpThreshold,Settings.smoothwidth, ...
        Settings.peakgroup,Settings.smoothtype);

    % Direct sound delay: maximum energy
    direct(ii) = peaks_ii(peaks_ii(:,3) == max(peaks_ii(:,3)),2);

    % Only delay information
    peaksCell{ii} = peaks_ii(:,2);
end
numPeaks = max(cellfun('size',peaksCell,1));    % Max number of peaks
peaks = reshape(cell2mat(cellfun(@(x) [x; zeros(numPeaks-numel(x),1)],peaksCell,'uni',0)),numPeaks,nArrays);    % Reshape to [nPeaks x nEM32]

end