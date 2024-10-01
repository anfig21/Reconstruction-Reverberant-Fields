function Structure = windowRIR(Data,Tini,Tfin,plotFlag)
%Structure = windowRIR(Data,Tini,Tfin,plotFlag) Applies Hanning window to
%RIR between time given by Tini and Tfin. Obtains the frequency response
%and the noise spectrum. Estimates the noise power after windowing.
%   Input:
%       - Data      : raw data. Structure
%       - Tini      : initial time. Scalar
%       - Tfin      : end time. Scalar
%       - plotFlag  : 'true'  - RIR in time domain
%                     'false' - Default value
%   Output:
%       - Structure : structure with windowed data
%
% Author: Antonio Figueroa Durán
% Date: August 2023

%% ERROR HANDLING
% plotFlag default value
if nargin < 4, plotFlag = false;
elseif nargin < 3, error('windowRIR Error: Not enough input parameters.'), end

%% MAIN CODE
Structure.T = [Tini Tfin];
Structure.N = floor(Data.Fs*Structure.T);

% Windowing
Structure.Mic.h = nan(size(Data.Mic.h));
for ii = 1:Data.Mic.nArrays
    w = createWindow(Data,Tini(ii),Tfin(ii));

    Structure.Mic.h(:,:,ii) = repmat(w,1,32).*Data.Mic.h(:,:,ii);
    Structure.Mic.H(:,:,ii) = fftUniBi(Structure.Mic.h(:,:,ii));
end

%% PLOT RIR
if plotFlag
    figure
    subplot(211), plot((0:Data.Nsamples-1)*1e3/Data.Fs,Data.Mic.h), grid on
    xlim([min(Tini)*1e3-1 max(Tfin)*1e3+1]), ylabel('$h(t)$')
    applyAxisProperties(gca)
    subplot(212), plot((0:Data.Nsamples-1)*1e3/Data.Fs,Structure.Mic.h), grid on
    xlim([min(Tini)*1e3-1 max(Tfin)*1e3+1]), ylabel('$h_w(t)$'), xlabel('$t$/ms')
    applyAxisProperties(gca)
end

end

