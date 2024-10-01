function Data = dataHandling(Data)
%Data = dataHandling(Data) Processes NBI dataset
%   Input & Output:
%       - Data      : data structure.
%
% Author: Antonio Figueroa Durán
% Date: November 2022

%% PRE-PROCESSING
Data.Ref.pos = vertcat(Data.Line1.pos,Data.Line2.pos,Data.Line3.pos);
Data.Ref.h = horzcat(Data.Line1.h,Data.Line2.h,Data.Line3.h);

Data = rmfield(Data,{'Line1','Line2','Line3'});

% Recwin up to T = 1s
Data.Mic.h = Data.Mic.h(1:Data.Nsamples,:,:);
Data.Ref.h = Data.Ref.h(1:Data.Nsamples,:);

%% SPHERICAL ARRAYS
N = 4;              % Ambisonics order
kmax = N/Data.Mic.rArray;

Data.Mic.fNyq = Data.c*kmax/(2*pi);
Data.Mic.M = size(Data.Mic.pos,1);

%% SCATTERING
[Data.Mic.h,~] = removeScattering(Data.Mic.h,Data.Mic.Omega,Data.f*2*pi/Data.c);

%% LOW-PASS FILTER
% Filter design
% Data.Mic.BW = [300 Data.Mic.fNyq];
Data.Mic.BW = [300 4.5e3];
Fc = Data.Mic.BW;

hpFilt = designfilt('highpassfir','PassbandFrequency',Fc(1)*1.3, ...
            'StopbandFrequency',Fc(1),'PassbandRipple',0.5, ...
            'StopbandAttenuation',65,'DesignMethod','kaiserwin',...
            'SampleRate',Data.Fs);

lpFilt = designfilt('lowpassfir','PassbandFrequency',Fc(2), ...
            'StopbandFrequency',Fc(2)*1.1,'PassbandRipple',0.5, ...
            'StopbandAttenuation',65,'DesignMethod','kaiserwin',...
            'SampleRate',Data.Fs);

% Filtering
Data.Ref.h = filtfilt(lpFilt,Data.Ref.h); Data.Ref.h = filtfilt(hpFilt,Data.Ref.h);
Data.Mic.h = filtfilt(lpFilt,Data.Mic.h); Data.Mic.h = filtfilt(hpFilt,Data.Mic.h);

%% DOWNSAMPLING
Q = 1;   % Downsampling factor (reduce time samples)
Data.Fs = Data.Fs/Q;

Data.Ref.h = resample(Data.Ref.h,1,Q);
Data.Mic.h = resample(Data.Mic.h,1,Q);

Data.t = 0:1/Data.Fs:Data.T-1/Data.Fs;
Data.Nsamples = Data.Fs*Data.T;
Data.f = Data.Fs/Data.Nsamples*(0:Data.Nsamples/2-1);

%% REFERENCE POLARITY
Data.Ref.h = -Data.Ref.h;

%% NORMALISATION
% Normalisation wrt max of each dataset (em32 & reference set)
maxValue = max(abs(Data.Mic.h),[],'all');
Data.Ref.h = Data.Ref.h/max(abs(Data.Ref.h),[],'all');
Data.Mic.h = Data.Mic.h/maxValue;

% Frequency domain
for ii = 1:Data.Mic.nArrays
    Data.Mic.H(:,:,ii) = fftUniBi(Data.Mic.h(:,:,ii));
end
[Data.Ref.H,~] = fftUniBi(Data.Ref.h);

end

