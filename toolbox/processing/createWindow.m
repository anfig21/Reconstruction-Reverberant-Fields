function w = createWindow(Data,Tini,Tfin)
%w = createWindow(Data,Tini,Tfin) Generate window using half Hanning window
%of total length 0.5 ms (0.25 ms on each side).
%   Input:
%       - Data      : raw data. Structure
%       - Tini      : initial time. Scalar
%       - Tfin      : end time. Scalar
%   Output:
%       - w         : window. N x 1
%
% Author: Antonio Figueroa-Duran
% Date: June 2023

%% ERROR HANDLING
if nargin < 3, error('createWindow Error: Not enough input parameters.'), end

%% MAIN CODE
tHann = 0.25e-3;
NHann = ceil(Data.Fs*tHann/2);
h = hann(2*NHann);
hHalf1 = h(1:end/2);
hHalf2 = h(end/2+1:end);

T = [Tini Tfin];
N = floor(Data.Fs*T);
N(N==0) = 1;
Nsamples = max(N(:,2)-N(:,1));

% Windowing - 1/2 Hanning (hann) window on each side
w = vertcat(zeros(N(1)-NHann,1),hHalf1,ones(Nsamples,1),hHalf2,...
    zeros(Data.Nsamples-N(2)-NHann,1));
w = w(1:Data.Nsamples);
end

