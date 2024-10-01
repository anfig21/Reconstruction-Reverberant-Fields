function [] = plotFreqResponse(Data)
%plotFreqResponse(Data) Plots the measurement setup, the magnitude of the
%transfer function and the frequency response of signal and noise.
%   Input:
%       - Data      : data structure.
%
% Author: Antonio Figueroa Durán
% Date: November 2022

mic = 1;

% Frequency response
figure, plot(Data.f*1e-3,20*log10(abs(Data.Mic.H(:,mic)))), grid on
xlabel('$f$/kHz'), ylabel('$|H(j\omega)|$/dB')
legend('Mic 1')
applyAxisProperties(gca)
applyLegendProperties(gcf)

end

