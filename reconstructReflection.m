function Rec = reconstructReflection(Data,Rec,r,plotFlag)
%Rec = reconstructReflection(Data,Rec,r,plotFlag) Reconstruct the RIR
%contained in Rec at points r using a point source propagation model.
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
elseif nargin < 3, error('reconstructReflection Error: Not enough input parameters.'), end

%% RECONSTRUCTION
% Dimensions
Nf = length(Rec.f);
R = size(r,1);

k = (2*pi*Rec.f)/Data.c;        % Propagation vector
d = pdist2(Rec.rs,r);           % mics x point sources

% Dictionary
Rec.P = zeros(R,length(Data.f));
for ii = 1:Nf
    H = (1./d).*exp(-1i*d*k(ii));       % Dictionary (point sources)
    Rec.P(:,Data.f==Rec.f(ii)) = H.'*Rec.x(:,Data.f==Rec.f(ii));
end

% Double-sided spectrum
P2 = [real(Rec.P(:,1)) Rec.P(:,2:end)/2];
P2 = [P2 flip(conj(P2),2)];
Rec.h = ifft(P2*Data.Nsamples,[],2,'symmetric').';

%% PLOT: REFERENCE RIR
if plotFlag
    T = [17 40]*1e-3;       % Source near field

    figure, hold on
    s = pcolor(r(:,1),Data.t*1e3,Rec.h);
    set(s,'edgecolor','none')
    xlabel('$x$ / m'), ylabel('Time / ms')
    axis([0 Data.D(1) T(1)*1e3 T(2)*1e3])
    DR = colormapDR(Rec.h);
    colorbarpwn(-DR,DR)
    c = colorbar;
    for ii = 1:Data.Mic.nArrays
        xline(min(Data.Mic.pos(:,1,ii))), xline(max(Data.Mic.pos(:,1,ii)))
    end
    applyColorbarProperties(c,'Normalised $h(t)$')
    applyAxisProperties(gca)
end

end

