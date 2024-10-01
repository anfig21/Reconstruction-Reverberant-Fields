function [pInc,pSc] = removeScattering(pTotal,Omega,k)
%[pInc,pSc] = removeScattering(pTotal,k) Removes the scattering of the
%Eigenmike em32 using SH decomposition.
%   Input:
%       - pTotal    : raw measurements in time domain. Nt x 32 x Narrays
%       - Omega     : microphone coordinates (theta,phi). 32 x 2 x Narrays
%       - k         : wavenumber. 1 x nF
%   Output:
%       - pInc      : incident field. Nt x 32 x Narrays
%       - pSc       : scattered field. Nt x 32 x Narrays
% 
% Author: Antonio Figueroa-Duran
% Date: June 2024

%% ERROR HANDLING
arguments
    pTotal (:,32,:) double {mustBeNonempty}
    Omega (32,2,:) double {mustBeNonempty}
    k (1,:) double {mustBeNonempty,mustBeNonnegative}
end

%% MAIN CODE
% Dimensions
[Np,nMic,nArrays] = size(pTotal);

% DC component
if k(1) == 0, k(1) = k(2); end

% Array processing
N = 4;
n = (0:N)';

% Eigenmike radius
aSph = 4.2e-3;

% Radial functions: Bessel
j_ka = radFunction(k,aSph,N+1,'open').';

% Radial functions: Hankel
h_ka = radFunction(k,aSph,N+1,'closed').';

% Derivative of the Bessel function
j_ka = [zeros(1,Np/2); j_ka];               % Add -1'th order
dj_ka = (n.*j_ka(1:end-2,:)-(n+1).*j_ka(3:end,:))./(2*n+1);
j_ka([1 end],:) = [];                       % Delete -1'th and 5'th order

% Derivative of the Hankel function
h_ka = [zeros(1,Np/2); h_ka];               % Add -1'th order
dh_ka = (n.*h_ka(1:end-2,:)-(n+1).*h_ka(3:end,:))./(2*n+1);
h_ka([1 end],:) = [];                       % Delete -1'th and 5'th order

% Repmat to size up to (N+1)^2
j_kaN = nan((N+1)^2,Np/2);
dj_kaN = nan((N+1)^2,Np/2);
h_kaN = nan((N+1)^2,Np/2);
dh_kaN = nan((N+1)^2,Np/2);
for ii = 0:N
    j_kaN(ii^2+(1:2*ii+1),:) = repmat(j_ka(ii+1,:),2*ii+1,1);
    dj_kaN(ii^2+(1:2*ii+1),:) = repmat(dj_ka(ii+1,:),2*ii+1,1);

    h_kaN(ii^2+(1:2*ii+1),:) = repmat(h_ka(ii+1,:),2*ii+1,1);
    dh_kaN(ii^2+(1:2*ii+1),:) = repmat(dh_ka(ii+1,:),2*ii+1,1);
end

% ---- Spherical Fourier Transform ---- %
pInc = nan(Np,nMic,nArrays);
pSc = nan(Np,nMic,nArrays);
for ii = 1:nArrays
    % Spherical harmonics
    Y = sh(Omega(:,:,ii),N);

    p_ii = squeeze(pTotal(:,:,ii));

    % Frequency domain
    [P_ii,~] = fftUniBi(p_ii);

    % Total field SH coefficients
    c_nm = (4*pi/nMic)*pinv(Y)*P_ii.';

    % Incident field
    A = 1j*(k*aSph).^2.*dh_kaN.*c_nm;
    PInc_ii = Y*(j_kaN.*A);
    PInc_ii = PInc_ii.';

    % Scattered field
    B = -A.*dj_kaN./dh_kaN;
    PSc_ii = Y*(h_kaN.*B);
    PSc_ii = PSc_ii.';

    % Inc: Time domain
    PInc2 = [real(PInc_ii(1,:)); PInc_ii(2:end,:,:)/2];
    PInc2 = [PInc2; flip(conj(PInc2))];
    pInc(:,:,ii) = ifft(PInc2*Np,[],1,'symmetric');

    % Sc: Time domain
    PSc2 = [real(PSc_ii(1,:)); PSc_ii(2:end,:,:)/2];
    PSc2 = [PSc2; flip(conj(PSc2))];
    pSc(:,:,ii) = ifft(PSc2*Np,[],1,'symmetric');
end

end