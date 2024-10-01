function DOA = dirDOA_SRP_PHAT(P,pos,f,Fs,c,N)
%DOA = dirDOA_SRP_PHAT(P,pos,f,Fs,c,N) Estimates the DOA for the given
%frequency-domain signal P using SRP-PHAT.
%   Input:
%       - P     : frequency response at the specific frequency bins. Nf x M
%       - pos   : microphone positions. M x 3
%       - f     : frequency vector. 1 x Nf
%       - Fs    : sampling frequency. Scalar
%       - c     : speed of sound in the medium. Scalar
%       - N     : number of looking directions on the 2D sphere. Scalar
%   Output:
%       - DOA   : DOA estimation via SRP-PHAT. 3 x 1
%
% Author: Antonio Figueroa Durán
% Date: March 2024

%% ERROR HANDLING
arguments
    P (:,:) double
    pos (:,3) double
    f (1,:) double
    Fs (1,1) double {mustBeNonnegative}
    c (1,1) double {mustBePositive}
    N (1,1) double {mustBePositive} = 1e3
end

%% MAIN CODE
% Looking directions
uk = fibonacciSampling(N)';
% resRad = deg2rad(res);
% 
% % Steering coordinates
% Ntheta = ceil(pi/resRad);   theta = linspace(0,pi,Ntheta);
% Nphi = ceil(2*pi/resRad);   phi = linspace(0,(1-1/Nphi)*2*pi,Nphi);
% 
% % Steering vector
% [Utheta,Uphi] = meshgrid(theta,phi);
% UkSph = [ones(size(Utheta(:),1),1) Utheta(:) Uphi(:)];
% [uk(:,1), uk(:,2), uk(:,3)] = sph2cart2(UkSph(:,1),UkSph(:,2),UkSph(:,3));

% Dimensions
[Nf,M] = size(P);
N = size(uk,1);
w = 2*pi*f;
pos = pos-mean(pos);

% Delay
distMat = permute(pdist2(uk,pos)',[1 3 2]);
tau_kl = (repmat(permute(distMat,[2 1 3]),M,1,1)-repmat(distMat,1,M,1))/c;
mask = tril(ones(M,M),-1);
tau_kl = mask.*tau_kl;

% GCC-PHAT
Phi = nan(Nf,M,M,N);
for ff = 1:Nf
    M_kl = tril(repmat(P(ff,:),M,1).*repmat(P(ff,:)',1,M),-1);  % GCC
    Phi_kl = M_kl.*exp(1j*w(ff)*tau_kl);    % Delay
    Phi(ff,:,:,:) = Phi_kl./abs(Phi_kl);    % Phase Transform
end
Phi(isnan(Phi)) = 0;
clear Phi_kl M_kl tau_kl mask distMat

% Reconstruct spectrum
fDelta = f(2)-f(1);
fTot = [0:fDelta:f(1) f(2:end-1) f(end):fDelta:Fs/2-fDelta];

srp = zeros(length(fTot)*2,N);
for k = 1:M
    for l = k+1:M
        PhiTot = zeros(length(fTot),N);
        PhiTot(ismember(fTot,f),:) = squeeze(Phi(:,l,k,:));
        P2 = cat(1,real(PhiTot(1,:)),PhiTot(2:end,:)/2);
        P2 = cat(1,P2,flip(conj(P2),1));
        srp = srp + ifft(P2*2*length(fTot),'symmetric');
    end
end

srp_map = sum(abs(srp).^2);

% Plot
% figure, surf(rad2deg(Utheta),rad2deg(Uphi),reshape(srp_map,size(Utheta))/max(srp_map,[],'all')), view(2)
% xlabel('theta'),ylabel('phi'), colorbar

[~,idx] = max(srp_map);
DOA = uk(idx,:);

end

